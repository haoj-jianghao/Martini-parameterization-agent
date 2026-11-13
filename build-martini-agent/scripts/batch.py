#!/usr/bin/env python3
"""Run a batch as one agent session per molecule (SKILL.md "Rules" and "Batch mode").

Each molecule gets its own folder and its own headless agent session, started with the prompt in
assets/agent_prompt.md and the environment variable MARTINI_MOLDIR set to that folder; design.py and cgbuild.py
refuse any other folder in that session. No session sees another molecule, so no decision can be shared or
hard-coded across molecules. Deterministic; makes no chemical decisions.

  batch.py run    --molecules FILE.tsv --batch NAME --agent claude|codex|<name> --model MODEL
                  [--jobs 4] [--gmx-threads 2] [--timeout 5400] [--budget-usd 10] [--dry-run] [--skip-preflight]
  batch.py prompt --molecules FILE.tsv --id ID [--batch-dir DIR]

FILE.tsv: tab-separated with header; columns id, name, smiles, aa_topology (path to the user's OPLS-AA
topology, relative to the TSV), optional aa_engine and aa_charge_model.
Agent commands: environment.json `agent_commands` {name: [argv template]} overrides the built-ins below.
Template fields: {prompt} {model} {moldir} {project} {budget}.

Before any agent session starts, every molecule folder is prepared with its agent_prompt.txt and
agent_env.sh. Then a preflight checks, without any agent or API call, that the environment
the sessions inherit can write the agent's state folder and reach its service (environment.json
`agent_preflight` {name: {"writable": [...], "connect": [...]}} overrides the built-ins; an entry may list
alternatives separated by "|", the first one whose variables are all set is used; "connect" takes host:port or
a URL). If the preflight fails, or a session ends without doing any molecule work (no input/request.json and
no design log), no further sessions are started, batch.json records launch status "failed", and FALLBACK.md
in the batch folder says how to finish: run the same command where the checks pass, or start one fresh
sub-agent per remaining molecule with its prepared prompt and environment (SKILL.md "Batch mode").
Exit code 3 means the launch failed.

Writes runs/<NAME>-<UTC>/: one folder per molecule, logs/<ID>.out/.err, batch.json (launch record),
FALLBACK.md if the launch failed, then `cgbuild.py summary`.
"""
import argparse, csv, datetime, hashlib, json, os, re, shlex, socket, subprocess, sys, tempfile, threading, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlparse

SCRIPT_DIR = Path(__file__).resolve().parent
SKILL = SCRIPT_DIR.parent
ENV_PATH = Path(os.environ.get('MARTINI_ENV', SKILL.parent / 'environment.json')).resolve()
PROJECT = ENV_PATH.parent
ENV = json.loads(ENV_PATH.read_text())

BUILTIN = {
    'claude': ['claude', '-p', '{prompt}', '--model', '{model}', '--output-format', 'json',
               '--allowedTools', 'Bash Read Write Edit Glob Grep', '--disallowedTools', 'Agent Task WebFetch WebSearch',
               '--add-dir', '{project}', '--max-budget-usd', '{budget}'],
    'codex': ['codex', 'exec', '--model', '{model}', '--sandbox', 'workspace-write', '--skip-git-repo-check',
              '-c', 'shell_environment_policy.inherit=all', '-C', '{moldir}', '{prompt}'],
}
# What each agent CLI needs from the environment it inherits: a writable state folder and its service.
BUILTIN_PREFLIGHT = {
    'claude': {'writable': ['$CLAUDE_CONFIG_DIR|$HOME/.claude'], 'connect': ['$ANTHROPIC_BASE_URL|api.anthropic.com:443']},
    'codex': {'writable': ['$CODEX_HOME|$HOME/.codex'], 'connect': ['$OPENAI_BASE_URL|api.openai.com:443']},
}
LAUNCH_FAILED_EXIT = 3


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')


def read_tsv(path):
    path = Path(path).resolve()
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    need = {'id', 'name', 'smiles', 'aa_topology'}
    if not rows or not need <= set(rows[0]):
        sys.exit(f'{path}: needs a header with columns {sorted(need)}')
    ids = [r['id'] for r in rows]
    if len(set(ids)) != len(ids):
        sys.exit(f'{path}: duplicate ids')
    for r in rows:
        top = Path(r['aa_topology'])
        r['aa_topology'] = str(top if top.is_absolute() else (path.parent / top).resolve())
        if not Path(r['aa_topology']).exists():
            sys.exit(f"{r['id']}: OPLS-AA topology {r['aa_topology']} not found (the user supplies one per molecule)")
    return rows


def prompt_for(row, moldir, batch_dir, gmx_threads):
    t = (SKILL / 'assets/agent_prompt.md').read_text()
    return t.format(id=row['id'], name=row['name'], smiles=row['smiles'], aa_topology=row['aa_topology'],
                    aa_engine=row.get('aa_engine') or 'not stated', aa_charge_model=row.get('aa_charge_model') or 'not stated',
                    moldir=moldir, skill_md=SKILL / 'SKILL.md', env=ENV_PATH, gmx_threads=gmx_threads,
                    python=ENV['python_command'][0], batch_dir=batch_dir)


def session_vars(moldir, batch_dir, gmx_threads, launch):
    return {'MARTINI_ENV': str(ENV_PATH), 'MARTINI_MOLDIR': str(moldir), 'MARTINI_BATCH': str(batch_dir),
            'MARTINI_GMX_THREADS': str(gmx_threads), 'MARTINI_LAUNCH': launch}


def child_env(moldir, batch_dir, agent, model, gmx_threads):
    e = dict(os.environ)
    e.pop('OMP_NUM_THREADS', None)
    e.update(session_vars(moldir, batch_dir, gmx_threads, 'batch.py'))
    e.update({'MARTINI_AGENT': agent, 'MARTINI_MODEL': model})
    return e


def env_script(moldir, batch_dir, gmx_threads):
    """agent_env.sh: the session variables for a fallback sub-agent (agent and model stay unset: unknown)."""
    lines = ['# Session variables for this molecule when it is run by a fallback sub-agent (SKILL.md "Batch mode").',
             '# Source it at the start of every shell command: . ' + shlex.quote(str(moldir / 'agent_env.sh'))]
    lines += [f'export {k}={shlex.quote(v)}' for k, v in session_vars(moldir, batch_dir, gmx_threads, 'fallback').items()]
    lines += ['unset OMP_NUM_THREADS MARTINI_AGENT MARTINI_MODEL', '']
    return '\n'.join(lines)


# ---------------------------------------------------------------- preflight (no agent, no API call)

def resolve(alternatives, env=None):
    """First alternative ("a|b|c") whose $VARS are all set, expanded; None if none resolves."""
    env = os.environ if env is None else env
    for alt in alternatives.split('|'):
        names = re.findall(r'\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?', alt)
        if all(env.get(n) for n in names):
            return re.sub(r'\$\{?([A-Za-z_][A-Za-z0-9_]*)\}?', lambda m: env[m.group(1)], alt).strip()
    return None


def host_port(target):
    if '://' in target:
        u = urlparse(target)
        return u.hostname, u.port or (443 if u.scheme == 'https' else 80)
    host, _, port = target.rpartition(':')
    return (host, int(port)) if host and port.isdigit() else (target, 443)


def check_writable(path):
    """Create and delete a file in `path` (or its nearest existing parent, which the agent would create it in)."""
    p = Path(path).expanduser()
    while not p.exists() and p != p.parent:
        p = p.parent
    try:
        with tempfile.NamedTemporaryFile(dir=p, prefix='.martini_preflight_'):
            pass
        return None
    except OSError as e:
        return f'{path} is not writable (tested in {p}): {e.strerror or e}'


def check_connect(target, timeout=5.0, env=None):
    env = os.environ if env is None else env
    proxy = next((env[k] for k in ('HTTPS_PROXY', 'https_proxy', 'ALL_PROXY', 'all_proxy') if env.get(k)), None)
    host, port = host_port(proxy or target)
    via = f' (through proxy {host}:{port})' if proxy else ''
    try:
        socket.create_connection((host, port), timeout=timeout).close()
        return None
    except OSError as e:
        return f'cannot reach {target}{via}: {e.strerror or e}'


def preflight(agent, env=None):
    env = os.environ if env is None else env
    spec = {**BUILTIN_PREFLIGHT, **ENV.get('agent_preflight', {})}.get(agent)
    out = {'spec': spec, 'checked': [], 'problems': [],
           'sandbox_markers': sorted(k for k in env if 'SANDBOX' in k.upper())}
    if spec is None:
        out['checked'].append(f'no preflight defined for agent {agent!r} (environment.json agent_preflight)')
        return out
    for kind, fn in (('writable', check_writable), ('connect', check_connect)):
        for item in spec.get(kind, []):
            target = resolve(item, env)
            if target is None:
                out['problems'].append(f'{kind}: none of {item!r} resolves (unset variables)')
                continue
            out['checked'].append(f'{kind} {target}')
            err = fn(target)
            if err:
                out['problems'].append(err)
    return out


def fallback_text(bd, ids, reason, cmdline, jobs):
    py = shlex.join(list(ENV['python_command']) + [str(SCRIPT_DIR / 'cgbuild.py'), 'summary', str(bd)])
    L = [f'# Batch {bd.name}: launch failed', '', f'Reason: {reason}', '',
         f'Molecules without a session: {", ".join(ids)}', '',
         '## Option A: rerun where the checks pass', '',
         'Run the same command from a terminal outside any agent sandbox, or from an agent session that can write',
         "the agent CLI's state folder and reach its service (README \"Batch mode\"). It creates a new batch folder:", '',
         '```', cmdline, '```', '',
         '## Option B: fallback, one fresh sub-agent per molecule', '',
         'For each molecule above, start one new sub-agent with no context from other molecules, and give it the text of',
         '`<ID>/agent_prompt.txt` as its task (full paths below). The prompt tells it to source `<ID>/agent_env.sh`, which binds',
         'its shell commands to that folder. Do not add, remove or summarise anything in the prompt, and never give a sub-agent more than',
         f'one molecule or two sub-agents the same molecule. Run at most {jobs} at a time.', '']
    for i in ids:
        L.append(f'- {i}: prompt `{bd / i / "agent_prompt.txt"}`, environment `{bd / i / "agent_env.sh"}`')
    L += ['', 'When all have ended, build the batch report from the records:', '', '```', py, '```', '',
          'The summary shows, for each molecule, how it was launched (batch.py or fallback, from its design log), whether',
          'all its decisions came from one session bound to its own folder, and how far it got. The agent and model of a',
          'fallback sub-agent are not recorded by the launcher; report them as stated by the dispatcher, or as unknown.', '']
    return '\n'.join(L)


# ---------------------------------------------------------------- run

def did_molecule_work(md):
    return (md / 'input/request.json').exists() or (md / 'decisions/design_log.jsonl').exists()


def cmd_run(a, argv_used):
    rows = read_tsv(a.molecules)
    templates = {**BUILTIN, **ENV.get('agent_commands', {})}
    if a.agent not in templates:
        sys.exit(f'unknown agent {a.agent}; built-ins {sorted(BUILTIN)}, or add it to environment.json agent_commands')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    bd = PROJECT / 'runs' / f'{a.batch}-{stamp}'
    (bd / 'logs').mkdir(parents=True)
    argv_abs = list(argv_used)
    for i, x in enumerate(argv_abs[:-1]):
        if x == '--molecules':
            argv_abs[i + 1] = str(Path(argv_abs[i + 1]).resolve())
    cmdline = shlex.join(list(ENV['python_command']) + [str(Path(__file__).resolve())] + argv_abs)
    rec = {'schema_version': 2, 'batch': bd.name, 'created_utc': now(), 'agent': a.agent, 'model': a.model,
           'command_template': templates[a.agent], 'command_line': cmdline, 'jobs': a.jobs,
           'gmx_threads_per_job': a.gmx_threads, 'timeout_s': a.timeout, 'molecules_tsv': str(Path(a.molecules).resolve()),
           'skill_md_sha256': hashlib.sha256((SKILL / 'SKILL.md').read_bytes()).hexdigest(),
           'rule': 'one agent session per molecule; MARTINI_MOLDIR binds each session to its folder',
           'molecules': [r['id'] for r in rows], 'runs': []}

    prompts = {}
    for row in rows:                       # every folder is prepared before any session starts
        md = bd / row['id']
        md.mkdir()
        prompts[row['id']] = prompt_for(row, md, bd, a.gmx_threads)
        (md / 'agent_prompt.txt').write_text(prompts[row['id']])
        (md / 'agent_env.sh').write_text(env_script(md, bd, a.gmx_threads))

    def finish(launch, code=0):
        rec['launch'] = launch
        rec['finished_utc'] = now()
        (bd / 'batch.json').write_text(json.dumps(rec, indent=1) + '\n')
        if launch['status'] == 'failed':
            todo = launch['molecules_without_session']
            (bd / 'FALLBACK.md').write_text(fallback_text(bd, todo, launch['reason'], cmdline, a.jobs))
            print(f"batch.py: launch failed: {launch['reason']}\n"
                  f"  {len(todo)} molecule(s) have no session; their folders and prompts are ready.\n"
                  f"  Next step: {bd / 'FALLBACK.md'} (rerun outside the sandbox, or one fresh sub-agent per molecule).",
                  file=sys.stderr)
        print(bd)
        sys.exit(code)

    if a.dry_run:
        rec['runs'] = [{'id': r['id'], 'moldir': str(bd / r['id']), 'dry_run': True} for r in rows]
        finish({'status': 'dry_run'})

    pf = preflight(a.agent) if not a.skip_preflight else {'skipped': True, 'problems': []}
    rec['preflight'] = pf
    if pf['problems']:
        marks = f"; sandbox variables set: {', '.join(pf['sandbox_markers'])}" if pf.get('sandbox_markers') else ''
        finish({'status': 'failed', 'stage': 'preflight',
                'reason': 'the agent sessions could not start here: ' + '; '.join(pf['problems']) + marks,
                'molecules_without_session': [r['id'] for r in rows]}, LAUNCH_FAILED_EXIT)

    stop = threading.Event()

    def one(row):
        md = bd / row['id']
        p = prompts[row['id']]
        argv = [x.format(prompt=p, model=a.model, moldir=md, project=PROJECT, budget=a.budget_usd) for x in templates[a.agent]]
        r = {'id': row['id'], 'moldir': str(md), 'argv_without_prompt': [x for x in argv if x != p]}
        if stop.is_set():
            r['launch'] = 'not_started'
            return r
        r['start_utc'] = now()
        t0 = time.time()
        with open(bd / 'logs' / f"{row['id']}.out", 'w') as out, open(bd / 'logs' / f"{row['id']}.err", 'w') as err:
            try:
                rc = subprocess.run(argv, cwd=md, env=child_env(md, bd, a.agent, a.model, a.gmx_threads),
                                    stdout=out, stderr=err, timeout=a.timeout).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
            except OSError as e:          # agent CLI missing or not executable
                err.write(f'batch.py: could not start {argv[0]}: {e}\n')
                rc = 'not_executable'
        worked = did_molecule_work(md)
        r.update({'exit': rc, 'wall_s': round(time.time() - t0, 1), 'end_utc': now(),
                  'launch': 'ok' if worked else 'failed_no_molecule_work',
                  'finalized': (md / 'validation.json').exists()})
        if not worked:
            stop.set()
        print(f"{row['id']}: exit {rc}, {r['wall_s']} s, molecule work {worked}, finalized {r['finalized']}", flush=True)
        return r

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rec['runs'] = list(ex.map(one, rows))

    e = dict(os.environ)
    e.pop('MARTINI_MOLDIR', None)
    e['MARTINI_ENV'] = str(ENV_PATH)
    todo = [r['id'] for r in rec['runs'] if r['launch'] != 'ok']
    if todo:
        first = next(r for r in rec['runs'] if r['launch'] == 'failed_no_molecule_work')
        tail = (bd / 'logs' / f"{first['id']}.err").read_text(errors='replace').strip().splitlines()[-3:]
        launch = {'status': 'failed', 'stage': 'session',
                  'reason': f"the session for {first['id']} ended (exit {first['exit']}) without doing any molecule work; "
                            f"last lines of logs/{first['id']}.err: " + ' | '.join(tail),
                  'molecules_without_session': todo}
        rec['launch'] = launch
        (bd / 'batch.json').write_text(json.dumps(rec, indent=1) + '\n')
        subprocess.run(list(ENV['python_command']) + [str(SCRIPT_DIR / 'cgbuild.py'), 'summary', str(bd)], env=e)
        finish(launch, LAUNCH_FAILED_EXIT)
    rec['launch'] = {'status': 'ok'}
    (bd / 'batch.json').write_text(json.dumps(rec, indent=1) + '\n')
    subprocess.run(list(ENV['python_command']) + [str(SCRIPT_DIR / 'cgbuild.py'), 'summary', str(bd)], env=e)
    finish(rec['launch'])


def cmd_prompt(a, _argv):
    row = next((r for r in read_tsv(a.molecules) if r['id'] == a.id), None)
    if row is None:
        sys.exit(f'no molecule {a.id}')
    bd = Path(a.batch_dir or PROJECT / 'runs' / 'BATCH')
    print(prompt_for(row, bd / row['id'], bd, 2))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = ap.add_subparsers(dest='cmd', required=True)
    r = sp.add_parser('run')
    r.add_argument('--molecules', required=True); r.add_argument('--batch', required=True)
    r.add_argument('--agent', required=True); r.add_argument('--model', required=True)
    r.add_argument('--jobs', type=int, default=4); r.add_argument('--gmx-threads', type=int, default=2)
    r.add_argument('--timeout', type=int, default=5400); r.add_argument('--budget-usd', default='10')
    r.add_argument('--dry-run', action='store_true')
    r.add_argument('--skip-preflight', action='store_true',
                   help='start the sessions without the state-folder and network checks (recorded in batch.json)')
    r.set_defaults(fn=cmd_run)
    q = sp.add_parser('prompt')
    q.add_argument('--molecules', required=True); q.add_argument('--id', required=True); q.add_argument('--batch-dir')
    q.set_defaults(fn=cmd_prompt)
    argv = sys.argv[1:]
    a = ap.parse_args(argv)
    if os.environ.get('MARTINI_MOLDIR'):
        sys.exit('batch.py: this is an agent session bound to one molecule; batches are started by the user, outside any agent session')
    a.fn(a, argv)


if __name__ == '__main__':
    main()
