#!/usr/bin/env python3
"""Run a batch as one agent session per molecule (SKILL.md "Rules" and "Batch mode").

Each molecule gets its own folder and its own headless agent session, started with the prompt in
assets/agent_prompt.md and the environment variable MARTINI_MOLDIR set to that folder; design.py and cgbuild.py
refuse any other folder in that session. No session sees another molecule, so no decision can be shared or
hard-coded across molecules. Deterministic; makes no chemical decisions.

  batch.py run    --molecules FILE.tsv --batch NAME --agent claude|codex|<name> --model MODEL
                  [--jobs 4] [--gmx-threads 2] [--timeout 5400] [--budget-usd 10] [--dry-run]
  batch.py prompt --molecules FILE.tsv --id ID [--batch-dir DIR]

FILE.tsv: tab-separated with header; columns id, name, smiles, aa_topology (path to the user's OPLS-AA
topology, relative to the TSV), optional aa_engine and aa_charge_model.
Agent commands: environment.json `agent_commands` {name: [argv template]} overrides the built-ins below.
Template fields: {prompt} {model} {moldir} {project} {budget}.
Writes runs/<NAME>-<UTC>/: one folder per molecule, logs/<ID>.out/.err, batch.json (launch record),
then `cgbuild.py summary`.
"""
import argparse, csv, datetime, json, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

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


def child_env(moldir, batch_dir, agent, model, gmx_threads):
    e = dict(os.environ)
    e.pop('OMP_NUM_THREADS', None)
    e.update({'MARTINI_ENV': str(ENV_PATH), 'MARTINI_MOLDIR': str(moldir), 'MARTINI_BATCH': str(batch_dir),
              'MARTINI_AGENT': agent, 'MARTINI_MODEL': model, 'MARTINI_GMX_THREADS': str(gmx_threads)})
    return e


def cmd_run(a):
    rows = read_tsv(a.molecules)
    templates = {**BUILTIN, **ENV.get('agent_commands', {})}
    if a.agent not in templates:
        sys.exit(f'unknown agent {a.agent}; built-ins {sorted(BUILTIN)}, or add it to environment.json agent_commands')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    bd = PROJECT / 'runs' / f'{a.batch}-{stamp}'
    (bd / 'logs').mkdir(parents=True)
    rec = {'schema_version': 1, 'batch': bd.name, 'created_utc': now(), 'agent': a.agent, 'model': a.model,
           'command_template': templates[a.agent], 'jobs': a.jobs, 'gmx_threads_per_job': a.gmx_threads,
           'timeout_s': a.timeout, 'molecules_tsv': str(Path(a.molecules).resolve()),
           'skill_md_sha256': __import__('hashlib').sha256((SKILL / 'SKILL.md').read_bytes()).hexdigest(),
           'rule': 'one agent session per molecule; MARTINI_MOLDIR binds each session to its folder', 'runs': []}

    def one(row):
        md = bd / row['id']
        md.mkdir()
        p = prompt_for(row, md, bd, a.gmx_threads)
        (md / 'agent_prompt.txt').write_text(p)
        argv = [x.format(prompt=p, model=a.model, moldir=md, project=PROJECT, budget=a.budget_usd) for x in templates[a.agent]]
        r = {'id': row['id'], 'moldir': str(md), 'start_utc': now(), 'argv_without_prompt': [x for x in argv if x != p]}
        if a.dry_run:
            r['dry_run'] = True
            return r
        t0 = time.time()
        with open(bd / 'logs' / f"{row['id']}.out", 'w') as out, open(bd / 'logs' / f"{row['id']}.err", 'w') as err:
            try:
                rc = subprocess.run(argv, cwd=md, env=child_env(md, bd, a.agent, a.model, a.gmx_threads),
                                    stdout=out, stderr=err, timeout=a.timeout).returncode
            except subprocess.TimeoutExpired:
                rc = 'timeout'
        r.update({'exit': rc, 'wall_s': round(time.time() - t0, 1), 'end_utc': now(),
                  'finalized': (md / 'validation.json').exists()})
        print(f"{row['id']}: exit {rc}, {r['wall_s']} s, finalized {r['finalized']}", flush=True)
        return r

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        rec['runs'] = list(ex.map(one, rows))
    rec['finished_utc'] = now()
    (bd / 'batch.json').write_text(json.dumps(rec, indent=1) + '\n')
    if not a.dry_run:
        e = dict(os.environ)
        e.pop('MARTINI_MOLDIR', None)
        e['MARTINI_ENV'] = str(ENV_PATH)
        subprocess.run(list(ENV['python_command']) + [str(SCRIPT_DIR / 'cgbuild.py'), 'summary', str(bd)], env=e)
    print(bd)


def cmd_prompt(a):
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
    r.set_defaults(fn=cmd_run)
    q = sp.add_parser('prompt')
    q.add_argument('--molecules', required=True); q.add_argument('--id', required=True); q.add_argument('--batch-dir')
    q.set_defaults(fn=cmd_prompt)
    a = ap.parse_args()
    if os.environ.get('MARTINI_MOLDIR'):
        sys.exit('batch.py: this is an agent session bound to one molecule; batches are started by the user, outside any agent session')
    a.fn(a)


if __name__ == '__main__':
    main()
