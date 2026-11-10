#!/usr/bin/env python3
"""Lock/verify reference bytes, not scientific correctness. No external packages."""
import argparse
import hashlib
import json
from pathlib import Path
import sys


def inventory(root):
    result = {}
    for path in sorted(root.rglob('*')):
        rel = path.relative_to(root)
        if '.git' in rel.parts or '__pycache__' in rel.parts:
            continue
        if path.is_symlink():
            raise ValueError(f'Symlinks are not allowed in locked references: {rel}')
        if path.is_file():
            digest = hashlib.sha256()
            with path.open('rb') as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b''):
                    digest.update(chunk)
            result[rel.as_posix()] = digest.hexdigest()
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=['create', 'verify'])
    parser.add_argument('root', type=Path)
    parser.add_argument('lock', type=Path)
    args = parser.parse_args()
    root, lock = args.root.resolve(), args.lock.resolve()
    if not root.is_dir():
        raise ValueError(f'Reference root missing: {root}')
    if lock.is_relative_to(root):
        raise ValueError('Store lock outside reference root to avoid self-hashing.')
    current = inventory(root)
    if not current:
        raise ValueError('Empty reference tree cannot be locked/verified.')
    if args.mode == 'create':
        lock.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation deliberately prevents silent replacement of a baseline.
        with lock.open('x', encoding='utf-8') as handle:
            json.dump({'schema_version': 1, 'algorithm': 'sha256', 'files': current}, handle, indent=2)
        print(f'Locked {len(current)} files. Completeness/content review still required.')
        return
    saved = json.loads(lock.read_text(encoding='utf-8'))['files']
    missing = sorted(set(saved) - set(current))
    added = sorted(set(current) - set(saved))
    changed = sorted(k for k in set(current) & set(saved) if current[k] != saved[k])
    print(json.dumps({'missing': missing, 'added': added, 'changed': changed}, indent=2))
    if missing or added or changed:
        sys.exit(1)
    print(f'PASS: {len(saved)} files match; scientific readiness not assessed.')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        sys.exit(str(exc))
