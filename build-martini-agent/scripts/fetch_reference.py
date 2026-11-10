#!/usr/bin/env python3
"""Fetch one HTTPS reference; validate basic content and optionally extract a ZIP.
Never executes downloaded code. Records URL, redirect target, UTC time and SHA-256.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import shutil
import stat
import sys
import urllib.request
import zipfile


def checked_members(data):
    archive = zipfile.ZipFile(io.BytesIO(data))
    total = 0
    for info in archive.infolist():
        name = PurePosixPath(info.filename)
        if name.is_absolute() or '..' in name.parts or '\\' in info.filename or ':' in info.filename:
            raise ValueError(f'Unsafe ZIP path: {info.filename}')
        if stat.S_ISLNK(info.external_attr >> 16):
            raise ValueError(f'ZIP symlink rejected: {info.filename}')
        total += info.file_size
    if total > 2_000_000_000:
        raise ValueError('Archive expands beyond 2 GB; inspect separately.')
    bad = archive.testzip()
    if bad:
        raise ValueError(f'Corrupt ZIP member: {bad}')
    return archive


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('url')
    p.add_argument('output', type=Path)
    p.add_argument('--kind', choices=['pdf', 'zip', 'html', 'text'], required=True)
    p.add_argument('--extract', type=Path)
    args = p.parse_args()
    if not args.url.startswith('https://'):
        raise ValueError('Use an HTTPS source.')
    meta = Path(str(args.output) + '.source.json')
    if args.output.exists() or meta.exists() or (args.extract and args.extract.exists()):
        raise ValueError('Destination exists. Use a new snapshot; do not overwrite.')
    if args.extract and args.kind != 'zip':
        raise ValueError('--extract requires --kind zip.')
    req = urllib.request.Request(args.url, headers={'User-Agent': 'MartiniReferenceHarness/1.0'})
    with urllib.request.urlopen(req, timeout=60) as response:
        final_url = response.geturl()
        if not final_url.startswith('https://'):
            raise ValueError('Insecure redirect rejected.')
        data = response.read(250_000_001)
    if len(data) > 250_000_000:
        raise ValueError('Download exceeds 250 MB; inspect separately.')
    if not data:
        raise ValueError('Empty response.')
    if args.kind == 'pdf' and not data.startswith(b'%PDF-'):
        raise ValueError('Not a PDF; may be a publisher access page.')
    archive = checked_members(data) if args.kind == 'zip' else None
    if args.kind == 'text' and (b'<html' in data[:2048].lower() or b'<!doctype html' in data[:2048].lower()):
        raise ValueError('Expected text; received HTML.')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open('xb') as handle:
        handle.write(data)
    meta.write_text(json.dumps({'url': args.url, 'resolved_url': final_url,
        'retrieved_utc': datetime.now(timezone.utc).isoformat(),
        'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}, indent=2))
    if args.extract:
        args.extract.mkdir(parents=True)
        # Paths and links checked above; write without honoring archive permissions.
        for member in archive.infolist():
            target = args.extract / member.filename
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
            else:
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(member) as src, target.open('xb') as dest:
                    shutil.copyfileobj(src, dest)
    print(f'Saved {args.output}. Inspect contents before locking references.')


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        sys.exit(str(exc))
