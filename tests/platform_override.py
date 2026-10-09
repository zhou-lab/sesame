#!/usr/bin/env python3
"""Exercise explicit platform overrides with self-contained IDAT fixtures."""
import gzip
import os
from pathlib import Path
import struct
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]


def idat(prefix, count):
    fields = [
        (1000, struct.pack('<i', count)),
        (102, struct.pack(f'<{count}I', *range(1, count + 1))),
        (103, struct.pack('<H', 1) * count),
        (104, struct.pack('<H', 100) * count),
        (107, bytes([5]) * count),
    ]
    offset = 16 + 10 * len(fields)
    table = bytearray()
    for code, data in fields:
        table += struct.pack('<Hq', code, offset)
        offset += len(data)
    data = b'IDAT' + struct.pack('<qi', 3, len(fields)) + table
    data += b''.join(payload for _, payload in fields)
    for channel in ('Grn', 'Red'):
        Path(f'{prefix}_{channel}.idat').write_bytes(data)


with tempfile.TemporaryDirectory() as temp:
    work = Path(temp)
    store = work / 'store'
    for platform in ('EPIC', 'EPICv2', 'HM450', 'MSA', 'NOPE'):
        directory = store / platform
        directory.mkdir(parents=True)
        with gzip.open(directory / f'{platform}.ordering.tsv.gz', 'wt') as f:
            f.write('Probe_ID\tM\tU\tcol\ncg0000001\tNA\t1\t2\n')
    small = work / 'small'
    known = work / 'known'
    idat(small, 2)  # Not a registered platform count.
    idat(known, 622399)  # HM450: permits inference, then batch validation.
    env = dict(os.environ, YAME_DATA_HOME=str(store), XDG_DATA_HOME=str(work / 'xdg'))

    def run(label, options, samples, success, message=None):
        result = subprocess.run(
            [str(ROOT / 'sesame'), 'preprocess', '--prep', '', '--output', 'beta',
             '--threads', '2', '--out', str(work / label), *options, *map(str, samples)],
            env=env, capture_output=True, text=True)
        assert (result.returncode == 0) == success, (label, result.stderr)
        if message:
            assert message in result.stderr, (label, result.stderr)
        if success:
            assert (work / label / 'beta.cg').stat().st_size > 0
        print(f'ok {label}')

    for platform in ('EPIC', 'EPICv2', 'HM450', 'MSA'):
        run(platform, ['--platform', platform], [small, known], True)
    ordering = store / 'EPIC' / 'EPIC.ordering.tsv.gz'
    run('explicit-index', ['--platform', 'EPIC', '--index', str(ordering)], [small], True)
    run('custom-path', ['--platform', str(ordering)], [small], True)
    run('unknown-name', ['--platform', 'NOPE'], [small], False, 'not the batch platform NOPE')
    run('inference-unknown', [], [small], False, 'cannot identify platform from 2 beads')
    run('inference-known', [], [known], True, 'detected HM450')
    run('inference-mixed', [], [known, small], False, 'not the batch platform HM450')
