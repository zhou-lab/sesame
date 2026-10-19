#!/usr/bin/env python3
"""Offline HM27 preprocessing regression with synthetic IDATs and a tiny store.

Use the HM27 bead count with Type I probes in both channels. Check the CLI's
store lookup, detection, recommended-mask selection, and default pipeline.
This exercises integration; it is not an R numerical comparison.
"""
import gzip
import math
import os
from pathlib import Path
import struct
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent.parent
BIN = ROOT / "sesame"
YAME = ROOT / "YAME" / "yame"
NBEADS = 55300
NPROBES = 256


def run(*args, **kwargs):
    result = subprocess.run([str(a) for a in args], capture_output=True,
                            text=True, **kwargs)
    if result.returncode:
        raise AssertionError(f"{args}: exit {result.returncode}\n{result.stderr}")
    return result


def write_idat(path, green):
    means = []
    for address in range(1, NBEADS + 1):
        probe = (address - 1) // 2
        in_band = probe < NPROBES and (probe < NPROBES // 2) == green
        signal = (1000 + probe % 100 * 10 if address % 2 else
                  3000 + probe % 100 * 20)
        means.append(signal if in_band else 30 + address % 20)
    fields = [(1000, struct.pack("<i", NBEADS)),
              (102, struct.pack(f"<{NBEADS}I", *range(1, NBEADS + 1))),
              (103, struct.pack(f"<{NBEADS}H", *([5] * NBEADS))),
              (104, struct.pack(f"<{NBEADS}H", *means)),
              (107, bytes([10]) * NBEADS)]
    offset = 16 + 10 * len(fields)
    table = b""
    for code, data in fields:
        table += struct.pack("<Hq", code, offset)
        offset += len(data)
    path.write_bytes(b"IDAT" + struct.pack("<qi", 3, len(fields)) + table +
                     b"".join(data for _, data in fields))


def main():
    if not YAME.is_file():
        subprocess.run(["make", "-C", str(ROOT / "YAME")], check=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        store = work / "store" / "HM27"
        store.mkdir(parents=True)
        env = {**os.environ, "YAME_DATA_HOME": str(store.parent)}
        prefix = work / "sample"
        for channel in ("Grn", "Red"):
            write_idat(work / f"sample_{channel}.idat", channel == "Grn")
        with gzip.open(store / "HM27.ordering.tsv.gz", "wt") as out:
            out.write("Probe_ID\tM\tU\tcol\n")
            for i in range(NPROBES):
                col = "G" if i < NPROBES // 2 else "R"
                out.write(f"cg{i:08d}\t{2*i+1}\t{2*i+2}\t{col}\n")
        # The four tracks actually published for HM27, plus an unselected
        # M_general track. Each recommended track masks a different probe.
        names = ["M_1baseSwitchSNPcommon_5pt", "M_mapping", "M_nonuniq",
                 "M_SNPcommon_5pt", "M_general"]
        mask = store / "HM27.hg38.mask.cm"
        with mask.open("wb") as out:
            for j in range(len(names)):
                track = work / f"track{j}.cm"
                bits = "".join(f"{int(i == j)}\n" for i in range(NPROBES))
                run(YAME, "pack", "-f", "b", "-", track, input=bits)
                out.write(track.read_bytes())
        tracks = work / "tracks.txt"
        tracks.write_text("\n".join(names) + "\n")
        run(YAME, "index", "-s", tracks, mask)

        def preprocess(label, *options):
            out = work / label
            result = run(BIN, "preprocess", "--out", out, *options, prefix, env=env)
            beta = run(YAME, "unpack", out / "beta.cg").stdout.splitlines()
            assert len(beta) == NPROBES, len(beta)
            return result, beta, out

        detected, raw, _ = preprocess("detected", "--prep", "", "--output", "beta")
        assert "detected HM27 (55300 beads)" in detected.stderr, detected.stderr
        named, explicit, _ = preprocess("named", "--platform", "HM27", "--prep", "",
                                        "--output", "beta")
        assert "warning" not in named.stderr, named.stderr
        assert raw == explicit
        for i, value in enumerate(raw):
            m, u = 1000 + i % 100 * 10, 3000 + i % 100 * 20
            # yame unpack prints beta values to three decimal places.
            assert math.isclose(float(value), m / (m + u), abs_tol=5e-4), (i, value, m, u)
        print("ok   HM27 named and detected: raw Type I betas in both channels")

        _, masked, _ = preprocess("quality", "--platform", "HM27", "--prep", "Q",
                                   "--output", "beta")
        assert masked[:4] == ["NA"] * 4, masked[:4]
        assert masked[4:] == raw[4:]
        print("ok   HM27 Q: four recommended tracks selected; M_general excluded")

        _, default, out = preprocess("default")
        assert default[:4] == ["NA"] * 4
        assert any(v != "NA" and math.isfinite(float(v)) for v in default[4:])
        for filename in ("beta.cg", "intensity.cg", "pval.cg", "qc.tsv"):
            assert (out / filename).stat().st_size > 0, filename
        print("ok   HM27 default QCDPB: beta, intensity, pval and QC outputs")


if __name__ == "__main__":
    main()
