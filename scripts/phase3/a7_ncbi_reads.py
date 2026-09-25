#!/usr/bin/env python
"""KMD PROMPT 3, A7. Copy the paired reads of the 45 usable samples into NCBI_SUBMISSION/SRA_reads,
gzipped, with an MD5 file. Originals in raw_reads are only read, never moved or altered.

Samples are taken from the three DADA2 import manifests actually used (phase2), so the deposit is
exactly the denoised set. Every file written is checked to lie inside the project root."""
import csv, gzip, hashlib, shutil, sys, time
from pathlib import Path

ROOT = Path(r"<KMD_ROOT>")
P2 = ROOT / "__reanalysis_2026-06" / "phase2"
stg = [p for p in ROOT.iterdir() if p.is_dir() and p.name.startswith("__restructured_")]
if len(stg) != 1:
    sys.exit(f"expected one staging folder, found {[p.name for p in stg]}")
RAW = stg[0] / "raw_reads"
OUT = ROOT / "NCBI_SUBMISSION" / "SRA_reads"
OUT.mkdir(parents=True, exist_ok=True)
assert ROOT in OUT.resolve().parents
EXCLUDE = {"P6", "P15", "P24", "YT-P", "30-C"}


def remap(p):
    """Manifest paths name a staging folder that has since been renamed, or phase2 itself."""
    p = p.replace("<KMD_ROOT>/", "")
    parts = Path(p).parts
    if parts[0].startswith("__restructured_"):
        return RAW / parts[-1]
    return ROOT / Path(p)


samples = {}
for man in (P2 / "groupB" / "manifest.tsv", P2 / "P2_0_bloem" / "manifest_bloem29.tsv", P2 / "otp" / "man.tsv"):
    for r in csv.DictReader(open(man, encoding="utf-8"), delimiter="\t"):
        sid = r["sample-id"]
        f1 = r["forward-absolute-filepath"]
        f2 = r.get("reverse-absolute-filepath") or f1.replace("_R1_", "_R2_")
        samples[sid] = (remap(f1), remap(f2))
assert not EXCLUDE & set(samples), EXCLUDE & set(samples)
assert len(samples) == 45, len(samples)

stats = {}
for r in csv.reader(open(ROOT / "__reanalysis_2026-06/phase3/outputs_27-07-2026/assembly/TableS4_dada2_stats.tsv", encoding="utf-8"), delimiter="\t"):
    if len(r) > 3 and r[1] in samples:
        stats[r[1]] = int(r[2])


def opener(p):
    return gzip.open(p, "rb") if p.suffix == ".gz" else open(p, "rb")


def copy_count(src, dst):
    """Write src to dst as gzip (re-using gzip bytes unchanged when the source is already gzip),
    and return (records, first header, last header)."""
    if src.suffix == ".gz":
        shutil.copyfile(src, dst)
    else:
        with open(src, "rb") as fi, gzip.open(dst, "wb", compresslevel=6) as fo:
            shutil.copyfileobj(fi, fo, 1 << 22)
    n = 0; first = last = b""
    with gzip.open(dst, "rb") as fh:          # count from the written file, which also tests it
        for i, line in enumerate(fh):
            if i % 4 == 0:
                n += 1; last = line.split()[0]
                if n == 1: first = last
    return n, first, last


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""): h.update(b)
    return h.hexdigest()


log = []
t0 = time.time()
for sid in sorted(samples):
    s1, s2 = samples[sid]
    for s in (s1, s2):
        if not s.exists(): sys.exit(f"MISSING source for {sid}: {s}")
    d1, d2 = OUT / f"{sid}_R1.fastq.gz", OUT / f"{sid}_R2.fastq.gz"
    for d in (d1, d2): assert ROOT in d.resolve().parents
    n1, a1, z1 = copy_count(s1, d1)
    n2, a2, z2 = copy_count(s2, d2)
    pair_ok = n1 == n2 and a1.split(b"/")[0] == a2.split(b"/")[0] and z1.split(b"/")[0] == z2.split(b"/")[0]
    exp = stats.get(sid)
    log.append((sid, s1.name, s2.name, n1, n2, pair_ok, exp, exp == n1 if exp else None))
    print(f"{sid:6} R1 {n1:>8} R2 {n2:>8} pair {'ok' if pair_ok else 'MISMATCH'}  dada2-input {exp} {'ok' if exp == n1 else 'DIFF'}  {time.time()-t0:6.0f}s", flush=True)

with open(OUT / "md5sums.txt", "w", encoding="utf-8", newline="\n") as fh:
    for p in sorted(OUT.glob("*.fastq.gz")):
        fh.write(f"{md5(p)}  {p.name}\n")
with open(OUT.parent / "SRA_reads_copy_log.tsv", "w", encoding="utf-8", newline="\n") as fh:
    fh.write("sample\tsource_R1\tsource_R2\treads_R1\treads_R2\tpairs_match\tdada2_input\tmatches_dada2_input\n")
    for r in log: fh.write("\t".join(str(x) for x in r) + "\n")
print("files:", len(list(OUT.glob('*.fastq.gz'))), "done in", round(time.time() - t0), "s")
