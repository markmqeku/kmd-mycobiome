#!/usr/bin/env python
"""KMD PROMPT 3, A7. NCBI submission metadata, GenBank files and validation. Submits nothing.

Values come only from metadata_v3, KMD_CHECKPOINT.md 2a (dates, sites, design) and the Methods
(coordinates, library design). A field with no documented value is left as MARK_TO_CHOOSE and listed,
never filled by guess. Every file written is checked to lie inside the project root."""
import csv, gzip, hashlib, re
from collections import defaultdict
from pathlib import Path

ROOT = Path(r"<KMD_ROOT>")
O = ROOT / "__reanalysis_2026-06" / "phase3" / "outputs_27-07-2026"
NS = ROOT / "NCBI_SUBMISSION"
META_DIR, GB_DIR, READS = NS / "SRA_metadata", NS / "GenBank", NS / "SRA_reads"
for d in (META_DIR, GB_DIR):
    d.mkdir(parents=True, exist_ok=True)
    assert ROOT in d.resolve().parents

# ---------------- sites: documented values ----------------
def dms(d, m, s, hemi):
    v = d + m / 60 + s / 3600
    return f"{v:.5f} {hemi}"

SITE = {
    "Bloemfontein": dict(geo="South Africa: Free State, Bloemfontein", date="26-Apr-2016",
                         lat_lon=f"{dms(29,6,38.3,'S')} {dms(26,10,51.7,'E')}", run="run 3, 2 x 250 bp"),
    "Christiana":   dict(geo="South Africa: North West, Christiana", date="10-May-2017",
                         lat_lon=f"{dms(27,54,24.9,'S')} {dms(25,9,48.9,'E')}", run="run 13, 2 x 301 bp"),
    "Pretoria":     dict(geo="South Africa: Gauteng, Pretoria", date="18-May-2017",
                         lat_lon=f"{dms(25,44,50.1,'S')} {dms(28,15,31.6,'E')}", run="run 13, 2 x 301 bp"),
}
PO = {"Leaves": "leaf [PO:0025034]", "Twigs": "stem [PO:0009047]", "Inflorescence": "inflorescence [PO:0009049]",
      "Seeds": "seed [PO:0009010]"}          # malformation has no documented PO term: flagged
TISSUE_WORD = {"Leaves": "leaf", "Twigs": "twig", "Inflorescence": "inflorescence", "Seeds": "seed",
               "Malformation": "malformed shoot"}
STERIL = "surface sterilised in sodium hypochlorite, distilled water, 96 percent ethanol and distilled water; lyophilised; DNA extracted from 40 mg with NucleoSpin Plant II Midi"

meta = [r for r in csv.DictReader(open(O / "clean" / "metadata_v3_27-07-2026.tsv", encoding="utf-8"), delimiter="\t")
        if not r["sample-id"].startswith("#")]
meta = {r["sample-id"]: r for r in meta}
reads = sorted(p.name for p in READS.glob("*.fastq.gz"))
samples = sorted({n.rsplit("_R", 1)[0] for n in reads})
assert len(samples) == 45 and len(reads) == 90, (len(samples), len(reads))
missing_meta = [s for s in samples if s not in meta]
assert not missing_meta, missing_meta

flags = defaultdict(list)
bios, sra = [], []
for s in samples:
    m = meta[s]; site = SITE[m["location"]]
    cond = m["condition_std"]
    tissue = m["tissue"]
    age = m["age"] if m["age"] and m["age"] != "not specified" else "not specified"
    plant_struc = PO.get(tissue, "MARK_TO_CHOOSE")
    if plant_struc == "MARK_TO_CHOOSE": flags["plant_struc (malformation tissue: no documented Plant Ontology term)"].append(s)
    host_disease = {"symptomatic": "karee malformation disease", "asymptomatic": "none visible (asymptomatic tree at an affected site)",
                    "reference_site": "none recorded (reference site)"}[cond]
    bios.append({
        "*sample_name": s, "sample_title": f"Searsia lancea {TISSUE_WORD[tissue]}, {m['location']}, {cond}",
        "bioproject_accession": "", "*organism": "plant metagenome", "*collection_date": site["date"],
        "*env_broad_scale": "MARK_TO_CHOOSE", "*env_local_scale": "MARK_TO_CHOOSE", "*env_medium": "plant matter [ENVO:01001121]",
        "*geo_loc_name": site["geo"], "*host": "Searsia lancea", "*lat_lon": site["lat_lon"],
        "host_taxid": "298678", "isolation_source": f"surface-sterilised {TISSUE_WORD[tissue]} tissue of Searsia lancea",
        "plant_struc": plant_struc, "host_disease": host_disease, "samp_mat_process": STERIL,
        "condition": cond, "age_class": age, "tissue": tissue})
    sra.append({
        "sample_name": s, "library_ID": f"KMD_{s}", "title": f"ITS2 amplicon of Searsia lancea {TISSUE_WORD[tissue]}, {m['location']}, {cond}",
        "library_strategy": "AMPLICON", "library_source": "METAGENOMIC", "library_selection": "PCR", "library_layout": "paired",
        "platform": "ILLUMINA", "instrument_model": "Illumina MiSeq",
        "design_description": f"ITS2 amplified with ITS3 and ITS4 carrying Nextera overhang adapters, KAPA HiFi HotStart ReadyMix, 25 cycles, "
                              f"Nextera XT indices over 8 cycles; Illumina MiSeq {site['run']} paired-end",
        "filetype": "fastq", "filename": f"{s}_R1.fastq.gz", "filename2": f"{s}_R2.fastq.gz"})
flags["env_broad_scale (biome choice per site; candidates in COWORK_FINDINGS/DEPOSITION_RULES.md 5.3)"] = samples
flags["env_local_scale (local-scale term per site; candidates in DEPOSITION_RULES.md 5.3)"] = samples


def write_tsv(path, rows):
    assert ROOT in path.resolve().parents
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

write_tsv(META_DIR / "BioSample_MIMARKS_survey_plant-associated_6.0.tsv", bios)
write_tsv(META_DIR / "SRA_metadata.tsv", sra)

# ---------------- GenBank: the 68 lineage ASVs that placed within Tremellomycetes ----------------
def fasta(p):
    d, k = {}, None
    for l in open(p, encoding="utf-8"):
        l = l.rstrip("\n")
        if l.startswith(">"): k = l[1:]; d[k] = ""
        elif k: d[k] += l.strip()
    return d

gb71 = fasta(O / "clean" / "genbank_submission" / "Tremellomycetes_ASVs.fasta")
yeast = fasta(O / "yeast_placement" / "yeast_asvs.fasta")
hit = {l.split("\t")[0] for l in open(O / "yeast_placement" / "yeast_vs_tremello.tsv", encoding="utf-8") if l.strip()}
seq2yid = {s.upper(): h.split()[0] for h, s in yeast.items()}
tax = {}
for r in csv.reader(open(O / "final_taxonomy.tsv", encoding="utf-8"), delimiter="\t"):
    if len(r) >= 2: tax[r[0]] = r[1]
lines = [l.rstrip("\n").split("\t") for l in open(O / "clean" / "feature-table-clean.tsv", encoding="utf-8") if not l.startswith("# ")]
cols = lines[0][1:]
counts = {row[0]: {c: float(v) for c, v in zip(cols, row[1:])} for row in lines[1:]}

kept, excluded, src_rows = {}, [], []
for head, seq in gb71.items():
    sid = head.split()[0]
    md5 = re.search(r"md5=([0-9a-f]{32})", head).group(1)
    assert hashlib.md5(seq.upper().encode()).hexdigest() == md5, f"md5 mismatch for {sid}"
    yid = seq2yid.get(seq.upper())
    assert yid, f"{sid} not in placement input"
    if yid not in hit:
        excluded.append((sid, md5, yid, tax.get(md5, "not in final_taxonomy.tsv"),
                         int(sum(counts.get(md5, {}).values()))))
        continue
    present = [c for c, v in counts.get(md5, {}).items() if v > 0]
    sites = sorted({meta[c]["location"] for c in present})
    tissues = sorted({TISSUE_WORD[meta[c]["tissue"]] for c in present})
    dates = sorted({SITE[s_]["date"] for s_ in sites}, key=lambda d: d[-4:] + {"Apr": "04", "May": "05"}[d[3:6]] + d[:2])
    geo = SITE[sites[0]]["geo"] if len(sites) == 1 else "South Africa"
    cdate = dates[0] if len(dates) == 1 else f"{dates[0]}/{dates[-1]}"
    iso = f"surface-sterilised Searsia lancea tissue ({', '.join(tissues)})"
    kept[sid] = seq
    src_rows.append({"Sequence_ID": sid, "organism": "uncultured Tremellomycetes", "clone": sid, "isolation-source": iso,
                     "host": "Searsia lancea", "geo_loc_name": geo, "collection_date": cdate,
                     "fwd_primer_name": "ITS3", "rev_primer_name": "ITS4",
                     "note": f"ITS2 amplicon sequence variant (DADA2) from Illumina MiSeq paired-end reads; placed within Tremellomycetes; md5 {md5}; sites {', '.join(sites)}"})

assert len(kept) == 68 and len(excluded) == 3, (len(kept), len(excluded))
with open(GB_DIR / "Tremellomycetes_68_ASVs.fasta", "w", encoding="utf-8", newline="\n") as fh:
    for r in src_rows:
        sid = r["Sequence_ID"]
        fh.write(f">{sid} [organism=uncultured Tremellomycetes] [clone={sid}] [host=Searsia lancea] "
                 f"[isolation-source={r['isolation-source']}] [geo_loc_name={r['geo_loc_name']}] [collection_date={r['collection_date']}]\n")
        s = kept[sid]
        for i in range(0, len(s), 70): fh.write(s[i:i + 70] + "\n")
write_tsv(GB_DIR / "source_table.tsv", src_rows)
with open(GB_DIR / "excluded_3_ASVs.tsv", "w", encoding="utf-8", newline="\n") as fh:
    fh.write("sequence_id\tmd5\tplacement_input_id\tUNITE_10_final_taxonomy\treads_in_fungal_table\tplacement\n")
    for sid, md5, yid, t, n in excluded:
        fh.write(f"{sid}\t{md5}\t{yid}\t{t}\t{n}\tno hit in the Tremellomycetes reference set (yeast_vs_tremello.tsv); no confident NCBI hit (FINAL_NUMBERS section 9)\n")

# ---------------- validation ----------------
problems = []
for r in sra:
    for f in (r["filename"], r["filename2"]):
        if not (READS / f).exists(): problems.append(f"metadata names a missing file: {f}")
named = {r["filename"] for r in sra} | {r["filename2"] for r in sra}
extra = set(reads) - named
if extra: problems.append(f"read files not named in metadata: {sorted(extra)}")
for bad in ("P6", "P15", "P24", "YT-P", "30-C"):
    if any(n.startswith(bad + "_") for n in reads): problems.append(f"excluded sample present: {bad}")
md5_ok = 0
for line in open(READS / "md5sums.txt", encoding="utf-8"):
    h, n = line.split()
    hh = hashlib.md5(); fh = open(READS / n, "rb")
    for b in iter(lambda: fh.read(1 << 22), b""): hh.update(b)
    fh.close()
    if hh.hexdigest() == h: md5_ok += 1
    else: problems.append(f"md5 mismatch {n}")
# gzip copies of uncompressed originals must decompress to the original bytes
stg = [p for p in ROOT.iterdir() if p.name.startswith("__restructured_")][0] / "raw_reads"
roundtrip = 0
for r in csv.DictReader(open(NS / "SRA_reads_copy_log.tsv", encoding="utf-8"), delimiter="\t"):
    for src, dst in ((r["source_R1"], f"{r['sample']}_R1.fastq.gz"), (r["source_R2"], f"{r['sample']}_R2.fastq.gz")):
        if src.endswith(".fastq"):
            a = hashlib.md5(); b = hashlib.md5()
            with open(stg / src, "rb") as f1:
                for x in iter(lambda: f1.read(1 << 22), b""): a.update(x)
            with gzip.open(READS / dst, "rb") as f2:
                for x in iter(lambda: f2.read(1 << 22), b""): b.update(x)
            if a.hexdigest() == b.hexdigest(): roundtrip += 1
            else: problems.append(f"gzip copy differs from original: {src}")

with open(META_DIR / "FIELDS_WITHOUT_DOCUMENTED_VALUE.md", "w", encoding="utf-8", newline="\n") as fh:
    fh.write("# Fields left as MARK_TO_CHOOSE\n\nNo documented value exists for these, so none was invented.\n\n")
    for k, v in flags.items():
        fh.write(f"- **{k}**: {len(v)} samples ({', '.join(v) if len(v) < 45 else 'all 45'})\n")
print(f"BioSample rows {len(bios)}, SRA rows {len(sra)}, read files {len(reads)}, md5 verified {md5_ok}/90, "
      f"gzip round-trips verified {roundtrip}")
print(f"GenBank: {len(kept)} kept, {len(excluded)} excluded:")
for e in excluded: print("   ", e)
print("flags:", {k[:40]: len(v) for k, v in flags.items()})
print("PROBLEMS:", problems if problems else "none")
