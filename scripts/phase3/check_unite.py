import zipfile, os
B = "<KMD_ROOT>/__reanalysis_2026-06"
targets = [
    ("UNITE 8.2  seqs", f"{B}/phase2/unite/unite82_seqs.qza"),
    ("UNITE 10.0 seqs", f"{B}/phase3/outputs_27-07-2026/unite_current/unite10_seqs.qza"),
    ("UNITE 8.2  tax",  f"{B}/phase2/unite/unite82_tax.qza"),
    ("UNITE 10.0 tax",  f"{B}/phase3/outputs_27-07-2026/unite_current/unite10_tax.qza"),
]
for lab, p in targets:
    if not os.path.exists(p):
        print(f"{lab}: MISSING {p}"); continue
    z = zipfile.ZipFile(p)
    data_files = [n for n in z.namelist() if "/data/" in n and not n.endswith("/")]
    tot = 0
    for n in data_files:
        raw = z.read(n)
        if n.endswith(".fasta"):
            tot = raw.count(b">")
        elif n.endswith(".tsv"):
            tot = raw.count(b"\n") - 1
    print(f"{lab}: {os.path.getsize(p):>10,} bytes | data={[os.path.basename(n) for n in data_files]} | records={tot:,}")
