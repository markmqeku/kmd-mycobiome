import zipfile, re, sys
from xml.etree import ElementTree as ET
A = '{http://schemas.openxmlformats.org/drawingml/2006/main}'
p = sys.argv[1]
z = zipfile.ZipFile(p)
slides = sorted([n for n in z.namelist() if re.match(r"ppt/slides/slide\d+\.xml$", n)],
                key=lambda x: int(re.search(r"(\d+)", x).group(1)))
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
print(f"# {p}  ({len(slides)} slides)")
for n in slides:
    root = ET.fromstring(z.read(n))
    texts = [t.text for t in root.iter(A + 't') if t.text and t.text.strip()]
    if texts:
        print(f"\n===== {n} =====")
        for t in texts:
            print(t.strip())
notes = sorted([n for n in z.namelist() if re.match(r"ppt/notesSlides/notesSlide\d+\.xml$", n)])
for n in notes:
    root = ET.fromstring(z.read(n))
    texts = [t.text for t in root.iter(A + 't') if t.text and t.text.strip()]
    if texts:
        print(f"\n===== NOTES {n} =====")
        for t in texts: print(t.strip())
