#!/usr/bin/env python
"""Check embedded figure aspect ratios, to see whether any single figure is inflating the
page count of the review document."""
import os, struct
K = r"<KMD_ROOT>\__reanalysis_2026-06\phase3\outputs_27-07-2026"
FILES = [(r"figures\Figure1_host_by_site.png", "Figure 1"),
         (r"clean\figures\Figure2_bloem_hostload_vs_ASVs.png", "Figure 2"),
         (r"figures\Figure3_tremellomycetes_placement.png", "Figure 3"),
         (r"clean\figures\Figure4_composition.png", "Figure 4"),
         (r"clean\figures\Figure5_hill_numbers.png", "Figure 5"),
         (r"clean\figures\Figure6_ordination.png", "Figure 6"),
         (r"clean\figures\Figure7_curvibasidium.png", "Figure 7"),
         (r"figures\Figure8_placement_by_neighbourhood.png", "Figure 8")]
PAGE_TEXT_HEIGHT = 29.7 - 5.0   # A4 less 2.5 cm top and bottom margins
for rel, name in FILES:
    p = os.path.join(K, rel)
    with open(p, "rb") as f:
        f.read(16); w, h = struct.unpack(">II", f.read(8))
    at155 = 15.5 * h / w
    flag = "  <-- taller than the text area" if at155 > PAGE_TEXT_HEIGHT else ""
    print(f"{name}: {w}x{h} px, ratio {h/w:.2f}, height at 15.5 cm wide = {at155:5.1f} cm{flag}")
print(f"\nA4 text area height with 2.5 cm margins: {PAGE_TEXT_HEIGHT:.1f} cm")
