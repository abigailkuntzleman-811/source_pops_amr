#!/usr/bin/env python3
"""Figure 5: PCA of European-ancestry haplotypes (chromosome 1), one panel per
admixed population projected onto the 1000 Genomes European populations.

Plotting code from Jazeps Medina-Tretmanis. It needs:
  - smartpca output for the haploid, ancestry-masked chromosome 1 data
    (data.evec); set PCA_DIR in config.local.sh to the folder holding it
    (default: chr1_haploid/pca under the repository);
  - Jazeps's plotting module plot_haploid_smartpca.py, importable from this
    folder or on PYTHONPATH.

TODO (Jazeps): add plot_haploid_smartpca.py and the steps that build the
masked haploid input and run smartpca.

Usage:
    python figure_code/fig5/01_plot_pca.py

Output: $FIGURE_DIR/european_pca_v2.png (the name used in the manuscript)
"""

import os
import sys
from pathlib import Path

_root = next(p for p in Path(__file__).resolve().parents if (p / "config.py").exists())
sys.path.insert(0, str(_root))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import FIGURE_DIR  # noqa: E402

try:
    import plot_haploid_smartpca as plotter
except ImportError:
    sys.exit("Needs plot_haploid_smartpca.py (from Jazeps) in figure_code/fig5/ or on PYTHONPATH.")

root = Path(os.path.expanduser(os.environ.get("PCA_DIR", str(_root / "chr1_haploid" / "pca"))))
evec = root / "data.evec"
if not evec.exists():
    sys.exit(f"Missing {evec}; set PCA_DIR to the folder with the smartpca output.")
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
output = FIGURE_DIR / "european_pca_v2.png"

plotter.POPULATION_COLORS.update({
    "CEU": "#1F3A5F",
    "FIN": "#E69F00",
    "GBR": "#8C510A",
    "IBS": "#CC79A7",
    "TSI": "#4D4D4D",
})

plotter.TARGET_COLORS.update({
    "CLM": "#A3D9FF",
    "MXL": "#DA3E52",
    "PEL": "#96E6B3",
    "PUR": "#7E6B8F",
})

points = []
with evec.open(encoding="ascii") as handle:
    for line in handle:
        fields = line.split()
        if fields and not line.startswith("#") and len(fields) >= 4:
            points.append((float(fields[1]), float(fields[2]), fields[-1]))

plotter.plot_panels(points, ["PUR", "CLM", "MXL", "PEL"], output)
print(f"wrote {output}")
