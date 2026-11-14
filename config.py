"""Central path configuration for this project.

Every script imports paths from here rather than hardcoding them, so the
pipeline can be run on any machine by editing this file or by exporting the
corresponding environment variables (the same ones config.sh defines).

    from config import HMMIX_DIR, FIGURE_DIR, HMMIX_MIN_PROB
"""

import os
from pathlib import Path

# Repository root, resolved relative to this file.
PROJECT_ROOT = Path(__file__).resolve().parent


def _load_local_config() -> None:
    """Apply settings from config.local.sh, if it exists.

    That file is shell, but it is a plain list of `export NAME=value` lines,
    so it can be read here too. Values already set in the environment win,
    so an explicit export on the command line still overrides it.
    """
    local = PROJECT_ROOT / "config.local.sh"
    if not local.exists():
        return

    for line in local.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):]
        if "=" not in line:
            continue
        name, _, value = line.partition("=")
        name, value = name.strip(), value.strip().strip('"').strip("'")
        if name.isidentifier() and name not in os.environ:
            os.environ[name] = value


_load_local_config()


def _path(env_var: str, default: Path) -> Path:
    """Return the path from an environment variable, or a default."""
    value = os.environ.get(env_var)
    return Path(value) if value else default


# --- Inputs -----------------------------------------------------------------
DATA_RAW = _path("DATA_RAW", PROJECT_ROOT / "data" / "raw")

# Skov (2024) hmmix archaic introgression calls, doi:10.5281/zenodo.14136628
HMMIX_DIR = _path("HMMIX_DIR", DATA_RAW / "skov_archaic")

# Local ancestry inference BED files.
# NOTE: not publicly available. Point this at wherever you have them.
LAI_DIR = _path("LAI_DIR", DATA_RAW / "local_ancestry")

# --- Outputs ----------------------------------------------------------------
DATA_PROCESSED = _path("DATA_PROCESSED", PROJECT_ROOT / "data" / "processed")
RESULTS_DIR = _path("RESULTS_DIR", PROJECT_ROOT / "results")
FIGURE_DIR = _path("FIGURE_DIR", RESULTS_DIR / "figures")
TABLE_DIR = _path("TABLE_DIR", RESULTS_DIR / "tables")

RECOVERED_DIR = _path("RECOVERED_DIR", DATA_PROCESSED / "recovered_archaic_by_admix_prop")

# Per-population, per-haplotype BED files (02_segments_to_bed.py).
BY_POP_DIR = _path("BY_POP_DIR", DATA_PROCESSED / "hmmix_by_pop")

# Pairwise sharing results (03_pairwise_sharing.sh).
PAIRWISE_DIR = _path("PAIRWISE_DIR", DATA_PROCESSED / "pairwise_comparisons")

# --- Analysis parameters ----------------------------------------------------
# Minimum hmmix mean posterior probability for a segment to be retained.
HMMIX_MIN_PROB = float(os.environ.get("HMMIX_MIN_PROB", 0.80))

# Admixed American populations analysed in this study.
AMR_POPS = ["CLM", "PUR", "PEL", "MXL"]

# Autosomes only.
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)]


def ensure_dirs() -> None:
    """Create output directories if they do not already exist."""
    for d in (DATA_PROCESSED, FIGURE_DIR, TABLE_DIR, RECOVERED_DIR):
        d.mkdir(parents=True, exist_ok=True)
