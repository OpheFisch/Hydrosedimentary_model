"""
Example entry point for Monte Carlo / GLUE simulations.

Run from the repository root, for example:
    python scripts/run_glue.py

The values below reproduce the structure of the original experiments, but
the default number of simulations is intentionally small for a test run.
Set N to 500_000 for the production runs used in the study.
"""

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from hydrosedimentary_model import monte_carlo_glue_debug


def summer_major_example(N=1000):
    """Example corresponding to the summer major-flood configuration."""
    criteria_bounds = [
        [0.1, 4],          # Qmax
        [0.1 * 0.7, 2 * 1.3],  # QSRmax
        [5, 210],          # V
        [2 * 0.7, 110 * 1.3],  # V_SR
        [20, 192],         # SSCmax
        [13, 1921],        # Vs
        [-0.5, -0.1],      # HI
        [-10, 10],         # phase lag (min)
        [300, 1700],       # deposited mass
        [0, 100],          # remobilized mass
    ]

    return monte_carlo_glue_debug(
        n_simulations=N,
        rainfall_intensity=11 / 15,
        cumulative_rainfall=26,
        month=6,
        initial_deposit_g=0,
        output_name="_ete_versants_depot_remob_majeur",
        hillslope=True,
        deposition=True,
        remobilization=True,
        criteria_bounds=criteria_bounds,
        output_dir=ROOT / "results" / "glue",
    )


if __name__ == "__main__":
    result = summer_major_example(N=1000)
    print(result.head())
