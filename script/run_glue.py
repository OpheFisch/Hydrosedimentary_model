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




'''
"""Example corresponding to the summer average-flood configuration."""

rainfall_intensity=4/15
cumulative_rainfall=16
month=6
depot_init=0
criteria_bounds = [
                [0.1, 2.8],     # Qmax
                [0.05, 1.5*1.3],     # QSRmax 
                [1,99],      # V
                [1,52*1.3],      # V_SR 
                [1, 111],      # SSCmax
                [1, 948],      # Vs 
                [ -0.5, -0.1 ],     # HI
                [-10 , 10], ## phase lag (min)
                [100 , 1650],  ## deposited mass (tons) 
                [0 , 100], ## remobilized mass (tons)
            ]






"""Example corresponding to the summer intermediate-flood configuration."""

rainfall_intensity=6/15
cumulative_rainfall=19
month=6
depot_init=0
criteria_bounds = [
                [0.1, 3],     # Qmax
                [0.0,1.5*1.3],     # QSRmax
                [3,144],      # V  
                [1,55*1.3],      # V_SR
                [0.5,120],      # SSCmax
                [0.9,1200],      # Vs 
                [ -0.5, -0.1 ],     # HI
                [-10 , 10], ## phase lag (min)
                [100 , 1700], ## deposited mass (tons) 
                [0 , 100], ## remobilized mass (tons)
]


              
"""Example corresponding to the summer minor-flood configuration."""

rainfall_intensity=0.9/15
cumulative_rainfall=4
month=6
depot_init=0
criteria_bounds = [
                [0.1, 1.5],     # Qmax
                [0.0,1],     # QSRmax
                [2,98],      # V 
                [0,50],      # V_SR
                [0.9,66],      # SSCmax
                [1,418],      # Vs 
                [ -0.5, -0.16 ],     # HI
                [-10 , 10], ## phase lag (min)
                [0 , 1700], ## deposited mass (tons) 
                [0 , 100], ## remobilized mass (tons)
]






"""Example corresponding to the winter average-flood configuration."""

rainfall_intensity=2/15
cumulative_rainfall=37
month=10
initial_deposit_g=4*10**17*0.3

criteria_bounds = [
                [1, 15],     # Qmax
                [0.5*0.7, 5.5*1.3],     # QSRmax 
                [10,900],      # V
                [5*0.7,250*1.3],      # V_SR  
                [2 ,24],      # SSCmax
                [2 ,3398],      # Vs 
                [0. ,0.4],     # HI
                [-100 , 100], ## phase lag (min)
                [0 , 100], ## deposited mass (tons) 
                [20 , 1500], ## remobilized mass (tons)
            ]



"""Example corresponding to the winter major-flood configuration."""

rainfall_intensity=4/15
cumulative_rainfall=71
month=10
initial_deposit_g=4*10**17*0.3

criteria_bounds = [
                [12, 28],     # Qmax
                [2*0.7,6*1.3],     # QSRmax
                [373,1527],      # V
                [100*0.7,350*1.3],      # V_SR
                [8,34],      # SSCmax
                [2000,6894],      # Vs 
                [ 0., 0.4 ],     # HI
                [-10 , 10], ## phase lag (min)
                [0 , 100], ## deposited mass (tons)
                [0 , 2000], ## remobilized mass (tons)
]



"""Example corresponding to the winter intermediate-flood configuration."""

rainfall_intensity=2/15
cumulative_rainfall=37
month=10
initial_deposit_g=4*10**17*0.3


criteria_bounds = [
                [1, 11],     # Qmax
                [0.02,4*1.3],     # QSRmax
                [28,896],      # V
                [30*0.7,150*1.3],      # V_SR
                [0.4,26],      # SSCmax
                [2,2478],      # Vs 
                [ 0.,0.42 ],     # HI
                [-10 , 10], ## phase lag (min)
                [0 , 100], ## deposited mass (tons)
                [0 , 1500], ## remobilized mass (tons)
]




"""Example corresponding to the winter minor-flood configuration."""

rainfall_intensity=0.8/15
cumulative_rainfall=13
month=10
initial_deposit_g=4*10**17*0.3


criteria_bounds = [
                [0.1,4],     # Qmax
                [0,3],     # QSRmax
                [2,320],      # V
                [0,150],      # V_SR
                [0.4,8],      # SSCmax
                [3,248],      # Vs 
                [ 0.0, 0.3 ],     # HI 
                [-10 , 10], ## phase lag (min)
                [0 , 1000], ## deposited mass (tons)
                [0 , 1000], ## remobilized mass (tons)
]

'''
