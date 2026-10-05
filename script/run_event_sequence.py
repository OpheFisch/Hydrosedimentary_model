"""
Simulate a succession of hydrosedimentary flood events.

The script draws accepted parameter sets from GLUE results, generates a
sequence of characteristic flood events, and propagates the sediment-deposit
stock from one event to the next.

Run from the repository root with:

    python scripts/run_event_sequence.py

User-defined settings are grouped at the beginning of the file.
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from routines_modele import generate_triangular_rainfall, simulate_hydrosedimentary_event


# ============================================================
# User configuration
# ============================================================

GLUE_DIR = ROOT / "results" / "glue"
OUTPUT_FILE = ROOT / "results" / "event_sequence.csv"

n_year = 100
n_simulations = 1000
groundwater_reduction = 0.0
future_simulation=True


PARAMETER_NAMES =  [ "alpha_SR","alpha_GW","T_Q","a_hillslope","T_hillslope","a_remob", "T_remob", "Qc_remob","a_depo", "Qc_depo"]
    

FLOOD_TYPES = {
    "M-S": {
        "season": 'summer',
        "glue_file": "glue_valid_results_ete_versants_depot_remob_mineur3.csv",
        "rainfall_intensity": 0.9 / 15,
        "rainfall_total": 4,
        "events_per_year": 2,
    },
    "I-S": {
        "season": 'summer',
        "glue_file": "glue_valid_results_ete_versants_depot_remob_inter3.csv",
        "rainfall_intensity": 6 / 15,
        "rainfall_total": 19,
        "events_per_year": 4,
    },
    "E-S": {
        "season": 'summer',
        "glue_file": "glue_valid_results_ete_versants_depot_remob_majeur3.csv",
        "rainfall_intensity": 11 / 15,
        "rainfall_total": 26,
        "events_per_year": 3,
    },
    "M-W": {
        "season": 'winter',
        "glue_file": "glue_valid_results_hiver_versants_depot_remob_mineur3.csv",
        "rainfall_intensity": 0.7 / 15,
        "rainfall_total": 12,
        "events_per_year": 2,
    },
    "I-W": {
        "season": 'winter',
        "glue_file": "glue_valid_results_hiver_versants_depot_remob_inter3.csv",
        "rainfall_intensity": 2 / 15,
        "rainfall_total": 37,
        "events_per_year": 5,
    },
    "E-W": {
        "season": 'winter',
        "glue_file": "glue_valid_results_hiver_versants_depot_remob_majeur3.csv",
        "rainfall_intensity": 4 / 15,
        "rainfall_total": 71,
        "events_per_year": 4,
    },
}


def apply_rcp_81_91_scenario(flood_types):
    """Apply the rainfall-frequency and intensity changes of the scenario."""
    scenario = {key: values.copy() for key, values in flood_types.items()}

    frequency_factors = {
        "E-S": 0.8,
        "I-S": 0.5,
        "M-S": 0.6,
        "E-W": 1.3,
        "I-W": 0.9,
        "M-W": 0.9,
    }

    intensity_factors = {
        "E-S": 1.1,
        "E-W": 1.2,
        "I-W": 1.1,
    }

    for flood_type, factor in frequency_factors.items():
        scenario[flood_type]["events_per_year"] *= factor

    for flood_type, factor in intensity_factors.items():
        scenario[flood_type]["rainfall_intensity"] *= factor

    return scenario


def load_glue_results(glue_dir, flood_types):
    """Load accepted GLUE parameter sets for each flood type."""
    glue_results = {}

    for flood_type, config in flood_types.items():
        file_path = glue_dir / config["glue_file"]

        if not file_path.exists():
            raise FileNotFoundError(
                f"Missing GLUE file for {flood_type}: {file_path}"
            )

        data = pd.read_csv(file_path)

        missing_columns = [
            name for name in PARAMETER_NAMES if name not in data.columns
        ]
        if missing_columns:
            raise ValueError(
                f"Missing parameter columns in {file_path}: {missing_columns}"
            )

        glue_results[flood_type] = data

    return glue_results
        
    
            

def build_yearly_event_sequence(flood_types, n_years):
    """Generate the flood-event sequence for the requested number of years."""
    all_events = []

    for year in range(1, n_years + 1):
        for season in ["summer", "winter"]:
            #select the type of flood per season
            season_types = {k:v for k,v in flood_types.items() if v["season"]==season}
            season_events = []
            # construct the list of events
            for flood_type, cfg in season_types.items():
                n_events = np.random.poisson(cfg["events_per_year"])
                for _ in range(n_events):
                    season_events.append({"year": year, "flood_type": flood_type})
    
            # Shuffle only within the season
            np.random.shuffle(season_events)

            # Add this season to the complete chronological sequence
            all_events.extend(season_events)

    return all_events


def simulate_event_sequence(
    event_sequence,
    glue_results,
    groundwater_reduction,flood_types
):
    """
    Simulate the complete flood sequence.

    The deposited sediment mass produced by one event becomes the initial
    sediment stock available for remobilisation during the following event.
    """
    results = []
    deposit_stock_g = 0.0

    for event_id, event in enumerate(event_sequence, start=1):
        rain_type = event["flood_type"]
        config = flood_types[rain_type]

        accepted_parameters = glue_results[rain_type]
        sampled_row = accepted_parameters.iloc[
            np.random.randint(0, len(accepted_parameters))
        ]

        parameters = sampled_row[PARAMETER_NAMES].to_numpy(dtype=float)

        rainfall = generate_triangular_rainfall(
            config["rainfall_intensity"],
            config["rainfall_total"],
        )["PRCP"].values


        erosion, final_deposit_t, criteria = simulate_hydrosedimentary_event(
            rainfall_mm_min=rainfall,
            initial_deposit_g=deposit_stock_g,
            parameters=parameters,
            deposition=True,
            remobilization=True,
            groundwater_reduction_fraction=groundwater_reduction)
        

        exported_mass_g = criteria[5]*1e6
        deposited_mass_g = criteria[7]*1e6
        remobilized_mass_g = criteria[8]*1e6

        hillslope_exported_mass_g = (
            exported_mass_g - remobilized_mass_g + deposited_mass_g
        )

        deposit_stock_g = final_deposit_t*1e6

        results.append(
            {
                "event_id": event_id,
                "year": event["year"],
                "flood_type": rain_type,
                "exported_mass_g": exported_mass_g,
                "deposited_mass_g": deposited_mass_g,
                "remobilized_mass_g": remobilized_mass_g,
                "hillslope_exported_mass_g": hillslope_exported_mass_g,
                "deposit_stock_after_event_g": deposit_stock_g,
            }
        )

    return pd.DataFrame(results)


def main():
    """Run the complete flood-sequence simulation."""
    if future_simulation :
        flood_types = apply_rcp_81_91_scenario(FLOOD_TYPES)
    else :
        flood_types =FLOOD_TYPES
    glue_results = load_glue_results(GLUE_DIR, flood_types)

    

    all_simulations = []

    for simulation_id in range(1, n_simulations + 1):
        event_sequence = build_yearly_event_sequence(
        flood_types,n_years,)
        
        simulation_results = simulate_event_sequence(
            event_sequence,
            glue_results,
            groundwater_reduction,flood_types
        )
        simulation_results["simulation_id"] = simulation_id
        all_simulations.append(simulation_results)

    results = pd.concat(all_simulations, ignore_index=True)

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    results.to_csv(OUTPUT_FILE, index=False)

    print(f"Results saved to: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
