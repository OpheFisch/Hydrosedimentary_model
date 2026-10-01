"""Model driver and Monte Carlo / GLUE calibration utilities."""

from pathlib import Path

import numpy as np
import pandas as pd

try:
    from .routines_modele import generate_triangular_rainfall, simulate_hydrosedimentary_event
except ImportError:  # Allows execution with ``src`` directly on PYTHONPATH.
    from routines_modele import generate_triangular_rainfall, simulate_hydrosedimentary_event


CRITERIA_NAMES = [
    "Qmax",
    "QSRmax",
    "V",
    "VSR",
    "SSCmax",
    "Vs",
    "HI",
    "phase_lag",
    "deposited_mass",
    "remobilized_mass",
]


def get_active_parameters(month, hillslope=True, remobilization=True, deposition=True):
    """Return calibration parameter names and bounds for a given month."""
    parameter_names = [
        "alpha_SR", "tps_SR", "alpha_GW", "tps_GW",
        "a_SR", "tps_sed_versants",
    ]

    if 5 <= month <= 9:
        parameter_bounds = [
            [0.004, 0.2], [18, 250], [0.01, 0.9], [18, 250], [17, 211], [0, 34]
        ]
    else:
        parameter_bounds = [
            [0.01, 0.35], [232, 3273], [0.03, 1], [232, 3273], [0, 13], [0, 34]
        ]

    if remobilization:
        parameter_names += ["a_lit", "tps_sed_lit", "Qc_remob"]
        parameter_bounds += [[0.01, 7], [-60, 0], [3, 6]]

    if deposition:
        parameter_names += ["a_depot", "Qc_depot"]
        parameter_bounds += [[4e-3, 1.7], [1, 5]]

    return parameter_names, parameter_bounds


def build_full_parameter_vector(active_parameters, hillslope=True, remobilization=True, deposition=True):
    """Reconstruct the complete 11-parameter model vector."""
    index = 0
    alpha_SR = active_parameters[index]; index += 1
    tps_SR = int(round(active_parameters[index])); index += 1
    alpha_GW = active_parameters[index]; index += 1
    tps_GW = active_parameters[index]; index += 1

    if hillslope:
        a_SR = active_parameters[index]; index += 1
        tps_sed_versants = int(round(active_parameters[index])); index += 1
    else:
        a_SR = 0.0
        tps_sed_versants = 0
        index += 2

    if remobilization:
        a_lit = active_parameters[index]; index += 1
        tps_sed_lit = int(round(active_parameters[index])); index += 1
        Qc_remob = active_parameters[index]; index += 1
    else:
        a_lit = tps_sed_lit = Qc_remob = 0.0

    if deposition:
        a_depot = active_parameters[index]; index += 1
        Qc_depot = active_parameters[index]
    else:
        a_depot = Qc_depot = 0.0

    return (
        alpha_SR, tps_SR, alpha_GW, tps_GW, a_SR, tps_sed_versants,
        a_lit, tps_sed_lit, Qc_remob, a_depot, Qc_depot,
    )


def criteria_are_valid(criteria, criteria_bounds):
    """Return True when every model criterion falls within its accepted range."""
    return all(low <= value <= high for value, (low, high) in zip(criteria, criteria_bounds))


def _default_criteria_bounds(month):
    """Return the original default GLUE acceptance ranges."""
    if 5 <= month <= 9:
        return [
            [0.1, 2.75], [0.17, 1.2], [1, 99], [2, 52], [1, 111],
            [1, 948], [-0.49, -0.2], [-10, 10], [100, 1000], [0, 10],
        ]
    return [
        [0.1, 14], [0.3, 4.6], [2, 917], [10, 329], [0.4, 24],
        [2, 3398], [0.1, 0.4], [-100, 100], [0, 10], [100, 1000],
    ]


def monte_carlo_glue(
    n_simulations, rainfall_intensity, cumulative_rainfall, month, initial_deposit_g,
    output_name, hillslope=True, deposition=True, remobilization=True,
    criteria_bounds=None, output_dir="results/glue",
):
    """Run a Monte Carlo / GLUE calibration and save accepted simulations."""
    if criteria_bounds is None:
        criteria_bounds = _default_criteria_bounds(month)

    parameter_names, parameter_bounds = get_active_parameters(
        month, hillslope, remobilization, deposition
    )
    rainfall = generate_triangular_rainfall(rainfall_intensity, cumulative_rainfall)["PRCP"].values

    valid_results = []

    for simulation_index in range(n_simulations):
        # Draw only parameters that are active for the selected configuration.
        active_parameters = [
            np.random.uniform(lower, upper) for lower, upper in parameter_bounds
        ]

        # The calibration imposes the same response time for surface runoff and
        # groundwater, as in the original experiments.
        active_parameters[3] = float(active_parameters[1])

        full_parameters = build_full_parameter_vector(
            active_parameters, hillslope, remobilization, deposition
        )
        erosion, _, criteria = simulate_hydrosedimentary_event(
            rainfall, initial_deposit_g, full_parameters, deposition, remobilization
        )

        if criteria_are_valid(criteria, criteria_bounds):
            hillslope_sediment_t = (
                np.nansum(erosion["SSC_versants"] * erosion["Q"]) * 60 * 1e-3
            )
            maximum_deposition_mass_t = (
                np.nansum(erosion["deposition_concentration"] * erosion["Q"])
                * 60
                * 1e-3
            )
            valid_results.append(
                list(active_parameters) + list(criteria)
                + [hillslope_sediment_t, maximum_deposition_mass_t]
            )

        if simulation_index % 100 == 0:
            print(f"Simulation {simulation_index}/{n_simulations}")

    columns = parameter_names + CRITERIA_NAMES + ["Vs_versants", "M_depot_max"]
    result = pd.DataFrame(valid_results, columns=columns)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path / f"glue_valid_results{output_name}.csv", index=False)

    return result


