"""Core numerical routines for the hydrosedimentary model.

The model represents an event-scale catchment response through four main
steps:

1. convert rainfall into surface-runoff and groundwater discharge;
2. convert surface runoff into hillslope sediment concentration;
3. remobilize sediment stored in the channel when discharge exceeds a
   threshold;
4. remove sediment from suspension through deposition.

The functions in this module are deliberately written with NumPy arrays so
that they can be called many times during Monte Carlo / GLUE calibration.
"""

import numpy as np
import pandas as pd


def generate_triangular_rainfall(max_intensity_mm_min, cumulative_rainfall_mm):
    """Generate a symmetric triangular rainfall event.

    Parameters
    ----------
    max_intensity_mm_min : float
        Peak rainfall intensity in mm/min.
    cumulative_rainfall_mm : float
        Total rainfall depth in mm.

    Returns
    -------
    pandas.DataFrame
        A table containing time ``t`` (min) and rainfall intensity ``PRCP``
        (mm/min).

    Notes
    -----
    The triangular event is constructed so that its integral is equal to the
    requested cumulative rainfall depth.
    """
    duration_min = int(2 * cumulative_rainfall_mm / max_intensity_mm_min)
    time_min = np.arange(duration_min)
    half_duration = int(duration_min / 2)

    rainfall_intensity = (
        [max_intensity_mm_min * i / half_duration for i in range(half_duration)]
        + [
            max_intensity_mm_min / half_duration * (duration_min - i)
            for i in range(half_duration, duration_min)
        ]
    )

    return pd.DataFrame({"t": time_min, "PRCP": rainfall_intensity})


def compute_hillslope_sediment_concentration(
    total_discharge_m3_s, surface_runoff_m3_s, sediment_transfer_time_min, hillslope_coefficient
):
    """Compute hillslope sediment concentration from surface runoff.

    A linear rating curve is used:

    ``SSC = hillslope_coefficient * surface_runoff``

    The resulting concentration is shifted in time to represent sediment
    transfer from the hillslopes to the channel. Suspended-sediment flux is
    then obtained by multiplying concentration by total discharge.

    Parameters
    ----------
    total_discharge_m3_s : array-like
        Total catchment discharge.
    surface_runoff_m3_s : array-like
        Surface-runoff component of discharge.
    sediment_transfer_time_min : int
        Sediment transfer time in minutes. Positive values delay the sediment
        response.
    hillslope_coefficient : float
        Coefficient of the linear sediment rating curve.

    Returns
    -------
    tuple of numpy.ndarray
        Sediment concentration (g/L) and sediment flux (g/s).
    """
    n_steps = len(total_discharge_m3_s)
    hillslope_concentration = hillslope_coefficient * np.maximum(surface_runoff_m3_s, 0)

    shifted_concentration = np.zeros(n_steps)
    if sediment_transfer_time_min >= 0:
        shifted_concentration[sediment_transfer_time_min:] = hillslope_concentration[
            : n_steps - sediment_transfer_time_min
        ]
    else:
        shifted_concentration[:sediment_transfer_time_min] = hillslope_concentration[
            -sediment_transfer_time_min :
        ]

    sediment_flux_g_s = shifted_concentration * total_discharge_m3_s * 1e3
    return shifted_concentration, sediment_flux_g_s


def compute_channel_remobilization(
    total_discharge_m3_s,
    initial_deposit_g,
    remobilization_coefficient,
    critical_discharge_m3_s,
    sediment_transfer_time_min,
):
    """Compute channel-sediment remobilization.

    Remobilization follows a linear threshold relationship. Sediment is
    available only while the initial channel deposit remains positive.

    Parameters
    ----------
    total_discharge_m3_s : array-like
        Total discharge at the outlet.
    initial_deposit_g : float
        Initial mass of sediment stored in the channel, in g.
    remobilization_coefficient : float
        Coefficient of the remobilization relationship.
    critical_discharge_m3_s : float
        Discharge threshold above which remobilization occurs.
    sediment_transfer_time_min : int
        Transfer time between remobilization and outlet concentration. 
        Negative values are allowed and represent sediment peak arraving before discharge peak

    Returns
    -------
    tuple of numpy.ndarray, numpy.ndarray, float
        Remobilized concentration (g/L), remobilized flux (g/s), and the
        remaining channel deposit (g).
    """
    n_steps = len(total_discharge_m3_s)

    remobilized_concentration = remobilization_coefficient * np.maximum(
        total_discharge_m3_s - critical_discharge_m3_s, 0
    )

    shifted_concentration = np.zeros(n_steps)
    if sediment_transfer_time_min >= 0:
        if sediment_transfer_time_min < n_steps:
            shifted_concentration[sediment_transfer_time_min:] = remobilized_concentration[
                : n_steps - sediment_transfer_time_min
            ]
    else:
        if -sediment_transfer_time_min < n_steps:
            shifted_concentration[:sediment_transfer_time_min] = remobilized_concentration[
                -sediment_transfer_time_min :
            ]

    remobilized_flux_g_s = shifted_concentration * total_discharge_m3_s * 1e3

    # Update the sediment stock through time. The time step is one minute,
    # hence the conversion from g/s to g uses a factor of 60.
    remaining_deposit_g = initial_deposit_g
    for step in range(n_steps): 
        available_g = remaining_deposit_g
        requested_g = remobilized_flux_g_s[step] * 60
        if available_g>=requested_g:
            remaining_deposit_g -= requested_g
        # Once the available deposit is exhausted, no further remobilization
        # can occur during this event.
        if available_g<requested_g:
            remaining_deposit_g =0
            remobilized_flux_g_s[step] = available_g / 60
            if step<n_steps-1:
                shifted_concentration[step+1:] = 0
                remobilized_flux_g_s[step+1:] = 0
            break

    return shifted_concentration, remobilized_flux_g_s, remaining_deposit_g


def simulate_reservoir_discharge(
    rainfall_mm_min, runoff_coefficient, response_time_min, catchment_area_km2,
    time_step_min, rainfall_threshold_mm_min
):
    """Simulate discharge using a linear reservoir.

    Rainfall above ``rainfall_threshold_mm_min`` recharges the reservoir.
    Reservoir outflow is controlled by the response time, producing a
    discharge time series in m3/s.
    """
    n_steps = len(rainfall_mm_min) + int(7 * response_time_min / time_step_min)
    reservoir_volume_m3 = np.zeros(n_steps)
    discharge_m3_s = np.zeros(n_steps)

    recession_coefficient = 1 / max(response_time_min * time_step_min, 1e-6)

    for step in range(1, n_steps):
        recharge_m3_min = 0.0
        if step < len(rainfall_mm_min) and rainfall_mm_min[step] > rainfall_threshold_mm_min:
            recharge_m3_min = (
                runoff_coefficient
                * rainfall_mm_min[step]
                * catchment_area_km2
                * 1e3
                / time_step_min
            )

        reservoir_volume_m3[step] = reservoir_volume_m3[step - 1] + time_step_min * (
            recharge_m3_min - recession_coefficient * reservoir_volume_m3[step - 1]
        )
        discharge_m3_s[step] = recession_coefficient * reservoir_volume_m3[step - 1] / 60

    return discharge_m3_s


def compute_hysteresis_index(discharge_m3_s, sediment_concentration_g_l, n_levels=100):
    """Compute the Lloyd hysteresis index between discharge and SSC.

    The rising and falling limbs of the hydrograph are normalized and compared
    at common discharge levels. Positive or negative values indicate the
    relative position of the sediment-concentration loop.
    """
    discharge_m3_s = np.asarray(discharge_m3_s)
    sediment_concentration_g_l = np.asarray(sediment_concentration_g_l)

    if len(discharge_m3_s) < 5 or np.nanmax(discharge_m3_s) == np.nanmin(discharge_m3_s):
        return np.nan

    discharge_normalized = (discharge_m3_s - np.nanmin(discharge_m3_s)) / (
        np.nanmax(discharge_m3_s) - np.nanmin(discharge_m3_s)
    )
    concentration_normalized = (
        sediment_concentration_g_l - np.nanmin(sediment_concentration_g_l)
    ) / (np.nanmax(sediment_concentration_g_l) - np.nanmin(sediment_concentration_g_l))

    peak_index = np.nanargmax(discharge_m3_s)
    discharge_rising = discharge_normalized[: peak_index + 1]
    concentration_rising = concentration_normalized[: peak_index + 1]
    discharge_falling = discharge_normalized[peak_index:]
    concentration_falling = concentration_normalized[peak_index:]

    minimum_common_discharge = max(discharge_rising[0], discharge_falling[-1])
    discharge_levels = np.linspace(minimum_common_discharge, 1.0, n_levels)

    hysteresis_values = []
    for discharge_level in discharge_levels:
        concentration_rise = np.interp(
            discharge_level, discharge_rising, concentration_rising
        )
        concentration_fall = np.interp(
            discharge_level, discharge_falling[::-1], concentration_falling[::-1]
        )
        hysteresis_values.append(concentration_rise - concentration_fall)

    return np.nanmean(hysteresis_values)



PARAMETER_NAMES = [
    "alpha_SR",
    "alpha_GW",
    "T_Q",
    "a_hillslope",
    "T_hillslope",
    "a_remob",
    "T_remob",
    "Qc_remob",
    "a_depo",
    "Qc_depo",
]


def simulate_hydrosedimentary_event(
    rainfall_mm_min,
    initial_deposit_g,
    parameters,
    deposition=True,
    remobilization=True,
    groundwater_reduction_fraction=0.0,
):
    """Simulate one complete hydrosedimentary event.

    Parameters
    ----------
    rainfall_mm_min : array-like
        Rainfall intensity time series (mm/min).
    initial_deposit_g : float
        Initial channel sediment stock (g).
    parameters : sequence of 10 floats
        Model parameters in the following order:
        PARAMETER_NAMES.
    deposition : bool, default=True
        Whether channel deposition is activated.
    remobilization : bool, default=True
        Whether channel sediment remobilization is activated.
    groundwater_reduction_fraction : float, default=0.0
        Fractional reduction in groundwater contribution, between 0 and 1.

    Returns
    -------
    erosion : dict
        Simulated discharge and sediment time series.
    float
        Final channel sediment stock (tonnes).
    list
        Nine model criteria used by GLUE, in the order defined by
        CRITERIA_NAMES.
    """
    import numpy as np

    # --------------------------------------------------------------
    # 1. Validate and unpack parameters
    # --------------------------------------------------------------
    parameters = np.asarray(parameters, dtype=float)

    if parameters.size != len(PARAMETER_NAMES):
        raise ValueError(
            f"Expected {len(PARAMETER_NAMES)} parameters "
            f"({PARAMETER_NAMES}), received {parameters.size}."
        )

    if not 0.0 <= groundwater_reduction_fraction <= 1.0:
        raise ValueError(
            "groundwater_reduction_fraction must be between 0 and 1."
        )

    (
        alpha_SR,
        alpha_GW,
        T_Q,
        a_hillslope,
        T_hillslope,
        a_remob,
        T_remob,
        Qc_remob,
        a_depo,
        Qc_depo,
    ) = parameters

    rainfall_mm_min = np.asarray(rainfall_mm_min, dtype=float)

    if rainfall_mm_min.ndim != 1 or rainfall_mm_min.size == 0:
        raise ValueError("rainfall_mm_min must be a non-empty 1D array.")

    if not np.all(np.isfinite(rainfall_mm_min)):
        raise ValueError("rainfall_mm_min contains non-finite values.")

    if np.any(rainfall_mm_min < 0):
        raise ValueError("Rainfall intensity cannot be negative.")

    if initial_deposit_g < 0:
        raise ValueError("initial_deposit_g cannot be negative.")

    T_Q = int(T_Q)
    T_hillslope = int(T_hillslope)
    T_remob = int(T_remob)

    # Preserve the original groundwater coefficient for the correction.
    alpha_GW_original = alpha_GW
    alpha_GW_effective = alpha_GW * (
        1.0 - groundwater_reduction_fraction
    )

    # --------------------------------------------------------------
    # 2. Hydrology
    # --------------------------------------------------------------
    surface_runoff_m3_s = simulate_reservoir_discharge(
        rainfall_mm_min,
        alpha_SR,
        T_Q,
        catchment_area_km2=20,
        time_step_min=1,
        rainfall_threshold_mm_min=0,
    )

    groundwater_m3_s = simulate_reservoir_discharge(
        rainfall_mm_min,
        alpha_GW_effective,
        T_Q,
        catchment_area_km2=20,
        time_step_min=1,
        rainfall_threshold_mm_min=0,
    )

    original_groundwater_m3_s = simulate_reservoir_discharge(
        rainfall_mm_min,
        alpha_GW_original,
        T_Q,
        catchment_area_km2=20,
        time_step_min=1,
        rainfall_threshold_mm_min=0,
    )

    n_steps = min(
        len(surface_runoff_m3_s),
        len(groundwater_m3_s),
        len(original_groundwater_m3_s),
    )

    surface_runoff_m3_s = surface_runoff_m3_s[:n_steps]
    groundwater_m3_s = groundwater_m3_s[:n_steps]
    original_groundwater_m3_s = original_groundwater_m3_s[:n_steps]

    total_discharge_m3_s = (
        surface_runoff_m3_s + groundwater_m3_s
    )
    original_total_discharge_m3_s = (
        surface_runoff_m3_s + original_groundwater_m3_s
    )

    # --------------------------------------------------------------
    # 3. Hillslope sediment production
    # --------------------------------------------------------------
    hillslope_concentration_g_l, _ = (
        compute_hillslope_sediment_concentration(
            total_discharge_m3_s,
            surface_runoff_m3_s,
            T_hillslope,
            a_hillslope,
        )
    )

    # Preserve the reference sediment flux when groundwater is reduced.
    groundwater_correction_factor = np.divide(
        original_total_discharge_m3_s,
        total_discharge_m3_s,
        out=np.ones_like(total_discharge_m3_s),
        where=total_discharge_m3_s > 0,
    )

    hillslope_concentration_g_l = (
        hillslope_concentration_g_l * groundwater_correction_factor
    )

    hillslope_flux_g_s = (
        hillslope_concentration_g_l * total_discharge_m3_s * 1e3
    )

    # --------------------------------------------------------------
    # 4. Channel sediment remobilization
    # --------------------------------------------------------------
    remaining_deposit_g = float(initial_deposit_g)
    remobilized_concentration_g_l = np.zeros(n_steps)
    remobilized_flux_g_s = np.zeros(n_steps)
    remobilized_mass_t = 0.0

    if remobilization and remaining_deposit_g > 0 and a_remob > 0:
        (
            remobilized_concentration_g_l,
            remobilized_flux_g_s,
            remaining_deposit_g,
        ) = compute_channel_remobilization(
            total_discharge_m3_s,
            remaining_deposit_g,
            a_remob,
            Qc_remob,
            T_remob,
        )

        remobilized_mass_t = (
            np.sum(remobilized_flux_g_s) * 60 * 1e-6
        )

    total_suspended_concentration_g_l = (
        hillslope_concentration_g_l + remobilized_concentration_g_l
    )

    # --------------------------------------------------------------
    # 5. Sediment deposition
    # --------------------------------------------------------------
    deposited_mass_during_event_g = 0.0
    deposition_concentration_g_l = np.zeros(n_steps)

    if deposition and a_depo > 0:
        deposition_rate_g_m2_s = (
            a_depo * total_suspended_concentration_g_l * 1e3
        )

        deposition_rate_g_m2_s *= np.maximum(
            Qc_depo - total_discharge_m3_s,
            0.0,
        )

        deposition_concentration_g_l = (
            deposition_rate_g_m2_s
            / (total_discharge_m3_s * 1e3 + 1e-9)
        )

        outlet_sediment_concentration_g_l = np.maximum(
            total_suspended_concentration_g_l
            - deposition_concentration_g_l,
            0.0,
        )

        deposited_mass_during_event_g = (
            np.sum(
                (
                    total_suspended_concentration_g_l
                    - outlet_sediment_concentration_g_l
                )
                * total_discharge_m3_s
                * 1e3
            )
            * 60
        )

        remaining_deposit_g += deposited_mass_during_event_g

    else:
        outlet_sediment_concentration_g_l = (
            total_suspended_concentration_g_l.copy()
        )

    # --------------------------------------------------------------
    # 6. Model criteria
    # --------------------------------------------------------------
    maximum_discharge_m3_s = np.max(total_discharge_m3_s)
    maximum_surface_runoff_m3_s = np.max(surface_runoff_m3_s)

    total_runoff_volume_m3 = (
        np.sum(total_discharge_m3_s) * 60 * 1e-3
    )
    surface_runoff_volume_m3 = (
        np.sum(surface_runoff_m3_s) * 60 * 1e-3
    )

    maximum_outlet_concentration_g_l = np.max(
        outlet_sediment_concentration_g_l
    )

    exported_sediment_t = (
        np.sum(
            outlet_sediment_concentration_g_l * total_discharge_m3_s
        )
        * 60
        * 1e-3
    )

    hysteresis_index = compute_hysteresis_index(
        total_discharge_m3_s,
        outlet_sediment_concentration_g_l,
    )

    if not deposition:
        deposited_mass_during_event_g = 0.0

    if not remobilization:
        remobilized_mass_t = 0.0

    criteria = [
        maximum_discharge_m3_s,
        maximum_surface_runoff_m3_s,
        total_runoff_volume_m3,
        surface_runoff_volume_m3,
        maximum_outlet_concentration_g_l,
        exported_sediment_t,
        hysteresis_index,
        deposited_mass_during_event_g * 1e-6,
        remobilized_mass_t,
    ]

    # --------------------------------------------------------------
    # 7. Collect time series and return results
    # --------------------------------------------------------------
    erosion = {
        "Q": total_discharge_m3_s,
        "Q_SR": surface_runoff_m3_s,
        "Q_GW": groundwater_m3_s,
        "SSC_versants": hillslope_concentration_g_l,
        "SSC_remob": remobilized_concentration_g_l,
        "SSC_out": outlet_sediment_concentration_g_l,
        "deposition_concentration": deposition_concentration_g_l,
    }

    final_deposit_t = remaining_deposit_g * 1e-6

    return erosion, final_deposit_t, criteria