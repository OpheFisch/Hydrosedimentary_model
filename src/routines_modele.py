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
        remaining_deposit_g -= remobilized_flux_g_s[step] * 60
        remaining_deposit_g = max(remaining_deposit_g, 0)

        # Once the available deposit is exhausted, no further remobilization
        # can occur during this event.
        if remaining_deposit_g == 0:
            shifted_concentration[step:] = 0
            remobilized_flux_g_s[step:] = 0
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


def simulate_hydrosedimentary_event(
    rainfall_mm_min,
    initial_deposit_g,
    parameters,
    deposition=True,
    remobilization=True,
    groundwater_reduction_fraction=0,
):
    """Run one complete hydrosedimentary event.

    The simulation is organized into four physical steps:

    1. **Hydrology:** surface runoff and groundwater discharge are simulated
       with two linear reservoirs.
    2. **Hillslope erosion:** sediment concentration is linked to surface runoff
       and delayed by the hillslope-to-channel transfer time.
    3. **Channel remobilization:** previously deposited sediment is released
       when discharge exceeds the remobilization threshold.
    4. **Deposition:** suspended sediment is removed according to the
       deposition relationship and the resulting outlet concentration is
       calculated.

    Parameters
    ----------
    rainfall_mm_min : array-like
        Rainfall intensity through time (mm/min).
    initial_deposit_g : float
        Initial mass of sediment stored in the channel (g).
    parameters : sequence of 11 floats
        Model parameters in the following order::

        alpha_SR, tps_SR, alpha_GW, tps_GW, a_SR,
        tps_sed_versants, a_lit, tps_sed_lit, Qc_remob,
        a_depot, Qc_depot
    deposition : bool, default=True
        Activate the deposition process.
    remobilization : bool, default=True
        Activate channel remobilization.
    groundwater_reduction_fraction : float, default=0
        Fractional reduction applied to groundwater runoff.

    Returns
    -------
    erosion : dict
        Simulated discharge and sediment time series.
    float
        Final channel sediment stock in tonnes.
    list
        Ten model criteria used by GLUE.
    """
    (
        runoff_coefficient_surface,
        response_time_surface,
        runoff_coefficient_groundwater,
        response_time_groundwater,
        hillslope_coefficient,
        hillslope_transfer_time,
        channel_coefficient,
        channel_transfer_time,
        critical_remobilization_discharge,
        deposition_coefficient,
        critical_deposition_discharge,
    ) = parameters

    # Groundwater reduction is applied only to the perturbed groundwater
    # contribution. The unmodified groundwater series is retained so that the
    # model can account for the effect of the reduction on sediment flux.
    original_groundwater_coefficient = runoff_coefficient_groundwater
    runoff_coefficient_groundwater *= 1 - groundwater_reduction_fraction

    response_time_surface = int(response_time_surface)
    response_time_groundwater = int(response_time_groundwater)
    hillslope_transfer_time = int(hillslope_transfer_time)
    channel_transfer_time = int(channel_transfer_time)

    # ------------------------------------------------------------------
    # 1. HYDROLOGY
    # ------------------------------------------------------------------
    surface_runoff_m3_s = simulate_reservoir_discharge(
        rainfall_mm_min,
        runoff_coefficient_surface,
        response_time_surface,
        catchment_area_km2=20,
        time_step_min=1,
        rainfall_threshold_mm_min=0,
    )
    groundwater_m3_s = simulate_reservoir_discharge(
        rainfall_mm_min,
        runoff_coefficient_groundwater,
        response_time_groundwater,
        catchment_area_km2=20,
        time_step_min=1,
        rainfall_threshold_mm_min=0,
    )
    original_groundwater_m3_s = simulate_reservoir_discharge(
        rainfall_mm_min,
        original_groundwater_coefficient,
        response_time_groundwater,
        catchment_area_km2=20,
        time_step_min=1,
        rainfall_threshold_mm_min=0,
    )

    n_steps = min(len(surface_runoff_m3_s), len(groundwater_m3_s))
    surface_runoff_m3_s = surface_runoff_m3_s[:n_steps]
    groundwater_m3_s = groundwater_m3_s[:n_steps]
    original_groundwater_m3_s = original_groundwater_m3_s[:n_steps]

    total_discharge_m3_s = surface_runoff_m3_s + groundwater_m3_s
    original_total_discharge_m3_s = surface_runoff_m3_s + original_groundwater_m3_s

    # ------------------------------------------------------------------
    # 2. HILLSLOPE SEDIMENT PRODUCTION
    # ------------------------------------------------------------------
    hillslope_concentration_g_l, _ = compute_hillslope_sediment_concentration(
        total_discharge_m3_s,
        surface_runoff_m3_s,
        hillslope_transfer_time,
        hillslope_coefficient,
    )

    # Correct the sediment concentration for the groundwater perturbation by
    # preserving the original sediment flux relative to the reference flow.
    with np.errstate(divide="ignore", invalid="ignore"):
        groundwater_correction_factor = np.divide(
            original_total_discharge_m3_s,
            total_discharge_m3_s,
            out=np.ones_like(total_discharge_m3_s),
            where=total_discharge_m3_s > 0,
        )
    hillslope_concentration_g_l *= groundwater_correction_factor
    hillslope_flux_g_s = hillslope_concentration_g_l * total_discharge_m3_s * 1e3

    # ------------------------------------------------------------------
    # 3. CHANNEL REMOBILIZATION
    # ------------------------------------------------------------------
    remaining_deposit_g = initial_deposit_g
    remobilized_concentration_g_l = np.zeros(n_steps)
    remobilized_flux_g_s = np.zeros(n_steps)
    remobilized_mass_t = 0.0

    if remobilization and initial_deposit_g > 0 and channel_coefficient > 0:
        (
            remobilized_concentration_g_l,
            remobilized_flux_g_s,
            remaining_deposit_g,
        ) = compute_channel_remobilization(
            total_discharge_m3_s,
            initial_deposit_g,
            channel_coefficient,
            critical_remobilization_discharge,
            channel_transfer_time,
        )
        remobilized_mass_t = np.sum(remobilized_flux_g_s) * 60 * 1e-6

    total_suspended_concentration_g_l = (
        hillslope_concentration_g_l + remobilized_concentration_g_l
    )

    # ------------------------------------------------------------------
    # 4. DEPOSITION
    # ------------------------------------------------------------------
    deposited_mass_during_event_g = 0.0

    if deposition and deposition_coefficient > 0:
        # Deposition is expressed as a mass flux per unit area (g/m2/s).
        deposition_rate_g_m2_s = (
            deposition_coefficient * total_suspended_concentration_g_l * 1e3
        )
        deposition_rate_g_m2_s[total_discharge_m3_s > critical_deposition_discharge] = 0
        deposition_rate_g_m2_s *= (
            critical_deposition_discharge - total_discharge_m3_s
        )
        deposition_rate_g_m2_s[deposition_rate_g_m2_s < 0] = 0

        # Convert the deposition rate into a concentration-equivalent quantity
        # before subtracting it from the suspended sediment concentration.
        deposition_concentration = deposition_rate_g_m2_s / (
            total_discharge_m3_s * 1e3 + 1e-9
        )

        outlet_sediment_concentration_g_l = np.maximum(
            total_suspended_concentration_g_l - deposition_concentration,
            0,
        )

        deposited_mass_during_event_g = (
            np.sum(
                (total_suspended_concentration_g_l - outlet_sediment_concentration_g_l)
                * total_discharge_m3_s
                * 1e3
            )
            * 60
        )
        remaining_deposit_g += deposited_mass_during_event_g
    else:
        outlet_sediment_concentration_g_l = total_suspended_concentration_g_l
        deposition_concentration = np.zeros(n_steps)

    # ------------------------------------------------------------------
    # 5. MODEL CRITERIA
    # ------------------------------------------------------------------
    maximum_discharge_m3_s = np.max(total_discharge_m3_s)
    maximum_surface_runoff_m3_s = np.max(surface_runoff_m3_s)
    total_runoff_volume_m3 = np.sum(total_discharge_m3_s) * 60 * 1e-3
    surface_runoff_volume_m3 = np.sum(surface_runoff_m3_s) * 60 * 1e-3
    maximum_outlet_concentration_g_l = np.max(outlet_sediment_concentration_g_l)
    exported_sediment_t = (
        np.sum(outlet_sediment_concentration_g_l * total_discharge_m3_s)
        * 60
        * 1e-3
    )
    hysteresis_index = compute_hysteresis_index(
        total_discharge_m3_s, outlet_sediment_concentration_g_l
    )

    # Phase lag is currently fixed to zero in the original model.
    phase_lag_min = 0

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
        phase_lag_min,
        deposited_mass_during_event_g * 1e-6,
        remobilized_mass_t,
    ]

    erosion = {
        "Q": total_discharge_m3_s,
        "Q_SR": surface_runoff_m3_s,
        "Q_GW": groundwater_m3_s,
        "SSC_versants": hillslope_concentration_g_l,
        "SSC_remob": remobilized_concentration_g_l,
        "SSC_out": outlet_sediment_concentration_g_l,
        "deposition_concentration": deposition_concentration,
    }

    return erosion, remaining_deposit_g * 1e-6, criteria
