# Hydrosedimentary model

Python implementation of the event-scale hydrosedimentary model used in the accompanying scientific study.

## Model workflow

For each rainfall event, the model follows four main steps:

1. **Hydrology** — rainfall is transformed into surface runoff and groundwater discharge using linear reservoirs.
2. **Hillslope sediment production** — sediment concentration is related to surface runoff and shifted according to the hillslope-to-channel transfer time.
3. **Channel remobilization** — previously deposited sediment is remobilized when discharge exceeds a critical threshold.
4. **Deposition** — suspended sediment is removed according to the deposition relationship, producing the outlet sediment concentration.

The model then computes discharge, sediment-export, deposition and hysteresis criteria used in Monte Carlo / GLUE calibration.

## Repository structure

```text
hydrosedimentary_model_public/
├── README.md
├── requirements.txt
├── scripts/
│   └── run_glue.py
└── src/
    ├── __init__.py
    ├── hydrosedimentary_model.py
    └── routines_modele.py
```

`routines_modele.py` contains the numerical model core. `hydrosedimentary_model.py` contains parameter handling and Monte Carlo / GLUE utilities. 

## Installation

```bash
python -m pip install -r requirements.txt
```

## Example

From the repository root:

```bash
python scripts/run_glue.py
```

The example uses a small number of simulations. For the production experiments, adjust `N` in `scripts/run_glue.py` to the number used in the study.

## Main model function

The complete event simulation is performed by:

```python
from routines_modele import simulate_hydrosedimentary_event
```

The main calculation is documented directly in the source code, including the physical meaning and units of the intermediate variables.

## Important quantities

The model uses explicit names for the main physical quantities, including:

- `total_discharge_m3_s`
- `surface_runoff_m3_s`
- `groundwater_m3_s`
- `hillslope_concentration_g_l`
- `remobilized_concentration_g_l`
- `deposition_rate_g_m2_s`
- `deposition_concentration`
- `outlet_sediment_concentration_g_l`

## Reproducibility

The repository is intended to contain the code required to reproduce the model calculations. Site-specific input data and calibration results should be archived separately and linked to the publication or DOI.
