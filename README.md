# Hydrosedimentary model

Python implementation of the hydrosedimentary model used in the associated
scientific study.

The repository contains the model core, the Monte Carlo/GLUE calibration
workflow, and a script for simulating a succession of hydrosedimentary flood
events.

## Repository contents

```text
hydrosedimentary-model/
├── src/
│   ├── __init__.py
│   ├── hydrosedimentary_model.py
│   └── routines_modele.py
├── scripts/
│   ├── run_glue.py
│   └── run_event_sequence.py
├── requirements.txt
└── README.md
```

### Main modules

- `src/hydrosedimentary_model.py`: model-level functions, construction of
  parameter vectors, criterion validation, and Monte Carlo/GLUE simulations.
- `src/routines_modele.py`: numerical routines for rainfall generation,
  hydrological response, hillslope sediment production, channel
  remobilisation, deposition, and hysteresis calculations.
- `scripts/run_glue.py`: example entry point for a Monte Carlo/GLUE
  simulation.
- `scripts/run_event_sequence.py`: simulation of successive flood events.
  The sediment-deposit stock is transferred from one event to the next.

The repository contains only the routines required by the hydrosedimentary
model. Auxiliary routines from the original research workspace are not
included because they are not dependencies of the model core.

## Installation

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows:

```powershell
.venv\Scripts\activate
```

## Running a GLUE simulation

A small example is provided for testing the repository:

```bash
python scripts/run_glue.py
```

The example uses `N=1000`. The production number of simulations can be
changed in the script according to the protocol used in the study.

## Running a succession of floods

`scripts/run_event_sequence.py` propagates the sediment-deposit stock through
a sequence of flood events. Six characteristic flood types are represented:

- M-S: minor summer flood
- I-S: intermediate summer flood
- E-S: major summer flood
- M-W: minor winter flood
- I-W: intermediate winter flood
- E-W: major winter flood

The script samples one accepted parameter vector from the corresponding GLUE
distribution for each event and then simulates the event with the sediment
stock left by the previous event.

The default configuration applies the rainfall-frequency and intensity
changes used in the RCP 8.1-9.1 scenario of the study.

The GLUE CSV files required by the succession experiment must be placed in
the directory supplied with `--glue-dir`. For example:

```bash
python scripts/run_event_sequence.py \
    --glue-dir results/glue \
    --output results/event_sequence.csv
```

For a quick test:

```bash
python scripts/run_event_sequence.py \
    --glue-dir results/glue \
    --output results/event_sequence_test.csv \
    --n-years 1 \
    --n-simulations 1 \
    --seed 42
```

The `--seed` option can be used to make the random sequence reproducible.

The required GLUE files are not included in this repository unless they are
explicitly made available with the software release.

## Results

Simulation outputs are written to the directory specified by `--output`.
The `results/` directory is generated at runtime and is not required for
importing the model.

## Reproducibility and provenance

The original research code contained personal absolute paths and experimental
runs executed directly when scripts were imported. These elements were
removed from the public version so that the model can be used from a clean
repository without access to the original author's filesystem.

The public version focuses on the model components and analysis scripts
required for the associated scientific study. Function and variable names
have been expanded where appropriate to make the workflow easier to read and
reuse.

## Citation

When this repository is associated with a published article, the exact
version used for the study should be archived with a version-specific DOI,
for example through Zenodo. The DOI should then be added to this section and
to the manuscript's software/data availability statement.
