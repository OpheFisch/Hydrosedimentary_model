# Hydrosedimentary model

Python code for a hydrosedimentary model used in a scientific study.

The repository contains:
- the model routines in `src/`;
- a script for Monte Carlo/GLUE calibration;
- a script for simulating a sequence of flood events.

## Repository structure

```text
hydrosedimentary_model_public/
├── src/
│   ├── hydrosedimentary_model.py
│   └── routines_modele.py
├── scripts/
│   ├── run_glue.py
│   └── run_event_sequence.py
├── requirements.txt
└── README.md
```

## Installation

From the repository root, create and activate a virtual environment, then install the dependencies:

```bash
python -m venv .venv
```

On macOS or Linux:

```bash
source .venv/bin/activate
```

On Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install the required packages:

```bash
python -m pip install -r requirements.txt
```

## Running the scripts

Run commands from the repository root.

### GLUE calibration

```bash
python scripts/run_glue.py
```

The number of simulations and the experiment settings are defined in `scripts/run_glue.py`. Review these settings before launching a long run.

### Sequence of flood events

```bash
python scripts/run_event_sequence.py
```

The settings for this experiment are defined near the beginning of `scripts/run_event_sequence.py`, including the number of years, the number of simulations, and the input/output directories.

This script expects GLUE result CSV files in `results/glue/`. These files are not included in the repository and must be generated or supplied separately. The expected filenames and parameter columns are specified in the script.

The scripts currently use settings defined in the Python files; they do not accept command-line options such as `--glue-dir`, `--output`, or `--seed`.

## Results

Results are written to the locations configured in the scripts. Check those settings before running the experiments. Simulation outputs and GLUE CSV files are not included by default.

## Citation

If this code is used in a publication, cite the associated scientific study. When a stable release is available, consider archiving that version (for example, with Zenodo) and adding its DOI here.
