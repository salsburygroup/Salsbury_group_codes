# Salsbury Analysis

This repository contains utilities for processing and analyzing molecular dynamics trajectories.
It bundles a small `analysis` package with reusable helpers for clustering.

## Installation

Install the package and its dependencies with `pip`:

```bash
pip install -e .
```

## Usage

The main entry points are two scripts:

- `cluster_HDBSCAN.py` – perform HDBSCAN clustering on a trajectory
- `cluster_traj.py` – extract the most populated clusters and save representative structures
- `process.py` – run a simple pipeline using a YAML configuration

Run `python cluster_HDBSCAN.py --help` for detailed CLI options. The pipeline
script accepts `--config` pointing to a YAML file (defaults to `config.yaml`).

## Development

A minimal test suite using `pytest` lives in the `tests/` directory. Run the
checks along with syntax validation via:

```bash
python -m py_compile $(git ls-files '*.py')
pytest -q
```

### Configuration

A small YAML file ``config.yaml`` controls default paths used by the pipeline:

```yaml
analysis_dir: ANALYSIS
conda_path: /home/USER/miniconda3/bin/
script_path: group_python
```

Override these values by passing ``--config myconfig.yaml`` to ``process.py``.
