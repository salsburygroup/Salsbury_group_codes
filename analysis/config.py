"""Configuration utilities for the analysis package."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict
import yaml

_DEFAULT_CONFIG = {
    "analysis_dir": "ANALYSIS",
    "conda_path": "/home/USER/miniconda3/bin/",
    "script_path": "group_python",
}


def load_config(path: str | Path | None = None) -> Dict[str, Any]:
    """Load configuration from *path* or return defaults."""
    if path is None:
        path = Path("config.yaml")
    path = Path(path)
    if path.exists():
        with open(path, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
    else:
        data = {}
    cfg = _DEFAULT_CONFIG.copy()
    cfg.update(data)
    return cfg
