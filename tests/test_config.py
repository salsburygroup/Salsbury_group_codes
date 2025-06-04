from pathlib import Path
from analysis.config import load_config

def test_load_config(tmp_path):
    cfg_file = tmp_path / "cfg.yaml"
    cfg_file.write_text("analysis_dir: TEST")
    cfg = load_config(cfg_file)
    assert cfg["analysis_dir"] == "TEST"
    assert "conda_path" in cfg

