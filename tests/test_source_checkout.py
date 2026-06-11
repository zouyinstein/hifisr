from __future__ import annotations

import importlib.util
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_source_version_is_0_6_0():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    import hifisr_functions

    assert hifisr_functions.__version__ == "0.6.0"


def test_analysis_scripts_bootstrap_source_checkout(monkeypatch):
    root = str(ROOT)
    monkeypatch.setattr(sys, "path", [path for path in sys.path if path != root])

    spec = importlib.util.spec_from_file_location(
        "_hifisr_bootstrap_test",
        ROOT / "analysis_scripts" / "_bootstrap.py",
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)

    assert sys.path[0] == root


def test_w3_5_2_config_uses_project_results_dir():
    config_text = (ROOT / "workflow" / "config" / "w3_5_2.yaml").read_text()

    assert 'project_dir: "{snakefile_dir}"' in config_text
    assert 'results_dir: "{project_dir}/results"' in config_text
    assert "dev_root" not in config_text
    assert "results_name" not in config_text
    assert "results_3" not in config_text


def test_snakefile_uses_release_target_names():
    snakefile_text = (ROOT / "Snakefile").read_text()

    assert "rule polish_alignment_variant:" in snakefile_text
    assert "rule polish_alignment_variant_review_inputs:" in snakefile_text
    assert "rule verify_corrected_genome:" in snakefile_text
    assert "rule work2:" not in snakefile_text
    assert "rule work3:" not in snakefile_text
    assert "manual_review_inputs" not in snakefile_text
