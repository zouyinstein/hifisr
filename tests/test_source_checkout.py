# HiFiSR module guide:
# - base: command, file, and soft_paths helpers; import hifisr_functions.base as hfbase
# - reads: read extraction, filtering, sampling, and correction; import hifisr_functions.reads as hfreads
# - references: reference rotation, assembly, polishing, and alignment; import hifisr_functions.references as hfref
# - variants: read-variant calling, grouping, and frequency analysis; import hifisr_functions.variants as hfvar
# - transfer: organelle/nuclear transfer-fragment analysis; import hifisr_functions.transfer as hftrans
# - annotations: annotation tables and feature-level summaries; import hifisr_functions.annotations as hfanno
# - reports: read statistics, plots, Excel tables, and report outputs; import hifisr_functions.reports as hfrps

from __future__ import annotations

import ast
import importlib
import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]


def test_source_version_is_0_6_1():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    import hifisr_functions

    assert hifisr_functions.__version__ == "0.6.1"


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
    config_text = (ROOT / "workflow" / "config" / "w3_5_2_macOS.yaml").read_text()

    assert 'project_dir: "{snakefile_dir}"' in config_text
    assert 'results_dir: "{project_dir}/results"' in config_text
    assert 'soft_paths: "{project_dir}/deps/soft_paths_macOS.txt"' in config_text
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


def test_soft_paths_validation_checks_executables(tmp_path):
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    import hifisr_functions.base as hfbase

    soft_paths = tmp_path / "soft_paths.txt"
    soft_paths.write_text(f"python\t{sys.executable}\n")
    loaded = hfbase.load_soft_paths(
        str(soft_paths),
        validate=True,
        required_tools=["python"],
        print_paths=False,
    )

    assert loaded["python"] == sys.executable

    bad_soft_paths = tmp_path / "bad_soft_paths.txt"
    bad_soft_paths.write_text("missing_tool\t/no/such/tool\n")
    with pytest.raises(RuntimeError, match="path does not exist"):
        hfbase.load_soft_paths(str(bad_soft_paths), validate=True, print_paths=False)


def test_python_package_check_uses_requirements(tmp_path):
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    import hifisr_functions.base as hfbase

    requirements = tmp_path / "requirements.txt"
    requirements.write_text("pytest\n")
    report = hfbase.check_python_packages(str(requirements), sys.executable)

    assert report["packages_checked"] == 1
    assert report["missing"] == []


def test_legacy_python_packages_are_in_runtime_requirements():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    import hifisr_functions.base as hfbase

    expected_packages = {
        "biopython",
        "pysam",
        "pandas",
        "numpy",
        "openpyxl",
        "xlsxwriter",
        "matplotlib",
        "polars",
        "fastexcel",
    }
    requirements = hfbase.parse_python_requirements(ROOT / "requirements-dev.txt")
    package_names = {requirement["name"].lower() for requirement in requirements}

    assert expected_packages <= package_names


def test_conda_installable_tools_are_in_environment_yml():
    expected_tools = {
        "minimap2",
        "samtools",
        "seqkit",
        "blast",
        "bcftools",
        "bamtools",
        "pigz",
        "bandage",
        "hifiasm",
        "flye",
        "canu",
    }
    environment_text = (ROOT / "environment.yml").read_text()

    for tool in expected_tools:
        assert f"  - {tool}\n" in environment_text
    assert "  - mecat\n" not in environment_text


def test_hifisr_functions_have_purity_metadata():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    valid_markers = {"pure", "impure"}
    functions_dir = ROOT / "hifisr_functions"

    for path in sorted(functions_dir.glob("*.py")):
        if path.name == "__init__.py":
            continue

        module = importlib.import_module(f"hifisr_functions.{path.stem}")
        purity = getattr(module, "FUNCTION_PURITY", None)
        assert isinstance(purity, dict) and purity, f"{path.name} missing FUNCTION_PURITY"
        assert set(purity.values()) <= valid_markers, f"{path.name} has invalid markers"

        tree = ast.parse(path.read_text())
        public_defs = {
            node.name
            for node in ast.iter_child_nodes(tree)
            if isinstance(node, (ast.FunctionDef, ast.ClassDef))
            and not node.name.startswith("_")
        }
        if path.name == "transfer.py":
            public_defs.update(getattr(module, "__all__", []))

        assert set(purity) == public_defs, f"{path.name} purity metadata is incomplete"
