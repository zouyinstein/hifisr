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


def test_new_architecture_imports_keep_legacy_modules_available():
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))

    architecture_modules = [
        "hifisr_functions.cli.main",
        "hifisr_functions.core.config",
        "hifisr_functions.core.io",
        "hifisr_functions.core.runner",
        "hifisr_functions.graph.coordinate",
        "hifisr_functions.graph.evidence_projection",
        "hifisr_functions.graph.gfa",
        "hifisr_functions.graph.path",
        "hifisr_functions.read_ops.binning",
        "hifisr_functions.read_ops.recruit",
        "hifisr_functions.read_ops.stats",
        "hifisr_functions.report.html",
        "hifisr_functions.report.plots",
        "hifisr_functions.report.tables",
        "hifisr_functions.validation.confidence",
        "hifisr_functions.validation.flag",
        "hifisr_functions.validation.model_selection",
        "hifisr_functions.validation.numt_nupt",
        "hifisr_functions.variant_ops.frequency",
        "hifisr_functions.variant_ops.junction",
        "hifisr_functions.variant_ops.read_linkage",
        "hifisr_functions.variant_ops.snv_indel",
        "hifisr_functions.variant_ops.sv",
        "hifisr_functions.workflow.checkpoints",
        "hifisr_functions.workflow.init",
        "hifisr_functions.workflow.recipes",
    ]
    for module_name in architecture_modules:
        importlib.import_module(module_name)

    import hifisr_functions.base as hfbase
    import hifisr_functions.reads as hfreads
    import hifisr_functions.references as hfref
    import hifisr_functions.reports as hfrps
    import hifisr_functions.transfer as hftrans
    import hifisr_functions.variants as hfvar

    from hifisr_functions.core import config, io, runner
    from hifisr_functions.graph import coordinate, evidence_projection, gfa, path
    from hifisr_functions.read_ops import binning, recruit, stats
    from hifisr_functions.report import plots, tables
    from hifisr_functions.validation import numt_nupt
    from hifisr_functions.variant_ops import junction, snv_indel, sv
    from hifisr_functions.workflow import checkpoints, recipes

    assert Path(hfreads.__file__).name == "reads.py"
    assert Path(hfvar.__file__).name == "variants.py"

    assert config.load_soft_paths is hfbase.load_soft_paths
    assert config.REQUIRED_WORKFLOW_TOOL_GROUPS is hfbase.REQUIRED_WORKFLOW_TOOL_GROUPS
    assert io.get_file_lines is hfbase.get_file_lines
    assert runner.run_checked is hfbase.run_checked
    assert coordinate.get_subseq is hfref.get_subseq
    assert evidence_projection.aln_to_ref is hfref.aln_to_ref
    assert gfa.get_gfa_blastn_png is hfrps.get_gfa_blastn_png
    assert gfa.get_gfa_reference_pdf is hfrps.get_gfa_reference_pdf
    assert path.rotate_ref_to_non_repeat_region is hfref.rotate_ref_to_non_repeat_region
    assert binning.random_sampling is hfreads.random_sampling
    assert recruit.split_mtpt_reads is hfreads.split_mtpt_reads
    assert stats.get_fastq_stats is hfrps.get_fastq_stats
    assert plots.get_gfa_reference_pdf is hfrps.get_gfa_reference_pdf
    assert plots.plot_coverage is hfrps.plot_coverage
    assert tables.convert_blastn_alignments_to_table is hfrps.convert_blastn_alignments_to_table
    assert numt_nupt.run_transfer_blastn is hftrans.run_transfer_blastn
    assert junction.get_next_groups is hfvar.get_next_groups
    assert snv_indel.snv_or_indel is hfvar.snv_or_indel
    assert sv.run_multi_threads_blastn is hfvar.run_multi_threads_blastn
    assert checkpoints.validate_soft_paths is hfbase.validate_soft_paths
    assert recipes.flye_assemble is hfref.flye_assemble
