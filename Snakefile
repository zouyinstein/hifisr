configfile: "workflow/config/w3_5_2_macOS.yaml"

from pathlib import Path
import sys

shell.executable("/bin/bash")


SNAKEFILE_DIR = Path(workflow.basedir).resolve()


def _fmt(value):
    if isinstance(value, str):
        return value.format(**FORMAT_VARS)
    return value


def _path(value):
    value = _fmt(value)
    path = Path(value)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return str(path)


SAMPLE = config["sample"]
BASE_FORMAT_VARS = {
    "sample": SAMPLE,
    "snakefile_dir": str(SNAKEFILE_DIR),
}

project_dir_value = config.get("project_dir", "{snakefile_dir}")
if isinstance(project_dir_value, str):
    project_dir_value = project_dir_value.format(**BASE_FORMAT_VARS)
PROJECT_ROOT = Path(project_dir_value)
if not PROJECT_ROOT.is_absolute():
    PROJECT_ROOT = (SNAKEFILE_DIR / PROJECT_ROOT).resolve()
else:
    PROJECT_ROOT = PROJECT_ROOT.resolve()

BASE_FORMAT_VARS.update({
    "project_dir": str(PROJECT_ROOT),
    "project_root": str(PROJECT_ROOT),
})
results_dir_value = config.get("results_dir", "{project_dir}/results")
if isinstance(results_dir_value, str):
    results_dir_value = results_dir_value.format(**BASE_FORMAT_VARS)
RESULTS_DIR = Path(results_dir_value)
if not RESULTS_DIR.is_absolute():
    RESULTS_DIR = (PROJECT_ROOT / RESULTS_DIR).resolve()
else:
    RESULTS_DIR = RESULTS_DIR.resolve()

FORMAT_VARS = {
    "snakefile_dir": str(SNAKEFILE_DIR),
    "project_dir": str(PROJECT_ROOT),
    "project_root": str(PROJECT_ROOT),
    "results_dir": str(RESULTS_DIR),
    "sample": SAMPLE,
}

SCRIPT_DIR = PROJECT_ROOT / "analysis_scripts"
SOFT_PATHS = _path(config.get("soft_paths", "{project_root}/deps/soft_paths.txt"))
REQUIREMENTS = str(PROJECT_ROOT / "requirements-dev.txt")
SNAKEMAKE_PYTHON = sys.executable


def load_soft_paths(path):
    soft_paths = {}
    path = Path(path)
    if not path.exists():
        return soft_paths
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=1)
        if len(fields) == 2:
            soft_paths[fields[0]] = fields[1]
    return soft_paths


SOFT_PATH_DICT = load_soft_paths(SOFT_PATHS)
PYTHON = _path(config.get("python", SOFT_PATH_DICT.get("python", "{project_root}/.venv/bin/python")))
READS = _path(config["reads"])

GENOMES = config.get("genomes", ["mito", "plastid"])
DRAFT_SIZE = {"mito": "500", "plastid": "150"}
RUN2 = config.get("run2", {}).get("name", "run_2")
RUN3 = config.get("run3", {}).get("name", "run_3")
RUN3_GENOMES = set(config.get("run3", {}).get("genomes", []))

THREADS = config.get("threads", {})
THREADS_READS = int(THREADS.get("reads", 8))
THREADS_DRAFT = int(THREADS.get("draft", 8))
THREADS_POLISH = int(THREADS.get("polish", 8))
THREADS_VARIANTS = int(THREADS.get("variants", 8))
DOWNSTREAM_READ_LIMIT = config.get("downstream_read_limit", {})
DOWNSTREAM_READ_LIMIT_MITO = int(DOWNSTREAM_READ_LIMIT.get("mito", 50000))
DOWNSTREAM_READ_LIMIT_PLASTID = int(DOWNSTREAM_READ_LIMIT.get("plastid", 50000))

STANDARD_DRAFT_ASSEMBLY_MODES = {"mecat_flye", "flye"}
MITO_SIMPLE_DRAFT_ASSEMBLY_MODES = {"ms", "mh", "mx"}
PLASTID_SIMPLE_DRAFT_ASSEMBLY_MODES = {"pl", "ps", "ph"}
SIMPLE_DRAFT_ASSEMBLY_MODES = MITO_SIMPLE_DRAFT_ASSEMBLY_MODES | PLASTID_SIMPLE_DRAFT_ASSEMBLY_MODES
PLASTID_DRAFT_MODE_MAP = {
    "ms": "ps",
    "mh": "ph",
}
DRAFT_MODE_ALIASES = {
    "mecat+flye": "mecat_flye",
    "mecat_flye": "mecat_flye",
    "mecat": "mecat_flye",
    "flye": "flye",
    "all": "flye",
    "ms": "ms",
    "mh": "mh",
    "mx": "mx",
    "pl": "pl",
    "ps": "ps",
    "ph": "ph",
}


def normalize_draft_assembly_modes(value):
    if isinstance(value, str):
        raw_modes = [mode.strip() for mode in value.split(",")]
    else:
        raw_modes = list(value)
    modes = []
    for raw_mode in raw_modes:
        mode = DRAFT_MODE_ALIASES.get(str(raw_mode).strip(), str(raw_mode).strip())
        if not mode:
            continue
        if mode not in STANDARD_DRAFT_ASSEMBLY_MODES | SIMPLE_DRAFT_ASSEMBLY_MODES:
            raise ValueError(f"Unsupported draft assembly mode: {mode}")
        if mode not in modes:
            modes.append(mode)
    if not modes:
        raise ValueError("draft_assembly.modes must contain at least one mode")
    return modes


DRAFT_ASSEMBLY_CFG = config.get("draft_assembly", {})
DRAFT_ASSEMBLY_MODES = normalize_draft_assembly_modes(
    DRAFT_ASSEMBLY_CFG.get("modes", ["ms", "mh", "mx"])
)

REF_CFG = config.get("references", {})
ROTATE_REFS = bool(REF_CFG.get("rotate", False))
DRAFT_EDITED = config.get("draft_edited", {})
POS_REF_ALT = config.get("pos_ref_alt", {})

LOG_DIR = RESULTS_DIR / SAMPLE / "logs" / "snakemake"
TMP_DIR = RESULTS_DIR / ".tmp"
MPLCONFIGDIR = RESULTS_DIR / ".matplotlib"
XDG_CACHE_HOME = RESULTS_DIR / ".cache"


def _path_parent(value):
    value = _fmt(value)
    path = Path(value)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    return str(path.resolve().parent)


def tool_path_dirs():
    dirs = []
    for executable in [PYTHON] + list(SOFT_PATH_DICT.values()):
        if not executable:
            continue
        parent = _path_parent(executable)
        if parent not in dirs:
            dirs.append(parent)
    for path_dir in config.get("path_dirs", []):
        parent = _path(path_dir)
        if parent not in dirs:
            dirs.append(parent)
    return dirs


TOOL_PATH = ":".join(tool_path_dirs())
RUNTIME_DIRS = f'mkdir -p "{LOG_DIR}" "{RESULTS_DIR}" "{MPLCONFIGDIR}" "{XDG_CACHE_HOME}" "{TMP_DIR}"'


def common_env():
    path_export = f'export PATH="{TOOL_PATH}:$PATH"; ' if TOOL_PATH else ""
    return (
        path_export +
        f'export MPLCONFIGDIR="{MPLCONFIGDIR}"; '
        f'export XDG_CACHE_HOME="{XDG_CACHE_HOME}"; '
        f'export HIFISR_TMPDIR="{TMP_DIR}"'
    )


def sample_dir():
    return RESULTS_DIR / SAMPLE


def reads_dir():
    return sample_dir() / "reads"


def draft_dir(genome):
    return sample_dir() / "draft_assembly" / genome


def genome_reads(genome):
    return str(reads_dir() / f"{genome}.fastq.gz")


def selected_ref(genome):
    return _path(REF_CFG[genome])


def rotated_ref(genome):
    return str(sample_dir() / "references" / f"{genome}_rotated.fasta")


def ref_for_workflow(genome):
    if ROTATE_REFS:
        return rotated_ref(genome)
    return selected_ref(genome)


def edited_fasta(genome):
    default_name = {
        "mito": "mito_checked_draft.fasta",
        "plastid": "plastid_checked_draft.fasta",
    }[genome]
    primary = Path(_path(DRAFT_EDITED.get(genome, f"{draft_dir(genome)}/{default_name}")))
    backup = draft_dir(genome) / "backup_info" / primary.name
    if not primary.exists() and backup.exists():
        return str(backup)
    return str(primary)


def pos_ref_alt(genome):
    return _path(POS_REF_ALT.get(genome, f"{draft_dir(genome)}/pos_ref_alt.txt"))


def draft_modes_for_genome(genome):
    modes = []
    for mode in DRAFT_ASSEMBLY_MODES:
        if mode in STANDARD_DRAFT_ASSEMBLY_MODES:
            genome_mode = mode
        elif genome == "mito" and mode in MITO_SIMPLE_DRAFT_ASSEMBLY_MODES:
            genome_mode = mode
        elif genome == "plastid" and mode in PLASTID_SIMPLE_DRAFT_ASSEMBLY_MODES:
            genome_mode = mode
        elif genome == "plastid" and mode in PLASTID_DRAFT_MODE_MAP:
            genome_mode = PLASTID_DRAFT_MODE_MAP[mode]
        else:
            continue
        if genome_mode not in modes:
            modes.append(genome_mode)
    return modes


def draft_modes_arg(genome):
    modes = draft_modes_for_genome(genome)
    if not modes:
        raise ValueError(f"No draft assembly modes are configured for {genome}")
    return ",".join(modes)


def draft_simple_modes_for_genome(genome):
    return [mode for mode in draft_modes_for_genome(genome) if mode in SIMPLE_DRAFT_ASSEMBLY_MODES]


def draft_flye_modes_for_genome(genome):
    return [mode for mode in draft_modes_for_genome(genome) if mode in STANDARD_DRAFT_ASSEMBLY_MODES]


def draft_simple_modes_arg(genome):
    return ",".join(draft_simple_modes_for_genome(genome))


def draft_flye_modes_arg(genome):
    return ",".join(draft_flye_modes_for_genome(genome))


def polish_alignment_variant_aligned(genome):
    return str(draft_dir(genome) / f"{genome}_flye_polish_aligned.fasta")


def polish_alignment_variant_table(genome):
    return str(sample_dir() / genome / RUN2 / "variants_anno_combined_depth_frq_filter.xlsx")


def polish_alignment_variant_done(genome):
    return str(sample_dir() / genome / RUN2 / ".snakemake.done")


def corrected_genome_fasta(genome):
    return str(draft_dir(genome) / f"{genome}_flye_polish_aligned_cor.fasta")


def verify_corrected_genome_table(genome):
    return str(sample_dir() / genome / RUN3 / "variants_anno_combined_depth_frq_filter.xlsx")


def verify_corrected_genome_done(genome):
    return str(sample_dir() / genome / RUN3 / ".snakemake.done")


def draft_outputs(genome):
    size = DRAFT_SIZE[genome]
    d = draft_dir(genome)
    modes = draft_modes_for_genome(genome)
    outputs = []
    if "mecat_flye" in modes:
        outputs.extend([
            str(d / f"mecat_{genome}_{size}K_before_rr.gfa"),
            str(d / f"mecat_{genome}_{size}K_after_rr.gfa"),
            str(d / f"mecat_{genome}_{size}K_before_rr.png"),
            str(d / f"mecat_{genome}_{size}K_after_rr.png"),
        ])
    if "flye" in modes:
        outputs.extend([
            str(d / f"all_{genome}_{size}K_before_rr.gfa"),
            str(d / f"all_{genome}_{size}K_after_rr.gfa"),
            str(d / f"all_{genome}_{size}K_before_rr.png"),
            str(d / f"all_{genome}_{size}K_after_rr.png"),
        ])
    for mode in modes:
        if mode in SIMPLE_DRAFT_ASSEMBLY_MODES:
            outputs.append(str(d / f"simple_draft_asm_{genome}_{size}K_{mode}.gfa"))
            outputs.append(str(d / f"simple_draft_asm_{genome}_{size}K_{mode}.png"))
    return outputs


def clean_intermediate_paths():
    paths = []
    for genome in GENOMES:
        paths.append(Path(genome_reads(genome)))

    variant_runs = [(genome, RUN2) for genome in GENOMES]
    variant_runs.extend((genome, RUN3) for genome in RUN3_GENOMES)
    for genome, run_name in variant_runs:
        run_dir = sample_dir() / genome / run_name
        for filename in [
            "reads.fasta",
            "new_reads.fasta",
            "FL.fasta",
            "partial.fasta",
            "variant_cov.fasta",
        ]:
            paths.append(run_dir / filename)
    return paths


def clean_intermediate_dirs():
    return [TMP_DIR, MPLCONFIGDIR]


POLISH_ALIGNMENT_VARIANT_TARGETS = [
    polish_alignment_variant_table(g) for g in GENOMES
] + [polish_alignment_variant_done(g) for g in GENOMES]
POLISH_ALIGNMENT_TARGETS = []
for genome in GENOMES:
    POLISH_ALIGNMENT_TARGETS.extend([
        str(draft_dir(genome) / f"{genome}_flye_polish_aligned.fasta"),
        str(draft_dir(genome) / f"{genome}_flye_polish_aligned_blastn_alignments.xlsx"),
        str(draft_dir(genome) / f"{genome}_flye_polish_variants.txt"),
    ])
VERIFY_CORRECTED_GENOME_TARGETS = [
    verify_corrected_genome_table(g) for g in RUN3_GENOMES
] + [verify_corrected_genome_done(g) for g in RUN3_GENOMES]
FINAL_TARGETS = [
    verify_corrected_genome_table(g) if g in RUN3_GENOMES else polish_alignment_variant_table(g)
    for g in GENOMES
]


rule all:
    input:
        FINAL_TARGETS


rule check_runtime_dependencies:
    input:
        soft_paths=SOFT_PATHS,
        requirements=REQUIREMENTS,
        script=str(SCRIPT_DIR / "check_runtime_dependencies.py")
    log:
        str(LOG_DIR / "check_runtime_dependencies.log")
    shell:
        r"""
        {RUNTIME_DIRS}
        "{SNAKEMAKE_PYTHON}" "{input.script}" "{input.soft_paths}" \
          --requirements "{input.requirements}" > "{log}" 2>&1
        cat "{log}"
        """


rule references_ready:
    input:
        [ref_for_workflow(g) for g in GENOMES]


rule reads_ready:
    input:
        expand(str(reads_dir() / "{genome}.fastq.gz"), genome=GENOMES)


rule draft_for_manual_edit:
    input:
        draft_outputs("mito") + draft_outputs("plastid")


rule polish_alignment_ready:
    input:
        POLISH_ALIGNMENT_TARGETS


rule polish_alignment_variant:
    input:
        POLISH_ALIGNMENT_VARIANT_TARGETS


rule polish_alignment_variant_review_inputs:
    input:
        POLISH_ALIGNMENT_VARIANT_TARGETS


rule verify_corrected_genome:
    input:
        VERIFY_CORRECTED_GENOME_TARGETS


rule final:
    input:
        FINAL_TARGETS


rule clean:
    message:
        "Remove large regenerable intermediates under {RESULTS_DIR}/{SAMPLE}"
    run:
        import shutil

        def within_results(path):
            path = Path(path).absolute()
            root = RESULTS_DIR.absolute()
            return path == root or root in path.parents

        def size_bytes(path):
            if path.is_dir():
                return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())
            if path.exists():
                return path.stat().st_size
            return 0

        removed = []
        skipped = []
        for path in clean_intermediate_paths():
            path = Path(path)
            if not path.exists():
                continue
            if not within_results(path):
                skipped.append(str(path))
                continue
            removed.append((str(path), size_bytes(path)))
            path.unlink()

        for path in clean_intermediate_dirs():
            path = Path(path)
            if not path.exists():
                continue
            if not within_results(path):
                skipped.append(str(path))
                continue
            removed.append((str(path), size_bytes(path)))
            shutil.rmtree(path)

        freed = sum(size for _, size in removed)
        for path, size in removed:
            print(f"removed\t{size}\t{path}")
        if skipped:
            for path in skipped:
                print(f"skipped_outside_results\t{path}")
        print(f"clean_removed_files_or_dirs={len(removed)}")
        print(f"clean_freed_bytes={freed}")


rule rotate_reference:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: selected_ref(wc.genome)
    output:
        rotated=str(sample_dir() / "references" / "{genome}_rotated.fasta")
    params:
        env=common_env(),
        workdir=lambda wc: str(sample_dir() / "references" / f".rotate_{wc.genome}"),
        script=str(SCRIPT_DIR / "adjust_ref_fasta.py")
    log:
        str(LOG_DIR / "rotate_reference.{genome}.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        mkdir -p "$(dirname "{output.rotated}")"
        rm -rf "{params.workdir}"
        mkdir -p "{params.workdir}"
        cd "{params.workdir}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{wildcards.genome}" "{input.ref}" \
          > "{log}" 2>&1
        rotated="$(ls "{wildcards.genome}"_rotated_*.fasta | head -n 1)"
        cp "$rotated" "{output.rotated}"
        """


rule extract_mtpt_reads:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        reads=READS,
        mito_ref=lambda wc: ref_for_workflow("mito"),
        plastid_ref=lambda wc: ref_for_workflow("plastid")
    output:
        mito_fastq=str(reads_dir() / "mito.fastq.gz"),
        plastid_fastq=str(reads_dir() / "plastid.fastq.gz"),
        mito_stats=str(reads_dir() / "mito_id_length_qual.txt"),
        plastid_stats=str(reads_dir() / "plastid_id_length_qual.txt")
    threads: THREADS_READS
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "get_mtpt_reads.py"),
        mito_read_limit=DOWNSTREAM_READ_LIMIT_MITO,
        plastid_read_limit=DOWNSTREAM_READ_LIMIT_PLASTID
    log:
        str(LOG_DIR / "extract_mtpt_reads.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" \
          "{input.mito_ref}" "{input.plastid_ref}" "{input.reads}" "{threads}" \
          "{params.mito_read_limit}" "{params.plastid_read_limit}" \
          > "{log}" 2>&1
        """


rule draft_assembly_mito:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: ref_for_workflow("mito"),
        reads=lambda wc: genome_reads("mito"),
        full_reads=lambda wc: genome_reads("mito")
    output:
        draft_outputs("mito")
    threads: THREADS_DRAFT
    params:
        env=common_env(),
        simple_modes=lambda wc: draft_simple_modes_arg("mito"),
        flye_modes=lambda wc: draft_flye_modes_arg("mito"),
        script=str(SCRIPT_DIR / "get_draft_assembly.py"),
        flye_script=str(SCRIPT_DIR / "get_draft_assembly_flye.py")
    log:
        str(LOG_DIR / "draft_assembly.mito.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        if [ -n "{params.simple_modes}" ]; then
          "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" mito \
            "{input.ref}" "{input.full_reads}" "{threads}" "{params.simple_modes}" "{input.full_reads}" \
            > "{log}" 2>&1
        fi
        if [ -n "{params.flye_modes}" ]; then
          "{input.python}" "{params.flye_script}" "{input.soft_paths}" "{SAMPLE}" mito \
            "{input.ref}" "{input.reads}" "{threads}" "{params.flye_modes}" \
            >> "{log}" 2>&1
        fi
        """


rule draft_assembly_plastid:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: ref_for_workflow("plastid"),
        reads=lambda wc: genome_reads("plastid"),
        full_reads=lambda wc: genome_reads("plastid")
    output:
        draft_outputs("plastid")
    threads: THREADS_DRAFT
    params:
        env=common_env(),
        simple_modes=lambda wc: draft_simple_modes_arg("plastid"),
        flye_modes=lambda wc: draft_flye_modes_arg("plastid"),
        script=str(SCRIPT_DIR / "get_draft_assembly.py"),
        flye_script=str(SCRIPT_DIR / "get_draft_assembly_flye.py")
    log:
        str(LOG_DIR / "draft_assembly.plastid.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        if [ -n "{params.simple_modes}" ]; then
          "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" plastid \
            "{input.ref}" "{input.full_reads}" "{threads}" "{params.simple_modes}" "{input.full_reads}" \
            > "{log}" 2>&1
        fi
        if [ -n "{params.flye_modes}" ]; then
          "{input.python}" "{params.flye_script}" "{input.soft_paths}" "{SAMPLE}" plastid \
            "{input.ref}" "{input.reads}" "{threads}" "{params.flye_modes}" \
            >> "{log}" 2>&1
        fi
        """


rule polish_alignment:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: ref_for_workflow(wc.genome),
        draft=lambda wc: edited_fasta(wc.genome),
        reads=lambda wc: genome_reads(wc.genome)
    output:
        aligned=str(sample_dir() / "draft_assembly" / "{genome}" / "{genome}_flye_polish_aligned.fasta"),
        table=str(sample_dir() / "draft_assembly" / "{genome}" / "{genome}_flye_polish_aligned_blastn_alignments.xlsx"),
        variants=str(sample_dir() / "draft_assembly" / "{genome}" / "{genome}_flye_polish_variants.txt")
    threads: THREADS_POLISH
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "get_polished_assembly.py")
    log:
        str(LOG_DIR / "polish_alignment.{genome}.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" "{wildcards.genome}" \
          "{input.ref}" "{input.draft}" "{input.reads}" "{threads}" \
          > "{log}" 2>&1
        """


rule polish_alignment_variant_calling:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: polish_alignment_variant_aligned(wc.genome),
        reads=lambda wc: genome_reads(wc.genome)
    output:
        table=str(sample_dir() / "{genome}" / RUN2 / "variants_anno_combined_depth_frq_filter.xlsx"),
        done=touch(str(sample_dir() / "{genome}" / RUN2 / ".snakemake.done"))
    threads: THREADS_VARIANTS
    params:
        env=common_env(),
        run_dir=lambda wc: str(sample_dir() / wc.genome / RUN2),
        script=str(SCRIPT_DIR / "get_variants_in_reads.py")
    log:
        str(LOG_DIR / "polish_alignment_variant_calling.{genome}.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" "{wildcards.genome}" "{RUN2}" \
          "{input.ref}" "{input.reads}" "{threads}" \
          > "{log}" 2>&1
        """


rule correct_reference_with_manual_edits:
    input:
        python=PYTHON,
        ref=lambda wc: polish_alignment_variant_aligned(wc.genome),
        pos_ref_alt=lambda wc: pos_ref_alt(wc.genome)
    output:
        corrected=str(sample_dir() / "draft_assembly" / "{genome}" / "{genome}_flye_polish_aligned_cor.fasta")
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "correct_erroneous.py")
    log:
        str(LOG_DIR / "correct_reference_with_manual_edits.{genome}.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        "{input.python}" "{params.script}" "{input.ref}" "{output.corrected}" "{input.pos_ref_alt}" \
          > "{log}" 2>&1
        """


rule verify_corrected_genome_variants:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: corrected_genome_fasta(wc.genome),
        reads=lambda wc: genome_reads(wc.genome)
    output:
        table=str(sample_dir() / "{genome}" / RUN3 / "variants_anno_combined_depth_frq_filter.xlsx"),
        done=touch(str(sample_dir() / "{genome}" / RUN3 / ".snakemake.done"))
    threads: THREADS_VARIANTS
    params:
        env=common_env(),
        run_dir=lambda wc: str(sample_dir() / wc.genome / RUN3),
        script=str(SCRIPT_DIR / "get_variants_in_reads.py")
    log:
        str(LOG_DIR / "verify_corrected_genome_variants.{genome}.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" "{wildcards.genome}" "{RUN3}" \
          "{input.ref}" "{input.reads}" "{threads}" \
          > "{log}" 2>&1
        """
