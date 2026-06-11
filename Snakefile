configfile: "workflow/config/w3_5_2.yaml"

from pathlib import Path

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
    return str(reads_dir() / f"{SAMPLE}_{genome}.fastq")


def sample_reads(genome):
    return str(reads_dir() / "sample_reads" / f"sample_4000_{genome}.fastq")


def filt_id_qual(genome):
    return str(reads_dir() / "filt_reads" / f"filt_L10K_{genome}_id_length_qual.txt")


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
        "mito": "all_mito_500K_after_rr.edited.fasta",
        "plastid": "all_plastid_150K_after_rr.edited.fasta",
    }[genome]
    return _path(DRAFT_EDITED.get(genome, f"{draft_dir(genome)}/{default_name}"))


def pos_ref_alt(genome):
    return _path(POS_REF_ALT.get(genome, f"{draft_dir(genome)}/pos_ref_alt.txt"))


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
    return [
        str(d / f"mecat_{genome}_{size}K_before_rr.gfa"),
        str(d / f"mecat_{genome}_{size}K_after_rr.gfa"),
        str(d / f"all_{genome}_{size}K_before_rr.gfa"),
        str(d / f"all_{genome}_{size}K_after_rr.gfa"),
        str(d / f"mecat_{genome}_{size}K_before_rr.png"),
        str(d / f"mecat_{genome}_{size}K_after_rr.png"),
        str(d / f"all_{genome}_{size}K_before_rr.png"),
        str(d / f"all_{genome}_{size}K_after_rr.png"),
    ]


def clean_intermediate_paths():
    paths = []
    for genome in GENOMES:
        paths.append(Path(genome_reads(genome)))
        paths.append(Path(sample_reads(genome)))

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


rule references_ready:
    input:
        [ref_for_workflow(g) for g in GENOMES]


rule reads_ready:
    input:
        expand(str(reads_dir() / "{sample}_{genome}.fastq"), sample=[SAMPLE], genome=GENOMES),
        expand(str(reads_dir() / "sample_reads" / "sample_4000_{genome}.fastq"), genome=GENOMES)


rule draft_for_manual_edit:
    input:
        draft_outputs("mito") + draft_outputs("plastid")


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
        mito_fastq=str(reads_dir() / f"{SAMPLE}_mito.fastq"),
        plastid_fastq=str(reads_dir() / f"{SAMPLE}_plastid.fastq"),
        mito_stats=str(reads_dir() / "mito_id_length_qual.txt"),
        plastid_stats=str(reads_dir() / "plastid_id_length_qual.txt")
    threads: THREADS_READS
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "get_mtpt_reads.py")
    log:
        str(LOG_DIR / "extract_mtpt_reads.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" \
          "{input.mito_ref}" "{input.plastid_ref}" "{input.reads}" "{threads}" \
          > "{log}" 2>&1
        """


rule filter_read_ids:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        mito_stats=str(reads_dir() / "mito_id_length_qual.txt"),
        plastid_stats=str(reads_dir() / "plastid_id_length_qual.txt")
    output:
        mito_id_qual=filt_id_qual("mito"),
        plastid_id_qual=filt_id_qual("plastid"),
        mito_ids=str(reads_dir() / "filt_reads" / "filt_L10K_mito_ids.txt"),
        plastid_ids=str(reads_dir() / "filt_reads" / "filt_L10K_plastid_ids.txt")
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "filt_read_ids.py")
    log:
        str(LOG_DIR / "filter_read_ids.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" \
          "{input.mito_stats}" "{input.plastid_stats}" 10000 0 \
          > "{log}" 2>&1
        """


rule sample_reads:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        id_qual=lambda wc: filt_id_qual(wc.genome),
        genome_fastq=lambda wc: genome_reads(wc.genome)
    output:
        fastq=str(reads_dir() / "sample_reads" / "sample_4000_{genome}.fastq"),
        id_qual=str(reads_dir() / "sample_reads" / "sample_4000_{genome}_id_length_qual.txt")
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "sample_reads.py")
    log:
        str(LOG_DIR / "sample_reads.{genome}.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" "{wildcards.genome}" \
          "{input.id_qual}" 4000 \
          > "{log}" 2>&1
        """


rule draft_assembly_mito:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: ref_for_workflow("mito"),
        reads=sample_reads("mito")
    output:
        draft_outputs("mito")
    threads: THREADS_DRAFT
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "get_draft_assembly.py")
    log:
        str(LOG_DIR / "draft_assembly.mito.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" mito \
          "{input.ref}" "{input.reads}" "{threads}" \
          > "{log}" 2>&1
        """


rule draft_assembly_plastid:
    input:
        soft_paths=SOFT_PATHS,
        python=PYTHON,
        ref=lambda wc: ref_for_workflow("plastid"),
        reads=sample_reads("plastid")
    output:
        draft_outputs("plastid")
    threads: THREADS_DRAFT
    params:
        env=common_env(),
        script=str(SCRIPT_DIR / "get_draft_assembly.py")
    log:
        str(LOG_DIR / "draft_assembly.plastid.log")
    shell:
        r"""
        {params.env}
        {RUNTIME_DIRS}
        cd "{RESULTS_DIR}"
        "{input.python}" "{params.script}" "{input.soft_paths}" "{SAMPLE}" plastid \
          "{input.ref}" "{input.reads}" "{threads}" \
          > "{log}" 2>&1
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
        reads=lambda wc: sample_reads(wc.genome)
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
        reads=lambda wc: sample_reads(wc.genome)
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
