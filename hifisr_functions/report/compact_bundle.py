#!/usr/bin/env python3
"""Generate compact hifisr report bundles with HTML project entry points."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import html
import json
import math
import os
import platform
from pathlib import Path
import shutil
import time
from typing import Any

from hifisr_functions import __version__ as HIFISR_VERSION


GENOMES = {
    "mito": {
        "run": "run_3",
        "verified_prefix": "verified_mito",
        "reference_pattern": "mito_rotated_*.fasta",
        "raw_gfa": "merged_mito_raw.gfa",
        "raw_pdf": "merged_mito_raw.pdf",
    },
    "plastid": {
        "run": "run_2",
        "verified_prefix": "verified_plastid",
        "reference_pattern": "plastid_rotated_*.fasta",
        "raw_gfa": "merged_plastid_raw.gfa",
        "raw_pdf": "merged_plastid_raw.pdf",
    },
}

REPORT_SUBDIRS = ["final", "metadata"]
DEFAULT_FINAL_COVERAGE_Y_RULE = "twice_mean"


def ensure_dirs(report_dir: Path) -> None:
    for subdir in REPORT_SUBDIRS:
        (report_dir / subdir).mkdir(parents=True, exist_ok=True)


def rel(path: Path, root: Path) -> str:
    try:
        return path.relative_to(root).as_posix()
    except ValueError:
        return path.as_posix()


def href(from_dir: Path, target: Path) -> str:
    return Path(os.path.relpath(target, from_dir)).as_posix()


def copy_file(src: Path | None, dst: Path, status_rows: list[dict[str, str]], label: str, base: Path) -> str:
    if src is None or not src.is_file():
        status_rows.append(
            {
                "item": label,
                "source": "" if src is None else rel(src, base),
                "target": rel(dst, base),
                "status": "missing_source",
                "notes": "",
            }
        )
        return "missing_source"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    status_rows.append(
        {
            "item": label,
            "source": rel(src, base),
            "target": rel(dst, base),
            "status": "copied",
            "notes": "",
        }
    )
    return "copied"


def first_match(root: Path, pattern: str) -> Path | None:
    matches = sorted(root.glob(pattern))
    return matches[0] if matches else None


def count_fasta(path: Path | None) -> int | None:
    if path is None or not path.is_file():
        return None
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def count_lines(path: Path | None) -> int | None:
    if path is None or not path.is_file():
        return None
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return sum(1 for line in handle if line.strip())


def format_int(value: int | None) -> str:
    if value is None:
        return "NA"
    return f"{value:,}"


def format_float(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "NA"
    return f"{value:.{digits}f}"


def read_fasta_stats(path: Path) -> dict[str, Any]:
    stats: dict[str, Any] = {"records": 0, "length": None, "gc_percent": None}
    if not path.is_file():
        return stats
    length = 0
    gc = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                stats["records"] += 1
                continue
            seq = line.upper()
            length += len(seq)
            gc += seq.count("G") + seq.count("C")
    stats["length"] = length
    stats["gc_percent"] = (gc / length * 100.0) if length else None
    return stats


def infer_gfa_topology(node_count: int, link_count: int, path_count: int, repeat_node_count: int, unique_edge_count: int) -> str:
    if node_count == 0:
        return "missing"
    if repeat_node_count or path_count:
        return f"repeat graph ({repeat_node_count} repeat nodes; {path_count} paths)"
    if node_count == 1 and link_count == 0:
        return "single contig"
    if unique_edge_count == node_count:
        return "circular graph"
    if unique_edge_count == max(node_count - 1, 0):
        return "linear graph"
    return "branched graph"


def read_gfa_stats(path: Path) -> dict[str, Any]:
    stats: dict[str, Any] = {
        "nodes": 0,
        "links": 0,
        "paths": 0,
        "repeat_nodes": 0,
        "topology": "missing",
    }
    if not path.is_file():
        return stats
    unique_edges: set[tuple[str, str]] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            record_type = parts[0]
            if record_type == "S" and len(parts) >= 3:
                stats["nodes"] += 1
                if any(tag == "NC:Z:repeat_node" for tag in parts[3:]):
                    stats["repeat_nodes"] += 1
            elif record_type == "L" and len(parts) >= 5:
                stats["links"] += 1
                unique_edges.add(tuple(sorted((parts[1], parts[3]))))
            elif record_type == "P":
                stats["paths"] += 1
    stats["topology"] = infer_gfa_topology(
        int(stats["nodes"]),
        int(stats["links"]),
        int(stats["paths"]),
        int(stats["repeat_nodes"]),
        len(unique_edges),
    )
    return stats


def read_sampling_summary(base: Path) -> dict[str, dict[str, str]]:
    summary = base / "reads" / "backup_info" / "downstream_read_sampling_summary.tsv"
    rows: dict[str, dict[str, str]] = {}
    if not summary.is_file():
        return rows
    with summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            genome = row.get("genome")
            if genome:
                rows[genome] = row
    return rows


def read_coverage_vector(path: Path, fallback_shape: Any = None):
    import numpy as np

    if path.is_file() and path.stat().st_size > 0:
        return np.loadtxt(path, dtype=float, usecols=1)
    if fallback_shape is None:
        return None
    return np.zeros(fallback_shape)


def coverage_y_max(cov_combine, y_axis_rule: str) -> tuple[float, str]:
    import numpy as np

    if cov_combine.size == 0:
        return 1.0, "empty"
    if y_axis_rule == "original":
        max_y = float(np.max(cov_combine))
        return max_y + 100.0, f"original=max(FL+partial)+100={max_y + 100.0:.6f}"
    if y_axis_rule == "twice_mean":
        mean_y = float(np.mean(cov_combine))
        y_max = 2.0 * mean_y
        if not math.isfinite(y_max) or y_max <= 0:
            max_y = float(np.max(cov_combine))
            y_max = max_y + 100.0
            return y_max, f"fallback=max(FL+partial)+100={y_max:.6f}"
        return y_max, f"twice_mean=2*mean(FL+partial)={y_max:.6f}"
    raise ValueError(f"Unsupported coverage y-axis rule: {y_axis_rule}")


def plot_coverage_like_original(
    fl_cov: Path,
    partial_cov: Path,
    variant_cov: Path,
    start: int,
    end: int,
    out_png: Path,
    out_pdf: Path,
    y_axis_rule: str = DEFAULT_FINAL_COVERAGE_Y_RULE,
) -> tuple[str, str]:
    """Redraw coverage with the original hifisr plot style and one y-axis rule change."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception as exc:  # pragma: no cover - environment-dependent fallback
        return "plot_failed", f"plot dependencies unavailable: {exc}"

    cov_1 = read_coverage_vector(fl_cov)
    if cov_1 is None:
        return "missing_source", "FL_cov.txt missing or empty"
    cov_2 = read_coverage_vector(partial_cov, cov_1.shape)
    cov_3 = read_coverage_vector(variant_cov, cov_1.shape)
    cov_combine = cov_1 + cov_2

    start = max(1, int(start))
    end = min(int(end), len(cov_combine))
    if end < start:
        return "plot_failed", "coverage interval is empty"

    sl = slice(start - 1, end)
    y_max, y_note = coverage_y_max(cov_combine[sl], y_axis_rule)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(12, 3), dpi=600)
    plt.plot(cov_combine[sl], color="#EAB13E", label="FL")
    plt.plot(cov_1[sl], color="#D1D1D1", label="partial")
    plt.fill_between(np.arange(start - 1, end), cov_combine[sl], 1, color="#EAB13E", alpha=1)
    plt.fill_between(np.arange(start - 1, end), cov_1[sl], 1, color="#D1D1D1", alpha=1)
    plt.plot(cov_3[sl], color="#5CAB38", label="variant", linewidth=1)
    plt.grid(True, alpha=0.5)
    ax = plt.gca()
    ax.set_xlim([start, end + 1])
    ax.set_ylim([0, y_max])
    fig.tight_layout()
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close()
    return "redrawn", f"same_plot_code_as_original_coverage; y_axis_rule={y_axis_rule}; {y_note}"


def count_fl_ids(ids_dir: Path) -> int:
    if not ids_dir.is_dir():
        return 0
    count = 0
    for path in ids_dir.glob("*_FL_ids.txt"):
        line_count = count_lines(path)
        count += int(line_count or 0)
    return count


def bubble_colour(mid_olp_1: float) -> str:
    if mid_olp_1 >= 1000:
        return "#FF00FF"
    if mid_olp_1 >= 300:
        return "#00FF00"
    if mid_olp_1 >= 200:
        return "#0016FF"
    if mid_olp_1 >= 100:
        return "#E8720C"
    if mid_olp_1 >= 50:
        return "#E6E600"
    return "#B0B0B0"


def plot_coverage_bubble_combined(
    fl_cov: Path,
    partial_cov: Path,
    variant_cov: Path,
    bubble_table: Path,
    ids_dir: Path,
    reference_fasta: Path | None,
    out_png: Path,
    out_pdf: Path,
    y_axis_rule: str = DEFAULT_FINAL_COVERAGE_Y_RULE,
) -> tuple[str, str]:
    """Draw coverage and type_2 repeat bubbles with a shared linear x-axis."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
    except Exception as exc:  # pragma: no cover - environment-dependent fallback
        return "plot_failed", f"plot dependencies unavailable: {exc}"

    cov_1 = read_coverage_vector(fl_cov)
    if cov_1 is None:
        return "missing_source", "FL_cov.txt missing or empty"
    cov_2 = read_coverage_vector(partial_cov, cov_1.shape)
    cov_3 = read_coverage_vector(variant_cov, cov_1.shape)
    cov_combine = cov_1 + cov_2

    ref_len = None
    if reference_fasta is not None:
        ref_len = read_fasta_stats(reference_fasta).get("length")
    if not ref_len:
        ref_len = len(cov_combine)
    ref_len = int(ref_len)
    end = min(ref_len, len(cov_combine))
    if end < 1:
        return "plot_failed", "coverage interval is empty"

    y_max, y_note = coverage_y_max(cov_combine[:end], y_axis_rule)
    x = np.arange(1, end + 1)

    bubble_note = "bubble=missing"
    bubble_df = None
    if bubble_table.is_file():
        try:
            bubble_df = pd.read_excel(bubble_table, sheet_name="Sheet1")
            coords = bubble_df["(se1, ss2)"].astype(str).str.extract(r"\(?\s*(\d+)\s*,\s*(\d+)\s*\)?")
            bubble_df["se1"] = pd.to_numeric(coords[0], errors="coerce")
            bubble_df["ss2"] = pd.to_numeric(coords[1], errors="coerce")
            bubble_df["mid_olp_1"] = pd.to_numeric(bubble_df["mid_olp_1"], errors="coerce").fillna(0)
            bubble_df["subgroup_count"] = pd.to_numeric(bubble_df["subgroup_count"], errors="coerce").fillna(0)
            bubble_df = bubble_df.dropna(subset=["se1", "ss2"]).copy()
            fl_count = count_fl_ids(ids_dir)
            if fl_count > 0:
                bubble_df["subgroup_count_norm"] = bubble_df["subgroup_count"] / fl_count * 10000.0
            else:
                bubble_df["subgroup_count_norm"] = bubble_df["subgroup_count"]
            bubble_df["color"] = bubble_df["mid_olp_1"].apply(bubble_colour)
            bubble_note = f"bubble_rows={len(bubble_df)}; FL_id_count={fl_count}"
        except Exception as exc:
            bubble_df = None
            bubble_note = f"bubble_parse_failed={exc}"

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax_cov, ax_bubble) = plt.subplots(
        2,
        1,
        figsize=(9.0, 7.0),
        dpi=600,
        sharex=True,
        gridspec_kw={"height_ratios": [1.0, 3.2], "hspace": 0.06},
    )

    ax_cov.plot(x, cov_combine[:end], color="#EAB13E", label="FL")
    ax_cov.plot(x, cov_1[:end], color="#D1D1D1", label="partial")
    ax_cov.fill_between(x, cov_combine[:end], 1, color="#EAB13E", alpha=1)
    ax_cov.fill_between(x, cov_1[:end], 1, color="#D1D1D1", alpha=1)
    ax_cov.plot(x, cov_3[:end], color="#5CAB38", label="variant", linewidth=1)
    ax_cov.grid(True, alpha=0.5)
    ax_cov.set_ylim([0, y_max])
    ax_cov.set_ylabel("Coverage")
    ax_cov.tick_params(axis="x", labelbottom=False)

    if bubble_df is not None and len(bubble_df) > 0:
        ax_bubble.scatter(
            bubble_df["se1"],
            bubble_df["ss2"],
            s=bubble_df["subgroup_count_norm"] * 10,
            c=bubble_df["color"],
            alpha=0.5,
            linewidths=0,
        )
    ax_bubble.grid(True, alpha=0.5)
    ax_bubble.set_xlim([1, ref_len])
    ax_bubble.set_ylim([1, ref_len])
    ax_bubble.set_aspect("equal", adjustable="box")
    ax_bubble.set_xlabel("Linear coordinate (bp)")
    ax_bubble.set_ylabel("Bubble coordinate (bp)")

    fig.tight_layout()
    fig.canvas.draw()
    bubble_pos = ax_bubble.get_position()
    cov_pos = ax_cov.get_position()
    ax_cov.set_position([bubble_pos.x0, cov_pos.y0, bubble_pos.width, cov_pos.height])
    plt.savefig(out_pdf)
    plt.savefig(out_png)
    plt.close(fig)
    return "redrawn", f"shared_x_axis=1..{ref_len}; bubble_aspect=equal; coverage_y_axis={y_axis_rule}; {y_note}; {bubble_note}"


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=False)
        handle.write("\n")


def html_page(title: str, body: str) -> str:
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            '<meta name="viewport" content="width=device-width, initial-scale=1">',
            f"<title>{html.escape(title)}</title>",
            "<style>",
            "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;margin:0;color:#222;background:#fafafa;}",
            "main{max-width:1120px;margin:0 auto;padding:32px 24px 48px;}",
            "h1{font-size:30px;margin:0 0 8px;} h2{font-size:20px;margin:30px 0 10px;} p{line-height:1.55;}",
            ".lead{color:#555;margin-top:0;}.grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:14px;}",
            ".meta{font-size:13px;color:#666;margin:36px 0 0;white-space:nowrap;overflow-x:auto;}",
            ".card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:14px 16px;}",
            ".report-hero{display:grid;grid-template-columns:minmax(0,.82fr) minmax(0,1.18fr);gap:14px;margin:18px 0 4px;align-items:start;}",
            ".hero-left{display:grid;grid-template-rows:auto auto;gap:14px;align-content:start;}",
            ".figure-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px;align-items:start;margin-top:16px;}",
            ".figure-card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:10px;}",
            ".figure-card h3{display:flex;justify-content:space-between;gap:12px;font-size:14px;margin:0 0 8px;white-space:nowrap;}",
            ".figure-links{font-weight:500;color:#666;}",
            ".figure-card img{width:100%;height:340px;object-fit:contain;border:1px solid #e3e3e3;border-radius:6px;background:#fff;}",
            ".figure-card.graph img{height:368px;}.figure-card.combined{align-self:start;}.figure-card.combined img{height:auto;display:block;}",
            ".final-strip{display:block;margin:0;}",
            ".final-card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:10px 12px;}",
            ".final-card h3{font-size:14px;margin:0 0 6px;white-space:nowrap;}.final-card p{font-size:13px;color:#555;line-height:1.35;margin:3px 0;}",
            ".final-card strong{color:#222;font-weight:650;}",
            ".evidence-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:12px;margin-top:24px;}",
            ".evidence-card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:11px 12px;}",
            ".evidence-card h3{font-size:14px;margin:0 0 6px;white-space:nowrap;}.evidence-card p{font-size:13px;color:#555;margin:0 0 8px;line-height:1.35;}",
            ".evidence-links{font-size:13px;line-height:1.55;}",
            ".summary-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px;}",
            ".summary-card{background:#fff;border:1px solid #ddd;border-radius:8px;padding:10px;}.summary-card h3{font-size:14px;margin:0 0 8px;white-space:nowrap;}",
            ".summary-card table{border:0;}.summary-card th,.summary-card td{padding:7px 8px;}",
            "table{width:100%;border-collapse:collapse;background:#fff;border:1px solid #ddd;} th,td{padding:9px 10px;border-bottom:1px solid #eee;text-align:left;vertical-align:top;white-space:nowrap;} th{background:#f2f2f2;}",
            "code{background:#f0f0f0;padding:1px 4px;border-radius:4px;} a{color:#1257a8;text-decoration:none;} a:hover{text-decoration:underline;}",
            "img{max-width:100%;border:1px solid #ddd;background:#fff;border-radius:6px;}",
            ".tag{display:inline-block;font-size:12px;background:#e9f2ff;color:#164b83;border-radius:999px;padding:2px 8px;margin-left:6px;}",
            "@media(max-width:980px){.evidence-grid{grid-template-columns:repeat(2,minmax(0,1fr));}}",
            "@media(max-width:820px){.report-hero,.hero-left,.figure-grid,.summary-grid,.evidence-grid{grid-template-columns:1fr;}.figure-card img,.figure-card.graph img,.figure-card.combined img{height:240px;}}",
            "</style>",
            "</head>",
            "<body><main>",
            body,
            "</main></body></html>",
            "",
        ]
    )


def link(from_dir: Path, target: Path, label: str) -> str:
    return f'<a href="{html.escape(href(from_dir, target))}">{html.escape(label)}</a>'


def figure_card(from_dir: Path, image: Path, links: list[tuple[Path, str]], title: str, extra_class: str = "") -> str:
    if image.is_file():
        preview = f'<img src="{html.escape(href(from_dir, image))}" alt="{html.escape(title)} preview">'
    else:
        preview = f'<p class="meta">Preview missing: <code>{html.escape(href(from_dir, image))}</code></p>'
    link_text = " · ".join(link(from_dir, target, label) for target, label in links)
    class_attr = "figure-card" + (f" {html.escape(extra_class)}" if extra_class else "")
    return (
        f'<section class="{class_attr}">'
        f"<h3><span>{html.escape(title)}</span>"
        f'<span class="figure-links">{link_text}</span></h3>'
        f"{preview}"
        "</section>"
    )


def evidence_card(from_dir: Path, title: str, description: str, folder: Path, key_files: list[tuple[Path, str]]) -> str:
    links = [link(from_dir, folder, "folder")]
    links.extend(link(from_dir, path, label) for path, label in key_files)
    return (
        '<section class="evidence-card">'
        f"<h3>{html.escape(title)}</h3>"
        f"<p>{html.escape(description)}</p>"
        f'<div class="evidence-links">{" · ".join(links)}</div>'
        "</section>"
    )


def table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(header)}</th>" for header in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in row) + "</tr>")
    return f"<table><thead><tr>{head}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def status_for(path: Path) -> str:
    return "present" if path.exists() else "missing"


def row_link(from_dir: Path, target: Path, label: str, description: str) -> list[str]:
    return [
        link(from_dir, target, label),
        html.escape(description),
        f"<code>{html.escape(href(from_dir, target))}</code>",
        html.escape(status_for(target)),
    ]


def genome_label(genome: str) -> str:
    labels = {"mito": "Mitochondrial genome", "plastid": "Plastid genome"}
    return labels.get(genome, genome)


def genome_reference(base: Path, genome: str, cfg: dict[str, str]) -> Path | None:
    return first_match(base / "draft_assembly" / genome / "backup_info", str(cfg["reference_pattern"]))


def genome_linear_fasta(base: Path, genome: str, cfg: dict[str, str]) -> Path:
    if cfg["run"] == "run_3":
        return base / "draft_assembly" / genome / f"{genome}_flye_polish_aligned_cor.fasta"
    return base / "draft_assembly" / genome / f"{genome}_flye_polish_aligned.fasta"


def genome_status(base: Path, genome: str, cfg: dict[str, str]) -> tuple[str, str]:
    prefix = cfg["verified_prefix"]
    report_dir = base / "reports" / genome
    support_dir = base / "draft_assembly" / genome / "verified_gfa_read_support"
    run_dir = base / genome / cfg["run"]
    required = [
        report_dir / "index.html",
        report_dir / "final" / "final.gfa",
        report_dir / "final" / "final.fasta",
        report_dir / "final" / "final_graph.pdf",
        report_dir / "final" / "final_graph.svg",
        report_dir / "final" / "final_coverage.pdf",
        report_dir / "final" / "final_coverage_bubble.pdf",
        support_dir / f"{prefix}.gfa",
        support_dir / f"{prefix}_linear_to_node_coordinate.tsv",
        support_dir / "variant_by_nodes" / "variants_anno_combined_depth_frq_filter.by_verified_node.xlsx",
        run_dir / ".snakemake.done",
    ]
    missing = [rel(path, base) for path in required if not path.exists()]
    if missing:
        return "incomplete", "; ".join(missing[:4])
    return "success", "all compact report and upstream key outputs are present"


def genome_project_rows(base: Path, reports_dir: Path, genome: str, cfg: dict[str, str]) -> list[list[str]]:
    run_dir = base / genome / cfg["run"]
    support_dir = base / "draft_assembly" / genome / "verified_gfa_read_support"
    rows = [
        row_link(reports_dir, reports_dir / genome / "index.html", f"{genome} report entry", "Compact final report entry page."),
        row_link(reports_dir, reports_dir / genome / "final", "final deliverables", "Small copied final files: graph, FASTA, graph PDF, coverage PDF/PNG, bubble PDF."),
        row_link(reports_dir, base / "reads", "reads workspace", "Extracted organelle reads and read summary plots; linked only, not copied."),
        row_link(reports_dir, base / "reads" / "backup_info" / "downstream_read_sampling_summary.tsv", "read sampling summary", "Read extraction and downstream sampling counts."),
        row_link(reports_dir, run_dir, f"{genome} {cfg['run']}", "Final variant/coverage run used for this report."),
        row_link(reports_dir, base / "draft_assembly" / genome, "draft assembly workspace", "Draft graphs, manual checked draft, polish files, and manual correction inputs."),
        row_link(reports_dir, support_dir, "verified GFA support", "Full verified graph, PDFs, node coordinate maps, and repeat-support evidence."),
        row_link(reports_dir, base / "logs" / "snakemake", "Snakemake logs", "Workflow logs and report generation log."),
    ]
    return rows


def genome_input_rows(base: Path, reports_dir: Path, genome: str, cfg: dict[str, str]) -> list[list[str]]:
    reference = genome_reference(base, genome, cfg)
    run_dir = base / genome / cfg["run"]
    rows = [
        row_link(reports_dir, base / "reads" / f"{genome}.fastq.gz", "organellar reads", "Reads extracted for this genome."),
        row_link(reports_dir, base / "draft_assembly" / genome / f"{genome}_checked_draft.gfa", "checked draft GFA", "Manual checked graph used as verified GFA input."),
        row_link(reports_dir, genome_linear_fasta(base, genome, cfg), "verified linear FASTA input", "Linear FASTA used by verified GFA construction and coordinate projection."),
        row_link(reports_dir, run_dir / "variants_anno_combined_depth_frq_filter.xlsx", "filtered variant table", "SNV/InDel review table used for final reporting."),
        row_link(reports_dir, run_dir / "FL_cov.txt", "FL coverage", "Coverage vector used to redraw final coverage."),
        row_link(reports_dir, run_dir / "partial_cov.txt", "partial coverage", "Coverage vector used to redraw final coverage."),
    ]
    if reference is not None:
        rows.append(row_link(reports_dir, reference, "rotated reference", "Reference FASTA used for graph image export and validation."))
    else:
        rows.append(["rotated reference", "Reference FASTA used for graph image export and validation.", "<code>not found</code>", "missing"])
    return rows


def genome_command_rows(base: Path, reports_dir: Path, genome: str, cfg: dict[str, str]) -> list[list[str]]:
    support_dir = base / "draft_assembly" / genome / "verified_gfa_read_support"
    run_dir = base / genome / cfg["run"]
    return [
        [
            "<code>snakemake --snakefile Snakefile --configfile workflow/config/w3_5_2_gfa.yaml draft_for_manual_edit</code>",
            "Generate draft GFA/PDF files for manual graph inspection.",
            "manual breakpoint",
        ],
        [
            "<code>snakemake --snakefile Snakefile --configfile workflow/config/w3_5_2_gfa.yaml verified_gfa_read_support</code>",
            "Build verified GFA outputs and node-coordinate projection tables.",
            html.escape(status_for(support_dir / f"{cfg['verified_prefix']}.gfa")),
        ],
        [
            "<code>snakemake --snakefile Snakefile --configfile workflow/config/w3_5_2_gfa.yaml compact_report_bundle</code>",
            "Generate the compact report and HTML browser.",
            html.escape(status_for(base / "reports" / genome / "index.html")),
        ],
        [
            "<code>snakemake --snakefile Snakefile --configfile workflow/config/w3_5_2_gfa.yaml final</code>",
            "Run the full final target, including compact reports.",
            html.escape(status_for(run_dir / ".snakemake.done")),
        ],
    ]


def genome_manual_breakpoint_rows(base: Path, reports_dir: Path, genome: str, cfg: dict[str, str]) -> list[list[str]]:
    draft_dir = base / "draft_assembly" / genome
    rows = [
        row_link(reports_dir, draft_dir / f"{genome}_checked_draft.gfa", "checked draft graph", "Human-reviewed graph expected before verified GFA support is built."),
        row_link(reports_dir, draft_dir / f"{genome}_checked_draft.pdf", "checked draft graph PDF", "Human-readable graph snapshot for manual review."),
    ]
    if cfg["run"] == "run_3":
        rows.append(row_link(reports_dir, draft_dir / "pos_ref_alt.txt", "manual correction table", "Manual sequence corrections applied before the final variant verification run."))
    else:
        rows.append(["manual correction table", "Not active for this genome in the current config.", "<code>not used</code>", "not_applicable"])
    return rows


def genome_output_rows(base: Path, reports_dir: Path, genome: str, cfg: dict[str, str]) -> list[list[str]]:
    prefix = cfg["verified_prefix"]
    report_dir = base / "reports" / genome
    support_dir = base / "draft_assembly" / genome / "verified_gfa_read_support"
    return [
        row_link(reports_dir, report_dir / "final" / "final.gfa", "final.gfa", "Copied verified final graph."),
        row_link(reports_dir, report_dir / "final" / "final.fasta", "final.fasta", "Copied verified final sequence."),
        row_link(reports_dir, report_dir / "final" / "final_graph.pdf", "final_graph.pdf", "Copied verified graph PDF."),
        row_link(reports_dir, report_dir / "final" / "final_graph.svg", "final_graph.svg", "Copied verified graph SVG for browser preview."),
        row_link(reports_dir, report_dir / "final" / "final_coverage.pdf", "final_coverage.pdf", "Redrawn coverage using the original hifisr coverage plot code with adjusted y-axis."),
        row_link(reports_dir, report_dir / "final" / "final_bubble.pdf", "final_bubble.pdf", "Copied final bubble plot when available."),
        row_link(reports_dir, support_dir / f"{prefix}_linear_to_node_coordinate.tsv", "linear-to-node coordinate map", "Full source coordinate projection table in the original project tree."),
        row_link(reports_dir, support_dir / "variant_by_nodes" / "variants_anno_combined_depth_frq_filter.by_verified_node.xlsx", "filtered variants by node", "Node-split filtered SNV/InDel table."),
    ]


def genome_status_rows(base: Path, reports_dir: Path) -> list[list[str]]:
    rows = []
    for genome, cfg in GENOMES.items():
        status, note = genome_status(base, genome, cfg)
        rows.append(
            [
                link(reports_dir, reports_dir / genome / "index.html", genome_label(genome)),
                html.escape(status),
                html.escape(note),
            ]
        )
    return rows


def write_project_index(base: Path, generated_at: datetime, elapsed_seconds: float) -> None:
    reports_dir = base / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    meta_line = (
        '<p class="meta">'
        f"hifisr {html.escape(HIFISR_VERSION)} · "
        f"generated {html.escape(generated_at.strftime('%Y-%m-%d %H:%M:%S %Z'))} · "
        f"report build {elapsed_seconds:.1f}s · "
        f"{html.escape(platform.platform())} · "
        f"Python {html.escape(platform.python_version())}"
        "</p>"
    )
    genome_rows: list[list[str]] = []
    for genome, cfg in GENOMES.items():
        status, note = genome_status(base, genome, cfg)
        report_dir = reports_dir / genome
        final_dir = report_dir / "final"
        genome_rows.append(
            [
                link(reports_dir, report_dir / "index.html", genome_label(genome)),
                f"{link(reports_dir, final_dir / 'final.gfa', 'GFA')} · "
                f"{link(reports_dir, final_dir / 'final.fasta', 'FASTA')}",
                html.escape(status),
                html.escape(note),
            ]
        )
    workflow_rows = [
        ["Input", f"{link(reports_dir, base / 'reads', 'reads')} · checked draft GFA · rotated reference"],
        ["Manual", "checked draft graph; mito pos_ref_alt corrections"],
        ["Run", "<code>snakemake --snakefile Snakefile --configfile workflow/config/w3_5_2_gfa.yaml final</code>"],
        ["Report", "<code>analysis_scripts/generate_compact_report_bundle.py</code> via Snakemake"],
    ]
    sections = [
        "<h1>W3-5-2_gfa reports</h1>",
        '<p class="lead">Open a genome report first. Large evidence files stay linked in the project tree.</p>',
        "<h2>At A Glance</h2>",
        table(["Genome", "Final files", "Status", "Note"], genome_rows),
        "<h2>Workflow</h2>",
        table(["Stage", "Summary"], workflow_rows),
        meta_line,
    ]
    body = "\n".join(sections)
    (reports_dir / "index.html").write_text(html_page("W3-5-2_gfa project browser", body), encoding="utf-8")
    write_json(
        reports_dir / "project_manifest.json",
        {
            "schema": "hifisr-project-report-index",
            "version": 1,
            "base": str(base),
            "entrypoint": "index.html",
            "hifisr_version": HIFISR_VERSION,
            "generated_at": generated_at.isoformat(),
            "report_build_seconds": round(elapsed_seconds, 3),
            "platform": platform.platform(),
            "python": platform.python_version(),
            "compact_reports": {genome: f"{genome}/index.html" for genome in GENOMES},
        },
    )


def source_row(report_dir: Path, path: Path, label: str, description: str) -> list[str]:
    status = "present" if path.exists() else "missing"
    return [
        link(report_dir, path, label),
        html.escape(description),
        f"<code>{html.escape(href(report_dir, path))}</code>",
        html.escape(status),
    ]


def write_genome_index(
    base: Path,
    genome: str,
    cfg: dict[str, str],
    report_dir: Path,
    coverage_status: str,
    coverage_note: str,
) -> None:
    prefix = cfg["verified_prefix"]
    run_dir = base / genome / cfg["run"]
    support_dir = base / "draft_assembly" / genome / "verified_gfa_read_support"
    final_dir = report_dir / "final"
    metadata_dir = report_dir / "metadata"

    gfa_stats = read_gfa_stats(final_dir / "final.gfa")
    fasta_stats = read_fasta_stats(final_dir / "final.fasta")
    final_stats_html = (
        '<section class="final-card">'
        "<h3>Final files</h3>"
        f"<p>{link(report_dir, final_dir / 'final.gfa', 'final.gfa')} · "
        f"<strong>{format_int(gfa_stats['nodes'])}</strong> nodes · "
        f"<strong>{format_int(gfa_stats['links'])}</strong> links</p>"
        f"<p><strong>Topology</strong> · {html.escape(str(gfa_stats['topology']))}</p>"
        f"<p>{link(report_dir, final_dir / 'final.fasta', 'final.fasta')} · "
        f"<strong>{format_int(fasta_stats['length'])}</strong> bp · "
        f"<strong>{format_float(fasta_stats['gc_percent'])}%</strong> GC</p>"
        "</section>"
    )
    hero_html = (
        '<div class="report-hero">'
        '<div class="hero-left">'
        + figure_card(
            report_dir,
            final_dir / "final_graph.svg",
            [(final_dir / "final_graph.pdf", "PDF"), (final_dir / "final_graph.svg", "SVG")],
            "final_graph",
            "graph",
        )
        + final_stats_html
        + "</div>"
        + figure_card(
            report_dir,
            final_dir / "final_coverage_bubble.png",
            [(final_dir / "final_coverage_bubble.pdf", "PDF"), (final_dir / "final_coverage_bubble.png", "PNG")],
            "coverage_bubble",
            "combined",
        )
        + "</div>"
    )
    draft_dir = base / "draft_assembly" / genome
    manual_graph_files = [
        (draft_dir / f"{genome}_checked_draft.gfa", "checked GFA"),
        (draft_dir / f"{genome}_checked_draft.pdf", "checked PDF"),
    ]
    correction_table = draft_dir / "pos_ref_alt.txt"
    if correction_table.exists():
        manual_graph_files.append((correction_table, "corrections"))
    evidence_cards = [
        evidence_card(
            report_dir,
            "Reads",
            "Extracted organelle reads and sampling summary.",
            base / "reads",
            [
                (base / "reads" / f"{genome}.fastq.gz", f"{genome}.fastq.gz"),
                (base / "reads" / "backup_info" / "downstream_read_sampling_summary.tsv", "sampling TSV"),
            ],
        ),
        evidence_card(
            report_dir,
            "Variant Run",
            f"Final {cfg['run']} variant and coverage evidence.",
            run_dir,
            [
                (run_dir / "variants_anno_combined_depth_frq_filter.xlsx", "filtered variants"),
                (run_dir / "FL_cov.txt", "FL_cov"),
                (run_dir / "partial_cov.txt", "partial_cov"),
            ],
        ),
        evidence_card(
            report_dir,
            "Verified GFA",
            "Graph, sequence, coordinate map, and node-split variants.",
            support_dir,
            [
                (support_dir / f"{prefix}.gfa", "GFA"),
                (support_dir / f"{prefix}.fasta", "FASTA"),
                (support_dir / f"{prefix}_linear_to_node_coordinate.tsv", "node map"),
                (support_dir / "variant_by_nodes" / "variants_anno_combined_depth_frq_filter.by_verified_node.xlsx", "variants by node"),
            ],
        ),
        evidence_card(
            report_dir,
            "Manual Graph",
            "Human-checked draft graph and optional correction table.",
            draft_dir,
            manual_graph_files,
        ),
        evidence_card(
            report_dir,
            "Logs",
            "Workflow and compact report run logs.",
            base / "logs" / "snakemake",
            [
                (base / "logs" / "snakemake" / "compact_report_bundle.log", "report log"),
            ],
        ),
        evidence_card(
            report_dir,
            "Workflow",
            "Inputs, manual checkpoint, final outputs, and report manifests.",
            metadata_dir,
            [
                (genome_linear_fasta(base, genome, cfg), "linear FASTA"),
                (draft_dir / f"{genome}_checked_draft.gfa", "checked draft"),
                (metadata_dir / f"{prefix}_linear_to_node_coordinate.tsv", "node map"),
                (metadata_dir / "report_manifest.json", "report manifest"),
            ],
        ),
    ]
    body = "\n".join(
        [
            f"<h1>{html.escape(genome_label(genome))}</h1>",
            f'<p class="lead"><span class="tag">{html.escape(coverage_status)}</span> Compact final report. Large evidence files are linked only.</p>',
            hero_html,
            '<div class="evidence-grid">' + "".join(evidence_cards) + "</div>",
            f'<p><a href="{html.escape(href(report_dir, base / "reports" / "index.html"))}">Back to project browser</a></p>',
        ]
    )
    (report_dir / "index.html").write_text(html_page(f"{genome} compact report", body), encoding="utf-8")


def generate_genome_report(
    base: Path,
    genome: str,
    sampling: dict[str, dict[str, str]],
    final_coverage_y_rule: str,
) -> None:
    cfg = GENOMES[genome]
    run_dir = base / genome / cfg["run"]
    support_dir = base / "draft_assembly" / genome / "verified_gfa_read_support"
    draft_backup = base / "draft_assembly" / genome / "backup_info"
    report_dir = base / "reports" / genome

    if report_dir.exists():
        shutil.rmtree(report_dir)
    ensure_dirs(report_dir)
    status_rows: list[dict[str, str]] = []

    prefix = cfg["verified_prefix"]
    verified_gfa = support_dir / f"{prefix}.gfa"
    verified_fasta = support_dir / f"{prefix}.fasta"
    verified_pdf = support_dir / f"{prefix}.pdf"
    verified_svg = support_dir / f"{prefix}.svg"
    reference_fasta = first_match(draft_backup, str(cfg["reference_pattern"]))
    bubble = first_match(run_dir, "bubble*.pdf")
    bubble_png = bubble.with_suffix(".png") if bubble is not None else first_match(run_dir, "bubble*.png")
    bubble_table = run_dir / "combined_excel" / "type_2_subtype_rep_NA_summary_anno.xlsx"
    ids_dir = run_dir / "backup_info" / "IDs"
    original_coverage_pdf = first_match(run_dir, "coverage_*.pdf")
    fl_cov = run_dir / "FL_cov.txt"
    partial_cov = run_dir / "partial_cov.txt"
    variant_cov = run_dir / "variant_cov.txt"

    # Final deliverables: copied intentionally small, user-facing files.
    copy_file(verified_gfa, report_dir / "final" / "final.gfa", status_rows, "final.gfa", base)
    copy_file(verified_fasta, report_dir / "final" / "final.fasta", status_rows, "final.fasta", base)
    copy_file(verified_pdf, report_dir / "final" / "final_graph.pdf", status_rows, "final_graph.pdf", base)
    graph_svg_status = copy_file(verified_svg, report_dir / "final" / "final_graph.svg", status_rows, "final_graph.svg", base)
    copy_file(bubble, report_dir / "final" / "final_bubble.pdf", status_rows, "final_bubble.pdf", base)
    bubble_png_status = copy_file(bubble_png, report_dir / "final" / "final_bubble.png", status_rows, "final_bubble.png", base)

    # Coverage is generated directly into final/; original coverage remains linked from the HTML.
    coverage_status, coverage_note = plot_coverage_like_original(
        fl_cov,
        partial_cov,
        variant_cov,
        1,
        10**18,
        report_dir / "final" / "final_coverage.png",
        report_dir / "final" / "final_coverage.pdf",
        y_axis_rule=final_coverage_y_rule,
    )
    status_rows.append(
        {
            "item": "final_coverage_redraw",
            "source": f"{rel(fl_cov, base)};{rel(partial_cov, base)};{rel(variant_cov, base)}",
            "target": rel(report_dir / "final" / "final_coverage.png", base),
            "status": coverage_status,
            "notes": coverage_note,
        }
    )
    combined_status, combined_note = plot_coverage_bubble_combined(
        fl_cov,
        partial_cov,
        variant_cov,
        bubble_table,
        ids_dir,
        reference_fasta,
        report_dir / "final" / "final_coverage_bubble.png",
        report_dir / "final" / "final_coverage_bubble.pdf",
        y_axis_rule=final_coverage_y_rule,
    )
    status_rows.append(
        {
            "item": "final_coverage_bubble_redraw",
            "source": f"{rel(fl_cov, base)};{rel(partial_cov, base)};{rel(variant_cov, base)};{rel(bubble_table, base)}",
            "target": rel(report_dir / "final" / "final_coverage_bubble.png", base),
            "status": combined_status,
            "notes": combined_note,
        }
    )
    final_gfa_stats = read_gfa_stats(report_dir / "final" / "final.gfa")
    final_fasta_stats = read_fasta_stats(report_dir / "final" / "final.fasta")
    write_tsv(
        report_dir / "metadata" / "final_summary.tsv",
        ["file", "metric", "value"],
        [
            {"file": "final.gfa", "metric": "nodes", "value": final_gfa_stats["nodes"]},
            {"file": "final.gfa", "metric": "links", "value": final_gfa_stats["links"]},
            {"file": "final.gfa", "metric": "paths", "value": final_gfa_stats["paths"]},
            {"file": "final.gfa", "metric": "repeat_nodes", "value": final_gfa_stats["repeat_nodes"]},
            {"file": "final.gfa", "metric": "topology", "value": final_gfa_stats["topology"]},
            {"file": "final.fasta", "metric": "records", "value": final_fasta_stats["records"]},
            {"file": "final.fasta", "metric": "length_bp", "value": final_fasta_stats["length"]},
            {"file": "final.fasta", "metric": "gc_percent", "value": format_float(final_fasta_stats["gc_percent"])},
        ],
    )

    # Compact metadata: small, high-value traceability files only.
    metadata_sources = [
        (support_dir / "verified_node_source_map.tsv", "verified_node_source_map.tsv"),
        (support_dir / f"{prefix}_linear_to_node_coordinate.tsv", f"{prefix}_linear_to_node_coordinate.tsv"),
        (support_dir / "coordinate_consistency.tsv", "coordinate_consistency.tsv"),
        (support_dir / f"{prefix}_by_nodes.fasta", f"{prefix}_by_nodes.fasta"),
        (support_dir / "variant_by_nodes" / "variants_anno_combined_depth_frq_filter.by_verified_node.xlsx", "variants_anno_combined_depth_frq_filter.by_verified_node.xlsx"),
        (run_dir / "variants_anno_combined_depth_frq_filter.xlsx", "variants_anno_combined_depth_frq_filter.xlsx"),
        (reference_fasta, "reference.fasta"),
    ]
    for src, name in metadata_sources:
        copy_file(src, report_dir / "metadata" / name, status_rows, f"metadata_{name}", base)

    final_rows = [
        {"final_item": "final_gfa", "source_file": rel(verified_gfa, base), "target_file": "final/final.gfa", "status": "copied" if verified_gfa.is_file() else "missing_source", "notes": "verified final graph"},
        {"final_item": "final_fasta", "source_file": rel(verified_fasta, base), "target_file": "final/final.fasta", "status": "copied" if verified_fasta.is_file() else "missing_source", "notes": "verified final sequence"},
        {"final_item": "final_graph_pdf", "source_file": rel(verified_pdf, base), "target_file": "final/final_graph.pdf", "status": "copied" if verified_pdf.is_file() else "missing_source", "notes": "verified graph image"},
        {"final_item": "final_graph_svg", "source_file": rel(verified_svg, base), "target_file": "final/final_graph.svg", "status": graph_svg_status, "notes": "verified graph SVG exported with GFA_Editor"},
        {"final_item": "final_coverage", "source_file": f"{rel(fl_cov, base)};{rel(partial_cov, base)};{rel(variant_cov, base)}", "target_file": "final/final_coverage.pdf", "status": coverage_status, "notes": coverage_note},
        {"final_item": "final_coverage_bubble", "source_file": f"{rel(fl_cov, base)};{rel(partial_cov, base)};{rel(variant_cov, base)};{rel(bubble_table, base)}", "target_file": "final/final_coverage_bubble.pdf;final/final_coverage_bubble.png", "status": combined_status, "notes": combined_note},
        {"final_item": "final_bubble", "source_file": rel(bubble, base) if bubble else "", "target_file": "final/final_bubble.pdf", "status": "copied" if bubble and bubble.is_file() else "missing_source", "notes": "bubble plot"},
        {"final_item": "final_bubble_png", "source_file": rel(bubble_png, base) if bubble_png else "", "target_file": "final/final_bubble.png", "status": bubble_png_status, "notes": "bubble plot PNG saved by matplotlib"},
    ]
    write_tsv(report_dir / "metadata" / "final_manifest.tsv", ["final_item", "source_file", "target_file", "status", "notes"], final_rows)
    write_tsv(
        report_dir / "metadata" / "figure_manifest.tsv",
        ["figure_id", "figure_type", "source_file", "target_file", "status", "redraw_rule", "notes"],
        [
            {
                "figure_id": "final_graph",
                "figure_type": "graph",
                "source_file": rel(verified_pdf, base) + ";" + rel(verified_svg, base),
                "target_file": "final/final_graph.pdf;final/final_graph.svg",
                "status": graph_svg_status,
                "redraw_rule": "copied_pdf; copied_svg_exported_with_gfa_editor",
                "notes": "",
            },
            {
                "figure_id": "original_coverage",
                "figure_type": "coverage",
                "source_file": rel(original_coverage_pdf, base) if original_coverage_pdf else "",
                "target_file": "",
                "status": "linked" if original_coverage_pdf and original_coverage_pdf.is_file() else "missing_source",
                "redraw_rule": "not_copied",
                "notes": "original coverage remains in the run directory and is linked from HTML",
            },
            {
                "figure_id": "final_coverage",
                "figure_type": "coverage",
                "source_file": f"{rel(fl_cov, base)};{rel(partial_cov, base)};{rel(variant_cov, base)}",
                "target_file": "final/final_coverage.pdf",
                "status": coverage_status,
                "redraw_rule": f"same as hifisr_functions.reports.plot_coverage; y_axis_rule={final_coverage_y_rule}",
                "notes": coverage_note,
            },
            {
                "figure_id": "final_coverage_bubble",
                "figure_type": "coverage_bubble",
                "source_file": f"{rel(fl_cov, base)};{rel(partial_cov, base)};{rel(variant_cov, base)};{rel(bubble_table, base)}",
                "target_file": "final/final_coverage_bubble.pdf;final/final_coverage_bubble.png",
                "status": combined_status,
                "redraw_rule": f"coverage and hifisr_functions.reports.plot_bubble_type_2_rep_raw logic redrawn on shared x-axis; y_axis_rule={final_coverage_y_rule}",
                "notes": combined_note,
            },
            {
                "figure_id": "final_bubble",
                "figure_type": "bubble",
                "source_file": rel(bubble, base) if bubble else "",
                "target_file": "final/final_bubble.pdf;final/final_bubble.png",
                "status": bubble_png_status,
                "redraw_rule": "copied_pdf; copied_png_saved_by_matplotlib",
                "notes": "",
            },
        ],
    )
    write_tsv(report_dir / "metadata" / "report_collection_status.tsv", ["item", "source", "target", "status", "notes"], status_rows)

    sample_row = sampling.get(genome, {})
    report_manifest = {
        "schema": "hifisr-compact-report-bundle",
        "version": 2,
        "genome": genome,
        "run": cfg["run"],
        "status": "generated",
        "coverage_redraw_status": coverage_status,
        "coverage_redraw_note": coverage_note,
        "coverage_bubble_redraw_status": combined_status,
        "coverage_bubble_redraw_note": combined_note,
        "final_coverage_y_rule": final_coverage_y_rule,
        "sections": {subdir: f"{subdir}/" for subdir in REPORT_SUBDIRS},
        "source_directories": {
            "run": rel(run_dir, base),
            "verified_gfa_read_support": rel(support_dir, base),
            "draft_assembly": rel(base / "draft_assembly" / genome, base),
        },
        "read_sampling": sample_row,
        "record_counts": {
            "reference_fasta": count_fasta(reference_fasta),
            "verified_fasta": count_fasta(verified_fasta),
            "filtered_variants_lines": count_lines(run_dir / "variants_anno_combined_depth_frq_filter.xlsx"),
        },
        "final_stats": {
            "gfa": final_gfa_stats,
            "fasta": final_fasta_stats,
        },
    }
    write_json(report_dir / "metadata" / "report_manifest.json", report_manifest)
    write_genome_index(base, genome, cfg, report_dir, coverage_status, coverage_note)


def main() -> int:
    started_at = datetime.now().astimezone()
    start_time = time.perf_counter()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base",
        default="/Users/zouyinstein-home/Documents/Codex/hifisr-dev/hifisr/results/W3-5-2_gfa",
        help="W3-5-2_gfa result directory.",
    )
    parser.add_argument(
        "--final-coverage-y-rule",
        choices=["twice_mean", "original"],
        default=DEFAULT_FINAL_COVERAGE_Y_RULE,
        help="Y-axis rule for final_coverage; plotting style otherwise matches the original coverage plot.",
    )
    args = parser.parse_args()

    base = Path(args.base).expanduser().resolve()
    if not base.is_dir():
        raise SystemExit(f"Base directory does not exist: {base}")
    sampling = read_sampling_summary(base)
    for genome in GENOMES:
        generate_genome_report(base, genome, sampling, args.final_coverage_y_rule)
    elapsed_seconds = time.perf_counter() - start_time
    write_project_index(base, started_at, elapsed_seconds)
    print(f"Generated compact reports under {base / 'reports'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
