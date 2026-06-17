"""Helpers that project linear evidence onto verified GFA graph nodes."""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
import ast
import csv
import hashlib
import re

from Bio import SeqIO
import pandas as pd

from hifisr_functions.graph.verified_gfa import parse_gfa
from hifisr_functions.references import aln_to_ref


FUNCTION_PURITY = {
    "aln_to_ref": "impure",
    "build_node_projection_intervals": "impure",
    "compare_fasta_coordinates": "impure",
    "project_linear_interval": "pure",
    "project_sv_summary_files": "impure",
    "project_variant_dataframe": "pure",
    "project_variants_by_node": "pure",
    "read_node_source_rows": "impure",
    "write_linear_to_node_coordinate_map": "impure",
    "write_method_description": "impure",
    "write_simple_variant_node_outputs": "impure",
    "write_verified_gfa_by_nodes_fasta": "impure",
    "write_projection_outputs": "impure",
}


@dataclass(frozen=True)
class _ProjectionInterval:
    node: str
    node_class: str
    source_start: int
    source_end: int
    source_strand: str
    node_at_source_start: int
    node_at_source_end: int
    node_length: int
    source_copy_index: int
    source_kind: str
    source_notes: str

    @property
    def source_interval(self):
        return str(self.source_start) + "-" + str(self.source_end)

    def overlaps(self, start, end):
        return not (end < self.source_start or start > self.source_end)

    def project_point(self, position):
        if self.source_strand == "-":
            return self.node_at_source_start - (position - self.source_start)
        return self.node_at_source_start + (position - self.source_start)


def _read_first_fasta(path):
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 1:
        raise ValueError("Expected exactly one FASTA record in " + str(path))
    return records[0]


def _sequence_md5(sequence):
    return hashlib.md5(str(sequence).upper().encode()).hexdigest()


def compare_fasta_coordinates(linear_fasta, verified_fasta):
    linear_record = _read_first_fasta(linear_fasta)
    verified_record = _read_first_fasta(verified_fasta)
    linear_sequence = str(linear_record.seq).upper()
    verified_sequence = str(verified_record.seq).upper()
    first_difference = "."
    if linear_sequence != verified_sequence:
        for index, (left, right) in enumerate(zip(linear_sequence, verified_sequence), start=1):
            if left != right:
                first_difference = str(index)
                break
        if first_difference == "." and len(linear_sequence) != len(verified_sequence):
            first_difference = str(min(len(linear_sequence), len(verified_sequence)) + 1)
    return {
        "linear_fasta": str(linear_fasta),
        "verified_fasta": str(verified_fasta),
        "linear_record": linear_record.id,
        "verified_record": verified_record.id,
        "linear_length": len(linear_sequence),
        "verified_length": len(verified_sequence),
        "linear_md5": _sequence_md5(linear_sequence),
        "verified_md5": _sequence_md5(verified_sequence),
        "identical_sequence": linear_sequence == verified_sequence,
        "first_difference": first_difference,
    }


def _reverse_complement(sequence):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return str(sequence).translate(table)[::-1]


def _parse_parts(parts_text):
    interval_text, strand = str(parts_text).rsplit(":", 1)
    parts = []
    for token in interval_text.split(","):
        start, end = token.split("-", 1)
        parts.append((int(start), int(end)))
    return parts, strand


def _parse_notes(notes_text):
    notes = {}
    for item in str(notes_text).split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        notes[key] = value
    return notes


def _parse_core_copies(notes_text):
    core_copies = _parse_notes(notes_text).get("core_copies")
    if not core_copies:
        return []
    copies = []
    for token in core_copies.split(","):
        interval_text, strand = token.rsplit(":", 1)
        start, end = interval_text.split("-", 1)
        copies.append({"start": int(start), "end": int(end), "strand": strand})
    return copies


def _copy_strand_against_node(reference_sequence, segment_sequence, start, end, fallback="+"):
    reference_part = reference_sequence[start - 1:end].upper()
    segment_sequence = str(segment_sequence).upper()
    if reference_part == segment_sequence:
        return "+"
    if _reverse_complement(reference_part).upper() == segment_sequence:
        return "-"
    return fallback


def _intervals_for_parts(row, parts, strand):
    intervals = []
    ordered_parts = list(parts)
    if strand == "-":
        ordered_parts = list(reversed(ordered_parts))
    offset = 0
    for source_copy_index, (start, end) in enumerate(ordered_parts, start=1):
        length = end - start + 1
        if strand == "-":
            node_at_source_start = offset + length
            node_at_source_end = offset + 1
        else:
            node_at_source_start = offset + 1
            node_at_source_end = offset + length
        intervals.append(_ProjectionInterval(
            node=row["node"],
            node_class=row["node_class"],
            source_start=start,
            source_end=end,
            source_strand=strand,
            node_at_source_start=node_at_source_start,
            node_at_source_end=node_at_source_end,
            node_length=int(row["length"]),
            source_copy_index=source_copy_index,
            source_kind=row["source"],
            source_notes=row.get("notes", "."),
        ))
        offset += length
    return intervals


def _intervals_for_repeat(row, segment_sequence, reference_sequence):
    copies = _parse_core_copies(row.get("notes", "."))
    if not copies:
        parts, fallback_strand = _parse_parts(row["parts"])
        copies = [
            {"start": start, "end": end, "strand": fallback_strand}
            for start, end in parts
        ]
    intervals = []
    for source_copy_index, copy in enumerate(copies, start=1):
        strand = _copy_strand_against_node(
            reference_sequence,
            segment_sequence,
            copy["start"],
            copy["end"],
            fallback=copy.get("strand", "+"),
        )
        length = copy["end"] - copy["start"] + 1
        if strand == "-":
            node_at_source_start = length
            node_at_source_end = 1
        else:
            node_at_source_start = 1
            node_at_source_end = length
        intervals.append(_ProjectionInterval(
            node=row["node"],
            node_class=row["node_class"],
            source_start=copy["start"],
            source_end=copy["end"],
            source_strand=strand,
            node_at_source_start=node_at_source_start,
            node_at_source_end=node_at_source_end,
            node_length=int(row["length"]),
            source_copy_index=source_copy_index,
            source_kind=row["source"],
            source_notes=row.get("notes", "."),
        ))
    return intervals


def build_node_projection_intervals(node_source_map, verified_gfa, verified_fasta):
    gfa = parse_gfa(verified_gfa)
    reference_record = _read_first_fasta(verified_fasta)
    reference_sequence = str(reference_record.seq).upper()
    rows = read_node_source_rows(node_source_map)

    intervals = []
    for row in rows:
        if row["node_class"] == "repeat_node":
            segment = gfa.segments.get(row["node"])
            if segment is None:
                continue
            intervals.extend(_intervals_for_repeat(row, segment.sequence, reference_sequence))
        else:
            parts, strand = _parse_parts(row["parts"])
            intervals.extend(_intervals_for_parts(row, parts, strand))
    return sorted(intervals, key=lambda item: (item.source_start, item.source_end, item.node))


def read_node_source_rows(node_source_map):
    rows = []
    with open(node_source_map) as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def write_verified_gfa_by_nodes_fasta(path, verified_gfa, node_source_map):
    gfa = parse_gfa(verified_gfa)
    rows = read_node_source_rows(node_source_map)
    written = set()
    with open(path, "wt") as fout:
        for row in rows:
            node = row["node"]
            segment = gfa.segments.get(node)
            if segment is None or node in written:
                continue
            header = (
                ">"
                + node
                + "|node_class="
                + row.get("node_class", ".")
                + "|source_parts="
                + row.get("parts", ".")
                + "|direction=verified_gfa_segment"
            )
            print(header, file=fout)
            print(segment.sequence, file=fout)
            written.add(node)


def write_linear_to_node_coordinate_map(path, intervals):
    header = [
        "linear_start",
        "linear_end",
        "node",
        "node_linear_start",
        "node_linear_end",
        "strand",
        "node_length",
        "node_class",
        "source_copy_index",
        "source_kind",
    ]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for interval in intervals:
            writer.writerow({
                "linear_start": interval.source_start,
                "linear_end": interval.source_end,
                "node": interval.node,
                "node_linear_start": interval.node_at_source_start,
                "node_linear_end": interval.node_at_source_end,
                "strand": interval.source_strand,
                "node_length": interval.node_length,
                "node_class": interval.node_class,
                "source_copy_index": interval.source_copy_index,
                "source_kind": interval.source_kind,
            })


def _variant_ref_span(row):
    ref = row.get("ref", "")
    if pd.isna(ref):
        ref = ""
    ref_length = max(1, len(str(ref)))
    start = int(row["pos"])
    return start, start + ref_length - 1


def _projection_status(projections, start, end):
    if not projections:
        return "unmapped"
    span = end - start + 1
    covered = sum(item["linear_overlap_length"] for item in projections)
    nodes = {item["gfa_node"] for item in projections}
    if covered < span:
        return "partial_overlap"
    if len(nodes) == 1:
        return "single_node"
    return "split_across_nodes"


def project_linear_interval(intervals, start, end):
    projections = []
    for interval in intervals:
        if not interval.overlaps(start, end):
            continue
        overlap_start = max(start, interval.source_start)
        overlap_end = min(end, interval.source_end)
        node_a = interval.project_point(overlap_start)
        node_b = interval.project_point(overlap_end)
        node_start = min(node_a, node_b)
        node_end = max(node_a, node_b)
        projections.append({
            "gfa_node": interval.node,
            "node_class": interval.node_class,
            "node_start": node_start,
            "node_end": node_end,
            "node_strand": interval.source_strand,
            "node_length": interval.node_length,
            "linear_overlap_start": overlap_start,
            "linear_overlap_end": overlap_end,
            "linear_overlap_length": overlap_end - overlap_start + 1,
            "source_interval": interval.source_interval,
            "source_copy_index": interval.source_copy_index,
            "source_kind": interval.source_kind,
            "graph_locus": interval.node + ":" + str(node_start) + "-" + str(node_end),
        })
    return projections


def project_variant_dataframe(df, intervals, source_table):
    output_rows = []
    for row_index, (_, row) in enumerate(df.iterrows(), start=1):
        linear_start, linear_end = _variant_ref_span(row)
        base = row.to_dict()
        base["source_table"] = source_table
        base["variant_id"] = source_table + ":row" + f"{row_index:06d}"
        base["linear_start"] = linear_start
        base["linear_end"] = linear_end
        base["linear_span_length"] = linear_end - linear_start + 1
        projections = project_linear_interval(intervals, linear_start, linear_end)
        status = _projection_status(projections, linear_start, linear_end)
        if not projections:
            item = dict(base)
            item.update({
                "projection_status": status,
                "gfa_node": ".",
                "node_class": ".",
                "node_start": ".",
                "node_end": ".",
                "node_strand": ".",
                "node_length": ".",
                "linear_overlap_start": ".",
                "linear_overlap_end": ".",
                "linear_overlap_length": 0,
                "source_interval": ".",
                "source_copy_index": ".",
                "source_kind": ".",
                "graph_locus": ".",
            })
            output_rows.append(item)
            continue
        for projection in projections:
            item = dict(base)
            item["projection_status"] = status
            item.update(projection)
            output_rows.append(item)
    return pd.DataFrame(output_rows)


def _simple_node_coordinate(interval, start, end):
    node_a = interval.project_point(start)
    node_b = interval.project_point(end)
    node_start = min(node_a, node_b)
    node_end = max(node_a, node_b)
    if node_start == node_end:
        return str(node_start)
    return str(node_start) + "-" + str(node_end)


def _simple_variant_node(row, intervals):
    start, end = _variant_ref_span(row)
    matches = [
        interval for interval in intervals
        if start >= interval.source_start and end <= interval.source_end
    ]
    if len(matches) != 1:
        return "other", "."
    interval = matches[0]
    if start == interval.source_start or end == interval.source_end:
        return "other", "."
    return interval.node, _simple_node_coordinate(interval, start, end)


def project_variants_by_node(df, intervals):
    output = df.copy()
    nodes = []
    coordinates = []
    for _, row in output.iterrows():
        node, coordinate = _simple_variant_node(row, intervals)
        nodes.append(node)
        coordinates.append(coordinate)
    output["node"] = nodes
    output["node_linear_coordinate"] = coordinates
    return output


def write_simple_variant_node_outputs(variant_tables, intervals, node_source_map, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    node_order = [row["node"] for row in read_node_source_rows(node_source_map)]
    node_order.append("other")
    outputs = OrderedDict()
    for variant_table in variant_tables:
        variant_table = Path(variant_table)
        projected = project_variants_by_node(pd.read_excel(variant_table), intervals)
        output_path = output_dir / (variant_table.stem + ".by_verified_node.xlsx")
        used = set()
        with pd.ExcelWriter(output_path, engine="xlsxwriter") as writer:
            for node in node_order:
                node_df = projected[projected["node"] == node]
                if node_df.empty and node != "other":
                    continue
                sheet = _safe_sheet_name(node, used)
                node_df.to_excel(writer, sheet_name=sheet, index=False)
        outputs[variant_table.stem] = output_path
    return outputs


def _variant_node_summary(projected):
    usable = projected[projected["gfa_node"] != "."].copy()
    if usable.empty:
        return pd.DataFrame()
    summary = (
        usable.groupby(["source_table", "gfa_node", "node_class"], dropna=False)
        .agg(
            projected_variant_count=("variant_id", "nunique"),
            projection_row_count=("variant_id", "size"),
            mean_frequency=("frequency", "mean"),
            max_frequency=("frequency", "max"),
            mean_depth=("depth", "mean"),
            max_total_count=("total_count", "max"),
        )
        .reset_index()
    )
    type_counts = (
        usable.pivot_table(
            index=["source_table", "gfa_node"],
            columns="type",
            values="variant_id",
            aggfunc=pd.Series.nunique,
            fill_value=0,
        )
        .reset_index()
    )
    type_counts.columns = [
        str(col) if col in {"source_table", "gfa_node"} else "count_" + str(col)
        for col in type_counts.columns
    ]
    return summary.merge(type_counts, on=["source_table", "gfa_node"], how="left")


def _safe_sheet_name(name, used):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", str(name))[:25]
    if not safe:
        safe = "sheet"
    candidate = safe[:31]
    index = 1
    while candidate in used:
        suffix = "_" + str(index)
        candidate = safe[:31 - len(suffix)] + suffix
        index += 1
    used.add(candidate)
    return candidate


def _write_variant_outputs(projected, output_prefix):
    output_prefix = Path(output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)
    projected_tsv = output_prefix.with_suffix(".projected_by_node.tsv")
    summary_tsv = output_prefix.with_suffix(".node_summary.tsv")
    projected.to_csv(projected_tsv, sep="\t", index=False)
    summary = _variant_node_summary(projected)
    summary.to_csv(summary_tsv, sep="\t", index=False)

    node_dir = output_prefix.parent / (output_prefix.name + ".by_node_tsv")
    node_dir.mkdir(parents=True, exist_ok=True)
    sheet_index = []
    used = set()
    xlsx_path = output_prefix.with_suffix(".by_node.xlsx")
    with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as writer:
        summary.to_excel(writer, sheet_name="node_summary", index=False)
        for node, node_df in projected[projected["gfa_node"] != "."].groupby("gfa_node", sort=False):
            sheet = _safe_sheet_name(node, used)
            node_df.to_excel(writer, sheet_name=sheet, index=False)
            node_path = node_dir / (sheet + ".tsv")
            node_df.to_csv(node_path, sep="\t", index=False)
            sheet_index.append({
                "gfa_node": node,
                "sheet": sheet,
                "tsv": str(node_path),
                "variant_count": node_df["variant_id"].nunique(),
                "projection_row_count": len(node_df),
            })
        pd.DataFrame(sheet_index).to_excel(writer, sheet_name="sheet_index", index=False)
    sheet_index_tsv = node_dir / "sheet_index.tsv"
    pd.DataFrame(sheet_index).to_csv(sheet_index_tsv, sep="\t", index=False)
    return {
        "projected_tsv": projected_tsv,
        "summary_tsv": summary_tsv,
        "xlsx": xlsx_path,
        "node_dir": node_dir,
        "sheet_index_tsv": sheet_index_tsv,
    }


def _coordinate_labels(tuple_column):
    labels = re.findall(r"[a-z]+[0-9]+", str(tuple_column))
    return labels


def _parse_coordinate_tuple(value):
    if isinstance(value, str):
        return list(ast.literal_eval(value))
    if isinstance(value, tuple):
        return list(value)
    if isinstance(value, list):
        return value
    return []


def project_sv_summary_files(summary_files, intervals):
    event_rows = []
    breakpoint_rows = []
    for summary_file in sorted(Path(path) for path in summary_files):
        df = pd.read_excel(summary_file)
        tuple_columns = [
            column for column in df.columns
            if str(column).startswith("(") and str(column).endswith(")")
        ]
        if not tuple_columns:
            continue
        tuple_column = tuple_columns[0]
        labels = _coordinate_labels(tuple_column)
        sv_signature = summary_file.name.replace("_summary_backup.xlsx", "")
        for row_index, (_, row) in enumerate(df.iterrows(), start=1):
            coords = _parse_coordinate_tuple(row[tuple_column])
            strand_tokens = str(row.get("strand_str", ".")).split(",")
            event_id = sv_signature + ":row" + f"{row_index:06d}"
            event_nodes = []
            event_loci = []
            for coord_index, position in enumerate(coords, start=1):
                label = labels[coord_index - 1] if coord_index <= len(labels) else "coord" + str(coord_index)
                strand = strand_tokens[coord_index - 1] if coord_index <= len(strand_tokens) else "."
                projections = project_linear_interval(intervals, int(position), int(position))
                status = _projection_status(projections, int(position), int(position))
                if not projections:
                    breakpoint_rows.append({
                        "sv_event_id": event_id,
                        "sv_signature": sv_signature,
                        "source_file": str(summary_file),
                        "old_index": row.get("old_index", "."),
                        "coord_label": label,
                        "linear_pos": int(position),
                        "breakpoint_strand": strand,
                        "projection_status": status,
                        "gfa_node": ".",
                        "node_pos": ".",
                        "node_class": ".",
                        "source_interval": ".",
                        "graph_locus": ".",
                    })
                    event_nodes.append(".")
                    event_loci.append(".")
                    continue
                for projection in projections:
                    breakpoint_rows.append({
                        "sv_event_id": event_id,
                        "sv_signature": sv_signature,
                        "source_file": str(summary_file),
                        "old_index": row.get("old_index", "."),
                        "coord_label": label,
                        "linear_pos": int(position),
                        "breakpoint_strand": strand,
                        "projection_status": status,
                        "gfa_node": projection["gfa_node"],
                        "node_pos": projection["node_start"],
                        "node_class": projection["node_class"],
                        "source_interval": projection["source_interval"],
                        "graph_locus": projection["gfa_node"] + ":" + str(projection["node_start"]),
                    })
                    event_nodes.append(projection["gfa_node"])
                    event_loci.append(projection["gfa_node"] + ":" + str(projection["node_start"]))
            base = row.to_dict()
            base.update({
                "sv_event_id": event_id,
                "sv_signature": sv_signature,
                "source_file": str(summary_file),
                "breakpoint_count": len(coords),
                "breakpoint_linear_positions": ",".join(str(int(item)) for item in coords),
                "breakpoint_graph_nodes": ">".join(event_nodes),
                "breakpoint_graph_loci": ">".join(event_loci),
                "sv_projection_level": "breakpoint_signature",
                "sv_projection_note": (
                    "SV length/topology is not inferred here; each candidate is "
                    "represented by projected breakpoint coordinates and read-group counts."
                ),
            })
            event_rows.append(base)
    return pd.DataFrame(event_rows), pd.DataFrame(breakpoint_rows)


def _write_sv_outputs(event_df, breakpoint_df, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    event_tsv = output_dir / "sv_event_projection.tsv"
    breakpoint_tsv = output_dir / "sv_breakpoint_projection.tsv"
    summary_tsv = output_dir / "sv_event_node_summary.tsv"
    xlsx_path = output_dir / "sv_breakpoint_projection.xlsx"
    event_df.to_csv(event_tsv, sep="\t", index=False)
    breakpoint_df.to_csv(breakpoint_tsv, sep="\t", index=False)
    if breakpoint_df.empty:
        summary = pd.DataFrame()
    else:
        usable = breakpoint_df[breakpoint_df["gfa_node"] != "."]
        summary = (
            usable.groupby(["sv_signature", "gfa_node", "node_class"], dropna=False)
            .agg(
                event_count=("sv_event_id", "nunique"),
                breakpoint_count=("sv_event_id", "size"),
            )
            .reset_index()
        )
    summary.to_csv(summary_tsv, sep="\t", index=False)
    with pd.ExcelWriter(xlsx_path, engine="xlsxwriter") as writer:
        event_df.to_excel(writer, sheet_name="sv_events", index=False)
        breakpoint_df.to_excel(writer, sheet_name="sv_breakpoints", index=False)
        summary.to_excel(writer, sheet_name="node_summary", index=False)
    return {
        "event_tsv": event_tsv,
        "breakpoint_tsv": breakpoint_tsv,
        "summary_tsv": summary_tsv,
        "xlsx": xlsx_path,
    }


def write_method_description(path, coordinate_check):
    text = f"""# Linear-to-Graph Coordinate Projection

The run-specific variant tables are reported in the linear coordinate system of `{Path(coordinate_check['linear_fasta']).name}`. For this run, that FASTA is coordinate-identical to `{Path(coordinate_check['verified_fasta']).name}`: both records are `{coordinate_check['linear_record']}`, both are {coordinate_check['linear_length']} bp long, and their sequence MD5 checksums are identical (`{coordinate_check['linear_md5']}`). Therefore a variant position from `run_3` can be projected directly onto the verified GFA coordinate model.

The projection uses `verified_node_source_map.tsv` as the authoritative bridge between the linear verified FASTA and the graph. Each single-copy node is represented by one or more source intervals. If a node was reverse-complemented when it was written to the GFA, the projected node coordinate is reversed within that interval. Repeat nodes are handled as collapsed graph nodes: all repeat-core copies recorded in the `core_copies` field are mapped back to the same GFA segment, and the copy orientation is inferred by comparing the GFA segment sequence with the corresponding reference interval and its reverse complement.

For a plus-strand source interval, node-local coordinates are computed as `node_pos = node_interval_start + (linear_pos - source_start)`. For a minus-strand source interval, they are computed as `node_pos = node_interval_start - (linear_pos - source_start)`, where `node_interval_start` is the node coordinate corresponding to `source_start`. For SNVs, the projected interval is the single 1-based linear position. For InDels, the projected interval is the 1-based inclusive span covered by the REF allele (`pos` to `pos + len(ref) - 1`), with insertions anchored on their REF base. A variant that overlaps one graph interval receives a `single_node` projection; a variant crossing graph-node boundaries is emitted as multiple rows with `split_across_nodes`; uncovered bases are marked `partial_overlap` or `unmapped`.

SV-like evidence from `type_*_summary_backup.xlsx` is projected at the breakpoint-signature level instead of forcing a possibly misleading SV length estimate. Each coordinate in the tuple field, such as `(se1, ss2)` or `(se1, ss2, se2, ss3)`, is projected independently to a node-local breakpoint. The event-level table keeps the ordered graph-node path of the breakpoints and the original subgroup read counts, while the breakpoint table provides the report-ready node-local loci.
"""
    Path(path).write_text(text)


def write_projection_outputs(
    variant_tables,
    sv_summary_files,
    intervals,
    coordinate_check,
    output_dir,
):
    output_dir = Path(output_dir)
    snv_indel_dir = output_dir / "snv_indel_by_node"
    sv_dir = output_dir / "sv_breakpoints_by_node"
    snv_indel_dir.mkdir(parents=True, exist_ok=True)
    outputs = OrderedDict()
    projected_tables = []
    for variant_table in variant_tables:
        variant_table = Path(variant_table)
        df = pd.read_excel(variant_table)
        projected = project_variant_dataframe(df, intervals, variant_table.stem)
        projected_tables.append(projected)
        outputs[variant_table.stem] = _write_variant_outputs(
            projected,
            snv_indel_dir / variant_table.stem,
        )
    if projected_tables:
        combined = pd.concat(projected_tables, ignore_index=True)
        outputs["combined_variant_node_summary"] = snv_indel_dir / "combined_variant_node_summary.tsv"
        _variant_node_summary(combined).to_csv(outputs["combined_variant_node_summary"], sep="\t", index=False)

    event_df, breakpoint_df = project_sv_summary_files(sv_summary_files, intervals)
    outputs["sv"] = _write_sv_outputs(event_df, breakpoint_df, sv_dir)

    interval_rows = [interval.__dict__ for interval in intervals]
    intervals_tsv = output_dir / "linear_to_graph_projection_intervals.tsv"
    pd.DataFrame(interval_rows).to_csv(intervals_tsv, sep="\t", index=False)
    outputs["projection_intervals"] = intervals_tsv

    coordinate_tsv = output_dir / "coordinate_consistency.tsv"
    pd.DataFrame([coordinate_check]).to_csv(coordinate_tsv, sep="\t", index=False)
    outputs["coordinate_consistency"] = coordinate_tsv

    method_md = output_dir / "linear_to_graph_coordinate_projection.md"
    write_method_description(method_md, coordinate_check)
    outputs["method_description"] = method_md
    return outputs


__all__ = [
    "aln_to_ref",
    "build_node_projection_intervals",
    "compare_fasta_coordinates",
    "project_linear_interval",
    "project_sv_summary_files",
    "project_variant_dataframe",
    "project_variants_by_node",
    "read_node_source_rows",
    "write_linear_to_node_coordinate_map",
    "write_method_description",
    "write_simple_variant_node_outputs",
    "write_verified_gfa_by_nodes_fasta",
    "write_projection_outputs",
]
