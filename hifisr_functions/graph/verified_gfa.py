"""Build a reviewable GFA against a run_3 verified FASTA.

This module never performs repeat resolution. It copies a curated merged GFA
when available, or only compacts clearly linear non-repeat paths as a fallback.
"""

from __future__ import annotations

from collections import OrderedDict, defaultdict
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import math
import os
import re
import shutil
import subprocess
import tempfile

from Bio import SeqIO
import pandas as pd
import hifisr_functions.reports as hfrps


FUNCTION_PURITY = {
    "annotate_merged_gfa_evidence": "impure",
    "build_verified_gfa": "impure",
    "merge_unambiguous_gfa": "pure",
    "parse_gfa": "impure",
    "write_repeat_path_support": "impure",
    "write_repeat_read_path_assignments": "impure",
    "write_repeat_read_path_support": "impure",
    "write_gfa": "impure",
    "write_node_remapped_coverage": "impure",
    "write_auto_repeat_check": "impure",
}


@dataclass
class Segment:
    name: str
    sequence: str
    tags: list[str]


@dataclass
class Link:
    from_name: str
    from_orient: str
    to_name: str
    to_orient: str
    overlap: str
    tags: list[str]


@dataclass
class Gfa:
    headers: list[list[str]]
    segments: OrderedDict[str, Segment]
    links: list[Link]
    other_lines: list[list[str]]


@dataclass
class PafHit:
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    matches: int
    block_length: int
    mapq: int

    @property
    def identity(self) -> float:
        if self.block_length == 0:
            return 0.0
        return self.matches / self.block_length

    @property
    def aligned_fraction(self) -> float:
        if self.query_length == 0:
            return 0.0
        return (self.query_end - self.query_start) / self.query_length

    @property
    def target_start_1based(self) -> int:
        return self.target_start + 1

    @property
    def target_end_1based(self) -> int:
        return self.target_end

    @property
    def query_start_1based(self) -> int:
        return self.query_start + 1

    @property
    def query_end_1based(self) -> int:
        return self.query_end


@dataclass
class ReadAlignment:
    read_id: str
    read_length: int
    query_start: int
    query_end: int
    target_id: str
    target_start: int
    target_end: int
    strand: str
    identity: float
    mapq: int
    alignment_role: str


def parse_gfa(path):
    headers = []
    segments = OrderedDict()
    links = []
    other_lines = []
    with open(path) as fin:
        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if fields[0] == "H":
                headers.append(fields)
            elif fields[0] == "S":
                segments[fields[1]] = Segment(fields[1], fields[2], fields[3:])
            elif fields[0] == "L":
                links.append(Link(fields[1], fields[2], fields[3], fields[4], fields[5], fields[6:]))
            else:
                other_lines.append(fields)
    return Gfa(headers, segments, links, other_lines)


def write_gfa(gfa, path):
    with open(path, "wt") as fout:
        for fields in gfa.headers:
            print("\t".join(fields), file=fout)
        for segment in gfa.segments.values():
            print("\t".join(["S", segment.name, segment.sequence] + segment.tags), file=fout)
        for link in gfa.links:
            print(
                "\t".join([
                    "L",
                    link.from_name,
                    link.from_orient,
                    link.to_name,
                    link.to_orient,
                    link.overlap,
                ] + link.tags),
                file=fout,
            )
        for fields in gfa.other_lines:
            print("\t".join(fields), file=fout)


def reverse_complement(sequence):
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return sequence.translate(table)[::-1]


def flip_orient(orient):
    return "-" if orient == "+" else "+"


def natural_key(value):
    return [int(part) if part.isdigit() else part for part in re.split(r"(\d+)", value)]


def tag_value(tags, tag_name):
    prefix = tag_name + ":"
    for tag in tags:
        if tag.startswith(prefix):
            fields = tag.split(":", 2)
            if len(fields) == 3:
                return fields[2]
    return None


def numeric_tag(tags, tag_name):
    value = tag_value(tags, tag_name)
    if value is None:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def integer_tag(tags, tag_name):
    value = tag_value(tags, tag_name)
    if value is None:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def replace_tag(tags, tag):
    tag_name = tag.split(":", 1)[0]
    return [item for item in tags if not item.startswith(tag_name + ":")] + [tag]


def replace_tags(tags, new_tags):
    updated = list(tags)
    for tag in new_tags:
        updated = replace_tag(updated, tag)
    return updated


def segment_depth(segment):
    depth = numeric_tag(segment.tags, "DP")
    if depth is None:
        depth = numeric_tag(segment.tags, "dp")
    return depth


def segment_raw_ids(name, raw_gfa):
    if name in raw_gfa.segments:
        return [name]
    fields = name.split("_")
    if all(field in raw_gfa.segments for field in fields):
        return fields
    return []


def coverage_tags_for_segment(name, segment, raw_gfa):
    raw_ids = segment_raw_ids(name, raw_gfa)
    raw_segments = [raw_gfa.segments[item] for item in raw_ids]
    if not raw_segments:
        depth = segment_depth(segment)
        if depth is None:
            return ["CM:Z:coverage_unavailable"]
        return [
            "DP:f:" + f"{depth:.6f}",
            "CM:Z:existing_segment_DP",
        ]

    total_len = sum(len(item.sequence) for item in raw_segments)
    depth_bases = 0.0
    depth_len = 0
    ab_sum = 0
    ac_sum = 0
    has_ab = False
    has_ac = False
    for item in raw_segments:
        depth = segment_depth(item)
        if depth is not None:
            depth_bases += depth * len(item.sequence)
            depth_len += len(item.sequence)
        ab_value = integer_tag(item.tags, "AB")
        ac_value = integer_tag(item.tags, "AC")
        if ab_value is not None:
            ab_sum += ab_value
            has_ab = True
        if ac_value is not None:
            ac_sum += ac_value
            has_ac = True

    tags = [
        "CM:Z:raw_node_DP_length_weighted",
        "RL:i:" + str(total_len),
    ]
    if depth_len > 0:
        tags.insert(0, "DP:f:" + f"{depth_bases / depth_len:.6f}")
        tags.append("DB:f:" + f"{depth_bases:.3f}")
    if has_ab:
        tags.append("AB:i:" + str(ab_sum))
    if has_ac:
        tags.append("AC:i:" + str(ac_sum))
    return tags


def write_segment_fasta(segments, path):
    with open(path, "wt") as fout:
        for segment in segments.values():
            print(">" + segment.name, file=fout)
            print(segment.sequence, file=fout)


def parse_paf_line(line):
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 12:
        return None
    return PafHit(
        query_name=fields[0],
        query_length=int(fields[1]),
        query_start=int(fields[2]),
        query_end=int(fields[3]),
        strand=fields[4],
        target_name=fields[5],
        target_length=int(fields[6]),
        target_start=int(fields[7]),
        target_end=int(fields[8]),
        matches=int(fields[9]),
        block_length=int(fields[10]),
        mapq=int(fields[11]),
    )


def run_minimap2_paf(query_fasta, target_fasta, soft_paths_dict, threads):
    minimap2 = soft_paths_dict.get("minimap2", "minimap2")
    command = [
        minimap2,
        "-x",
        "asm5",
        "-c",
        "--eqx",
        "-N",
        "20",
        "-t",
        str(threads),
        str(target_fasta),
        str(query_fasta),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            "minimap2 failed with exit code "
            + str(completed.returncode)
            + "\n"
            + completed.stderr
        )
    hits_by_node = defaultdict(list)
    for line in completed.stdout.splitlines():
        hit = parse_paf_line(line)
        if hit is not None:
            hits_by_node[hit.query_name].append(hit)
    return hits_by_node, " ".join(command), completed.stderr


def best_hit(hits):
    if not hits:
        return None
    return sorted(
        hits,
        key=lambda hit: (
            hit.aligned_fraction,
            hit.identity,
            hit.mapq,
            hit.block_length,
        ),
        reverse=True,
    )[0]


def confidence_for_hits(hits):
    hit = best_hit(hits)
    if hit is None:
        return "unmapped", 0.0, "no_alignment"
    good_hits = [
        item
        for item in hits
        if item.identity >= 0.995 and item.aligned_fraction >= 0.80
    ]
    unique = len(good_hits) <= 1
    score = hit.identity * min(hit.aligned_fraction, 1.0)
    if hit.identity >= 0.995 and hit.aligned_fraction >= 0.98 and unique:
        return "high", score, "near_full_length_unique_alignment"
    if hit.identity >= 0.990 and hit.aligned_fraction >= 0.90 and unique:
        return "medium", score, "partial_or_terminal_uncertainty"
    if len(good_hits) > 1:
        return "repeat_or_multimap", score, "multiple_high_identity_alignments"
    return "low", score, "low_identity_or_partial_alignment"


def graph_degrees(gfa):
    degrees = defaultdict(int)
    for link in gfa.links:
        degrees[link.from_name] += 1
        degrees[link.to_name] += 1
    return degrees


def median_depth_for_gfa(gfa):
    depths = [
        segment_depth(segment)
        for segment in gfa.segments.values()
        if segment_depth(segment) is not None
    ]
    if not depths:
        return None
    depths_sorted = sorted(depths)
    middle = len(depths_sorted) // 2
    if len(depths_sorted) % 2 == 1:
        return depths_sorted[middle]
    return (depths_sorted[middle - 1] + depths_sorted[middle]) / 2


def repeat_like_nodes(gfa):
    degrees = graph_degrees(gfa)
    branch_degrees = [degrees[name] for name in gfa.segments if degrees[name] >= 3]
    if not branch_degrees:
        return set()
    branch_floor = min(branch_degrees)
    has_higher_order_branch = any(
        degrees[name] > branch_floor for name in gfa.segments
    )
    if has_higher_order_branch:
        return {
            name
            for name in gfa.segments
            if degrees[name] > branch_floor
        }
    return {
        name
        for name in gfa.segments
        if degrees[name] >= 3
    }


def link_support(link):
    for tag_name in ["RC", "SK", "PA"]:
        value = numeric_tag(link.tags, tag_name)
        if value is not None:
            return value, tag_name
    return 0.0, "missing"


def link_endpoint_side(link, node_name):
    if link.from_name == node_name:
        return "-" if link.from_orient == "-" else "+"
    if link.to_name == node_name:
        return "+" if link.to_orient == "-" else "-"
    return None


def other_link_endpoint(link, node_name):
    if link.from_name == node_name:
        return link.to_name + link.to_orient
    if link.to_name == node_name:
        return link.from_name + link.from_orient
    return "."


def repeat_path_support_rows(gfa, node_classes):
    rows = []
    for node_name, node_class in node_classes.items():
        if node_class != "repeat_node":
            continue
        side_links = {"-": [], "+": []}
        for link_index, link in enumerate(gfa.links, start=1):
            side = link_endpoint_side(link, node_name)
            if side in side_links:
                support, support_tag = link_support(link)
                side_links[side].append({
                    "index": link_index,
                    "link": link,
                    "side": side,
                    "support": support,
                    "support_tag": support_tag,
                    "endpoint": other_link_endpoint(link, node_name),
                })
        if len(side_links["-"]) == 2 and len(side_links["+"]) == 2:
            repeat_status = "unresolved"
        elif len(side_links["-"]) <= 1 and len(side_links["+"]) <= 1:
            repeat_status = "resolved"
        else:
            repeat_status = "ambiguous"
        combinations = []
        for left in side_links["-"]:
            for right in side_links["+"]:
                support = min(left["support"], right["support"])
                combinations.append({
                    "repeat_node": node_name,
                    "repeat_status": repeat_status,
                    "support_method": "edge_rc_min_proxy",
                    "left_edge_index": left["index"],
                    "right_edge_index": right["index"],
                    "left_endpoint": left["endpoint"],
                    "right_endpoint": right["endpoint"],
                    "left_support": left["support"],
                    "right_support": right["support"],
                    "path_support": support,
                    "path_ratio": 0.0,
                })
        total_support = sum(item["path_support"] for item in combinations)
        for item in combinations:
            if total_support > 0:
                item["path_ratio"] = item["path_support"] / total_support
            rows.append(item)
        if not combinations:
            rows.append({
                "repeat_node": node_name,
                "repeat_status": repeat_status,
                "support_method": "edge_rc_min_proxy",
                "left_edge_index": ".",
                "right_edge_index": ".",
                "left_endpoint": ".",
                "right_endpoint": ".",
                "left_support": ".",
                "right_support": ".",
                "path_support": ".",
                "path_ratio": ".",
            })
    return rows


def repeat_path_tags(node_name, rows):
    node_rows = [item for item in rows if item["repeat_node"] == node_name]
    if not node_rows:
        return []
    status = node_rows[0]["repeat_status"]
    method = node_rows[0]["support_method"]
    path_items = []
    for index, item in enumerate(node_rows, start=1):
        path_name = (
            str(item["left_endpoint"])
            + ">"
            + node_name
            + ">"
            + str(item["right_endpoint"])
        )
        support_value = item.get("path_count", item.get("path_support", "."))
        ratio_value = item.get("path_ratio", ".")
        if support_value == "." or ratio_value == ".":
            continue
        path_items.append(
            "p"
            + str(index)
            + "="
            + path_name
            + ":count="
            + f"{float(support_value):.3f}"
            + ":ratio="
            + f"{float(ratio_value):.6f}"
        )
    tags = [
        "RT:Z:" + status,
        "RM:Z:graph_degree_2in2out" if status == "unresolved" else "RM:Z:graph_degree",
        "PM:Z:" + method,
    ]
    if path_items:
        tags.append("PC:Z:" + ",".join(path_items))
    return tags


def _is_missing_value(value):
    return value is None or value == "" or value == "." or (
        isinstance(value, float) and math.isnan(value)
    )


def _sanitize_gfa_name(value):
    sanitized = re.sub(r"[^A-Za-z0-9_.:-]+", "_", str(value))
    return sanitized if sanitized else "unknown"


def _split_gfa_endpoint(endpoint):
    endpoint = str(endpoint)
    if len(endpoint) < 2 or endpoint[-1] not in "+-":
        return None
    return endpoint[:-1], endpoint[-1]


def _flip_gfa_endpoint(endpoint):
    parsed = _split_gfa_endpoint(endpoint)
    if parsed is None:
        return "."
    name, orient = parsed
    return name + flip_orient(orient)


def _repeat_path_numeric_tag(tag_name, tag_type, value, digits=6):
    if _is_missing_value(value):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if tag_type == "i":
        return tag_name + ":i:" + str(int(round(number)))
    return tag_name + ":f:" + f"{number:.{digits}f}"


def _repeat_path_text_tag(tag_name, value):
    if _is_missing_value(value):
        return None
    safe_value = str(value).replace("\t", "_").replace("\n", "_")
    return tag_name + ":Z:" + safe_value


def repeat_path_gfa_lines(repeat_rows):
    path_lines = []
    for index, row in enumerate(repeat_rows, start=1):
        repeat_node = row.get("repeat_node")
        if _is_missing_value(repeat_node):
            continue
        left_endpoint = row.get("left_endpoint", row.get("left_gfa_endpoint"))
        right_endpoint = row.get("right_endpoint", row.get("right_gfa_endpoint"))
        left_segment = _split_gfa_endpoint(left_endpoint)
        right_segment = _split_gfa_endpoint(_flip_gfa_endpoint(right_endpoint))
        if left_segment is None or right_segment is None:
            continue
        path_id = row.get("path_id", "p" + str(index))
        path_name = (
            "repeat_"
            + _sanitize_gfa_name(repeat_node)
            + "_"
            + _sanitize_gfa_name(path_id)
        )
        segment_names = [
            left_segment[0] + left_segment[1],
            str(repeat_node) + "+",
            right_segment[0] + right_segment[1],
        ]
        tags = [
            "PT:Z:repeat_path_support",
            _repeat_path_text_tag("RN", repeat_node),
            _repeat_path_text_tag("PI", path_id),
            _repeat_path_text_tag("RS", row.get("repeat_status")),
            _repeat_path_text_tag("PM", row.get("support_method")),
            _repeat_path_numeric_tag(
                "RC",
                "f",
                row.get("path_count", row.get("path_support")),
                digits=3,
            ),
            _repeat_path_numeric_tag("PR", "f", row.get("path_ratio")),
            _repeat_path_numeric_tag("UR", "i", row.get("unique_assigned_read_count")),
            _repeat_path_numeric_tag("AR", "i", row.get("ambiguous_read_count")),
            _repeat_path_numeric_tag("CR", "i", row.get("candidate_read_count")),
            _repeat_path_numeric_tag("FR", "i", row.get("total_fl_read_count")),
            _repeat_path_text_tag("LE", left_endpoint),
            _repeat_path_text_tag("RE", right_endpoint),
        ]
        path_lines.append([
            "P",
            path_name,
            ",".join(segment_names),
            "*,*",
        ] + [tag for tag in tags if tag is not None])
    return path_lines


def add_repeat_path_gfa_lines(gfa, repeat_rows):
    other_lines = [
        fields for fields in gfa.other_lines
        if not (
            fields
            and fields[0] == "P"
            and "PT:Z:repeat_path_support" in fields[4:]
        )
    ]
    other_lines.extend(repeat_path_gfa_lines(repeat_rows))
    return Gfa(gfa.headers, gfa.segments, gfa.links, other_lines)


def interval_overlap(start_a, end_a, start_b, end_b):
    start = max(int(start_a), int(start_b))
    end = min(int(end_a), int(end_b))
    if end < start:
        return 0
    return end - start + 1


def interval_overlap_fraction(start_a, end_a, start_b, end_b):
    overlap = interval_overlap(start_a, end_a, start_b, end_b)
    shortest = min(end_a - start_a + 1, end_b - start_b + 1)
    if shortest <= 0:
        return 0.0
    return overlap / shortest


def select_distinct_hits(hits, max_hits=None, target_overlap_cutoff=0.80, query_overlap_cutoff=0.80):
    selected = []
    sorted_hits = sorted(
        hits,
        key=lambda hit: (
            hit.identity,
            hit.aligned_fraction,
            hit.mapq,
            hit.block_length,
        ),
        reverse=True,
    )
    for hit in sorted_hits:
        redundant = False
        for old_hit in selected:
            target_overlap = interval_overlap_fraction(
                hit.target_start_1based,
                hit.target_end_1based,
                old_hit.target_start_1based,
                old_hit.target_end_1based,
            )
            query_overlap = interval_overlap_fraction(
                hit.query_start_1based,
                hit.query_end_1based,
                old_hit.query_start_1based,
                old_hit.query_end_1based,
            )
            if target_overlap >= target_overlap_cutoff and query_overlap >= query_overlap_cutoff:
                redundant = True
                break
        if not redundant:
            selected.append(hit)
        if max_hits is not None and len(selected) >= max_hits:
            break
    return sorted(selected, key=lambda hit: (hit.query_start, hit.target_start))


def selected_mapping_hits_for_node(segment, hits, node_class, max_repeat_copies=2):
    if not hits:
        return []
    if node_class == "repeat_node":
        candidates = [
            hit
            for hit in hits
            if hit.identity >= 0.990 and hit.aligned_fraction >= 0.75
        ]
        if not candidates:
            hit = best_hit(hits)
            candidates = [] if hit is None else [hit]
        return select_distinct_hits(candidates, max_hits=max_repeat_copies)

    candidates = [
        hit
        for hit in hits
        if hit.identity >= 0.970 and hit.aligned_fraction >= 0.05
    ]
    if not candidates:
        hit = best_hit(hits)
        candidates = [] if hit is None else [hit]
    return select_distinct_hits(candidates)


def map_query_interval_to_target(hit, query_start, query_end):
    query_start = max(int(query_start), hit.query_start_1based)
    query_end = min(int(query_end), hit.query_end_1based)
    if query_end < query_start:
        return None

    query_span = max(1, hit.query_end_1based - hit.query_start_1based + 1)
    target_span = max(1, hit.target_end_1based - hit.target_start_1based + 1)
    scale = target_span / query_span
    if hit.strand == "+":
        target_a = hit.target_start_1based + (query_start - hit.query_start_1based) * scale
        target_b = hit.target_start_1based + (query_end - hit.query_start_1based) * scale
    else:
        target_a = hit.target_end_1based - (query_start - hit.query_start_1based) * scale
        target_b = hit.target_end_1based - (query_end - hit.query_start_1based) * scale
    target_start = max(1, int(round(min(target_a, target_b))))
    target_end = min(hit.target_length, int(round(max(target_a, target_b))))
    if target_end < target_start:
        return None
    return target_start, target_end


def endpoint_anchor_intervals(
    node_name,
    segment,
    hits,
    endpoint_side,
    anchor_length,
    min_anchor_overlap,
    flank_anchor_offset=1000,
):
    node_length = len(segment.sequence)
    anchor_length = min(anchor_length, node_length)
    flank_anchor_offset = min(flank_anchor_offset, max(0, node_length - anchor_length))
    if endpoint_side == "-":
        query_start = flank_anchor_offset + 1
        query_end = flank_anchor_offset + anchor_length
    else:
        query_start = max(1, node_length - flank_anchor_offset - anchor_length + 1)
        query_end = node_length - flank_anchor_offset

    anchors = []
    for index, hit in enumerate(hits, start=1):
        interval = map_query_interval_to_target(hit, query_start, query_end)
        if interval is None:
            continue
        start, end = interval
        if end - start + 1 < min_anchor_overlap:
            continue
        anchors.append({
            "anchor_id": node_name + endpoint_side + ":copy" + str(index),
            "node": node_name,
            "side": endpoint_side,
            "target_name": hit.target_name,
            "start": start,
            "end": end,
            "strand": hit.strand,
            "hit_index": index,
        })
    return anchors


def repeat_anchor_intervals(node_name, segment, hits, anchor_length, min_anchor_overlap):
    node_length = len(segment.sequence)
    anchor_length = min(anchor_length, node_length)
    center = (node_length + 1) // 2
    query_start = max(1, center - anchor_length // 2)
    query_end = min(node_length, query_start + anchor_length - 1)
    anchors = []
    for index, hit in enumerate(hits, start=1):
        interval = map_query_interval_to_target(hit, query_start, query_end)
        if interval is None:
            continue
        start, end = interval
        if end - start + 1 < min_anchor_overlap:
            continue
        anchors.append({
            "anchor_id": node_name + ":copy" + str(index),
            "node": node_name,
            "side": "core",
            "target_name": hit.target_name,
            "start": start,
            "end": end,
            "strand": hit.strand,
            "hit_index": index,
        })
    return anchors


def flank_endpoint_for_repeat_link(link, repeat_node):
    if link.from_name == repeat_node:
        flank_node = link.to_name
    elif link.to_name == repeat_node:
        flank_node = link.from_name
    else:
        return None
    endpoint_side = link_endpoint_side(link, flank_node)
    if endpoint_side is None:
        return None
    return {
        "node": flank_node,
        "endpoint_side": endpoint_side,
        "gfa_endpoint": other_link_endpoint(link, repeat_node),
    }


def build_repeat_path_candidates(
    gfa,
    node_classes,
    selected_hits_by_node,
    target_name,
    flank_anchor_length=500,
    flank_anchor_offset=1000,
    repeat_anchor_length=500,
    min_anchor_overlap=200,
):
    candidates = []
    links_by_repeat = defaultdict(lambda: {"-": [], "+": []})
    for link_index, link in enumerate(gfa.links, start=1):
        for node_name in [link.from_name, link.to_name]:
            if node_classes.get(node_name) != "repeat_node":
                continue
            side = link_endpoint_side(link, node_name)
            if side in links_by_repeat[node_name]:
                links_by_repeat[node_name][side].append({
                    "index": link_index,
                    "link": link,
                    "repeat_side": side,
                })

    for repeat_node in sorted(links_by_repeat, key=natural_key):
        side_links = links_by_repeat[repeat_node]
        if len(side_links["-"]) == 2 and len(side_links["+"]) == 2:
            repeat_status = "unresolved"
        elif len(side_links["-"]) <= 1 and len(side_links["+"]) <= 1:
            repeat_status = "resolved"
        else:
            repeat_status = "ambiguous"

        repeat_segment = gfa.segments[repeat_node]
        repeat_anchors = repeat_anchor_intervals(
            repeat_node,
            repeat_segment,
            selected_hits_by_node.get(repeat_node, []),
            repeat_anchor_length,
            min_anchor_overlap,
        )
        repeat_copy_count = len(repeat_anchors)
        path_index = 0
        for left in side_links["-"]:
            left_flank = flank_endpoint_for_repeat_link(left["link"], repeat_node)
            if left_flank is None:
                continue
            left_segment = gfa.segments[left_flank["node"]]
            left_anchors = endpoint_anchor_intervals(
                left_flank["node"],
                left_segment,
                selected_hits_by_node.get(left_flank["node"], []),
                left_flank["endpoint_side"],
                flank_anchor_length,
                min_anchor_overlap,
                flank_anchor_offset=flank_anchor_offset,
            )
            for right in side_links["+"]:
                right_flank = flank_endpoint_for_repeat_link(right["link"], repeat_node)
                if right_flank is None:
                    continue
                right_segment = gfa.segments[right_flank["node"]]
                right_anchors = endpoint_anchor_intervals(
                    right_flank["node"],
                    right_segment,
                    selected_hits_by_node.get(right_flank["node"], []),
                    right_flank["endpoint_side"],
                    flank_anchor_length,
                    min_anchor_overlap,
                    flank_anchor_offset=flank_anchor_offset,
                )
                path_index += 1
                candidates.append({
                    "repeat_node": repeat_node,
                    "repeat_status": repeat_status,
                    "path_id": "p" + str(path_index),
                    "left_edge_index": left["index"],
                    "right_edge_index": right["index"],
                    "left_endpoint": left_flank["node"] + left_flank["endpoint_side"],
                    "right_endpoint": right_flank["node"] + right_flank["endpoint_side"],
                    "left_gfa_endpoint": left_flank["gfa_endpoint"],
                    "right_gfa_endpoint": right_flank["gfa_endpoint"],
                    "left_anchors": left_anchors,
                    "repeat_anchors": repeat_anchors,
                    "right_anchors": right_anchors,
                    "repeat_copy_count": repeat_copy_count,
                    "target_name": target_name,
                })
    return candidates


def load_fl_read_ids(run_dir):
    path = Path(run_dir) / "FL_ids.txt"
    if not path.exists():
        return set(), path
    with open(path) as fin:
        return {line.strip() for line in fin if line.strip()}, path


def parse_read_alignment_row(row):
    try:
        return ReadAlignment(
            read_id=row["read_id"],
            read_length=int(float(row["read_length"])),
            query_start=int(float(row["query_start"])),
            query_end=int(float(row["query_end"])),
            target_id=row["target_id"],
            target_start=int(float(row["target_start"])),
            target_end=int(float(row["target_end"])),
            strand=row["strand"],
            identity=float(row["identity"]),
            mapq=int(float(row["mapq"])),
            alignment_role=row.get("alignment_role", "."),
        )
    except (KeyError, TypeError, ValueError):
        return None


def load_fl_whole_read_evidence(run_dir):
    run_dir = Path(run_dir)
    fl_ids, fl_ids_path = load_fl_read_ids(run_dir)
    evidence_path = run_dir / "backup_info" / "whole_read_evidence.tsv"
    alignments = defaultdict(list)
    total_rows = 0
    kept_rows = 0
    if not evidence_path.exists():
        return alignments, {
            "fl_ids_path": fl_ids_path,
            "whole_read_evidence_path": evidence_path,
            "fl_read_count": len(fl_ids),
            "total_evidence_rows": 0,
            "kept_fl_evidence_rows": 0,
        }
    with open(evidence_path) as fin:
        reader = csv.DictReader(fin, delimiter="\t")
        for row in reader:
            total_rows += 1
            read_id = row.get("read_id", "")
            if fl_ids and read_id not in fl_ids:
                continue
            alignment = parse_read_alignment_row(row)
            if alignment is None:
                continue
            alignments[read_id].append(alignment)
            kept_rows += 1
    return alignments, {
        "fl_ids_path": fl_ids_path,
        "whole_read_evidence_path": evidence_path,
        "fl_read_count": len(fl_ids),
        "total_evidence_rows": total_rows,
        "kept_fl_evidence_rows": kept_rows,
    }


def orient_matches(left, right):
    return left in "+-" and right in "+-" and left == right


def orient_relative_to_anchor(alignment_strand, anchor_strand):
    if alignment_strand not in "+-" or anchor_strand not in "+-":
        return "."
    return "+" if alignment_strand == anchor_strand else "-"


def orient_relative_to_path(alignment_strand, anchor_strand, path_orient):
    # Raw read-vs-reference strand is not comparable across inverted repeat copies.
    node_orient = orient_relative_to_anchor(alignment_strand, anchor_strand)
    if node_orient == "." or path_orient not in "+-":
        return "."
    return "+" if node_orient == path_orient else "-"


def anchor_observations_for_read(
    alignments,
    anchors,
    min_anchor_overlap,
    path_orient=None,
):
    observations = []
    for anchor in anchors:
        for alignment in alignments:
            if alignment.target_id != anchor["target_name"]:
                continue
            overlap_start = max(alignment.target_start, anchor["start"])
            overlap_end = min(alignment.target_end, anchor["end"])
            overlap = overlap_end - overlap_start + 1
            if overlap < min_anchor_overlap:
                continue
            target_span = max(1, alignment.target_end - alignment.target_start + 1)
            query_span = max(1, alignment.query_end - alignment.query_start + 1)
            scale = query_span / target_span
            if alignment.strand == "+":
                query_a = alignment.query_start + (overlap_start - alignment.target_start) * scale
                query_b = alignment.query_start + (overlap_end - alignment.target_start) * scale
            else:
                query_a = alignment.query_end - (overlap_end - alignment.target_start) * scale
                query_b = alignment.query_end - (overlap_start - alignment.target_start) * scale
            observations.append({
                "anchor_id": anchor["anchor_id"],
                "strand": alignment.strand,
                "anchor_strand": anchor.get("strand", "."),
                "path_strand": (
                    orient_relative_to_path(
                        alignment.strand,
                        anchor.get("strand", "."),
                        path_orient,
                    )
                    if path_orient is not None
                    else "."
                ),
                "query_mid": (query_a + query_b) / 2,
                "query_start": min(query_a, query_b),
                "query_end": max(query_a, query_b),
                "overlap": overlap,
                "alignment_role": alignment.alignment_role,
                "mapq": alignment.mapq,
            })
    return observations


def find_candidate_read_support(alignments, candidate, min_anchor_overlap):
    left_endpoint = _split_gfa_endpoint(candidate["left_endpoint"])
    right_endpoint = _split_gfa_endpoint(candidate["right_endpoint"])
    left_path_orient = left_endpoint[1] if left_endpoint is not None else None
    right_path_orient = (
        flip_orient(right_endpoint[1]) if right_endpoint is not None else None
    )
    left_obs = anchor_observations_for_read(
        alignments,
        candidate["left_anchors"],
        min_anchor_overlap,
        path_orient=left_path_orient,
    )
    repeat_obs = anchor_observations_for_read(
        alignments,
        candidate["repeat_anchors"],
        min_anchor_overlap,
        path_orient="+",
    )
    right_obs = anchor_observations_for_read(
        alignments,
        candidate["right_anchors"],
        min_anchor_overlap,
        path_orient=right_path_orient,
    )
    best = None
    for left in left_obs:
        for repeat in repeat_obs:
            if not orient_matches(left["path_strand"], repeat["path_strand"]):
                continue
            for right in right_obs:
                if not orient_matches(repeat["path_strand"], right["path_strand"]):
                    continue
                if left["query_mid"] < repeat["query_mid"] < right["query_mid"]:
                    if repeat["path_strand"] != "+":
                        continue
                    orientation = "left_to_right"
                elif right["query_mid"] < repeat["query_mid"] < left["query_mid"]:
                    if repeat["path_strand"] != "-":
                        continue
                    orientation = "right_to_left"
                else:
                    continue
                score = left["overlap"] + repeat["overlap"] + right["overlap"]
                query_span = abs(right["query_mid"] - left["query_mid"])
                item = {
                    "orientation": orientation,
                    "score": score,
                    "query_span": query_span,
                    "left_query_mid": left["query_mid"],
                    "repeat_query_mid": repeat["query_mid"],
                    "right_query_mid": right["query_mid"],
                    "left_anchor": left["anchor_id"],
                    "repeat_anchor": repeat["anchor_id"],
                    "right_anchor": right["anchor_id"],
                    "strand": repeat["strand"],
                    "path_strand": repeat["path_strand"],
                }
                if best is None or (item["score"], item["query_span"]) > (best["score"], best["query_span"]):
                    best = item
    return best


def compute_repeat_read_path_support(
    gfa,
    node_classes,
    hits_by_node,
    run_dir,
    verified_record,
    selected_hits_by_node=None,
    flank_anchor_length=500,
    flank_anchor_offset=1000,
    repeat_anchor_length=500,
    min_anchor_overlap=200,
):
    if selected_hits_by_node is None:
        selected_hits_by_node = {
            name: selected_mapping_hits_for_node(
                segment,
                hits_by_node.get(name, []),
                node_classes.get(name, "unknown"),
            )
            for name, segment in gfa.segments.items()
        }
    candidates = build_repeat_path_candidates(
        gfa,
        node_classes,
        selected_hits_by_node,
        str(verified_record.id),
        flank_anchor_length=flank_anchor_length,
        flank_anchor_offset=flank_anchor_offset,
        repeat_anchor_length=repeat_anchor_length,
        min_anchor_overlap=min_anchor_overlap,
    )
    alignments_by_read, evidence_info = load_fl_whole_read_evidence(run_dir)
    candidate_support_counts = defaultdict(int)
    candidate_best_support = {}
    assigned_read_count_by_repeat = defaultdict(int)
    ambiguous_read_count_by_repeat = defaultdict(int)
    candidate_read_count_by_repeat = defaultdict(int)
    assignment_rows = []

    candidates_by_repeat = defaultdict(list)
    for candidate in candidates:
        candidates_by_repeat[candidate["repeat_node"]].append(candidate)

    for repeat_node, repeat_candidates in candidates_by_repeat.items():
        for read_id, alignments in alignments_by_read.items():
            supported = []
            for candidate in repeat_candidates:
                support = find_candidate_read_support(alignments, candidate, min_anchor_overlap)
                if support is not None:
                    supported.append((candidate, support))
            if not supported:
                continue
            candidate_read_count_by_repeat[repeat_node] += 1
            if len(supported) == 1:
                candidate, support = supported[0]
                key = (repeat_node, candidate["path_id"])
                candidate_support_counts[key] += 1
                assigned_read_count_by_repeat[repeat_node] += 1
                previous = candidate_best_support.get(key)
                if previous is None or support["score"] > previous["score"]:
                    candidate_best_support[key] = support
                assignment_rows.append({
                    "repeat_node": repeat_node,
                    "read_id": read_id,
                    "assignment_status": "assigned",
                    "path_id": candidate["path_id"],
                    "supported_path_count": 1,
                    "support_orientation": support["orientation"],
                    "strand": support["strand"],
                    "left_query_mid": support["left_query_mid"],
                    "repeat_query_mid": support["repeat_query_mid"],
                    "right_query_mid": support["right_query_mid"],
                    "left_anchor": support["left_anchor"],
                    "repeat_anchor": support["repeat_anchor"],
                    "right_anchor": support["right_anchor"],
                    "support_score": support["score"],
                })
            else:
                ambiguous_read_count_by_repeat[repeat_node] += 1
                supported_path_ids = sorted({candidate["path_id"] for candidate, _ in supported})
                assignment_rows.append({
                    "repeat_node": repeat_node,
                    "read_id": read_id,
                    "assignment_status": "ambiguous",
                    "path_id": ",".join(supported_path_ids),
                    "supported_path_count": len(supported_path_ids),
                    "support_orientation": ".",
                    "strand": ".",
                    "left_query_mid": ".",
                    "repeat_query_mid": ".",
                    "right_query_mid": ".",
                    "left_anchor": ".",
                    "repeat_anchor": ".",
                    "right_anchor": ".",
                    "support_score": ".",
                })

    support_rows = []
    for repeat_node, repeat_candidates in candidates_by_repeat.items():
        unique_total = sum(
            candidate_support_counts[(repeat_node, candidate["path_id"])]
            for candidate in repeat_candidates
        )
        for candidate in repeat_candidates:
            key = (repeat_node, candidate["path_id"])
            count = candidate_support_counts[key]
            ratio = count / unique_total if unique_total else 0.0
            left_anchor_count = len(candidate["left_anchors"])
            right_anchor_count = len(candidate["right_anchors"])
            notes = "FL_ids_only;whole_read_evidence.tsv"
            if candidate["repeat_copy_count"] != 2:
                notes += ";ambiguous_repeat_copy_count"
            if left_anchor_count == 0 or right_anchor_count == 0 or candidate["repeat_copy_count"] == 0:
                notes += ";missing_anchor"
            support_rows.append({
                "repeat_node": repeat_node,
                "repeat_status": candidate["repeat_status"],
                "support_method": "remapped_FL_whole_read_evidence",
                "path_id": candidate["path_id"],
                "left_edge_index": candidate["left_edge_index"],
                "right_edge_index": candidate["right_edge_index"],
                "left_endpoint": candidate["left_endpoint"],
                "right_endpoint": candidate["right_endpoint"],
                "left_gfa_endpoint": candidate["left_gfa_endpoint"],
                "right_gfa_endpoint": candidate["right_gfa_endpoint"],
                "left_anchor_count": left_anchor_count,
                "repeat_copy_count": candidate["repeat_copy_count"],
                "right_anchor_count": right_anchor_count,
                "path_count": count,
                "path_ratio": ratio,
                "unique_assigned_read_count": unique_total,
                "ambiguous_read_count": ambiguous_read_count_by_repeat[repeat_node],
                "candidate_read_count": candidate_read_count_by_repeat[repeat_node],
                "total_fl_read_count": evidence_info["fl_read_count"],
                "kept_fl_evidence_rows": evidence_info["kept_fl_evidence_rows"],
                "notes": notes,
            })
    return support_rows, assignment_rows, evidence_info, selected_hits_by_node


def read_coverage_depths(path):
    path = Path(path)
    if not path.exists():
        return []
    depths = [0.0]
    with open(path) as fin:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            try:
                depths.append(float(fields[-1]))
            except ValueError:
                depths.append(0.0)
    return depths


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if end < start:
            continue
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def coverage_for_intervals(depths, intervals):
    if not depths or not intervals:
        return None
    max_position = len(depths) - 1
    total_depth = 0.0
    total_bases = 0
    for start, end in merge_intervals(intervals):
        start = max(1, start)
        end = min(max_position, end)
        if end < start:
            continue
        for value in depths[start:end + 1]:
            total_depth += value
        total_bases += end - start + 1
    if total_bases == 0:
        return None
    return {
        "depth_mean": total_depth / total_bases,
        "depth_bases": total_depth,
        "covered_bases": total_bases,
        "intervals": merge_intervals(intervals),
    }


def compute_node_remapped_coverage(gfa, node_classes, selected_hits_by_node, run_dir):
    cov_path = Path(run_dir) / "FL_cov.txt"
    depths = read_coverage_depths(cov_path)
    rows = []
    row_by_node = {}
    for name, segment in gfa.segments.items():
        hits = selected_hits_by_node.get(name, [])
        intervals = [(hit.target_start_1based, hit.target_end_1based) for hit in hits]
        coverage = coverage_for_intervals(depths, intervals)
        if coverage is None:
            row = {
                "node": name,
                "node_class": node_classes.get(name, "unknown"),
                "coverage_source": "remapped_FL_to_verified_fasta_unavailable",
                "depth_mean": ".",
                "depth_bases": ".",
                "covered_bases": ".",
                "interval_count": 0,
                "intervals": ".",
                "FL_cov_path": str(cov_path),
            }
        else:
            interval_text = ",".join(str(start) + "-" + str(end) for start, end in coverage["intervals"])
            row = {
                "node": name,
                "node_class": node_classes.get(name, "unknown"),
                "coverage_source": "remapped_FL_to_verified_fasta",
                "depth_mean": coverage["depth_mean"],
                "depth_bases": coverage["depth_bases"],
                "covered_bases": coverage["covered_bases"],
                "interval_count": len(coverage["intervals"]),
                "intervals": interval_text,
                "FL_cov_path": str(cov_path),
            }
        rows.append(row)
        row_by_node[name] = row
    return rows, row_by_node


def coverage_tags_from_remapped_row(row):
    if row is None or row.get("depth_mean") == ".":
        return ["CM:Z:remapped_FL_to_verified_fasta_unavailable"]
    return [
        "DP:f:" + f"{float(row['depth_mean']):.6f}",
        "CM:Z:remapped_FL_to_verified_fasta",
        "RI:i:" + str(row["interval_count"]),
        "RB:f:" + f"{float(row['depth_bases']):.3f}",
    ]


def make_synthetic_hit(node_name, node_length, target_name, target_length, query_start, query_end, target_start, target_end, strand):
    return PafHit(
        query_name=node_name,
        query_length=node_length,
        query_start=query_start - 1,
        query_end=query_end,
        strand=strand,
        target_name=target_name,
        target_length=target_length,
        target_start=target_start - 1,
        target_end=target_end,
        matches=query_end - query_start + 1,
        block_length=query_end - query_start + 1,
        mapq=60,
    )


def blast_interval_sequence(record, start, end):
    if start <= end:
        return extract_record_interval(record, start, end, "+"), "+", start, end
    return extract_record_interval(record, end, start, "-"), "-", end, start


def map_oriented_offset_to_record(start, end, strand, offset, length):
    if strand == "+":
        core_start = start + offset
        core_end = core_start + length - 1
        return core_start, core_end, "+"
    core_end = end - offset
    core_start = core_end - length + 1
    return core_start, core_end, "-"


def oriented_base_at(record, position, strand):
    base = str(record.seq[position - 1])
    if strand == "-":
        return reverse_complement(base)
    return base


def repeat_copy_extension_position(copy, side, reference_length):
    if side == "left":
        if copy["strand"] == "+":
            return copy["start"] - 1 if copy["start"] > 1 else None
        return copy["end"] + 1 if copy["end"] < reference_length else None
    if copy["strand"] == "+":
        return copy["end"] + 1 if copy["end"] < reference_length else None
    return copy["start"] - 1 if copy["start"] > 1 else None


def extended_repeat_copy(copy, position, side):
    updated = dict(copy)
    if side == "left":
        if updated["strand"] == "+":
            updated["start"] = position
        else:
            updated["end"] = position
    else:
        if updated["strand"] == "+":
            updated["end"] = position
        else:
            updated["start"] = position
    return updated


def repeat_copies_overlap(copy_a, copy_b):
    return intervals_overlap_size(
        copy_a["start"],
        copy_a["end"],
        copy_b["start"],
        copy_b["end"],
    ) > 0


def try_extend_repeat_core_pair(core_pair, verified_record, side):
    reference_length = len(verified_record.seq)
    copy_a, copy_b = core_pair["copies"]
    pos_a = repeat_copy_extension_position(copy_a, side, reference_length)
    pos_b = repeat_copy_extension_position(copy_b, side, reference_length)
    if pos_a is None or pos_b is None:
        return False
    base_a = oriented_base_at(verified_record, pos_a, copy_a["strand"])
    base_b = oriented_base_at(verified_record, pos_b, copy_b["strand"])
    if base_a.upper() != base_b.upper():
        return False
    updated_a = extended_repeat_copy(copy_a, pos_a, side)
    updated_b = extended_repeat_copy(copy_b, pos_b, side)
    if repeat_copies_overlap(updated_a, updated_b):
        return False
    core_pair["copies"] = [updated_a, updated_b]
    core_pair["length"] += 1
    return True


def extend_repeat_core_pair_to_max_exact(core_pair, verified_record):
    left_extensions = 0
    right_extensions = 0
    while try_extend_repeat_core_pair(core_pair, verified_record, "left"):
        left_extensions += 1
    while try_extend_repeat_core_pair(core_pair, verified_record, "right"):
        right_extensions += 1
    core_pair["boundary_extension_left"] = left_extensions
    core_pair["boundary_extension_right"] = right_extensions
    core_pair["boundary_extension_method"] = "terminal_exact_base_extension"
    return core_pair


def longest_exact_common_substring(seq_a, seq_b):
    previous = [0] * (len(seq_b) + 1)
    best_start_a = 0
    best_start_b = 0
    best_length = 0
    for index_a, char_a in enumerate(seq_a, start=1):
        current = [0] * (len(seq_b) + 1)
        for index_b, char_b in enumerate(seq_b, start=1):
            if char_a == char_b:
                current[index_b] = previous[index_b - 1] + 1
                if current[index_b] > best_length:
                    best_length = current[index_b]
                    best_start_a = index_a - best_length
                    best_start_b = index_b - best_length
        previous = current
    return best_start_a, best_start_b, best_length


def run_self_blastn_repeat_hits(verified_fasta, soft_paths_dict):
    blastn = soft_paths_dict.get("blastn", "blastn")
    command = [
        blastn,
        "-query",
        str(verified_fasta),
        "-subject",
        str(verified_fasta),
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-dust",
        "no",
        "-soft_masking",
        "false",
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(
            "self blastn failed with exit code "
            + str(completed.returncode)
            + "\n"
            + completed.stderr
        )
    rows = []
    for line in completed.stdout.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 12:
            continue
        row = {
            "qseqid": fields[0],
            "sseqid": fields[1],
            "pident": float(fields[2]),
            "length": int(fields[3]),
            "mismatch": int(fields[4]),
            "gapopen": int(fields[5]),
            "qstart": int(fields[6]),
            "qend": int(fields[7]),
            "sstart": int(fields[8]),
            "send": int(fields[9]),
            "evalue": fields[10],
            "bitscore": float(fields[11]),
        }
        q_low, q_high = sorted([row["qstart"], row["qend"]])
        s_low, s_high = sorted([row["sstart"], row["send"]])
        if q_low == s_low and q_high == s_high:
            continue
        rows.append(row)
    return rows, " ".join(command), completed.stderr


def core_pair_key(row):
    q_low, q_high = sorted([row["qstart"], row["qend"]])
    s_low, s_high = sorted([row["sstart"], row["send"]])
    return tuple(sorted([(q_low, q_high), (s_low, s_high)]))


def repeat_core_pair_from_blast_row(row, verified_record):
    q_seq, q_strand, q_low, q_high = blast_interval_sequence(
        verified_record,
        row["qstart"],
        row["qend"],
    )
    s_seq, s_strand, s_low, s_high = blast_interval_sequence(
        verified_record,
        row["sstart"],
        row["send"],
    )
    q_offset, s_offset, core_length = longest_exact_common_substring(q_seq, s_seq)
    q_core_start, q_core_end, q_core_strand = map_oriented_offset_to_record(
        q_low,
        q_high,
        q_strand,
        q_offset,
        core_length,
    )
    s_core_start, s_core_end, s_core_strand = map_oriented_offset_to_record(
        s_low,
        s_high,
        s_strand,
        s_offset,
        core_length,
    )
    core_pair = {
        "length": core_length,
        "source_blast_length": row["length"],
        "source_pident": row["pident"],
        "source_gapopen": row["gapopen"],
        "source_mismatch": row["mismatch"],
        "copies": [
            {
                "start": q_core_start,
                "end": q_core_end,
                "strand": q_core_strand,
            },
            {
                "start": s_core_start,
                "end": s_core_end,
                "strand": s_core_strand,
            },
        ],
    }
    return extend_repeat_core_pair_to_max_exact(core_pair, verified_record)


def intervals_overlap_size(start_a, end_a, start_b, end_b):
    return interval_overlap(start_a, end_a, start_b, end_b)


def hit_interval_overlap(hit, interval):
    return intervals_overlap_size(
        hit.target_start_1based,
        hit.target_end_1based,
        interval["start"],
        interval["end"],
    )


def infer_repeat_core_pairs(verified_fasta, verified_record, soft_paths_dict, repeat_count):
    rows, command, stderr = run_self_blastn_repeat_hits(verified_fasta, soft_paths_dict)
    rows = [
        row
        for row in rows
        if row["length"] >= 1000 and row["pident"] >= 95.0
    ]
    rows.sort(key=lambda row: (row["length"], row["pident"], row["bitscore"]), reverse=True)
    selected = []
    seen = set()
    used_intervals = []
    for row in rows:
        key = core_pair_key(row)
        if key in seen:
            continue
        seen.add(key)
        q_low, q_high = sorted([row["qstart"], row["qend"]])
        s_low, s_high = sorted([row["sstart"], row["send"]])
        if any(
            intervals_overlap_size(q_low, q_high, old_start, old_end) > 100
            or intervals_overlap_size(s_low, s_high, old_start, old_end) > 100
            for old_start, old_end in used_intervals
        ):
            continue
        core_pair = repeat_core_pair_from_blast_row(row, verified_record)
        if core_pair["length"] == 0:
            continue
        selected.append(core_pair)
        used_intervals.append((core_pair["copies"][0]["start"], core_pair["copies"][0]["end"]))
        used_intervals.append((core_pair["copies"][1]["start"], core_pair["copies"][1]["end"]))
        if len(selected) >= repeat_count:
            break
    return selected, command, stderr


def assign_repeat_nodes_to_core_pairs(gfa, node_classes, hits_by_node, core_pairs):
    assignments = {}
    used_pairs = set()
    for node_name, node_class in node_classes.items():
        if node_class != "repeat_node":
            continue
        best_pair_index = None
        best_score = -1
        for pair_index, pair in enumerate(core_pairs):
            if pair_index in used_pairs:
                continue
            score = 0
            for hit in hits_by_node.get(node_name, []):
                for copy in pair["copies"]:
                    score += hit_interval_overlap(hit, copy)
            if score > best_score:
                best_score = score
                best_pair_index = pair_index
        if best_pair_index is None:
            continue
        assignments[node_name] = core_pairs[best_pair_index]
        used_pairs.add(best_pair_index)
    return assignments


def circular_sequence_parts(record, parts):
    return "".join(extract_record_interval(record, start, end, "+") for start, end in parts)


def repeat_copy_sequence(record, copy, strand):
    return extract_record_interval(record, copy["start"], copy["end"], strand)


def choose_repeat_copy_for_node(node_name, repeat_pair, hits_by_node):
    best_copy = repeat_pair["copies"][0]
    best_strand = repeat_pair["copies"][0]["strand"]
    best_score = -1
    for hit in hits_by_node.get(node_name, []):
        for copy in repeat_pair["copies"]:
            score = hit_interval_overlap(hit, copy)
            if score > best_score:
                best_score = score
                best_copy = copy
                best_strand = hit.strand
    return best_copy, best_strand


def sorted_repeat_copies(repeat_assignments):
    copies = []
    for repeat_node, pair in repeat_assignments.items():
        for copy_index, copy in enumerate(pair["copies"], start=1):
            item = dict(copy)
            item["repeat_node"] = repeat_node
            item["copy_index"] = copy_index
            copies.append(item)
    return sorted(copies, key=lambda copy: copy["start"])


def repeat_gap_parts(left_copy, right_copy, reference_length):
    start = left_copy["end"] + 1
    end = right_copy["start"] - 1
    if start <= end:
        return [(start, end)]
    parts = []
    if start <= reference_length:
        parts.append((start, reference_length))
    if end >= 1:
        parts.append((1, end))
    return parts


def gap_length(parts):
    return sum(end - start + 1 for start, end in parts)


def build_repeat_core_gaps(repeat_assignments, reference_length):
    copies = sorted_repeat_copies(repeat_assignments)
    gaps = []
    if not copies:
        return gaps
    for index, left_copy in enumerate(copies):
        right_copy = copies[(index + 1) % len(copies)]
        parts = repeat_gap_parts(left_copy, right_copy, reference_length)
        gaps.append({
            "left_repeat_node": left_copy["repeat_node"],
            "left_copy_index": left_copy["copy_index"],
            "right_repeat_node": right_copy["repeat_node"],
            "right_copy_index": right_copy["copy_index"],
            "parts": parts,
            "length": gap_length(parts),
        })
    return gaps


def hit_overlap_with_parts(hit, parts):
    return sum(
        intervals_overlap_size(
            hit.target_start_1based,
            hit.target_end_1based,
            start,
            end,
        )
        for start, end in parts
    )


def assign_single_copy_nodes_to_gaps(gfa, node_classes, hits_by_node, gaps):
    assignments = {}
    used_gaps = set()
    for node_name, node_class in node_classes.items():
        if node_class != "single_copy_node":
            continue
        best_gap_index = None
        best_score = -1
        for gap_index, gap in enumerate(gaps):
            if gap_index in used_gaps:
                continue
            score = sum(hit_overlap_with_parts(hit, gap["parts"]) for hit in hits_by_node.get(node_name, []))
            if score > best_score:
                best_score = score
                best_gap_index = gap_index
        if best_gap_index is None:
            continue
        assignments[node_name] = gaps[best_gap_index]
        used_gaps.add(best_gap_index)
    return assignments


def best_orientation_for_parts(node_name, parts, hits_by_node):
    best_strand = "+"
    best_score = -1
    for hit in hits_by_node.get(node_name, []):
        score = hit_overlap_with_parts(hit, parts)
        if score > best_score:
            best_score = score
            best_strand = hit.strand
    return best_strand


def synthetic_hits_for_parts(node_name, node_length, record, parts, strand):
    hits = []
    query_position = 1
    ordered_parts = list(parts)
    if strand == "-":
        ordered_parts = list(reversed(ordered_parts))
    for start, end in ordered_parts:
        part_length = end - start + 1
        hits.append(
            make_synthetic_hit(
                node_name,
                node_length,
                record.id,
                len(record.seq),
                query_position,
                query_position + part_length - 1,
                start,
                end,
                strand,
            )
        )
        query_position += part_length
    return hits


def verified_sequence_from_parts(record, parts, strand):
    sequence = circular_sequence_parts(record, parts)
    if strand == "-":
        sequence = reverse_complement(sequence)
    return sequence


def synthetic_hits_for_repeat_copies(node_name, node_sequence, record, copies):
    hits = []
    node_length = len(node_sequence)
    for copy in copies:
        plus_sequence = extract_record_interval(record, copy["start"], copy["end"], "+")
        minus_sequence = reverse_complement(plus_sequence)
        if node_sequence == plus_sequence:
            strand = "+"
        elif node_sequence == minus_sequence:
            strand = "-"
        else:
            strand = copy.get("strand", "+")
        hits.append(
            make_synthetic_hit(
                node_name,
                node_length,
                record.id,
                len(record.seq),
                1,
                node_length,
                copy["start"],
                copy["end"],
                strand,
            )
        )
    return hits


def infer_verified_sequences_from_repeat_cores(merged_gfa, verified_record, verified_fasta, hits_by_node, node_classes, soft_paths_dict):
    repeat_node_count = sum(1 for node_class in node_classes.values() if node_class == "repeat_node")
    if repeat_node_count == 0:
        return None
    core_pairs, self_blastn_command, self_blastn_stderr = infer_repeat_core_pairs(
        verified_fasta,
        verified_record,
        soft_paths_dict,
        repeat_node_count,
    )
    if len(core_pairs) < repeat_node_count:
        return None
    repeat_assignments = assign_repeat_nodes_to_core_pairs(
        merged_gfa,
        node_classes,
        hits_by_node,
        core_pairs,
    )
    if len(repeat_assignments) < repeat_node_count:
        return None
    gaps = build_repeat_core_gaps(repeat_assignments, len(verified_record.seq))
    single_assignments = assign_single_copy_nodes_to_gaps(
        merged_gfa,
        node_classes,
        hits_by_node,
        gaps,
    )

    sequences = {}
    synthetic_hits = defaultdict(list)
    rows = []
    for node_name, repeat_pair in repeat_assignments.items():
        copy, strand = choose_repeat_copy_for_node(node_name, repeat_pair, hits_by_node)
        sequence = repeat_copy_sequence(verified_record, copy, strand)
        sequences[node_name] = sequence
        synthetic_hits[node_name] = synthetic_hits_for_repeat_copies(
            node_name,
            sequence,
            verified_record,
            repeat_pair["copies"],
        )
        rows.append({
            "node": node_name,
            "node_class": "repeat_node",
            "source": "verified_fasta_self_blast_repeat_core",
            "length": len(sequence),
            "parts": str(copy["start"]) + "-" + str(copy["end"]) + ":" + strand,
            "notes": "core_copies="
            + ",".join(
                str(item["start"]) + "-" + str(item["end"]) + ":" + item["strand"]
                for item in repeat_pair["copies"]
            )
            + ";blast_length="
            + str(repeat_pair["source_blast_length"])
            + ";blast_pident="
            + f"{repeat_pair['source_pident']:.3f}"
            + ";boundary_extension_left="
            + str(repeat_pair.get("boundary_extension_left", 0))
            + ";boundary_extension_right="
            + str(repeat_pair.get("boundary_extension_right", 0)),
        })

    for node_name, gap in single_assignments.items():
        parts = gap["parts"]
        strand = best_orientation_for_parts(node_name, parts, hits_by_node)
        sequence = verified_sequence_from_parts(verified_record, parts, strand)
        sequences[node_name] = sequence
        synthetic_hits[node_name] = synthetic_hits_for_parts(
            node_name,
            len(sequence),
            verified_record,
            parts,
            strand,
        )
        rows.append({
            "node": node_name,
            "node_class": "single_copy_node",
            "source": "verified_fasta_between_repeat_cores",
            "length": len(sequence),
            "parts": ",".join(str(start) + "-" + str(end) for start, end in parts) + ":" + strand,
            "notes": "left="
            + gap["left_repeat_node"]
            + ":copy"
            + str(gap["left_copy_index"])
            + ";right="
            + gap["right_repeat_node"]
            + ":copy"
            + str(gap["right_copy_index"]),
        })

    for name, segment in merged_gfa.segments.items():
        if name in sequences:
            continue
        hit = best_hit(hits_by_node.get(name, []))
        if hit is None:
            sequences[name] = segment.sequence
            continue
        sequence = extract_verified_sequence(verified_record, hit)
        sequences[name] = sequence
        synthetic_hits[name] = [hit]
        rows.append({
            "node": name,
            "node_class": node_classes.get(name, "unknown"),
            "source": "verified_fasta_minimap2_fallback",
            "length": len(sequence),
            "parts": str(hit.target_start_1based) + "-" + str(hit.target_end_1based) + ":" + hit.strand,
            "notes": "fallback",
        })

    info = {
        "self_blastn_command": self_blastn_command,
        "self_blastn_stderr": self_blastn_stderr,
        "repeat_core_pairs": core_pairs,
    }
    return sequences, synthetic_hits, rows, info


def write_verified_node_source_map(path, rows):
    header = ["node", "node_class", "source", "length", "parts", "notes"]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in sorted(rows, key=lambda item: natural_key(item["node"])):
            writer.writerow(row)


def make_verified_gfa_from_sequences(merged_gfa, verified_sequences_by_node, selected_hits_by_node, node_classes):
    verified_segments = OrderedDict()
    for name, segment in merged_gfa.segments.items():
        sequence = verified_sequences_by_node.get(name, segment.sequence)
        hits = selected_hits_by_node.get(name, [])
        interval_text = ",".join(
            str(hit.target_start_1based) + "-" + str(hit.target_end_1based) + ":" + hit.strand
            for hit in hits
        )
        tags = replace_tags(segment.tags, [
            "VC:Z:verified_fasta_repeat_core_guided",
            "VZ:Z:verified_fasta",
            "VS:Z:" + (interval_text if interval_text else "."),
            "NC:Z:" + node_classes.get(name, "unknown"),
            "LN:i:" + str(len(sequence)),
            "OL:i:" + str(len(segment.sequence)),
            "VL:i:" + str(len(sequence)),
        ])
        verified_segments[name] = Segment(name, sequence, tags)
    return Gfa(
        headers=merged_gfa.headers,
        segments=verified_segments,
        links=merged_gfa.links,
        other_lines=merged_gfa.other_lines,
    )


def annotate_gfa_with_read_evidence(gfa, node_classes, coverage_by_node, repeat_rows):
    annotated_segments = OrderedDict()
    for name, segment in gfa.segments.items():
        tags = replace_tags(segment.tags, coverage_tags_from_remapped_row(coverage_by_node.get(name)))
        tags = replace_tags(tags, ["NC:Z:" + node_classes.get(name, "unknown")])
        if node_classes.get(name) == "repeat_node":
            tags = replace_tags(tags, repeat_path_tags(name, repeat_rows))
        annotated_segments[name] = Segment(name, segment.sequence, tags)
    return Gfa(gfa.headers, annotated_segments, gfa.links, gfa.other_lines)


def classify_nodes(gfa, hits_by_node):
    median_depth = median_depth_for_gfa(gfa)
    repeats = repeat_like_nodes(gfa)
    node_classes = {}
    for name in gfa.segments:
        if name in repeats:
            node_classes[name] = "repeat_node"
        else:
            node_classes[name] = "single_copy_node"
    return node_classes, median_depth


def load_first_fasta_record(path):
    return next(SeqIO.parse(path, "fasta"))


def extract_verified_sequence(verified_record, hit):
    sequence = str(verified_record.seq[hit.target_start:hit.target_end])
    if hit.strand == "-":
        sequence = reverse_complement(sequence)
    return sequence


def extract_record_interval(record, start, end, strand="+"):
    sequence = str(record.seq[start - 1:end])
    if strand == "-":
        sequence = reverse_complement(sequence)
    return sequence


def make_verified_gfa(merged_gfa, verified_record, hits_by_node, node_classes):
    verified_segments = OrderedDict()
    for name, segment in merged_gfa.segments.items():
        hits = hits_by_node.get(name, [])
        hit = best_hit(hits)
        confidence, score, reason = confidence_for_hits(hits)
        if hit is None:
            sequence = segment.sequence
            tags = replace_tags(segment.tags, [
                "VC:Z:unmapped",
                "VR:Z:no_alignment",
                "NC:Z:" + node_classes.get(name, "unknown"),
                "LN:i:" + str(len(sequence)),
                "OL:i:" + str(len(segment.sequence)),
                "VL:i:" + str(len(sequence)),
            ])
        else:
            if confidence == "low":
                sequence = segment.sequence
                sequence_source = "raw_due_to_low_confidence"
            else:
                sequence = extract_verified_sequence(verified_record, hit)
                sequence_source = "verified_fasta"
            tags = replace_tags(segment.tags, [
                "VC:Z:" + confidence,
                "VS:i:" + str(hit.target_start_1based),
                "VE:i:" + str(hit.target_end_1based),
                "VI:f:" + f"{hit.identity:.6f}",
                "VF:f:" + f"{hit.aligned_fraction:.6f}",
                "VQ:i:" + str(hit.mapq),
                "VR:Z:" + reason,
                "VZ:Z:" + sequence_source,
                "NC:Z:" + node_classes.get(name, "unknown"),
                "LN:i:" + str(len(sequence)),
                "OL:i:" + str(len(segment.sequence)),
                "VL:i:" + str(len(sequence)),
            ])
        verified_segments[name] = Segment(name, sequence, tags)
    return Gfa(
        headers=merged_gfa.headers,
        segments=verified_segments,
        links=merged_gfa.links,
        other_lines=merged_gfa.other_lines,
    )


def annotate_merged_gfa_evidence(merged_gfa, raw_gfa, node_classes):
    rows = repeat_path_support_rows(merged_gfa, node_classes)
    annotated_segments = OrderedDict()
    for name, segment in merged_gfa.segments.items():
        tags = replace_tags(segment.tags, coverage_tags_for_segment(name, segment, raw_gfa))
        tags = replace_tags(tags, ["NC:Z:" + node_classes.get(name, "unknown")])
        if node_classes.get(name) == "repeat_node":
            tags = replace_tags(tags, repeat_path_tags(name, rows))
        annotated_segments[name] = Segment(name, segment.sequence, tags)
    return Gfa(
        merged_gfa.headers,
        annotated_segments,
        merged_gfa.links,
        merged_gfa.other_lines,
    ), rows


def write_repeat_path_support(path, rows):
    header = [
        "repeat_node",
        "repeat_status",
        "support_method",
        "left_edge_index",
        "right_edge_index",
        "left_endpoint",
        "right_endpoint",
        "left_support",
        "right_support",
        "path_support",
        "path_ratio",
    ]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            for key in ["left_support", "right_support", "path_support", "path_ratio"]:
                if isinstance(formatted.get(key), float):
                    formatted[key] = f"{formatted[key]:.6f}"
            writer.writerow(formatted)


def write_repeat_read_path_support(path, rows):
    header = [
        "repeat_node",
        "repeat_status",
        "support_method",
        "path_id",
        "left_edge_index",
        "right_edge_index",
        "left_endpoint",
        "right_endpoint",
        "left_gfa_endpoint",
        "right_gfa_endpoint",
        "left_anchor_count",
        "repeat_copy_count",
        "right_anchor_count",
        "path_count",
        "path_ratio",
        "unique_assigned_read_count",
        "ambiguous_read_count",
        "candidate_read_count",
        "total_fl_read_count",
        "kept_fl_evidence_rows",
        "notes",
    ]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            if isinstance(formatted.get("path_ratio"), float):
                formatted["path_ratio"] = f"{formatted['path_ratio']:.6f}"
            writer.writerow(formatted)


def write_repeat_read_path_assignments(path, rows):
    header = [
        "repeat_node",
        "read_id",
        "assignment_status",
        "path_id",
        "supported_path_count",
        "support_orientation",
        "strand",
        "left_query_mid",
        "repeat_query_mid",
        "right_query_mid",
        "left_anchor",
        "repeat_anchor",
        "right_anchor",
        "support_score",
    ]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            for key in ["left_query_mid", "repeat_query_mid", "right_query_mid", "support_score"]:
                if isinstance(formatted.get(key), float):
                    formatted[key] = f"{formatted[key]:.3f}"
            writer.writerow(formatted)


def write_node_remapped_coverage(path, rows):
    header = [
        "node",
        "node_class",
        "coverage_source",
        "depth_mean",
        "depth_bases",
        "covered_bases",
        "interval_count",
        "intervals",
        "FL_cov_path",
    ]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        for row in rows:
            formatted = dict(row)
            for key in ["depth_mean", "depth_bases"]:
                if isinstance(formatted.get(key), float):
                    formatted[key] = f"{formatted[key]:.6f}"
            writer.writerow(formatted)


def merge_unambiguous_gfa(raw_gfa):
    """Compact linear, non-branching paths without resolving repeat branches.

    This is a conservative fallback. When a curated merged GFA is provided, use
    that instead; the fallback removes branching nodes from the graph, merges
    only the remaining linear components, and keeps branching nodes separate.
    """
    protected_nodes = repeat_like_nodes(raw_gfa)
    adjacency = defaultdict(list)
    for link in raw_gfa.links:
        if link.from_name in protected_nodes or link.to_name in protected_nodes:
            continue
        adjacency[link.from_name].append(link.to_name)
        adjacency[link.to_name].append(link.from_name)

    visited = set()
    components = []
    for name in raw_gfa.segments:
        if name in visited or name in protected_nodes:
            continue
        stack = [name]
        component = []
        visited.add(name)
        while stack:
            current = stack.pop()
            component.append(current)
            for nxt in adjacency[current]:
                if nxt not in visited:
                    visited.add(nxt)
                    stack.append(nxt)
        components.append(component)

    def ordered_linear_component(component):
        component_set = set(component)
        sub_degrees = {
            name: sum(1 for item in adjacency[name] if item in component_set)
            for name in component
        }
        if any(value > 2 for value in sub_degrees.values()):
            return None
        endpoints = [name for name, value in sub_degrees.items() if value <= 1]
        if len(component) > 1 and len(endpoints) != 2:
            return None
        start = sorted(endpoints or component, key=natural_key)[0]
        order = []
        previous = None
        current = start
        while current is not None:
            order.append(current)
            candidates = [
                item
                for item in adjacency[current]
                if item in component_set and item != previous
            ]
            if not candidates:
                break
            previous, current = current, sorted(candidates, key=natural_key)[0]
        if len(order) != len(component):
            return None
        return order

    links_by_pair = defaultdict(list)
    for link in raw_gfa.links:
        key = frozenset([link.from_name, link.to_name])
        links_by_pair[key].append(link)

    def next_orient(current, current_orient, nxt):
        for link in links_by_pair[frozenset([current, nxt])]:
            if (
                link.from_name == current
                and link.to_name == nxt
                and link.from_orient == current_orient
            ):
                return link.to_orient
            if (
                link.to_name == current
                and link.from_name == nxt
                and flip_orient(link.to_orient) == current_orient
            ):
                return flip_orient(link.from_orient)
        return None

    def orient_path(order):
        for start_orient in ["+", "-"]:
            orientations = [start_orient]
            ok = True
            for current, nxt in zip(order, order[1:]):
                orient = next_orient(current, orientations[-1], nxt)
                if orient is None:
                    ok = False
                    break
                orientations.append(orient)
            if ok:
                return orientations
        return ["+"] * len(order)

    path_records = []
    used_nodes = set()
    for component in components:
        order = ordered_linear_component(component)
        if order is None:
            for name in sorted(component, key=natural_key):
                path_records.append(([name], ["+"]))
                used_nodes.add(name)
            continue
        orientations = orient_path(order)
        path_records.append((order, orientations))
        used_nodes.update(order)
    for name in raw_gfa.segments:
        if name not in used_nodes:
            path_records.append(([name], ["+"]))
            used_nodes.add(name)

    old_to_new = {}
    old_to_path_info = {}
    merged_segments = OrderedDict()
    mergeable_path_count = 0
    for order, orientations in path_records:
        if len(order) > 1:
            mergeable_path_count += 1
        name = "_".join(order)
        sequence_parts = []
        for old_name, orient in zip(order, orientations):
            raw_sequence = raw_gfa.segments[old_name].sequence
            if orient == "-":
                raw_sequence = reverse_complement(raw_sequence)
            sequence_parts.append(raw_sequence)
            old_to_new[old_name] = name
            old_to_path_info[old_name] = {
                "path": order,
                "orientations": orientations,
                "index": len(sequence_parts) - 1,
            }
        sequence = "".join(sequence_parts)
        tags = [
            "LN:i:" + str(len(sequence)),
            "SC:Z:linear_compaction" if len(order) > 1 else "SC:Z:preserved_node",
            "RR:Z:disabled",
        ]
        merged_segments[name] = Segment(name, sequence, tags)

    def convert_endpoint(old_name, orient):
        info = old_to_path_info[old_name]
        path = info["path"]
        orientations = info["orientations"]
        index = info["index"]
        path_orient = orientations[index]
        if len(path) == 1:
            return "+" if orient == path_orient else "-"
        if index == 0 and orient == flip_orient(path_orient):
            return "-"
        if index == len(path) - 1 and orient == path_orient:
            return "+"
        return None

    def convert_link(link):
        options = [
            (link.from_name, link.from_orient, link.to_name, link.to_orient),
            (link.to_name, flip_orient(link.to_orient), link.from_name, flip_orient(link.from_orient)),
        ]
        for from_name, from_orient, to_name, to_orient in options:
            merged_from = old_to_new[from_name]
            merged_to = old_to_new[to_name]
            if merged_from == merged_to:
                continue
            converted_from = convert_endpoint(from_name, from_orient)
            converted_to = convert_endpoint(to_name, to_orient)
            if converted_from is not None and converted_to is not None:
                return Link(merged_from, converted_from, merged_to, converted_to, link.overlap, link.tags)
        return None

    merged_links = []
    seen_links = set()
    for link in raw_gfa.links:
        converted = convert_link(link)
        if converted is None:
            continue
        key = (
            converted.from_name,
            converted.from_orient,
            converted.to_name,
            converted.to_orient,
            converted.overlap,
        )
        if key in seen_links:
            continue
        seen_links.add(key)
        merged_links.append(converted)

    merge_mode = "auto_non_repeat_linear_compaction"
    if mergeable_path_count == 0:
        merge_mode = "input_already_merged_or_no_linear_compaction"
    return Gfa(raw_gfa.headers, merged_segments, merged_links, raw_gfa.other_lines), merge_mode


def read_soft_paths(path):
    soft_paths = {}
    with open(path) as fin:
        for line in fin:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t", 1)
            if len(fields) == 2:
                soft_paths[fields[0]] = fields[1]
    return soft_paths


def count_lines(path):
    if path is None or not Path(path).exists():
        return 0
    if Path(path).suffix in {".xlsx", ".xls"}:
        return len(pd.read_excel(path))
    with open(path) as fin:
        return sum(1 for _ in fin)


def fasta_length(path):
    if path is None or not Path(path).exists():
        return 0
    return sum(len(record.seq) for record in SeqIO.parse(path, "fasta"))


def coverage_stats(path):
    if path is None or not Path(path).exists():
        return {}
    values = []
    with open(path) as fin:
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            try:
                values.append(float(fields[-1]))
            except ValueError:
                continue
    if not values:
        return {}
    values_sorted = sorted(values)
    middle = len(values_sorted) // 2
    median = values_sorted[middle]
    if len(values_sorted) % 2 == 0:
        median = (values_sorted[middle - 1] + values_sorted[middle]) / 2
    return {
        "mean": sum(values) / len(values),
        "median": median,
        "min": min(values),
        "max": max(values),
    }


def write_coordinate_map(path, merged_gfa, hits_by_node, node_classes, max_link_gap=1000):
    degrees = graph_degrees(merged_gfa)
    best_by_node = {name: best_hit(hits_by_node.get(name, [])) for name in merged_gfa.segments}
    header = [
        "row_type",
        "verified_fasta",
        "verified_start",
        "verified_end",
        "gfa_node",
        "node_start",
        "node_end",
        "strand",
        "node_length",
        "node_class",
        "edge_link_candidate",
        "mapping_confidence",
        "confidence_score",
        "identity",
        "aligned_fraction",
        "mapq",
        "degree",
        "notes",
    ]
    rows = [header]
    for name, segment in merged_gfa.segments.items():
        hits = hits_by_node.get(name, [])
        hit = best_hit(hits)
        confidence, score, reason = confidence_for_hits(hits)
        if hit is None:
            rows.append([
                "segment",
                ".",
                ".",
                ".",
                name,
                ".",
                ".",
                ".",
                str(len(segment.sequence)),
                node_classes.get(name, "unknown"),
                "node",
                confidence,
                f"{score:.6f}",
                ".",
                ".",
                ".",
                str(degrees[name]),
                reason,
            ])
        else:
            rows.append([
                "segment",
                hit.target_name,
                str(hit.target_start_1based),
                str(hit.target_end_1based),
                name,
                str(hit.query_start_1based),
                str(hit.query_end_1based),
                hit.strand,
                str(len(segment.sequence)),
                node_classes.get(name, "unknown"),
                "node",
                confidence,
                f"{score:.6f}",
                f"{hit.identity:.6f}",
                f"{hit.aligned_fraction:.6f}",
                str(hit.mapq),
                str(degrees[name]),
                reason,
            ])

    for index, link in enumerate(merged_gfa.links, start=1):
        left = best_by_node.get(link.from_name)
        right = best_by_node.get(link.to_name)
        link_name = (
            link.from_name
            + link.from_orient
            + "->"
            + link.to_name
            + link.to_orient
        )
        if left is None or right is None or left.target_name != right.target_name:
            confidence = "low"
            score = 0.0
            verified_start = "."
            verified_end = "."
            notes = "missing_or_different_contig_node_mapping"
            target_name = "."
        else:
            gap = max(0, max(left.target_start, right.target_start) - min(left.target_end, right.target_end))
            verified_start = str(min(left.target_start_1based, right.target_start_1based))
            verified_end = str(max(left.target_end_1based, right.target_end_1based))
            target_name = left.target_name
            if gap <= max_link_gap:
                confidence = "medium"
                score = 1.0 / (1.0 + gap)
                notes = "coordinate_adjacent_or_overlapping;gap=" + str(gap)
            else:
                confidence = "low"
                score = 0.0
                notes = "not_coordinate_adjacent;gap=" + str(gap)
        rows.append([
            "link_candidate",
            target_name,
            verified_start,
            verified_end,
            link_name,
            ".",
            ".",
            ".",
            ".",
            "link",
            "link_candidate",
            confidence,
            f"{score:.6f}",
            ".",
            ".",
            ".",
            ".",
            "link_index=" + str(index) + ";" + notes,
        ])

    with open(path, "wt", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerows(rows)


def write_before_after_statistics(path, raw_gfa, merged_gfa, verified_gfa, verified_fasta, run_dir):
    run_dir = Path(run_dir) if run_dir is not None else None
    stats = OrderedDict()
    stats["raw_segment_count"] = len(raw_gfa.segments)
    stats["raw_link_count"] = len(raw_gfa.links)
    stats["raw_segment_total_length"] = sum(len(segment.sequence) for segment in raw_gfa.segments.values())
    stats["merged_segment_count"] = len(merged_gfa.segments)
    stats["merged_link_count"] = len(merged_gfa.links)
    stats["merged_segment_total_length"] = sum(len(segment.sequence) for segment in merged_gfa.segments.values())
    stats["verified_gfa_segment_count"] = len(verified_gfa.segments)
    stats["verified_gfa_link_count"] = len(verified_gfa.links)
    stats["verified_gfa_segment_total_length"] = sum(len(segment.sequence) for segment in verified_gfa.segments.values())
    stats["verified_fasta_total_length"] = fasta_length(verified_fasta)
    if run_dir is not None:
        stats["run_filtered_variant_count"] = count_lines(run_dir / "variants_anno_combined_depth_frq_filter.xlsx")
        stats["run_all_bcftools_call_count"] = count_lines(run_dir / "all_bcftools_calls.txt")
        stats["fl_read_count"] = count_lines(run_dir / "FL_ids.txt")
        stats["partial_read_count"] = count_lines(run_dir / "partial_ids.txt")
        stats["variant_cov_read_count"] = count_lines(run_dir / "variant_cov_ids.txt")
        cov = coverage_stats(run_dir / "variant_cov.txt")
        for key, value in cov.items():
            stats["variant_cov_" + key] = value
    with open(path, "wt", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["metric", "value"])
        for key, value in stats.items():
            if isinstance(value, float):
                writer.writerow([key, f"{value:.6f}"])
            else:
                writer.writerow([key, value])


def write_node_sequences(path, verified_gfa):
    with open(path, "wt") as fout:
        for segment in verified_gfa.segments.values():
            print(">" + segment.name, file=fout)
            print(segment.sequence, file=fout)


def write_build_manifest(path, values):
    with open(path, "wt", newline="") as fout:
        writer = csv.writer(fout, delimiter="\t")
        writer.writerow(["key", "value"])
        for key, value in values.items():
            writer.writerow([key, value])


def resolve_gfa_editor_cli(soft_paths_dict, explicit_path=None):
    candidates = []
    if explicit_path is not None:
        candidates.append(explicit_path)
    soft_path = soft_paths_dict.get("gfa_editor_cli")
    if soft_path:
        candidates.append(soft_path)
    candidates.append("/Users/zouyinstein-m4max/Documents/Codex/GFA_Editor/scripts/gfa_editor_cli.py")
    for candidate in candidates:
        if candidate is None:
            continue
        candidate_path = Path(candidate)
        if candidate_path.exists():
            return str(candidate_path)
    return None


def command_for_python_script(path, soft_paths_dict):
    path = str(path)
    if path.endswith(".py"):
        script_path = Path(path)
        local_python = script_path.resolve().parents[1] / ".venv" / "bin" / "python"
        if local_python.exists():
            return [str(local_python), path]
        if not os.access(path, os.X_OK):
            return [soft_paths_dict.get("python", "python"), path]
    return [path]


def export_verified_gfa_pdfs(output_dir, image_reference_fasta, soft_paths_dict):
    output_dir = Path(output_dir)
    image_reference_fasta = Path(image_reference_fasta)
    previous_cwd = Path.cwd()
    try:
        os.chdir(output_dir)
        pdf_outputs = hfrps.get_gfa_reference_images(
            str(image_reference_fasta),
            soft_paths_dict,
            output_format="pdf",
            renderer="gfa_editor",
        )
        svg_outputs = hfrps.get_gfa_reference_images(
            str(image_reference_fasta),
            soft_paths_dict,
            output_format="svg",
            renderer="gfa_editor",
            append_protocol=True,
        )
        return pdf_outputs + svg_outputs
    finally:
        os.chdir(previous_cwd)


def gfa_segment_total_length(path):
    if not Path(path).exists():
        return None
    gfa = parse_gfa(path)
    return sum(len(segment.sequence) for segment in gfa.segments.values())


def write_auto_repeat_check(path, row):
    header = [
        "status",
        "message",
        "reference_fasta",
        "verified_gfa",
        "resolved_gfa",
        "summary_json",
        "history_json",
        "gfa_editor_cli",
        "returncode",
        "selected_candidate",
        "best_candidate",
        "method",
        "score",
        "orientation",
        "length_delta",
        "continuous_bp",
        "reference_length",
        "resolved_segment_total_length",
        "command",
        "stdout",
        "stderr",
    ]
    with open(path, "wt", newline="") as fout:
        writer = csv.DictWriter(fout, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerow({key: row.get(key, ".") for key in header})


def run_gfa_editor_auto_repeat_check(
    genome,
    verified_gfa_path,
    verified_fasta,
    output_dir,
    soft_paths_dict,
    gfa_editor_cli=None,
):
    resolved_gfa = output_dir / ("verified_" + genome + ".auto_repeat_resolved.gfa")
    summary_json = output_dir / ("verified_" + genome + ".auto_repeat_summary.json")
    history_json = output_dir / ("verified_" + genome + ".auto_repeat_history.json")
    check_tsv = output_dir / ("verified_" + genome + ".auto_repeat_check.tsv")
    cli = resolve_gfa_editor_cli(soft_paths_dict, gfa_editor_cli)
    reference_records = list(SeqIO.parse(verified_fasta, "fasta"))
    reference_length = len(reference_records[0].seq) if len(reference_records) == 1 else None

    base_row = {
        "reference_fasta": str(verified_fasta),
        "verified_gfa": str(verified_gfa_path),
        "resolved_gfa": str(resolved_gfa),
        "summary_json": str(summary_json),
        "history_json": str(history_json),
        "gfa_editor_cli": "." if cli is None else str(cli),
        "reference_length": "." if reference_length is None else str(reference_length),
    }
    if cli is None:
        row = dict(base_row)
        row.update({
            "status": "SKIP",
            "message": "gfa_editor_cli_not_found",
            "returncode": ".",
        })
        write_auto_repeat_check(check_tsv, row)
        return {
            "auto_repeat_resolved_gfa": resolved_gfa,
            "auto_repeat_summary_json": summary_json,
            "auto_repeat_history_json": history_json,
            "auto_repeat_check": check_tsv,
            "auto_repeat_check_status": "SKIP",
            "auto_repeat_check_message": "gfa_editor_cli_not_found",
        }

    command = command_for_python_script(cli, soft_paths_dict) + [
        "auto-repeat",
        str(verified_gfa_path),
        str(resolved_gfa),
        "--reference-fasta",
        str(verified_fasta),
        "--summary-json",
        str(summary_json),
        "--history-json",
        str(history_json),
    ]
    completed = subprocess.run(command, capture_output=True, text=True)
    summary = {}
    if summary_json.exists():
        with open(summary_json) as fin:
            summary = json.load(fin)
    reference_selection = summary.get("reference_selection", {})
    best = reference_selection.get("best", {}) if isinstance(reference_selection, dict) else {}
    selected = summary.get("selected", {}) if isinstance(summary.get("selected", {}), dict) else {}
    resolved_length = gfa_segment_total_length(resolved_gfa)
    method = best.get("method", ".")
    score = best.get("score", ".")
    length_delta = best.get("lengthDelta", ".")
    continuous_bp = best.get("continuousBp", ".")
    selected_candidate = selected.get("id", ".")
    best_candidate = best.get("candidate", ".")
    try:
        score_ok = abs(float(score) - 1.0) < 1e-9
    except (TypeError, ValueError):
        score_ok = False
    try:
        length_delta_ok = int(length_delta) == 0
    except (TypeError, ValueError):
        length_delta_ok = False
    try:
        continuous_ok = reference_length is not None and int(continuous_bp) == reference_length
    except (TypeError, ValueError):
        continuous_ok = False
    length_ok = reference_length is not None and resolved_length == reference_length
    selected_ok = selected_candidate == best_candidate and selected_candidate != "."
    status = "PASS" if (
        completed.returncode == 0
        and method == "sequence-exact-circular"
        and score_ok
        and length_delta_ok
        and continuous_ok
        and length_ok
        and selected_ok
    ) else "FAIL"
    message = "resolved_gfa_matches_verified_fasta_circular_exact" if status == "PASS" else "resolved_gfa_reference_check_failed"
    row = dict(base_row)
    row.update({
        "status": status,
        "message": message,
        "returncode": str(completed.returncode),
        "selected_candidate": str(selected_candidate),
        "best_candidate": str(best_candidate),
        "method": str(method),
        "score": str(score),
        "orientation": str(best.get("orientation", ".")),
        "length_delta": str(length_delta),
        "continuous_bp": str(continuous_bp),
        "resolved_segment_total_length": "." if resolved_length is None else str(resolved_length),
        "command": " ".join(command),
        "stdout": completed.stdout.replace("\n", "\\n"),
        "stderr": completed.stderr.replace("\n", "\\n"),
    })
    write_auto_repeat_check(check_tsv, row)
    if status == "FAIL":
        raise RuntimeError(
            "verified GFA auto-repeat consistency check failed; see "
            + str(check_tsv)
        )
    return {
        "auto_repeat_resolved_gfa": resolved_gfa,
        "auto_repeat_summary_json": summary_json,
        "auto_repeat_history_json": history_json,
        "auto_repeat_check": check_tsv,
        "auto_repeat_check_status": status,
        "auto_repeat_check_message": message,
    }


def build_verified_gfa(
    genome,
    raw_gfa_path,
    verified_fasta_path,
    run_dir,
    output_dir,
    soft_paths_dict,
    threads=1,
    merged_gfa_template=None,
    max_link_gap=1000,
    flank_anchor_length=500,
    flank_anchor_offset=1000,
    repeat_anchor_length=500,
    min_anchor_overlap=200,
    gfa_editor_cli=None,
    image_reference_fasta=None,
):
    genome = str(genome)
    raw_gfa_path = Path(raw_gfa_path)
    verified_fasta_path = Path(verified_fasta_path)
    if image_reference_fasta is None:
        image_reference_fasta = verified_fasta_path
    else:
        image_reference_fasta = Path(image_reference_fasta)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    unmerged_raw = output_dir / ("unmerged_" + genome + "_raw.gfa")
    merged_raw = output_dir / ("merged_" + genome + "_raw.gfa")
    verified_fasta = output_dir / ("verified_" + genome + ".fasta")
    graph_coordinate_map = output_dir / "graph_coordinate_map.tsv"
    before_after_stats = output_dir / "before_after_graph_statistics.tsv"
    repeat_path_support = output_dir / "repeat_path_support.tsv"
    repeat_path_support_proxy = output_dir / "repeat_path_support_proxy.tsv"
    repeat_path_assignments = output_dir / "repeat_path_read_assignments.tsv"
    node_remapped_coverage = output_dir / "node_coverage_remapped.tsv"
    verified_node_source_map = output_dir / "verified_node_source_map.tsv"
    verified_gfa_path = output_dir / ("verified_" + genome + ".gfa")
    verified_nodes_fasta = output_dir / ("verified_" + genome + "_nodes.fasta")
    manifest = output_dir / "verified_gfa_build_manifest.tsv"

    shutil.copyfile(raw_gfa_path, unmerged_raw)
    shutil.copyfile(verified_fasta_path, verified_fasta)

    raw_gfa = parse_gfa(unmerged_raw)
    if merged_gfa_template is not None:
        shutil.copyfile(merged_gfa_template, merged_raw)
        merge_mode = "template"
    else:
        merged_gfa_auto, merge_mode = merge_unambiguous_gfa(raw_gfa)
        write_gfa(merged_gfa_auto, merged_raw)

    merged_gfa = parse_gfa(merged_raw)
    pre_mapping_node_classes, _ = classify_nodes(merged_gfa, {})
    proxy_repeat_rows = repeat_path_support_rows(merged_gfa, pre_mapping_node_classes)
    write_repeat_path_support(repeat_path_support_proxy, proxy_repeat_rows)
    verified_record = load_first_fasta_record(verified_fasta)

    with tempfile.TemporaryDirectory(prefix="hifisr_verified_gfa_") as tmpdir:
        node_fasta = Path(tmpdir) / "merged_nodes.fasta"
        write_segment_fasta(merged_gfa.segments, node_fasta)
        hits_by_node, minimap2_command, minimap2_stderr = run_minimap2_paf(
            node_fasta,
            verified_fasta,
            soft_paths_dict,
            threads,
        )

    node_classes, median_depth = classify_nodes(merged_gfa, hits_by_node)
    repeat_core_info = infer_verified_sequences_from_repeat_cores(
        merged_gfa,
        verified_record,
        verified_fasta,
        hits_by_node,
        node_classes,
        soft_paths_dict,
    )
    if repeat_core_info is None:
        verified_sequences_by_node = {
            name: (
                extract_verified_sequence(verified_record, best_hit(hits_by_node.get(name, [])))
                if best_hit(hits_by_node.get(name, [])) is not None
                else segment.sequence
            )
            for name, segment in merged_gfa.segments.items()
        }
        selected_hits_by_node = {
            name: selected_mapping_hits_for_node(
                segment,
                hits_by_node.get(name, []),
                node_classes.get(name, "unknown"),
            )
            for name, segment in merged_gfa.segments.items()
        }
        verified_node_rows = []
        self_blastn_command = "."
        self_blastn_stderr = "."
        sequence_projection_method = "minimap2_fallback"
    else:
        (
            verified_sequences_by_node,
            selected_hits_by_node,
            verified_node_rows,
            self_blastn_info,
        ) = repeat_core_info
        self_blastn_command = self_blastn_info["self_blastn_command"]
        self_blastn_stderr = self_blastn_info["self_blastn_stderr"]
        sequence_projection_method = "verified_fasta_self_blast_repeat_core"
    write_verified_node_source_map(verified_node_source_map, verified_node_rows)

    sequence_gfa = make_verified_gfa_from_sequences(
        merged_gfa,
        verified_sequences_by_node,
        selected_hits_by_node,
        node_classes,
    )
    repeat_rows, assignment_rows, evidence_info, selected_hits_by_node = compute_repeat_read_path_support(
        sequence_gfa,
        node_classes,
        hits_by_node,
        run_dir,
        verified_record,
        selected_hits_by_node=selected_hits_by_node,
        flank_anchor_length=flank_anchor_length,
        flank_anchor_offset=flank_anchor_offset,
        repeat_anchor_length=repeat_anchor_length,
        min_anchor_overlap=min_anchor_overlap,
    )
    coverage_rows, coverage_by_node = compute_node_remapped_coverage(
        sequence_gfa,
        node_classes,
        selected_hits_by_node,
        run_dir,
    )
    merged_gfa = annotate_gfa_with_read_evidence(
        merged_gfa,
        node_classes,
        coverage_by_node,
        repeat_rows,
    )
    merged_gfa = add_repeat_path_gfa_lines(merged_gfa, repeat_rows)
    write_gfa(merged_gfa, merged_raw)
    write_repeat_read_path_support(repeat_path_support, repeat_rows)
    write_repeat_read_path_assignments(repeat_path_assignments, assignment_rows)
    write_node_remapped_coverage(node_remapped_coverage, coverage_rows)
    verified_gfa = annotate_gfa_with_read_evidence(
        sequence_gfa,
        node_classes,
        coverage_by_node,
        repeat_rows,
    )
    verified_gfa = add_repeat_path_gfa_lines(verified_gfa, repeat_rows)
    write_gfa(verified_gfa, verified_gfa_path)
    write_node_sequences(verified_nodes_fasta, verified_gfa)
    write_coordinate_map(graph_coordinate_map, verified_gfa, selected_hits_by_node, node_classes, max_link_gap=max_link_gap)
    write_before_after_statistics(before_after_stats, raw_gfa, merged_gfa, verified_gfa, verified_fasta, run_dir)
    auto_repeat_check_info = run_gfa_editor_auto_repeat_check(
        genome,
        verified_gfa_path,
        verified_fasta,
        output_dir,
        soft_paths_dict,
        gfa_editor_cli=gfa_editor_cli,
    )
    gfa_image_outputs = export_verified_gfa_pdfs(
        output_dir,
        image_reference_fasta,
        soft_paths_dict,
    )
    gfa_image_protocol = output_dir / "gfa_image_export_protocol.tsv"
    write_build_manifest(
        manifest,
        OrderedDict([
            ("genome", genome),
            ("raw_gfa", str(raw_gfa_path)),
            ("unmerged_raw_gfa", str(unmerged_raw)),
            ("merged_raw_gfa", str(merged_raw)),
            ("verified_fasta_input", str(verified_fasta_path)),
            ("gfa_image_reference_fasta", str(image_reference_fasta)),
            ("verified_fasta_copy", str(verified_fasta)),
            ("run_dir", str(run_dir)),
            ("merge_mode", merge_mode),
            ("repeat_resolution", "disabled"),
            ("sequence_projection_method", sequence_projection_method),
            ("coverage_source", "remapped_FL_to_verified_fasta"),
            ("repeat_path_support_method", "remapped_FL_whole_read_evidence"),
            ("repeat_path_support_FL_ids", str(evidence_info["fl_ids_path"])),
            ("repeat_path_support_whole_read_evidence", str(evidence_info["whole_read_evidence_path"])),
            ("repeat_path_support_total_FL_reads", str(evidence_info["fl_read_count"])),
            ("repeat_path_support_kept_FL_evidence_rows", str(evidence_info["kept_fl_evidence_rows"])),
            ("flank_anchor_length", str(flank_anchor_length)),
            ("flank_anchor_offset", str(flank_anchor_offset)),
            ("repeat_anchor_length", str(repeat_anchor_length)),
            ("min_anchor_overlap", str(min_anchor_overlap)),
            ("median_segment_depth", "." if median_depth is None else f"{median_depth:.6f}"),
            ("minimap2_command", minimap2_command),
            ("minimap2_stderr", minimap2_stderr.replace("\n", "\\n")),
            ("self_blastn_command", self_blastn_command),
            ("self_blastn_stderr", self_blastn_stderr.replace("\n", "\\n")),
            ("gfa_editor_auto_repeat_check_status", auto_repeat_check_info["auto_repeat_check_status"]),
            ("gfa_editor_auto_repeat_check_message", auto_repeat_check_info["auto_repeat_check_message"]),
            ("gfa_editor_auto_repeat_check", str(auto_repeat_check_info["auto_repeat_check"])),
            ("gfa_editor_auto_repeat_resolved_gfa", str(auto_repeat_check_info["auto_repeat_resolved_gfa"])),
            ("gfa_editor_auto_repeat_summary_json", str(auto_repeat_check_info["auto_repeat_summary_json"])),
            ("gfa_editor_auto_repeat_history_json", str(auto_repeat_check_info["auto_repeat_history_json"])),
            ("gfa_image_export_protocol", str(gfa_image_protocol)),
            ("gfa_image_outputs", ",".join(str(Path(output_dir) / output) for output in gfa_image_outputs)),
        ]),
    )
    return {
        "unmerged_raw_gfa": unmerged_raw,
        "merged_raw_gfa": merged_raw,
        "verified_fasta": verified_fasta,
        "verified_gfa": verified_gfa_path,
        "verified_nodes_fasta": verified_nodes_fasta,
        "graph_coordinate_map": graph_coordinate_map,
        "before_after_graph_statistics": before_after_stats,
        "repeat_path_support": repeat_path_support,
        "repeat_path_support_proxy": repeat_path_support_proxy,
        "repeat_path_read_assignments": repeat_path_assignments,
        "node_coverage_remapped": node_remapped_coverage,
        "verified_node_source_map": verified_node_source_map,
        "auto_repeat_resolved_gfa": auto_repeat_check_info["auto_repeat_resolved_gfa"],
        "auto_repeat_summary_json": auto_repeat_check_info["auto_repeat_summary_json"],
        "auto_repeat_history_json": auto_repeat_check_info["auto_repeat_history_json"],
        "auto_repeat_check": auto_repeat_check_info["auto_repeat_check"],
        "gfa_image_export_protocol": gfa_image_protocol,
        "gfa_image_outputs": ",".join(str(Path(output_dir) / output) for output in gfa_image_outputs),
        "manifest": manifest,
    }
