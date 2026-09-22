#!/usr/bin/env python3
"""Check DAS/Rucio disk staging for CMS datasets."""

import argparse
import csv
import datetime as dt
import json
import os
import subprocess
import sys
from collections import defaultdict


DEFAULT_DATASETS = [
    "/Cosmics/pORun2025-PromptReco-v1/AOD",
    "/Cosmics/OORun2025-PromptReco-v1/AOD",
    "/Cosmics/NeNeRun2025-PromptReco-v1/AOD",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Query DAS and Rucio to summarize which RSEs have disk replicas "
            "for each dataset. Run from an lxplus/CMSSW environment with a CMS proxy."
        )
    )
    parser.add_argument("datasets", nargs="*", default=DEFAULT_DATASETS)
    parser.add_argument(
        "-o",
        "--output",
        default="CosmicMuography/data/dataset_staging_report",
        help="Output prefix; .json and .md will be written.",
    )
    parser.add_argument(
        "--top-sites",
        type=int,
        default=12,
        help="Number of disk RSEs to show per dataset in the Markdown summary.",
    )
    parser.add_argument(
        "--timeout",
        type=int,
        default=90,
        help="Timeout in seconds for each DAS/Rucio subprocess call.",
    )
    return parser.parse_args()


def run_command(command, timeout):
    try:
        completed = subprocess.run(
            command,
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired as exc:
        return {
            "ok": False,
            "returncode": None,
            "stdout": exc.stdout or "",
            "stderr": f"timeout after {timeout}s",
            "command": command,
        }
    return {
        "ok": completed.returncode == 0,
        "returncode": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
        "command": command,
    }


def das_query(query, timeout, json_mode=False, limit=0):
    command = ["dasgoclient", f"--query={query}", f"--limit={limit}"]
    if json_mode:
        command.append("--json")
    result = run_command(command, timeout)
    if not result["ok"]:
        raise RuntimeError(
            "DAS query failed: {query}\n{stderr}".format(
                query=query,
                stderr=result["stderr"].strip() or result["stdout"].strip(),
            )
        )
    return result["stdout"]


def rucio_dataset_csv(did, timeout):
    result = run_command(["rucio", "replica", "list", "dataset", "--csv", did], timeout)
    if not result["ok"]:
        raise RuntimeError(
            "Rucio replica query failed: {did}\n{stderr}".format(
                did=did,
                stderr=result["stderr"].strip() or result["stdout"].strip(),
            )
        )
    return result["stdout"]


def load_summary(dataset, timeout):
    text = das_query(f"summary dataset={dataset} instance=prod/global", timeout, True)
    payload = json.loads(text)
    if not payload:
        return {}
    summaries = payload[0].get("summary", [])
    return summaries[0] if summaries else {}


def load_block_summary(block, timeout):
    text = das_query(f"summary block={block} instance=prod/global", timeout, True)
    payload = json.loads(text)
    if not payload:
        return {}
    summaries = payload[0].get("summary", [])
    return summaries[0] if summaries else {}


def load_lines(query, timeout):
    text = das_query(query, timeout, False)
    return [line.strip() for line in text.splitlines() if line.strip()]


def parse_rucio_csv(text):
    rows = []
    for row in csv.reader(text.splitlines()):
        if len(row) < 3:
            continue
        try:
            found = int(row[1])
            total = int(row[2])
        except ValueError:
            continue
        rows.append({"rse": row[0], "found": found, "total": total})
    return rows


def is_tape_rse(rse):
    return "tape" in rse.lower()


def summarize_dataset(dataset, timeout):
    print(f"dataset {dataset}", flush=True)
    summary = load_summary(dataset, timeout)
    blocks = load_lines(f"block dataset={dataset} instance=prod/global", timeout)
    das_sites = load_lines(f"site dataset={dataset} instance=prod/global", timeout)

    rse_totals = defaultdict(
        lambda: {
            "found_files": 0,
            "reported_total_files": 0,
            "blocks_with_any": 0,
            "complete_blocks": 0,
            "rows": 0,
            "is_tape": False,
        }
    )
    block_records = []

    for iblock, block in enumerate(blocks, start=1):
        print(f"  block {iblock}/{len(blocks)}", flush=True)
        block_summary = load_block_summary(block, timeout)
        block_nfiles = int(block_summary.get("nfiles", block_summary.get("num_file", 0)) or 0)
        rows = parse_rucio_csv(rucio_dataset_csv(f"cms:{block}", timeout))
        positive_disk_rses = []
        positive_tape_rses = []
        complete_disk_rses = []
        disk_replica_file_copies = 0
        for row in rows:
            rse = row["rse"]
            total = rse_totals[rse]
            total["found_files"] += row["found"]
            total["reported_total_files"] += row["total"]
            total["rows"] += 1
            total["is_tape"] = total["is_tape"] or is_tape_rse(rse)
            if row["found"] > 0:
                total["blocks_with_any"] += 1
                if row["found"] >= block_nfiles and block_nfiles > 0:
                    total["complete_blocks"] += 1
                if is_tape_rse(rse):
                    positive_tape_rses.append(rse)
                else:
                    disk_replica_file_copies += row["found"]
                    positive_disk_rses.append(rse)
                    if row["found"] >= block_nfiles and block_nfiles > 0:
                        complete_disk_rses.append(rse)
        block_records.append(
            {
                "name": block,
                "nfiles": block_nfiles,
                "positive_disk_rses": sorted(positive_disk_rses),
                "positive_tape_rses": sorted(positive_tape_rses),
                "complete_disk_rses": sorted(complete_disk_rses),
                "disk_replica_file_copies": disk_replica_file_copies,
            }
        )

    nfiles = int(summary.get("nfiles", summary.get("num_file", 0)) or 0)
    for rse, info in rse_totals.items():
        info["coverage_fraction"] = (float(info["found_files"]) / nfiles) if nfiles else 0.0

    disk_rses = {
        rse: info
        for rse, info in rse_totals.items()
        if not info["is_tape"] and info["found_files"] > 0
    }
    tape_rses = {
        rse: info
        for rse, info in rse_totals.items()
        if info["is_tape"] and info["found_files"] > 0
    }
    full_disk_rses = sorted(
        rse for rse, info in disk_rses.items() if nfiles and info["found_files"] >= nfiles
    )
    blocks_with_disk = sum(1 for block in block_records if block["positive_disk_rses"])
    blocks_complete_on_disk = sum(1 for block in block_records if block["complete_disk_rses"])
    files_in_blocks_complete_on_disk = sum(
        block["nfiles"] for block in block_records if block["complete_disk_rses"]
    )
    disk_replica_file_copies = sum(info["found_files"] for info in disk_rses.values())
    best_disk_rse = ""
    best_disk_files = 0
    if disk_rses:
        best_disk_rse, best_disk_info = max(
            disk_rses.items(), key=lambda item: (item[1]["found_files"], item[0])
        )
        best_disk_files = best_disk_info["found_files"]

    return {
        "dataset": dataset,
        "summary": summary,
        "das_sites": das_sites,
        "das_disk_sites": [site for site in das_sites if "Tape" not in site],
        "das_tape_sites": [site for site in das_sites if "Tape" in site],
        "blocks": block_records,
        "nblocks": len(block_records),
        "nfiles": nfiles,
        "nevents": int(summary.get("nevents", summary.get("num_event", 0)) or 0),
        "file_size_bytes": int(summary.get("file_size", 0) or 0),
        "rse_totals": dict(sorted(rse_totals.items())),
        "disk_rses": dict(sorted(disk_rses.items())),
        "tape_rses": dict(sorted(tape_rses.items())),
        "full_disk_rses": full_disk_rses,
        "blocks_with_disk": blocks_with_disk,
        "blocks_complete_on_disk": blocks_complete_on_disk,
        "files_in_blocks_complete_on_disk": files_in_blocks_complete_on_disk,
        "disk_replica_file_copies": disk_replica_file_copies,
        "best_disk_rse": best_disk_rse,
        "best_disk_files": best_disk_files,
    }


def fmt_bytes(value):
    value = float(value or 0.0)
    for unit in ("B", "KB", "MB", "GB", "TB", "PB"):
        if abs(value) < 1024.0 or unit == "PB":
            return f"{value:.2f} {unit}"
        value /= 1024.0
    return f"{value:.2f} PB"


def top_rses(items, limit):
    return sorted(
        items.items(),
        key=lambda item: (item[1]["found_files"], item[1]["complete_blocks"], item[0]),
        reverse=True,
    )[:limit]


def write_markdown(path, report, top_sites):
    lines = []
    lines.append("# CMS Cosmic AOD Dataset Staging Report")
    lines.append("")
    lines.append(f"Generated: `{report['generated_at']}`")
    lines.append("")
    lines.append(
        "Disk staging is derived from per-block Rucio replica counts. "
        "`found_files / files` is summed over blocks for each RSE; tape RSEs are shown separately."
    )
    lines.append("")
    for item in report["datasets"]:
        lines.append(f"## `{item['dataset']}`")
        lines.append("")
        lines.append(
            "- DAS summary: {files} files, {blocks} blocks, {events} events, {size}".format(
                files=item["nfiles"],
                blocks=item["nblocks"],
                events=item["nevents"],
                size=fmt_bytes(item["file_size_bytes"]),
            )
        )
        lines.append(
            "- Blocks with at least one disk replica: {blocks}/{total}".format(
                blocks=item["blocks_with_disk"],
                total=item["nblocks"],
            )
        )
        lines.append(
            "- Blocks complete on at least one disk RSE: {blocks}/{total} ({files} files in complete disk blocks)".format(
                blocks=item["blocks_complete_on_disk"],
                total=item["nblocks"],
                files=item["files_in_blocks_complete_on_disk"],
            )
        )
        if item["best_disk_rse"]:
            lines.append(
                "- Disk file-copies observed: {copies}; best single disk RSE: `{rse}` has {files}/{total} files ({coverage:.1%})".format(
                    copies=item["disk_replica_file_copies"],
                    rse=item["best_disk_rse"],
                    files=item["best_disk_files"],
                    total=item["nfiles"],
                    coverage=(float(item["best_disk_files"]) / item["nfiles"]) if item["nfiles"] else 0.0,
                )
            )
        else:
            lines.append("- Disk file-copies observed: 0")
        if item["full_disk_rses"]:
            lines.append("- Full dataset on disk at: " + ", ".join(f"`{rse}`" for rse in item["full_disk_rses"]))
        else:
            lines.append("- Full dataset on a single disk RSE: none found")
        lines.append("- DAS disk sites: " + ", ".join(f"`{site}`" for site in item["das_disk_sites"]))
        lines.append("")
        lines.append("| Disk RSE | found files | coverage | blocks with any | complete blocks |")
        lines.append("|---|---:|---:|---:|---:|")
        disk_rows = top_rses(item["disk_rses"], top_sites)
        if disk_rows:
            for rse, info in disk_rows:
                lines.append(
                    "| `{rse}` | {found} | {coverage:.1%} | {any_blocks} | {complete_blocks} |".format(
                        rse=rse,
                        found=info["found_files"],
                        coverage=info["coverage_fraction"],
                        any_blocks=info["blocks_with_any"],
                        complete_blocks=info["complete_blocks"],
                    )
                )
        else:
            lines.append("| none | 0 | 0.0% | 0 | 0 |")
        lines.append("")
        if item["tape_rses"]:
            lines.append("| Tape RSE | found files | coverage |")
            lines.append("|---|---:|---:|")
            for rse, info in top_rses(item["tape_rses"], top_sites):
                lines.append(
                    "| `{rse}` | {found} | {coverage:.1%} |".format(
                        rse=rse,
                        found=info["found_files"],
                        coverage=info["coverage_fraction"],
                    )
                )
            lines.append("")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def main():
    args = parse_args()
    report = {
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "datasets": [],
    }
    for dataset in args.datasets:
        report["datasets"].append(summarize_dataset(dataset, args.timeout))

    output_prefix = args.output
    output_dir = os.path.dirname(os.path.abspath(output_prefix))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    json_path = output_prefix + ".json"
    md_path = output_prefix + ".md"
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    write_markdown(md_path, report, args.top_sites)
    print(f"wrote {json_path}")
    print(f"wrote {md_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
