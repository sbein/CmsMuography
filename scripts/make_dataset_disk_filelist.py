#!/usr/bin/env python3
"""Build file lists for dataset files with available disk replicas."""

import argparse
import csv
import datetime as dt
import json
import os
import subprocess
import sys
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Use DAS to list dataset files and Rucio to find available disk "
            "replicas with ROOT PFNs. Writes one-PFN-per-file and all-PFN lists."
        )
    )
    parser.add_argument(
        "dataset",
        nargs="?",
        default="/Cosmics/OORun2025-PromptReco-v1/AOD",
        help="CMS dataset name.",
    )
    parser.add_argument(
        "-o",
        "--output-prefix",
        default="CosmicMuography/data/OORun2025_AOD_disk_root",
        help="Output prefix for .json/.txt/.tsv products.",
    )
    parser.add_argument(
        "--protocols",
        default="root",
        help="Rucio protocols filter, e.g. root or davs. Default is ROOT/XRootD.",
    )
    parser.add_argument(
        "--domain",
        choices=("wan", "lan", "all"),
        default="wan",
        help="Rucio PFN domain.",
    )
    parser.add_argument(
        "--prefer-rses",
        default="T2_CH_CERN,T2_CH_CSCS,T1_ES_PIC_Disk,T2_PL_Cyfronet,T2_FI_HIP,T2_BE_UCL,T2_RU_JINR",
        help="Comma-separated RSE preference order for the one-PFN-per-file list.",
    )
    parser.add_argument("--max-files", type=int, default=-1, help="Limit DAS files for tests.")
    parser.add_argument("--timeout", type=int, default=90, help="Subprocess timeout in seconds.")
    parser.add_argument("--progress-every", type=int, default=10)
    return parser.parse_args()


def run_command(command, timeout):
    completed = subprocess.run(
        command,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        timeout=timeout,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            "command failed ({code}): {cmd}\n{err}".format(
                code=completed.returncode,
                cmd=" ".join(command),
                err=completed.stderr.strip() or completed.stdout.strip(),
            )
        )
    return completed.stdout


def das_files(dataset, timeout, max_files):
    limit = max_files if max_files and max_files > 0 else 0
    output = run_command(
        ["dasgoclient", f"--query=file dataset={dataset} instance=prod/global", f"--limit={limit}"],
        timeout,
    )
    return [line.strip() for line in output.splitlines() if line.startswith("/store/")]


def parse_replica_table(text):
    replicas = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line.startswith("|"):
            continue
        columns = [column.strip() for column in line.strip("|").split("|")]
        if len(columns) < 5 or columns[0] == "SCOPE" or columns[0].startswith("-"):
            continue
        rse_replica = columns[4]
        if ": " not in rse_replica:
            continue
        rse, pfn = rse_replica.split(": ", 1)
        replicas.append({"rse": rse.strip(), "pfn": pfn.strip()})
    return replicas


def rucio_replicas(lfn, protocols, domain, timeout):
    did = "cms:" + lfn
    output = run_command(
        ["rucio", "replica", "list", "file", "--protocols", protocols, "--domain", domain, did],
        timeout,
    )
    return parse_replica_table(output)


def is_tape_rse(rse):
    return "tape" in rse.lower()


def choose_replica(replicas, preferred_rses):
    if not replicas:
        return None
    by_rse = {replica["rse"]: replica for replica in replicas}
    for rse in preferred_rses:
        if rse in by_rse:
            return by_rse[rse]
    return replicas[0]


def write_lines(path, lines):
    with open(path, "w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


def main():
    args = parse_args()
    preferred_rses = [rse.strip() for rse in args.prefer_rses.split(",") if rse.strip()]

    files = das_files(args.dataset, args.timeout, args.max_files)
    records = []
    chosen_pfns = []
    staged_lfns = []
    all_replica_rows = []
    rse_counts = Counter()

    for index, lfn in enumerate(files, start=1):
        if args.progress_every > 0 and (index == 1 or index % args.progress_every == 0):
            print(f"file {index}/{len(files)}", flush=True)
        replicas = rucio_replicas(lfn, args.protocols, args.domain, args.timeout)
        disk_replicas = [replica for replica in replicas if not is_tape_rse(replica["rse"])]
        chosen = choose_replica(disk_replicas, preferred_rses)
        for replica in disk_replicas:
            rse_counts[replica["rse"]] += 1
            all_replica_rows.append(
                {
                    "lfn": lfn,
                    "rse": replica["rse"],
                    "pfn": replica["pfn"],
                }
            )
        if chosen:
            staged_lfns.append(lfn)
            chosen_pfns.append(chosen["pfn"])
        records.append(
            {
                "lfn": lfn,
                "disk_replicas": disk_replicas,
                "chosen_replica": chosen,
            }
        )

    output_prefix = args.output_prefix
    output_dir = os.path.dirname(os.path.abspath(output_prefix))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    report = {
        "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
        "dataset": args.dataset,
        "protocols": args.protocols,
        "domain": args.domain,
        "n_dataset_files": len(files),
        "n_files_with_disk_replicas": len(staged_lfns),
        "n_disk_replica_pfns": len(all_replica_rows),
        "preferred_rses": preferred_rses,
        "rse_file_counts": dict(sorted(rse_counts.items())),
        "records": records,
    }

    json_path = output_prefix + ".json"
    one_pfn_path = output_prefix + "_one_pfn_per_file.txt"
    lfn_path = output_prefix + "_lfns.txt"
    tsv_path = output_prefix + "_all_pfns.tsv"
    summary_path = output_prefix + "_summary.txt"

    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    write_lines(one_pfn_path, chosen_pfns)
    write_lines(lfn_path, staged_lfns)
    with open(tsv_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["lfn", "rse", "pfn"], delimiter="\t")
        writer.writeheader()
        writer.writerows(all_replica_rows)
    with open(summary_path, "w", encoding="utf-8") as handle:
        handle.write(f"dataset: {args.dataset}\n")
        handle.write(f"dataset files: {len(files)}\n")
        handle.write(f"files with disk {args.protocols} replicas: {len(staged_lfns)}\n")
        handle.write(f"disk replica PFNs: {len(all_replica_rows)}\n")
        handle.write("RSE counts:\n")
        for rse, count in rse_counts.most_common():
            handle.write(f"  {rse}: {count}\n")

    print(f"dataset files: {len(files)}")
    print(f"files with disk {args.protocols} replicas: {len(staged_lfns)}")
    print(f"disk replica PFNs: {len(all_replica_rows)}")
    print(f"wrote {json_path}")
    print(f"wrote {one_pfn_path}")
    print(f"wrote {lfn_path}")
    print(f"wrote {tsv_path}")
    print(f"wrote {summary_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
