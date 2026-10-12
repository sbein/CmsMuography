#!/usr/bin/env python3
"""Validate ROOT/XRootD PFNs with xrdfs stat."""

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
from urllib.parse import urlparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run xrdfs stat over a text file of root:// PFNs and write the readable subset."
    )
    parser.add_argument("input", help="Text file containing one root:// PFN per line.")
    parser.add_argument(
        "-o",
        "--output-prefix",
        default=None,
        help="Output prefix; defaults to input path without .txt plus _xrdfs.",
    )
    parser.add_argument("--timeout", type=int, default=20, help="timeout seconds per PFN.")
    parser.add_argument("--progress-every", type=int, default=10)
    return parser.parse_args()


def read_pfns(path):
    pfns = []
    with open(path, encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line and not line.startswith("#"):
                pfns.append(line)
    return pfns


def split_root_url(url):
    parsed = urlparse(url)
    if parsed.scheme != "root" or not parsed.netloc:
        raise ValueError(f"not a root:// PFN: {url}")
    return f"root://{parsed.netloc}", parsed.path


def xrdfs_stat(url, timeout):
    endpoint, path = split_root_url(url)
    completed = subprocess.run(
        ["timeout", str(timeout), "xrdfs", endpoint, "stat", path],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return {
        "ok": completed.returncode == 0,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
        "endpoint": endpoint,
        "path": path,
    }


def write_lines(path, lines):
    with open(path, "w", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line + "\n")


def main():
    args = parse_args()
    output_prefix = args.output_prefix
    if output_prefix is None:
        base = args.input[:-4] if args.input.endswith(".txt") else args.input
        output_prefix = base + "_xrdfs"

    pfns = read_pfns(args.input)
    readable = []
    records = []
    for index, pfn in enumerate(pfns, start=1):
        if args.progress_every > 0 and (index == 1 or index % args.progress_every == 0):
            print(f"pfn {index}/{len(pfns)}", flush=True)
        try:
            result = xrdfs_stat(pfn, args.timeout)
        except Exception as exc:
            result = {
                "ok": False,
                "returncode": None,
                "stdout": "",
                "stderr": str(exc),
                "endpoint": "",
                "path": "",
            }
        records.append({"pfn": pfn, **result})
        if result["ok"]:
            readable.append(pfn)

    output_dir = os.path.dirname(os.path.abspath(output_prefix))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    readable_path = output_prefix + "_readable.txt"
    failures_path = output_prefix + "_failures.tsv"
    json_path = output_prefix + ".json"

    write_lines(readable_path, readable)
    with open(failures_path, "w", encoding="utf-8") as handle:
        handle.write("pfn\treturncode\terror\n")
        for record in records:
            if not record["ok"]:
                error = (record["stderr"] or record["stdout"]).replace("\t", " ").replace("\n", " ")
                handle.write(f"{record['pfn']}\t{record['returncode']}\t{error}\n")
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "generated_at": dt.datetime.now(dt.timezone.utc).isoformat(),
                "input": args.input,
                "n_input_pfns": len(pfns),
                "n_readable_pfns": len(readable),
                "timeout": args.timeout,
                "records": records,
            },
            handle,
            indent=2,
            sort_keys=True,
        )

    print(f"input PFNs: {len(pfns)}")
    print(f"xrdfs-readable PFNs: {len(readable)}")
    print(f"wrote {readable_path}")
    print(f"wrote {failures_path}")
    print(f"wrote {json_path}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
