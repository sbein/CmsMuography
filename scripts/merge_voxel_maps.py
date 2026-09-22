#!/usr/bin/env python3
"""Merge per-file cosmic voxel-map ROOT outputs."""

import argparse
import json
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge TH3/cutflow/tracks outputs from cosmic_voxel_map.py."
    )
    parser.add_argument("inputs", nargs="*", help="Input ROOT files.")
    parser.add_argument("--input-list", help="Text file containing input ROOT files.")
    parser.add_argument("-o", "--output", required=True, help="Merged output ROOT file.")
    parser.add_argument("--hist-name", default="cavern_muon_voxels")
    parser.add_argument("--tree-name", default="tracks")
    parser.add_argument("--skip-bad-files", action="store_true")
    return parser.parse_args()


def expand_inputs(args):
    inputs = list(args.inputs)
    if args.input_list:
        with open(args.input_list, encoding="utf-8") as handle:
            inputs.extend(
                line.strip()
                for line in handle
                if line.strip() and not line.lstrip().startswith("#")
            )
    return inputs


def main():
    args = parse_args()
    inputs = expand_inputs(args)
    if not inputs:
        raise SystemExit("No input files were provided.")

    import ROOT

    good_inputs = []
    merged_hist = None
    merged_cutflow = None
    configs = []
    total_hist_entries = 0.0

    for path in inputs:
        root_file = ROOT.TFile.Open(path)
        if not root_file or root_file.IsZombie():
            if args.skip_bad_files:
                print(f"Skipping bad ROOT file: {path}")
                continue
            raise RuntimeError(f"Could not open {path}")
        hist = root_file.Get(args.hist_name)
        if not hist:
            if args.skip_bad_files:
                print(f"Skipping ROOT file without {args.hist_name}: {path}")
                root_file.Close()
                continue
            raise RuntimeError(f"Could not find {args.hist_name} in {path}")
        cutflow = root_file.Get("cutflow")
        if merged_hist is None:
            merged_hist = hist.Clone(args.hist_name)
            merged_hist.SetDirectory(0)
            merged_hist.Reset("ICES")
        merged_hist.Add(hist)
        total_hist_entries += float(hist.GetEntries())
        if cutflow:
            if merged_cutflow is None:
                merged_cutflow = cutflow.Clone("cutflow")
                merged_cutflow.SetDirectory(0)
                merged_cutflow.Reset("ICES")
            merged_cutflow.Add(cutflow)
        config = root_file.Get("config_json")
        if config:
            try:
                configs.append(json.loads(str(config.GetString().Data())))
            except Exception:
                configs.append({"unparsed_config_from": path})
        good_inputs.append(path)
        root_file.Close()

    if merged_hist is None:
        raise RuntimeError("No mergeable histograms were found.")

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    output.cd()

    merged_hist.SetEntries(total_hist_entries)
    merged_hist.Write()
    if merged_cutflow:
        merged_cutflow.Write()

    chain = ROOT.TChain(args.tree_name)
    for path in good_inputs:
        chain.Add(path)
    if chain.GetEntries() > 0:
        tracks = chain.CloneTree(-1, "fast")
        tracks.Write(args.tree_name)

    merged_config = {
        "merge_inputs": good_inputs,
        "n_merge_inputs": len(good_inputs),
        "source_configs": configs,
        "note": "Merged from per-file cosmic_voxel_map.py outputs.",
    }
    output.WriteObject(ROOT.TObjString(json.dumps(merged_config, indent=2, sort_keys=True)), "config_json")
    output.Close()

    print(
        "Merged {n} files into {output} with {entries:g} histogram entries.".format(
            n=len(good_inputs),
            output=args.output,
            entries=total_hist_entries,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
