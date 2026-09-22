#!/usr/bin/env python3
"""Rebuild a clean TH3 voxel map from a saved cosmic-muon tracks tree."""

import argparse
import json
import math
import os

from cosmic_voxel_map import (
    DEFAULT_VOLUMES,
    axis_titles,
    bin_counts,
    effective_bounds,
    fill_voxels,
    make_cutflow,
    ray_box_intersection,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read the tracks tree written by cosmic_voxel_map.py and refill the "
            "cavern_muon_voxels TH3. This is useful for recovering a clean map "
            "from an interrupted ROOT file whose tree was recoverable."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input ROOT file with a tracks tree.")
    parser.add_argument("-o", "--output", required=True, help="Output clean ROOT file.")
    parser.add_argument("--tree-name", default="tracks", help="Input tree name.")
    parser.add_argument(
        "--coordinate-system",
        choices=("cavern", "cms"),
        default="cavern",
        help="Coordinate system of vx/vy/vz and ux/uy/uz in the tracks tree.",
    )
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--ymin", type=float, default=None)
    parser.add_argument("--ymax", type=float, default=None)
    parser.add_argument("--zmin", type=float, default=None)
    parser.add_argument("--zmax", type=float, default=None)
    parser.add_argument("--voxel-size", type=float, default=100.0)
    parser.add_argument("--nx", type=int, default=0)
    parser.add_argument("--ny", type=int, default=0)
    parser.add_argument("--nz", type=int, default=0)
    parser.add_argument(
        "--trace-axis",
        choices=("x", "y", "z"),
        default=None,
        help="Trace axis used for half-volume clipping; defaults match cosmic_voxel_map.py.",
    )
    parser.add_argument(
        "--trace-sign",
        choices=("positive", "negative"),
        default="positive",
        help="Trace sign used for half-volume clipping.",
    )
    parser.add_argument(
        "--allow-negative-trace-side",
        action="store_true",
        help="Do not clip to the selected positive/negative trace half-space.",
    )
    parser.add_argument("--max-tracks", type=int, default=-1, help="Maximum tracks to rebuild.")
    parser.add_argument("--max-voxels-per-track", type=int, default=20000)
    parser.add_argument("--weight", type=float, default=1.0)
    parser.add_argument("--progress-every", type=int, default=100000)
    parser.add_argument(
        "--copy-tracks",
        action="store_true",
        help="Write a clean copy of the source tracks tree to the output file.",
    )
    parser.add_argument(
        "--reuse-stored-intersections",
        action="store_true",
        help=(
            "Use t_enter/t_exit from the input tracks tree. By default intersections "
            "are recomputed against the requested output bounds, which is needed when "
            "changing the voxel volume."
        ),
    )
    return resolve_args(parser.parse_args())


def resolve_args(args):
    defaults = DEFAULT_VOLUMES[args.coordinate_system]
    for key, value in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    if args.trace_axis is None:
        args.trace_axis = "z" if args.coordinate_system == "cavern" else "y"
    return args


def finite(value):
    try:
        value = float(value)
    except Exception:
        return False
    return math.isfinite(value)


def main():
    args = parse_args()

    import ROOT

    source = ROOT.TFile.Open(args.input)
    if not source or source.IsZombie():
        raise RuntimeError(f"Could not open input file: {args.input}")
    tracks = source.Get(args.tree_name)
    if not tracks:
        raise RuntimeError(f"Could not find tree {args.tree_name!r} in {args.input}")

    bounds = effective_bounds(args)
    bins = bin_counts(args, bounds)
    nx, ny, nz = bins
    x_title, y_title, z_title = axis_titles(args.coordinate_system)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")

    voxel_counts = ROOT.TH3F(
        "cavern_muon_voxels",
        (
            f"cosmic-muon voxel crossings rebuilt from tracks along {args.trace_sign} "
            f"{args.coordinate_system} {args.trace_axis};{x_title};{y_title};{z_title}"
        ),
        nx,
        bounds["xmin"],
        bounds["xmax"],
        ny,
        bounds["ymin"],
        bounds["ymax"],
        nz,
        bounds["zmin"],
        bounds["zmax"],
    )
    cutflow, cut_bins = make_cutflow(ROOT)

    max_tracks = tracks.GetEntries() if args.max_tracks < 0 else min(args.max_tracks, tracks.GetEntries())
    total_voxel_fills = 0
    rebuilt_tracks = 0
    truncated_tracks = 0
    accepted_event_keys = set()

    for entry in range(max_tracks):
        tracks.GetEntry(entry)
        point = (float(tracks.vx), float(tracks.vy), float(tracks.vz))
        direction = (float(tracks.ux), float(tracks.uy), float(tracks.uz))
        if not all(finite(value) for value in point + direction):
            continue
        if args.reuse_stored_intersections:
            t_enter = float(tracks.t_enter)
            t_exit = float(tracks.t_exit)
            if not all(finite(value) for value in (t_enter, t_exit)):
                continue
            if t_exit < t_enter:
                continue
        else:
            intersection = ray_box_intersection(point, direction, bounds)
            if intersection is None:
                continue
            t_enter, t_exit = intersection

        n_voxels, truncated = fill_voxels(
            voxel_counts,
            point,
            direction,
            t_enter,
            t_exit,
            bounds,
            bins,
            args.weight,
            args.max_voxels_per_track,
        )
        if n_voxels <= 0:
            continue

        rebuilt_tracks += 1
        total_voxel_fills += n_voxels
        accepted_event_keys.add((int(tracks.run), int(tracks.lumi), int(tracks.event)))
        if truncated:
            truncated_tracks += 1
        if args.progress_every > 0 and rebuilt_tracks % args.progress_every == 0:
            print(
                f"rebuilt {rebuilt_tracks}/{max_tracks} tracks "
                f"with {total_voxel_fills} voxel fills"
            )

    cutflow.SetBinContent(cut_bins["events"], len(accepted_event_keys))
    cutflow.SetBinContent(cut_bins["accepted_tracks"], rebuilt_tracks)
    cutflow.SetBinContent(cut_bins["voxel_fills"], total_voxel_fills)
    cutflow.SetBinContent(cut_bins["truncated_tracks"], truncated_tracks)
    voxel_counts.SetEntries(total_voxel_fills)

    config = {
        "source_file": args.input,
        "source_tree": args.tree_name,
        "source_tree_entries": int(tracks.GetEntries()),
        "rebuilt_tracks": rebuilt_tracks,
        "accepted_events_with_tracks": len(accepted_event_keys),
        "coordinate_system": args.coordinate_system,
        "bounds_cm": bounds,
        "bins": {"nx": nx, "ny": ny, "nz": nz},
        "voxel_size_request_cm": args.voxel_size,
        "trace_axis": args.trace_axis,
        "trace_sign": args.trace_sign,
        "allow_negative_trace_side": args.allow_negative_trace_side,
        "intersection_source": "stored_tree" if args.reuse_stored_intersections else "recomputed_from_requested_bounds",
        "charge_handling": "positive and negative reconstructed muons are combined in the same histograms",
        "note": (
            "This file was rebuilt from the saved accepted-track tree, so the cutflow "
            "contains accepted tracks and accepted events with tracks, not the full "
            "preselection event/candidate accounting from the original AOD loop."
        ),
    }

    output.cd()
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "config_json")
    voxel_counts.Write()
    cutflow.Write()
    if args.copy_tracks:
        tracks.SetBranchStatus("*", 1)
        clean_tracks = tracks.CloneTree(max_tracks, "fast")
        clean_tracks.Write(args.tree_name)
    output.Close()
    source.Close()

    print(
        "Wrote {output} from {tracks} tracks and {voxels} voxel fills.".format(
            output=args.output,
            tracks=rebuilt_tracks,
            voxels=total_voxel_fills,
        )
    )


if __name__ == "__main__":
    main()
