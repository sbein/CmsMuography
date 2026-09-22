#!/usr/bin/env python3
"""Build a first-pass solid-ground expected voxel map from EXPACS/PARMA flux tables."""

import argparse
import csv
import json
import math
import os

from cosmic_voxel_map import fill_voxels, ray_box_intersection
from suppress_voxel_pedestals import (
    DEFAULT_GROUND_COEFFICIENTS_CM,
    ground_z_cavern_cm,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Use the accepted-track geometry as a detector-acceptance kernel, weight each "
            "ray by an EXPACS/PARMA surface muon flux integral above the solid-ground "
            "range threshold, and fill an expected voxel occupancy."
        )
    )
    parser.add_argument("--data-root", required=True, help="ROOT file containing the observed TH3.")
    parser.add_argument("--data-hist", default="cavern_muon_voxels")
    parser.add_argument("--tracks-root", required=True, help="ROOT file containing the accepted tracks tree.")
    parser.add_argument("--tracks-tree", default="tracks")
    parser.add_argument("--flux-table", required=True, help="CSV from expacs_muon_grid.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--nchunks", type=int, default=1)
    parser.add_argument("--max-tracks", type=int, default=-1)
    parser.add_argument("--progress-every", type=int, default=100000)
    parser.add_argument("--density-g-cm3", type=float, default=2.30)
    parser.add_argument("--dedx-a-mev-gcm2", type=float, default=2.0)
    parser.add_argument("--dedx-b-per-gcm2", type=float, default=4.0e-6)
    parser.add_argument(
        "--ground-coefficients-cm",
        default=",".join(str(value) for value in DEFAULT_GROUND_COEFFICIENTS_CM),
        help="Same six-coefficient cavern ground model used by suppress_voxel_pedestals.py.",
    )
    parser.add_argument("--max-voxels-per-track", type=int, default=30000)
    parser.add_argument(
        "--no-cos-projection",
        action="store_true",
        help="Do not multiply the surface angular flux by cos(theta).",
    )
    parser.add_argument(
        "--weight-floor",
        type=float,
        default=0.0,
        help="Minimum expected weight needed before filling a ray.",
    )
    return parser.parse_args()


class FluxTable:
    def __init__(self, path):
        rows = []
        with open(path, encoding="utf-8") as handle:
            reader = csv.DictReader(line for line in handle if not line.startswith("#"))
            for row in reader:
                rows.append({key: float(value) for key, value in row.items()})
        if not rows:
            raise RuntimeError(f"No flux rows found in {path}")

        self.energy_bins = sorted(
            {(row["energy_low_mev"], row["energy_high_mev"], row["energy_mid_mev"]) for row in rows}
        )
        self.cos_bins = sorted(
            {(row["cos_low"], row["cos_high"], row["cos_mid"]) for row in rows}
        )
        self.ne = len(self.energy_bins)
        self.ncos = len(self.cos_bins)
        self.flux = [[0.0 for _ in range(self.ne)] for _ in range(self.ncos)]

        energy_index = {item: index for index, item in enumerate(self.energy_bins)}
        cos_index = {item: index for index, item in enumerate(self.cos_bins)}
        for row in rows:
            ie = energy_index[(row["energy_low_mev"], row["energy_high_mev"], row["energy_mid_mev"])]
            ic = cos_index[(row["cos_low"], row["cos_high"], row["cos_mid"])]
            self.flux[ic][ie] = row["flux"]

    def cos_index(self, cosz):
        cosz = min(max(cosz, 0.0), 0.999999999)
        best = 0
        best_distance = 999.0
        for index, (_low, _high, mid) in enumerate(self.cos_bins):
            distance = abs(cosz - mid)
            if distance < best_distance:
                best = index
                best_distance = distance
        return best

    def integral_above(self, emin_mev, cosz):
        ic = self.cos_index(cosz)
        total = 0.0
        for ie, (elow, ehigh, _emid) in enumerate(self.energy_bins):
            if ehigh <= emin_mev:
                continue
            width = ehigh - max(elow, emin_mev)
            if width > 0.0:
                total += self.flux[ic][ie] * width
        return total


class GroundArgs:
    def __init__(self, coefficients):
        self.ground_coefficients_cm = coefficients


def hist_bounds_and_bins(hist):
    return {
        "xmin": hist.GetXaxis().GetXmin(),
        "xmax": hist.GetXaxis().GetXmax(),
        "ymin": hist.GetYaxis().GetXmin(),
        "ymax": hist.GetYaxis().GetXmax(),
        "zmin": hist.GetZaxis().GetXmin(),
        "zmax": hist.GetZaxis().GetXmax(),
    }, (hist.GetNbinsX(), hist.GetNbinsY(), hist.GetNbinsZ())


def clone_empty_th3(root, source, name, title):
    hist = root.TH3F(
        name,
        title,
        source.GetNbinsX(),
        source.GetXaxis().GetXmin(),
        source.GetXaxis().GetXmax(),
        source.GetNbinsY(),
        source.GetYaxis().GetXmin(),
        source.GetYaxis().GetXmax(),
        source.GetNbinsZ(),
        source.GetZaxis().GetXmin(),
        source.GetZaxis().GetXmax(),
    )
    hist.GetXaxis().SetTitle(source.GetXaxis().GetTitle())
    hist.GetYaxis().SetTitle(source.GetYaxis().GetTitle())
    hist.GetZaxis().SetTitle(source.GetZaxis().GetTitle())
    return hist


def solid_ground_range_emin_mev(depth_g_cm2, a_mev_gcm2, b_per_gcm2):
    depth_g_cm2 = max(depth_g_cm2, 0.0)
    if b_per_gcm2 <= 0.0:
        return a_mev_gcm2 * depth_g_cm2
    exponent = min(b_per_gcm2 * depth_g_cm2, 50.0)
    return (a_mev_gcm2 / b_per_gcm2) * (math.exp(exponent) - 1.0)


def ground_intersection_t(point, direction, ground_args):
    if direction[2] <= 0.0:
        return None

    def residual(tvalue):
        x = point[0] + tvalue * direction[0]
        y = point[1] + tvalue * direction[1]
        z = point[2] + tvalue * direction[2]
        return z - ground_z_cavern_cm(ground_args, x, y)

    low = 0.0
    flow = residual(low)
    if flow >= 0.0:
        return 0.0
    high = 1000.0
    for _attempt in range(40):
        if residual(high) >= 0.0:
            break
        high *= 1.5
    else:
        return None

    for _iteration in range(80):
        mid = 0.5 * (low + high)
        if residual(mid) >= 0.0:
            high = mid
        else:
            low = mid
    return high


def main():
    args = parse_args()
    if args.nchunks <= 0:
        raise RuntimeError("--nchunks must be positive")
    if not (0 <= args.chunk_index < args.nchunks):
        raise RuntimeError("--chunk-index must be in [0, nchunks)")

    import ROOT

    data_file = ROOT.TFile.Open(args.data_root)
    if not data_file or data_file.IsZombie():
        raise RuntimeError(f"Could not open {args.data_root}")
    data_hist = data_file.Get(args.data_hist)
    if not data_hist:
        raise RuntimeError(f"Could not find {args.data_hist} in {args.data_root}")
    data_hist.SetDirectory(0)
    data_file.Close()

    tracks_file = ROOT.TFile.Open(args.tracks_root)
    if not tracks_file or tracks_file.IsZombie():
        raise RuntimeError(f"Could not open {args.tracks_root}")
    tracks = tracks_file.Get(args.tracks_tree)
    if not tracks:
        raise RuntimeError(f"Could not find {args.tracks_tree} in {args.tracks_root}")

    flux_table = FluxTable(args.flux_table)
    ground_args = GroundArgs(args.ground_coefficients_cm)
    bounds, bins = hist_bounds_and_bins(data_hist)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    expected = clone_empty_th3(
        ROOT,
        data_hist,
        "solid_ground_expected_voxels",
        "EXPACS solid-ground expected voxel crossings;{};{};{}".format(
            data_hist.GetXaxis().GetTitle(),
            data_hist.GetYaxis().GetTitle(),
            data_hist.GetZaxis().GetTitle(),
        ),
    )

    max_entries = tracks.GetEntries() if args.max_tracks < 0 else min(args.max_tracks, tracks.GetEntries())
    processed = 0
    filled = 0
    total_raw_weight = 0.0
    total_voxel_weight = 0.0
    skipped_no_ground = 0
    skipped_no_flux = 0

    for entry in range(args.chunk_index, max_entries, args.nchunks):
        tracks.GetEntry(entry)
        point = (float(tracks.vx), float(tracks.vy), float(tracks.vz))
        direction = (float(tracks.ux), float(tracks.uy), float(tracks.uz))
        norm = math.sqrt(sum(value * value for value in direction))
        if norm <= 0.0:
            continue
        direction = tuple(value / norm for value in direction)
        if direction[2] <= 0.0:
            continue

        processed += 1
        t_ground = ground_intersection_t(point, direction, ground_args)
        if t_ground is None:
            skipped_no_ground += 1
            continue
        depth = args.density_g_cm3 * t_ground
        emin = solid_ground_range_emin_mev(depth, args.dedx_a_mev_gcm2, args.dedx_b_per_gcm2)
        cosz = max(direction[2], 0.0)
        weight = flux_table.integral_above(emin, cosz)
        if not args.no_cos_projection:
            weight *= cosz
        if weight <= args.weight_floor:
            skipped_no_flux += 1
            continue

        segment = ray_box_intersection(point, direction, bounds)
        if segment is None:
            continue
        n_voxels, _truncated = fill_voxels(
            expected,
            point,
            direction,
            segment[0],
            segment[1],
            bounds,
            bins,
            weight,
            args.max_voxels_per_track,
        )
        if n_voxels <= 0:
            continue
        filled += 1
        total_raw_weight += weight
        total_voxel_weight += weight * n_voxels
        if args.progress_every > 0 and processed % args.progress_every == 0:
            print(
                f"chunk {args.chunk_index}/{args.nchunks}: processed={processed} "
                f"filled={filled} voxel_weight={total_voxel_weight:g}"
            )

    expected.SetEntries(total_voxel_weight)
    output.cd()
    expected.Write()
    config = {
        "data_root": args.data_root,
        "data_hist": args.data_hist,
        "tracks_root": args.tracks_root,
        "tracks_tree": args.tracks_tree,
        "flux_table": args.flux_table,
        "flux_model": {
            "name": "EXPACS/PARMA",
            "url": "https://phits.jaea.go.jp/expacs/",
            "citations": [
                "T. Sato, PLOS ONE 10(12): e0144679 (2015), doi:10.1371/journal.pone.0144679",
                "T. Sato, PLOS ONE 11(8): e0160390 (2016), doi:10.1371/journal.pone.0160390",
            ],
        },
        "chunk_index": args.chunk_index,
        "nchunks": args.nchunks,
        "max_tracks": int(max_entries),
        "processed_chunk_tracks": processed,
        "filled_tracks": filled,
        "raw_weight_sum": total_raw_weight,
        "voxel_weight_sum": total_voxel_weight,
        "skipped_no_ground_intersection": skipped_no_ground,
        "skipped_no_flux": skipped_no_flux,
        "density_g_cm3": args.density_g_cm3,
        "dedx_a_mev_gcm2": args.dedx_a_mev_gcm2,
        "dedx_b_per_gcm2": args.dedx_b_per_gcm2,
        "ground_coefficients_cm": args.ground_coefficients_cm,
        "note": (
            "First-pass null model: accepted data track directions are used as the detector "
            "acceptance kernel; EXPACS/PARMA supplies surface muon flux; solid ground is "
            "approximated by a continuous-slowing-down muon range threshold."
        ),
    }
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "config_json")
    output.Close()
    tracks_file.Close()
    print(
        "Wrote {output}: processed={processed} filled={filled} voxel_weight={weight:g}".format(
            output=args.output,
            processed=processed,
            filled=filled,
            weight=total_voxel_weight,
        )
    )


if __name__ == "__main__":
    main()
