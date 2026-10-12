#!/usr/bin/env python3
"""Create tiny ROOT/CSV fixtures for the end-to-end smoke pipeline."""

import argparse
import csv
import math
import os
from array import array


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--output-dir", required=True)
    parser.add_argument("--ntracks", type=int, default=24)
    return parser.parse_args()


def write_flux_table(path):
    rows = []
    energy_edges = [30000.0, 60000.0, 120000.0, 250000.0]
    cos_edges = [0.20, 0.45, 0.70, 1.00]
    for i in range(len(energy_edges) - 1):
        for j in range(len(cos_edges) - 1):
            e_low = energy_edges[i]
            e_high = energy_edges[i + 1]
            cos_low = cos_edges[j]
            cos_high = cos_edges[j + 1]
            e_mid = math.sqrt(e_low * e_high)
            cos_mid = 0.5 * (cos_low + cos_high)
            flux = 1.0e-7 * (e_mid / 30000.0) ** -2.7 * (0.3 + cos_mid)
            rows.append(
                {
                    "energy_low_mev": e_low,
                    "energy_high_mev": e_high,
                    "energy_mid_mev": e_mid,
                    "cos_low": cos_low,
                    "cos_high": cos_high,
                    "cos_mid": cos_mid,
                    "flux": flux,
                }
            )
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "energy_low_mev",
                "energy_high_mev",
                "energy_mid_mev",
                "cos_low",
                "cos_high",
                "cos_mid",
                "flux",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def write_tracks_root(path, ntracks):
    import ROOT

    output = ROOT.TFile.Open(path, "RECREATE")
    tree = ROOT.TTree("tracks", "synthetic accepted cosmic-muon tracks for smoke tests")

    ints = {name: array("i", [0]) for name in ("run", "lumi", "event")}
    floats = {
        name: array("f", [0.0])
        for name in (
            "vx",
            "vy",
            "vz",
            "ux",
            "uy",
            "uz",
            "t_enter",
            "t_exit",
        )
    }
    for name, value in ints.items():
        tree.Branch(name, value, f"{name}/I")
    for name, value in floats.items():
        tree.Branch(name, value, f"{name}/F")

    for index in range(ntracks):
        phi = 2.0 * math.pi * index / max(ntracks, 1)
        radius = 250.0 + 15.0 * (index % 5)
        ints["run"][0] = 1
        ints["lumi"][0] = 1 + index // 8
        ints["event"][0] = 1000 + index
        floats["vx"][0] = radius * math.cos(phi)
        floats["vy"][0] = 0.65 * radius * math.sin(phi)
        floats["vz"][0] = 100.0 + 20.0 * (index % 4)
        ux = 0.10 * math.cos(phi + 0.4)
        uy = 0.06 * math.sin(phi - 0.2)
        uz = math.sqrt(max(1.0 - ux * ux - uy * uy, 0.0))
        norm = math.sqrt(ux * ux + uy * uy + uz * uz)
        floats["ux"][0] = ux / norm
        floats["uy"][0] = uy / norm
        floats["uz"][0] = uz / norm
        floats["t_enter"][0] = 0.0
        floats["t_exit"][0] = 4500.0
        tree.Fill()

    output.cd()
    tree.Write()
    output.Close()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    flux_path = os.path.join(args.output_dir, "smoke_expacs_flux.csv")
    tracks_path = os.path.join(args.output_dir, "smoke_tracks.root")
    write_flux_table(flux_path)
    write_tracks_root(tracks_path, args.ntracks)
    print(f"wrote {flux_path}")
    print(f"wrote {tracks_path}")


if __name__ == "__main__":
    main()
