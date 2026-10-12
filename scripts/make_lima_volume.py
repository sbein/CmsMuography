#!/usr/bin/env python3
"""Compute per-voxel excess/deficit significance from observed and expected TH3 maps."""

import argparse
import json
import math
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build Li-Ma/profile-likelihood significance TH3 volumes from observed and expected maps."
    )
    parser.add_argument("--observed-root", required=True)
    parser.add_argument("--observed-hist", default="cavern_muon_voxels")
    parser.add_argument("--expected-root", required=True)
    parser.add_argument("--expected-hist", default="generator_expected_voxels")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument(
        "--normalize",
        choices=("none", "global"),
        default="global",
        help="Scale the expected integral to the observed integral before computing significances.",
    )
    parser.add_argument(
        "--method",
        choices=("poisson-known", "lima"),
        default="lima",
        help="poisson-known treats expected as fixed; lima uses an off region with --alpha.",
    )
    parser.add_argument(
        "--alpha",
        type=float,
        default=0.05,
        help="Li-Ma exposure ratio N_expected = alpha * N_off when --method lima.",
    )
    parser.add_argument("--min-expected", type=float, default=1.0e-9)
    return parser.parse_args()


def require_th3(root_file, name):
    hist = root_file.Get(name)
    if not hist:
        raise RuntimeError(f"Could not find {name} in {root_file.GetName()}")
    if not hist.InheritsFrom("TH3"):
        raise RuntimeError(f"{name} is not a TH3")
    hist.SetDirectory(0)
    return hist


def clone_empty(root, source, name, title):
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


def poisson_known_z(observed, expected):
    if expected <= 0.0:
        return 0.0
    if observed <= 0.0:
        q = 2.0 * expected
    else:
        q = 2.0 * (observed * math.log(observed / expected) - (observed - expected))
    q = max(q, 0.0)
    sign = 1.0 if observed >= expected else -1.0
    return sign * math.sqrt(q)


def lima_z(observed, expected, alpha):
    if alpha <= 0.0 or expected <= 0.0:
        return 0.0
    off = expected / alpha
    total = observed + off
    if total <= 0.0:
        return 0.0
    term = 0.0
    if observed > 0.0:
        term += observed * math.log(((1.0 + alpha) / alpha) * observed / total)
    if off > 0.0:
        term += off * math.log((1.0 + alpha) * off / total)
    q = max(2.0 * term, 0.0)
    sign = 1.0 if observed >= expected else -1.0
    return sign * math.sqrt(q)


def check_compatible(observed, expected):
    if (
        observed.GetNbinsX() != expected.GetNbinsX()
        or observed.GetNbinsY() != expected.GetNbinsY()
        or observed.GetNbinsZ() != expected.GetNbinsZ()
    ):
        raise RuntimeError("Observed and expected histograms have different bin counts")


def main():
    args = parse_args()

    import ROOT

    observed_file = ROOT.TFile.Open(args.observed_root)
    expected_file = ROOT.TFile.Open(args.expected_root)
    if not observed_file or observed_file.IsZombie():
        raise RuntimeError(f"Could not open {args.observed_root}")
    if not expected_file or expected_file.IsZombie():
        raise RuntimeError(f"Could not open {args.expected_root}")
    observed = require_th3(observed_file, args.observed_hist)
    expected = require_th3(expected_file, args.expected_hist)
    observed_file.Close()
    expected_file.Close()
    check_compatible(observed, expected)

    obs_integral = observed.Integral()
    exp_integral = expected.Integral()
    scale = 1.0
    if args.normalize == "global":
        if exp_integral <= 0.0:
            raise RuntimeError("Cannot globally normalize an expected histogram with zero integral")
        scale = obs_integral / exp_integral
        expected.Scale(scale)

    excess = clone_empty(
        ROOT,
        observed,
        "cavity_excess_significance",
        "cavity-candidate excess significance;{};{};{}".format(
            observed.GetXaxis().GetTitle(),
            observed.GetYaxis().GetTitle(),
            observed.GetZaxis().GetTitle(),
        ),
    )
    signed = clone_empty(
        ROOT,
        observed,
        "signed_significance",
        "signed observed-minus-expected significance;{};{};{}".format(
            observed.GetXaxis().GetTitle(),
            observed.GetYaxis().GetTitle(),
            observed.GetZaxis().GetTitle(),
        ),
    )
    ratio = clone_empty(
        ROOT,
        observed,
        "observed_over_expected",
        "observed / expected;{};{};{}".format(
            observed.GetXaxis().GetTitle(),
            observed.GetYaxis().GetTitle(),
            observed.GetZaxis().GetTitle(),
        ),
    )

    zsum = 0.0
    nz = 0
    max_z = -999.0
    min_z = 999.0
    for ix in range(1, observed.GetNbinsX() + 1):
        for iy in range(1, observed.GetNbinsY() + 1):
            for iz in range(1, observed.GetNbinsZ() + 1):
                obs = observed.GetBinContent(ix, iy, iz)
                exp = expected.GetBinContent(ix, iy, iz)
                if exp < args.min_expected:
                    continue
                if args.method == "lima":
                    z = lima_z(obs, exp, args.alpha)
                else:
                    z = poisson_known_z(obs, exp)
                signed.SetBinContent(ix, iy, iz, z)
                excess.SetBinContent(ix, iy, iz, max(z, 0.0))
                ratio.SetBinContent(ix, iy, iz, obs / exp if exp > 0.0 else 0.0)
                zsum += z
                nz += 1
                max_z = max(max_z, z)
                min_z = min(min_z, z)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    output.cd()
    observed.Write("observed_voxels")
    expected.Write("expected_voxels_scaled")
    signed.Write()
    excess.Write()
    ratio.Write()
    config = {
        "observed_root": args.observed_root,
        "observed_hist": args.observed_hist,
        "expected_root": args.expected_root,
        "expected_hist": args.expected_hist,
        "normalization": args.normalize,
        "expected_scale": scale,
        "observed_integral": obs_integral,
        "expected_integral_raw": exp_integral,
        "expected_integral_scaled": expected.Integral(),
        "method": args.method,
        "alpha": args.alpha,
        "min_expected": args.min_expected,
        "filled_voxels": nz,
        "mean_signed_z": zsum / nz if nz else 0.0,
        "max_signed_z": max_z if nz else 0.0,
        "min_signed_z": min_z if nz else 0.0,
    }
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "config_json")
    output.Close()
    print(
        "Wrote {output}: scale={scale:g} voxels={voxels} z_range=[{zmin:g},{zmax:g}]".format(
            output=args.output,
            scale=scale,
            voxels=nz,
            zmin=config["min_signed_z"],
            zmax=config["max_signed_z"],
        )
    )


if __name__ == "__main__":
    main()
