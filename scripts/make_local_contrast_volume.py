#!/usr/bin/env python3
"""Build a broad-trend-subtracted transmission contrast volume."""

import argparse
import json
import math
import os

from calibrate_lima_expectation import (
    axis_centers,
    clone_empty,
    fill_hist_from_array,
    hist_to_array,
    require_th3,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Subtract a very broad 3D Gaussian baseline from a log-transmission "
            "TH3. Positive residuals are a diagnostic cavity/shaft-like contrast."
        )
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--log-transmission-hist", default="log_transmission_excess")
    parser.add_argument("--expected-hist", default="calibrated_expected_voxels")
    parser.add_argument("--observed-hist", default="observed_voxels")
    parser.add_argument("--min-expected", type=float, default=1.0e-6)
    parser.add_argument("--x-sigma-cm", type=float, default=4000.0)
    parser.add_argument("--y-sigma-cm", type=float, default=4000.0)
    parser.add_argument("--z-sigma-cm", type=float, default=5000.0)
    return parser.parse_args()


def gaussian_kernel(sigma_bins):
    import numpy as np

    if sigma_bins <= 0.0:
        return np.array([1.0], dtype="float64")
    half_width = max(1, int(math.ceil(4.0 * sigma_bins)))
    offsets = np.arange(-half_width, half_width + 1, dtype="float64")
    kernel = np.exp(-0.5 * (offsets / sigma_bins) ** 2)
    kernel /= kernel.sum()
    return kernel


def convolve_along_axis(values, kernel, axis_index):
    import numpy as np

    moved = np.ascontiguousarray(np.moveaxis(values, axis_index, 0))
    result = np.zeros_like(moved)
    half_width = len(kernel) // 2
    flat_in = moved.reshape((moved.shape[0], -1))
    flat_out = result.reshape((result.shape[0], -1))
    for icol in range(flat_in.shape[1]):
        padded = np.pad(flat_in[:, icol], (half_width, half_width), mode="constant")
        flat_out[:, icol] = np.convolve(padded, kernel, mode="valid")
    return np.moveaxis(result, 0, axis_index)


def weighted_gaussian_smooth(values, weights, sigma_bins):
    import numpy as np

    numerator = values * weights
    denominator = weights.copy()
    for axis_index, sigma in enumerate(sigma_bins):
        kernel = gaussian_kernel(sigma)
        numerator = convolve_along_axis(numerator, kernel, axis_index)
        denominator = convolve_along_axis(denominator, kernel, axis_index)
    smooth = np.zeros_like(values)
    valid = denominator > 0.0
    smooth[valid] = numerator[valid] / denominator[valid]
    return smooth


def main():
    args = parse_args()

    import numpy as np
    import ROOT

    source = ROOT.TFile.Open(args.input)
    if not source or source.IsZombie():
        raise RuntimeError(f"Could not open {args.input}")
    log_hist = require_th3(source, args.log_transmission_hist)
    expected_hist = require_th3(source, args.expected_hist)
    observed_hist = require_th3(source, args.observed_hist)
    source.Close()

    log_values = hist_to_array(log_hist)
    expected = hist_to_array(expected_hist)
    weights = (expected > args.min_expected).astype("float64")

    sigma_bins = (
        args.x_sigma_cm / log_hist.GetXaxis().GetBinWidth(1),
        args.y_sigma_cm / log_hist.GetYaxis().GetBinWidth(1),
        args.z_sigma_cm / log_hist.GetZaxis().GetBinWidth(1),
    )
    smooth = weighted_gaussian_smooth(log_values, weights, sigma_bins)
    contrast = np.zeros_like(log_values)
    valid = weights > 0.0
    contrast[valid] = log_values[valid] - smooth[valid]
    cavity_contrast = np.maximum(contrast, 0.0)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create {args.output}")
    output.cd()
    observed_hist.Write("observed_voxels")
    expected_hist.Write("calibrated_expected_voxels")
    log_hist.Write("log_transmission_excess")

    products = {
        "smooth_log_transmission_baseline": (
            smooth,
            "broad Gaussian baseline of log transmission residual",
        ),
        "local_log_transmission_contrast": (
            contrast,
            "local log-transmission contrast after broad-baseline subtraction",
        ),
        "local_cavity_contrast": (
            cavity_contrast,
            "positive local transmission contrast; diagnostic cavity-like score",
        ),
    }
    for name, (values, title) in products.items():
        hist = clone_empty(
            ROOT,
            log_hist,
            name,
            f"{title};{log_hist.GetXaxis().GetTitle()};{log_hist.GetYaxis().GetTitle()};{log_hist.GetZaxis().GetTitle()}",
        )
        fill_hist_from_array(hist, values)
        hist.SetEntries(float(np.count_nonzero(values)))
        hist.Write()

    config = {
        "input": args.input,
        "log_transmission_hist": args.log_transmission_hist,
        "expected_hist": args.expected_hist,
        "observed_hist": args.observed_hist,
        "min_expected": args.min_expected,
        "sigma_cm": {
            "x": args.x_sigma_cm,
            "y": args.y_sigma_cm,
            "z": args.z_sigma_cm,
        },
        "sigma_bins": {
            "x": sigma_bins[0],
            "y": sigma_bins[1],
            "z": sigma_bins[2],
        },
        "valid_voxels": int(np.count_nonzero(valid)),
        "contrast_range": [float(contrast.min()), float(contrast.max())],
        "note": (
            "Diagnostic only: this subtracts a broad residual trend from the "
            "log-transmission map. Positive values are cavity-like transmission "
            "excesses relative to nearby large-scale behavior, not standalone "
            "Li-Ma significances."
        ),
    }
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "local_contrast_config_json")
    output.Close()
    print(
        "Wrote {output}: valid_voxels={valid} contrast_range=[{lo:g},{hi:g}]".format(
            output=args.output,
            valid=config["valid_voxels"],
            lo=config["contrast_range"][0],
            hi=config["contrast_range"][1],
        )
    )


if __name__ == "__main__":
    main()
