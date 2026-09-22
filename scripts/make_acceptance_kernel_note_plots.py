#!/usr/bin/env python3
"""Make compact figures for the data-seeded EXPACS acceptance-kernel note."""

import argparse
import math
import os

import ROOT

from suppress_voxel_pedestals import clone_hist_with_scaled_axes, make_slice, view_axis_order


DEFAULT_INPUT = (
    "test_output/"
    "OORun2025_AOD_disk55_no_kipt_zonly_smoothcalib_lima_alpha0p05_120m_rhcavern.root"
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", default=DEFAULT_INPUT)
    parser.add_argument("-o", "--output-dir", default="docs/figures")
    parser.add_argument(
        "--slice-z-cm",
        type=float,
        default=250.0,
        help="Cavern z position for the horizontal XY slice, in cm.",
    )
    parser.add_argument("--width", type=int, default=1200)
    parser.add_argument("--height", type=int, default=950)
    return parser.parse_args()


def finite_values(hist2, positive=False):
    values = []
    for ix in range(1, hist2.GetNbinsX() + 1):
        for iy in range(1, hist2.GetNbinsY() + 1):
            value = float(hist2.GetBinContent(ix, iy))
            if not math.isfinite(value):
                continue
            if positive and value <= 0.0:
                continue
            values.append(value)
    values.sort()
    return values


def quantile(values, fraction):
    if not values:
        return 0.0
    if len(values) == 1:
        return values[0]
    position = min(max(fraction, 0.0), 1.0) * (len(values) - 1)
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return values[low]
    weight = position - low
    return values[low] * (1.0 - weight) + values[high] * weight


def set_robust_range(hist2, logz=False, symmetric=False, lower=0.01, upper=0.995):
    values = finite_values(hist2, positive=logz)
    if not values:
        return
    if logz:
        hist2.SetMinimum(max(quantile(values, lower), 1.0e-4))
        hist2.SetMaximum(max(quantile(values, upper), hist2.GetMinimum() * 10.0))
        return
    if symmetric:
        span = max(abs(quantile(values, lower)), abs(quantile(values, upper)))
        if span > 0:
            hist2.SetMinimum(-span)
            hist2.SetMaximum(span)
        return
    zmin = quantile(values, lower)
    zmax = quantile(values, upper)
    if zmax > zmin:
        hist2.SetMinimum(zmin)
        hist2.SetMaximum(zmax)


def xy_slice(root, hist3_cm, name, slice_z_cm):
    hist3_m = clone_hist_with_scaled_axes(root, hist3_cm, f"{name}_m", "m")
    order = view_axis_order(hist3_m)
    slice_z_m = slice_z_cm / 100.0
    hist2 = make_slice(root, hist3_m, *order["xy_slice_z"], slice_z_m, f"{name}_xy_z")
    hist2.SetDirectory(0)
    return hist2


def draw_plot(root, hist2, output_path, title, logz=False, symmetric=False):
    canvas = root.TCanvas(os.path.splitext(os.path.basename(output_path))[0], title, 1200, 950)
    pad = canvas.cd()
    pad.SetRightMargin(0.16)
    pad.SetLeftMargin(0.12)
    pad.SetBottomMargin(0.12)
    pad.SetTopMargin(0.08)
    pad.SetLogz(bool(logz))
    hist2.SetTitle(title)
    hist2.SetStats(False)
    set_robust_range(hist2, logz=logz, symmetric=symmetric)
    hist2.Draw("COLZ")
    canvas.SaveAs(output_path)
    return canvas


def main():
    args = parse_args()
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)

    infile = ROOT.TFile.Open(args.input)
    if not infile or infile.IsZombie():
        raise RuntimeError(f"Could not open {args.input}")

    os.makedirs(args.output_dir, exist_ok=True)
    z_m = args.slice_z_cm / 100.0
    z_tag = f"z{z_m:.1f}m".replace(".", "p")

    plot_specs = [
        (
            "observed_voxels",
            f"data-seeded ray occupancy kernel, xy slice at z = {z_m:.1f} m",
            f"acceptance_kernel_observed_xy_{z_tag}.png",
            True,
            False,
        ),
        (
            "expected_voxels_global_scaled",
            f"EXPACS solid-ground expectation on same rays, xy slice at z = {z_m:.1f} m",
            f"acceptance_kernel_expected_global_xy_{z_tag}.png",
            True,
            False,
        ),
        (
            "smooth_calibration_factor",
            f"smooth null calibration field, xy slice at z = {z_m:.1f} m",
            f"acceptance_kernel_smooth_calibration_xy_{z_tag}.png",
            False,
            False,
        ),
        (
            "cavity_excess_significance_calibrated",
            f"cavity-oriented Li-Ma residual, xy slice at z = {z_m:.1f} m",
            f"acceptance_kernel_cavity_lima_xy_{z_tag}.png",
            False,
            False,
        ),
    ]

    canvases = []
    for hist_name, title, filename, logz, symmetric in plot_specs:
        hist3 = infile.Get(hist_name)
        if not hist3:
            raise RuntimeError(f"Missing histogram {hist_name} in {args.input}")
        hist2 = xy_slice(ROOT, hist3, hist_name, args.slice_z_cm)
        output_path = os.path.join(args.output_dir, filename)
        canvases.append(draw_plot(ROOT, hist2, output_path, title, logz=logz, symmetric=symmetric))
        print(output_path)

    infile.Close()
    return canvases


if __name__ == "__main__":
    main()
