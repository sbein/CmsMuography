#!/usr/bin/env python3
"""Make slice-stack displays for Li-Ma/significance voxel volumes."""

import argparse
import json
import math
import os
import re
from array import array

from suppress_voxel_pedestals import (
    CAVERN_HALL_ARCH_SPRING_Z_CM,
    CAVERN_HALL_BOX_CM,
    CAVERN_SHAFT_CM,
    DEFAULT_GROUND_COEFFICIENTS_CM,
    axis,
    cavern_cm_to_hist_units_point,
    cm_to_hist_units,
    cms_x_bearing_degrees,
    clone_hist_with_scaled_axes,
    copy_axis_style,
    draw_cavern_outline_2d,
    draw_cavern_outline_3d,
    draw_cms_outline_2d,
    draw_cms_outline_3d,
    draw_compass_rose_3d,
    draw_ground_surface_3d,
    draw_label,
    draw_orientation_overlay_2d,
    ground_color,
    ground_z_cavern_cm,
    hist_axis_unit,
    hist_coordinate_system,
    hist_units_point_to_cavern_cm,
    make_slice,
    require_th3,
    vertical_axis,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Copy Li-Ma input/output TH3s into a display ROOT file and build "
            "slice-stack canvases instead of threshold/projection comparison grids."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input Li-Ma ROOT file.")
    parser.add_argument("-o", "--output", required=True, help="Output display ROOT file.")
    parser.add_argument(
        "--hist-names",
        default=(
            "observed_voxels,expected_voxels_global_scaled,calibrated_expected_voxels,"
            "cavity_excess_significance_calibrated,signed_significance_calibrated,"
            "log_transmission_excess"
        ),
        help="Comma-separated TH3s to draw as slice stacks.",
    )
    parser.add_argument(
        "--copy-hists",
        default=(
            "observed_voxels,expected_voxels_global_scaled,calibrated_expected_voxels,"
            "smooth_calibration_factor,cavity_excess_significance_calibrated,"
            "signed_significance_calibrated,log_transmission_excess"
        ),
        help="Comma-separated TH3s to copy into the output inputs directory.",
    )
    parser.add_argument(
        "--z-slices-cm",
        default="250,750,1250,1750,2500,3500,5000,7500,10250",
        help="Comma-separated cavern z/up slice values for XY panels.",
    )
    parser.add_argument(
        "--y-slices-cm",
        default="-2500,-1750,-1000,-500,0,500,1000,1750,2500",
        help="Comma-separated cavern y/beam slice values for XZ panels.",
    )
    parser.add_argument("--canvas-columns", type=int, default=3)
    parser.add_argument("--canvas-rows", type=int, default=3)
    parser.add_argument("--canvas-width", type=int, default=3600)
    parser.add_argument("--canvas-height", type=int, default=3600)
    parser.add_argument(
        "--plot-units",
        choices=("m", "cm"),
        default="m",
        help="Length units used for copied display histograms and canvases.",
    )
    parser.add_argument(
        "--display-rebin",
        type=int,
        default=2,
        help=(
            "Coarsen display TH3s by this integer factor along x, y, and z before "
            "writing plots. Counts are summed; intensive maps are averaged."
        ),
    )
    parser.add_argument("--linear-color", action="store_true", help="Disable log color axes.")
    parser.add_argument(
        "--linear-hists",
        default=(
            "signed_significance,signed_significance_calibrated,log_transmission_excess,"
            "smooth_calibration_factor,unet_cavity_probability,unet_cavity_score,"
            "unet_density_multiplier,unet_density_deficit,unet_clustered_density_deficit,"
            "unet_density_deficit_cluster_mask,unet_density_deficit_cluster_id"
        ),
        help="Comma-separated histogram names that should use linear color even when log color is otherwise enabled.",
    )
    parser.add_argument(
        "--symmetric-hists",
        default="signed_significance,signed_significance_calibrated,log_transmission_excess,unet_density_residual",
        help="Comma-separated histogram names whose 2D slice colors should be symmetric around zero.",
    )
    parser.add_argument(
        "--log-3d-color",
        dest="log_3d_color",
        action="store_true",
        default=True,
        help="Use log color on positive 3D BOX2Z canvases (default).",
    )
    parser.add_argument(
        "--linear-3d-color",
        dest="log_3d_color",
        action="store_false",
        help="Use linear color on all 3D BOX2Z canvases.",
    )
    parser.add_argument("--log-min", type=float, default=1.0e-3)
    parser.add_argument("--no-3d", action="store_true", help="Do not draw any 3D canvases.")
    parser.add_argument(
        "--no-geometry-baseline",
        action="store_true",
        help="Do not draw the empty analysis-volume/CMS/cavern reference canvas.",
    )
    parser.add_argument("--draw-3d-no-cut", action="store_true", help="Also draw the heavy no-cut 3D canvas.")
    parser.add_argument("--no-3d-grid", action="store_true", help="Do not draw the 3D threshold-grid canvas.")
    parser.add_argument(
        "--three-d-top-fractions",
        default="0.30,0.20,0.12,0.08,0.05,0.03,0.02,0.01,0.005",
        help="Comma-separated positive-bin fractions retained in the 3D threshold grid.",
    )
    parser.add_argument(
        "--comparison-hists",
        default="observed_voxels,expected_voxels_global_scaled,calibrated_expected_voxels",
        help=(
            "Comma-separated TH3 names to draw side by side in one 3D comparison "
            "canvas; for example observed_voxels,expected_voxels_global_scaled,"
            "calibrated_expected_voxels."
        ),
    )
    parser.add_argument(
        "--comparison-top-fraction",
        type=float,
        default=0.05,
        help="Positive-voxel fraction retained independently in each 3D comparison panel.",
    )
    parser.add_argument("--no-cms-outline", action="store_true")
    parser.add_argument("--no-cavern-outline", action="store_true")
    parser.add_argument(
        "--no-cms-outline-2d",
        action="store_true",
        help="Keep CMS wireframe on optional 3D canvases but omit it from 2D slice-stack panels.",
    )
    parser.add_argument("--no-orientation-overlay", action="store_true")
    parser.add_argument("--cms-x-bearing-deg", type=float, default=None)
    parser.add_argument("--no-ground-overlay", action="store_true")
    parser.add_argument(
        "--ground-coefficients-cm",
        default=",".join(str(value) for value in DEFAULT_GROUND_COEFFICIENTS_CM),
    )
    parser.add_argument("--ground-grid-step-cm", type=float, default=500.0)
    parser.add_argument(
        "--png-dir",
        help="Optional directory where canvases are also exported as PNG.",
    )
    return parser.parse_args()


def split_csv(text):
    return [item.strip() for item in text.split(",") if item.strip()]


def split_float_csv(text):
    return [float(item) for item in split_csv(text)]


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_")


def value_token(value):
    token = f"{value:g}".replace("-", "m").replace(".", "p")
    return token


QUANTITY_INFO = {
    "observed_voxels": {
        "short": "Observed muon crossings",
        "detail": "observed track-painted voxel crossings",
        "symbol": "Nobs",
    },
    "expected_voxels_scaled": {
        "short": "Uniform-ground null crossings",
        "detail": "uniform-ground expected voxel crossings, globally scaled",
        "symbol": "Nexp",
    },
    "expected_voxels_global_scaled": {
        "short": "Uniform-ground null crossings",
        "detail": "uniform-ground expected voxel crossings, globally scaled",
        "symbol": "Nexp",
    },
    "calibrated_expected_voxels": {
        "short": "Smooth-null crossings",
        "detail": "smooth-calibrated expected voxel crossings",
        "symbol": "Nnull",
    },
    "signed_significance": {
        "short": "Signed Li-Ma Z",
        "detail": "signed data-minus-null Li-Ma significance",
        "symbol": "Zsigned",
    },
    "signed_significance_calibrated": {
        "short": "Signed Li-Ma Z",
        "detail": "signed data-minus-smooth-null Li-Ma significance",
        "symbol": "Zsigned",
    },
    "cavity_excess_significance": {
        "short": "Cavity Li-Ma Z",
        "detail": "positive data excess over the null; cavity-like Li-Ma significance",
        "symbol": "Zcav",
    },
    "cavity_excess_significance_calibrated": {
        "short": "Cavity Li-Ma Z",
        "detail": "positive data excess over the smooth null; cavity-like Li-Ma significance",
        "symbol": "Zcav",
    },
    "observed_over_expected": {
        "short": "Data/null ratio",
        "detail": "observed crossings divided by expected crossings",
        "symbol": "Nobs/Nexp",
    },
    "observed_over_expected_calibrated": {
        "short": "Data/smooth-null ratio",
        "detail": "observed crossings divided by smooth-calibrated expected crossings",
        "symbol": "Nobs/Nnull",
    },
    "log_transmission_excess": {
        "short": "Log transmission excess",
        "detail": "ln((observed+pseudo)/(smooth null+pseudo)); positive is cavity-like",
        "symbol": "lnT",
    },
    "smooth_calibration_factor": {
        "short": "Smooth calibration factor",
        "detail": "smooth multiplicative null correction C(z)*Cxy(x,y)",
        "symbol": "Cnull",
    },
    "unet_cavity_probability": {
        "short": "U-Net cavity probability",
        "detail": "trained U-Net cavity-probability score applied to real-data feature volume",
        "symbol": "P(cavity)",
    },
    "unet_cavity_score": {
        "short": "U-Net cavity-deficit score",
        "detail": "derived U-Net cavity-deficit score 1-rho/rho0",
        "symbol": "1-rho/rho0",
    },
    "unet_density_deficit": {
        "short": "U-Net density deficit",
        "detail": "trained U-Net inferred density deficit score 1-rho/rho0",
        "symbol": "1-rho/rho0",
    },
    "unet_clustered_density_deficit": {
        "short": "Clustered U-Net density deficit",
        "detail": "U-Net density-deficit voxels retained in loose-quantile connected clusters",
        "symbol": "clustered 1-rho/rho0",
    },
    "unet_density_multiplier": {
        "short": "U-Net density multiplier",
        "detail": "trained U-Net inferred relative density multiplier rho/rho0",
        "symbol": "rho/rho0",
    },
    "unet_density_residual": {
        "short": "U-Net density residual",
        "detail": "trained U-Net inferred density residual rho/rho0 - 1",
        "symbol": "rho/rho0-1",
    },
    "unet_density_deficit_cluster_mask": {
        "short": "U-Net deficit cluster mask",
        "detail": "mask of density-deficit voxels retained in connected clusters",
        "symbol": "cluster mask",
    },
    "unet_density_deficit_cluster_id": {
        "short": "U-Net deficit cluster ID",
        "detail": "connected-component label for retained density-deficit clusters",
        "symbol": "cluster id",
    },
}


def quantity_info(hist_name):
    base_name = hist_name.split("/")[-1]
    fallback = base_name.replace("_", " ")
    return QUANTITY_INFO.get(
        base_name,
        {
            "short": fallback,
            "detail": fallback,
            "symbol": base_name,
        },
    )


def quantity_short(hist_name):
    return quantity_info(hist_name)["short"]


def quantity_detail(hist_name):
    return quantity_info(hist_name)["detail"]


def quantity_symbol(hist_name):
    return quantity_info(hist_name)["symbol"]


def root_title_text(text):
    return text.replace(";", ",")


def axis_range(axis_obj):
    return axis_obj.GetXmin(), axis_obj.GetXmax()


def clipped_values(values, axis_obj):
    low, high = axis_range(axis_obj)
    return [min(max(value, low), high) for value in values]


def display_rebin_mode(hist_name):
    base_name = hist_name.split("/")[-1].lower()
    summed_names = {
        "observed_voxels",
        "expected_voxels_scaled",
        "expected_voxels_global_scaled",
        "calibrated_expected_voxels",
        "solid_ground_expected_voxels",
        "generator_expected_voxels",
    }
    if base_name in summed_names or base_name.endswith("_voxels"):
        return "sum"
    return "mean"


def rebin_display_th3(root, hist, hist_name, factor):
    factor = int(factor)
    if factor <= 1:
        return hist
    if not hist.InheritsFrom("TH3"):
        return hist
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    nz = hist.GetNbinsZ()
    if nx % factor or ny % factor or nz % factor:
        raise RuntimeError(
            f"Cannot rebin {hist_name} by {factor}: dimensions are {nx}x{ny}x{nz}"
        )

    rebinned = root.TH3F(
        f"{hist.GetName()}_rebin{factor}",
        hist.GetTitle(),
        nx // factor,
        hist.GetXaxis().GetXmin(),
        hist.GetXaxis().GetXmax(),
        ny // factor,
        hist.GetYaxis().GetXmin(),
        hist.GetYaxis().GetXmax(),
        nz // factor,
        hist.GetZaxis().GetXmin(),
        hist.GetZaxis().GetXmax(),
    )
    rebinned.SetDirectory(0)
    copy_axis_style(hist.GetXaxis(), rebinned.GetXaxis(), hist_axis_unit(hist))
    copy_axis_style(hist.GetYaxis(), rebinned.GetYaxis(), hist_axis_unit(hist))
    copy_axis_style(hist.GetZaxis(), rebinned.GetZaxis(), hist_axis_unit(hist))

    mode = display_rebin_mode(hist_name)
    divisor = float(factor * factor * factor)
    if hist.GetSumw2N() > 0:
        rebinned.Sumw2()
    total_entries = 0.0
    for ix in range(1, rebinned.GetNbinsX() + 1):
        source_x0 = (ix - 1) * factor + 1
        for iy in range(1, rebinned.GetNbinsY() + 1):
            source_y0 = (iy - 1) * factor + 1
            for iz in range(1, rebinned.GetNbinsZ() + 1):
                source_z0 = (iz - 1) * factor + 1
                value_sum = 0.0
                error2_sum = 0.0
                for dx in range(factor):
                    for dy in range(factor):
                        for dz in range(factor):
                            source_bin = hist.GetBin(source_x0 + dx, source_y0 + dy, source_z0 + dz)
                            value_sum += float(hist.GetBinContent(source_bin))
                            if hist.GetSumw2N() > 0:
                                error = float(hist.GetBinError(source_bin))
                                error2_sum += error * error
                if mode == "mean":
                    value = value_sum / divisor
                    error = math.sqrt(error2_sum) / divisor if hist.GetSumw2N() > 0 else 0.0
                else:
                    value = value_sum
                    error = math.sqrt(error2_sum) if hist.GetSumw2N() > 0 else 0.0
                if value != 0.0:
                    rebinned.SetBinContent(ix, iy, iz, value)
                    total_entries += value
                if hist.GetSumw2N() > 0 and error > 0.0:
                    rebinned.SetBinError(ix, iy, iz, error)
    rebinned.SetEntries(total_entries if mode == "sum" else hist.GetEntries())
    return rebinned


def positive_minmax_2d(hist):
    minimum = None
    maximum = None
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            value = float(hist.GetBinContent(ix, iy))
            if value <= 0.0:
                continue
            minimum = value if minimum is None else min(minimum, value)
            maximum = value if maximum is None else max(maximum, value)
    return minimum, maximum


def positive_minmax_3d(hist):
    minimum = None
    maximum = None
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            for iz in range(1, hist.GetNbinsZ() + 1):
                value = float(hist.GetBinContent(ix, iy, iz))
                if value <= 0.0:
                    continue
                minimum = value if minimum is None else min(minimum, value)
                maximum = value if maximum is None else max(maximum, value)
    return minimum, maximum


def minmax_3d(hist):
    minimum = None
    maximum = None
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            for iz in range(1, hist.GetNbinsZ() + 1):
                value = float(hist.GetBinContent(ix, iy, iz))
                minimum = value if minimum is None else min(minimum, value)
                maximum = value if maximum is None else max(maximum, value)
    return minimum, maximum


def minmax_2d(hist):
    minimum = None
    maximum = None
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            value = float(hist.GetBinContent(ix, iy))
            minimum = value if minimum is None else min(minimum, value)
            maximum = value if maximum is None else max(maximum, value)
    return minimum, maximum


def should_use_log_color(hist_name, args):
    base_name = hist_name.split("/")[-1]
    return not args.linear_color and base_name not in set(split_csv(args.linear_hists))


def configure_color_axis(pad, hist, hist_name, args, is_3d=False):
    symmetric_names = set(split_csv(args.symmetric_hists))
    base_name = hist_name.split("/")[-1]
    if base_name in symmetric_names and not is_3d:
        minimum, maximum = minmax_2d(hist)
        if minimum is not None and maximum is not None:
            bound = max(abs(minimum), abs(maximum))
            if bound > 0.0:
                hist.SetMinimum(-bound)
                hist.SetMaximum(bound)
    if not should_use_log_color(hist_name, args):
        return False
    positive_min, positive_max = positive_minmax_3d(hist) if is_3d else positive_minmax_2d(hist)
    if positive_min is None or positive_max is None or positive_max <= 0.0:
        return False
    hist.SetMinimum(max(args.log_min, positive_min))
    if is_3d:
        return False
    pad.SetLogz(True)
    return True


def color_range_3d(hist, hist_name, args):
    if should_use_log_color(hist_name, args):
        positive_min, positive_max = positive_minmax_3d(hist)
        if positive_min is None or positive_max is None or positive_max <= 0.0:
            return None
        return max(args.log_min, positive_min), positive_max

    minimum, maximum = minmax_3d(hist)
    if minimum is None or maximum is None:
        return None
    base_name = hist_name.split("/")[-1]
    if base_name in set(split_csv(args.symmetric_hists)):
        bound = max(abs(minimum), abs(maximum))
        if bound > 0.0:
            return -bound, bound
    if maximum <= minimum:
        padding = max(0.05, 0.05 * max(abs(minimum), 1.0))
        return minimum - padding, maximum + padding
    return minimum, maximum


def combined_color_range_3d(hists, hist_names, args):
    log_flags = [should_use_log_color(name, args) for name in hist_names]
    if not log_flags or any(flag != log_flags[0] for flag in log_flags):
        return None
    ranges = [color_range_3d(hist, name, args) for hist, name in zip(hists, hist_names)]
    ranges = [item for item in ranges if item is not None]
    if not ranges:
        return None
    if log_flags[0]:
        return min(item[0] for item in ranges), max(item[1] for item in ranges)
    if any(name.split("/")[-1] in set(split_csv(args.symmetric_hists)) for name in hist_names):
        bound = max(max(abs(item[0]), abs(item[1])) for item in ranges)
        return -bound, bound
    return min(item[0] for item in ranges), max(item[1] for item in ranges)


def apply_log10_color_transform_3d(hist, args, original_range=None):
    if original_range is None:
        original_range = color_range_3d(hist, hist.GetName(), args)
    if original_range is None:
        return None
    positive_min, positive_max = original_range
    if positive_min is None or positive_max is None or positive_max <= 0.0:
        return None
    floor = max(args.log_min, positive_min)
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            for iz in range(1, hist.GetNbinsZ() + 1):
                value = float(hist.GetBinContent(ix, iy, iz))
                if value > 0.0:
                    hist.SetBinContent(ix, iy, iz, math.log10(max(value, floor)))
    log_minimum = math.log10(floor)
    log_maximum = math.log10(max(positive_max, floor))
    if log_maximum <= log_minimum:
        # ROOT draws no BOX2Z cells when a sparse threshold leaves a single
        # positive value and the displayed color range collapses to a point.
        padding = max(0.05, 0.05 * max(abs(log_minimum), 1.0))
        log_minimum -= padding
        log_maximum += padding
    hist.SetMinimum(log_minimum)
    hist.SetMaximum(log_maximum)
    return floor, max(positive_max, floor)


def maybe_apply_3d_log_color_transform(hist, hist_name, args, original_range=None):
    if not args.log_3d_color or not should_use_log_color(hist_name, args):
        return None
    return apply_log10_color_transform_3d(hist, args, original_range)


def draw_3d_log_palette_axis(root, pad, hist, original_range):
    if not original_range:
        return []
    low, high = original_range
    if low <= 0.0 or high <= low:
        return []
    pad.Update()
    palette = hist.GetListOfFunctions().FindObject("palette")
    if not palette:
        return []
    palette.SetLabelSize(0.0)
    palette.SetTickLength(0.0)
    palette.GetAxis().SetLabelSize(0.0)
    palette.GetAxis().SetTickLength(0.0)

    x_axis = max(palette.GetX1NDC() - 0.004, 0.0)
    axis_obj = root.TGaxis(
        x_axis,
        palette.GetY1NDC(),
        x_axis,
        palette.GetY2NDC(),
        low,
        high,
        506,
        "-G",
    )
    axis_obj.SetNDC(True)
    axis_obj.SetLabelFont(42)
    axis_obj.SetLabelSize(0.022)
    axis_obj.SetLabelOffset(0.004)
    axis_obj.SetTickSize(0.010)
    axis_obj.SetMoreLogLabels(True)
    axis_obj.Draw()
    return [axis_obj]


def draw_3d_panel_labels(root, hist_name, selection_label, log_color_range):
    color_mode = "log-scale color" if log_color_range else "linear color"
    symbol = quantity_symbol(hist_name)
    return [
        draw_label(
            root,
            f"{quantity_short(hist_name)}; {color_mode}: {symbol}",
            size=0.036,
            x=0.12,
            y=0.935,
        ),
        draw_label(
            root,
            selection_label,
            size=0.027,
            x=0.12,
            y=0.895,
            font=42,
        ),
    ]


def draw_3d_threshold_label(root, hist_name, threshold, fraction):
    return draw_label(
        root,
        f"Top {100.0 * fraction:g}%: {quantity_symbol(hist_name)} > {threshold:.3g}",
        size=0.036,
        x=0.12,
        y=0.935,
    )


def coordinate_value_to_cavern_horizontal(source_hist, coordinate_system, axis_index, value, fixed_axis, fixed_value):
    point = [0.0, 0.0, 0.0]
    point[axis_index] = value
    point[fixed_axis] = fixed_value
    cavern_x, cavern_y, _cavern_z = hist_units_point_to_cavern_cm(
        source_hist, coordinate_system, tuple(point)
    )
    return cavern_x, cavern_y


def draw_ground_vertical_slice(root, hist2, source_hist, axis_a, axis_b, axis_c, fixed_value, args):
    coordinate_system = hist_coordinate_system(source_hist)
    up_axis = vertical_axis(coordinate_system)
    if up_axis not in (axis_a, axis_b):
        return []
    horizontal_axis = axis_b if axis_a == up_axis else axis_a
    color = ground_color(root)
    horizontal_range = axis(source_hist, horizontal_axis)
    npoints = 181
    graph = root.TGraph(npoints)
    for ipoint in range(npoints):
        value = horizontal_range.GetXmin() + (
            horizontal_range.GetXmax() - horizontal_range.GetXmin()
        ) * ipoint / (npoints - 1)
        cavern_x, cavern_y = coordinate_value_to_cavern_horizontal(
            source_hist, coordinate_system, horizontal_axis, value, axis_c, fixed_value
        )
        hist_point = cavern_cm_to_hist_units_point(
            source_hist,
            coordinate_system,
            cavern_x,
            cavern_y,
            ground_z_cavern_cm(args, cavern_x, cavern_y),
        )
        graph.SetPoint(ipoint, hist_point[axis_a], hist_point[axis_b])
    graph.SetLineColor(color)
    graph.SetLineWidth(3)
    graph.SetLineStyle(2)
    graph.Draw("L SAME")
    return [graph]


def draw_ground_horizontal_slice(root, source_hist, axis_a, axis_b, axis_c, fixed_value, args):
    coordinate_system = hist_coordinate_system(source_hist)
    up_axis = vertical_axis(coordinate_system)
    if up_axis in (axis_a, axis_b) or axis_c != up_axis:
        return []
    color = ground_color(root)
    ax = axis(source_hist, axis_a)
    ay = axis(source_hist, axis_b)
    contours = root.TH2F(
        f"ground_surface_at_z_{value_token(fixed_value)}_{axis_a}_{axis_b}",
        "",
        50,
        ax.GetXmin(),
        ax.GetXmax(),
        50,
        ay.GetXmin(),
        ay.GetXmax(),
    )
    low = None
    high = None
    for ix in range(1, contours.GetNbinsX() + 1):
        x_value = contours.GetXaxis().GetBinCenter(ix)
        for iy in range(1, contours.GetNbinsY() + 1):
            y_value = contours.GetYaxis().GetBinCenter(iy)
            point = [0.0, 0.0, 0.0]
            point[axis_a] = x_value
            point[axis_b] = y_value
            cavern_x, cavern_y, _cavern_z = hist_units_point_to_cavern_cm(
                source_hist, coordinate_system, tuple(point)
            )
            ground_z = cm_to_hist_units(source_hist, ground_z_cavern_cm(args, cavern_x, cavern_y))
            contours.SetBinContent(ix, iy, ground_z)
            low = ground_z if low is None else min(low, ground_z)
            high = ground_z if high is None else max(high, ground_z)
    if low is None or high is None or fixed_value < low or fixed_value > high:
        return []
    levels = array("d", [fixed_value])
    contours.SetContour(1, levels)
    contours.SetLineColor(color)
    contours.SetLineWidth(3)
    contours.SetLineStyle(2)
    contours.SetStats(False)
    contours.Draw("CONT3 SAME")
    return [contours]


def draw_ground_for_slice(root, hist2, source_hist, axis_a, axis_b, axis_c, fixed_value, args):
    if args.no_ground_overlay:
        return []
    objects = draw_ground_vertical_slice(root, hist2, source_hist, axis_a, axis_b, axis_c, fixed_value, args)
    if objects:
        return objects
    return draw_ground_horizontal_slice(root, source_hist, axis_a, axis_b, axis_c, fixed_value, args)


def draw_slice_panel(root, pad, hist2, source_hist, hist_name, axis_a, axis_b, axis_c, fixed_value, label, args):
    pad.SetRightMargin(0.16)
    pad.SetLeftMargin(0.12)
    pad.SetBottomMargin(0.11)
    pad.SetTopMargin(0.08)
    configure_color_axis(pad, hist2, hist_name, args, is_3d=False)
    hist2.SetStats(False)
    hist2.Draw("COLZ")
    keepalive = []
    coordinate_system = hist_coordinate_system(source_hist)
    if not args.no_cavern_outline:
        keepalive.extend(draw_cavern_outline_2d(root, source_hist, axis_a, axis_b, axis_c, fixed_value))
    if not args.no_cms_outline and not args.no_cms_outline_2d:
        keepalive.extend(draw_cms_outline_2d(root, source_hist, axis_a, axis_b, axis_c, fixed_value))
    keepalive.extend(draw_ground_for_slice(root, hist2, source_hist, axis_a, axis_b, axis_c, fixed_value, args))
    if not args.no_orientation_overlay:
        keepalive.extend(
            draw_orientation_overlay_2d(
                root,
                coordinate_system,
                axis_a,
                axis_b,
                cms_x_bearing_degrees(args),
            )
        )
    keepalive.append(draw_label(root, label, size=0.045))
    return keepalive


def draw_3d_canvas(root, hist3, hist_name, args):
    canvas = root.TCanvas(
        f"canvas_{safe_name(hist_name)}_3d_no_cut",
        f"{quantity_short(hist_name)} 3D no cut",
        args.canvas_width,
        args.canvas_height,
    )
    pad = canvas.cd()
    pad.SetRightMargin(0.08)
    pad.SetLeftMargin(0.10)
    pad.SetBottomMargin(0.10)
    pad.SetTopMargin(0.08)
    pad.SetTheta(25)
    pad.SetPhi(35)
    draw_hist = hist3.Clone(f"{safe_name(hist_name)}_3d_draw")
    draw_hist.SetDirectory(0)
    draw_hist.SetTitle("")
    for iaxis in range(3):
        ax = axis(draw_hist, iaxis)
        ax.SetLabelSize(0.018)
        ax.SetTitleSize(0.022)
        ax.SetTitleOffset(2.0)
    log_color_range = maybe_apply_3d_log_color_transform(draw_hist, hist_name, args)
    draw_hist.SetStats(False)
    draw_hist.Draw("BOX2Z")
    coordinate_system = hist_coordinate_system(draw_hist)
    keepalive = []
    keepalive.append(draw_hist)
    keepalive.extend(draw_3d_log_palette_axis(root, pad, draw_hist, log_color_range))
    keepalive.extend(draw_ground_surface_3d(root, draw_hist, coordinate_system, args))
    if not args.no_cavern_outline:
        keepalive.extend(draw_cavern_outline_3d(root, draw_hist))
    if not args.no_cms_outline:
        keepalive.extend(draw_cms_outline_3d(root, draw_hist))
    if not args.no_orientation_overlay:
        keepalive.extend(draw_compass_rose_3d(root, draw_hist, cms_x_bearing_degrees(args)))
    keepalive.extend(draw_3d_panel_labels(root, hist_name, "no threshold", log_color_range))
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def draw_geometry_baseline_canvas(root, hist3, args):
    canvas = root.TCanvas(
        "canvas_geometry_baseline_3d",
        "CMS cavern muography analysis volume and reference geometry",
        args.canvas_width,
        args.canvas_height,
    )
    pad = canvas.cd()
    pad.SetRightMargin(0.06)
    pad.SetLeftMargin(0.10)
    pad.SetBottomMargin(0.10)
    pad.SetTopMargin(0.08)
    pad.SetTheta(25)
    pad.SetPhi(35)

    frame_hist = hist3.Clone("geometry_baseline_frame")
    frame_hist.SetDirectory(0)
    frame_hist.Reset("ICES")
    frame_hist.SetTitle("")
    configure_3d_axes(frame_hist)
    frame_hist.SetStats(False)
    frame_hist.SetMinimum(0.0)
    frame_hist.SetMaximum(1.0)
    frame_hist.Draw("BOX")

    coordinate_system = hist_coordinate_system(frame_hist)
    keepalive = [frame_hist]
    keepalive.extend(draw_ground_surface_3d(root, frame_hist, coordinate_system, args))
    if not args.no_cavern_outline:
        keepalive.extend(draw_cavern_outline_3d(root, frame_hist))
    if not args.no_cms_outline:
        keepalive.extend(draw_cms_outline_3d(root, frame_hist))
    if not args.no_orientation_overlay:
        keepalive.extend(draw_compass_rose_3d(root, frame_hist, cms_x_bearing_degrees(args)))
    keepalive.append(
        draw_label(
            root,
            "CMS cavern muography reference geometry",
            size=0.036,
            x=0.12,
            y=0.935,
        )
    )
    keepalive.append(
        draw_label(
            root,
            "empty 3D voxel volume; detector, cavern, shaft, ground, and compass overlays",
            size=0.027,
            x=0.12,
            y=0.895,
            font=42,
        )
    )
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def positive_values_3d(hist):
    values = []
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            for iz in range(1, hist.GetNbinsZ() + 1):
                value = float(hist.GetBinContent(ix, iy, iz))
                if value > 0.0:
                    values.append(value)
    values.sort()
    return values


def thresholds_for_top_fractions(values, fractions):
    thresholds = []
    if not values:
        return thresholds
    for fraction in fractions:
        fraction = min(max(fraction, 1.0e-12), 1.0)
        target_count = max(1, int(math.ceil(fraction * len(values))))
        thresholds.append(values[-target_count])
    return thresholds


def make_thresholded_3d_hist(hist3, name, threshold):
    clone = hist3.Clone(name)
    clone.SetDirectory(0)
    clone.Reset("ICES")
    total = 0.0
    for ix in range(1, hist3.GetNbinsX() + 1):
        for iy in range(1, hist3.GetNbinsY() + 1):
            for iz in range(1, hist3.GetNbinsZ() + 1):
                value = float(hist3.GetBinContent(ix, iy, iz))
                if value >= threshold and value > 0.0:
                    clone.SetBinContent(ix, iy, iz, value)
                    total += value
    clone.SetEntries(total)
    return clone


def configure_3d_axes(hist):
    for iaxis in range(3):
        ax = axis(hist, iaxis)
        ax.SetLabelSize(0.018)
        ax.SetTitleSize(0.022)
        ax.SetTitleOffset(2.0)


def draw_3d_threshold_panel(
    root,
    pad,
    hist3,
    hist_name,
    threshold,
    fraction,
    index,
    args,
    shared_color_range=None,
):
    pad.SetRightMargin(0.08)
    pad.SetLeftMargin(0.10)
    pad.SetBottomMargin(0.10)
    pad.SetTopMargin(0.08)
    pad.SetTheta(25)
    pad.SetPhi(35)
    draw_hist = make_thresholded_3d_hist(hist3, f"{safe_name(hist_name)}_3d_cut_{index}", threshold)
    draw_hist.SetTitle("")
    configure_3d_axes(draw_hist)
    log_color_range = maybe_apply_3d_log_color_transform(
        draw_hist,
        hist_name,
        args,
        shared_color_range,
    )
    if log_color_range is None and shared_color_range is not None:
        draw_hist.SetMinimum(shared_color_range[0])
        draw_hist.SetMaximum(shared_color_range[1])
    draw_hist.SetStats(False)
    draw_hist.Draw("BOX2Z")
    coordinate_system = hist_coordinate_system(draw_hist)
    keepalive = [draw_hist]
    keepalive.extend(draw_3d_log_palette_axis(root, pad, draw_hist, log_color_range))
    keepalive.extend(draw_ground_surface_3d(root, draw_hist, coordinate_system, args))
    if not args.no_cavern_outline:
        keepalive.extend(draw_cavern_outline_3d(root, draw_hist))
    if not args.no_cms_outline:
        keepalive.extend(draw_cms_outline_3d(root, draw_hist))
    if not args.no_orientation_overlay:
        keepalive.extend(draw_compass_rose_3d(root, draw_hist, cms_x_bearing_degrees(args)))
    keepalive.append(draw_3d_threshold_label(root, hist_name, threshold, fraction))
    return keepalive


def draw_3d_threshold_grid_canvas(root, hist3, hist_name, args):
    fractions = split_float_csv(args.three_d_top_fractions)
    values = positive_values_3d(hist3)
    thresholds = thresholds_for_top_fractions(values, fractions)
    canvas = root.TCanvas(
        f"canvas_{safe_name(hist_name)}_3d_threshold_grid",
        f"{quantity_short(hist_name)} 3D threshold grid",
        args.canvas_width,
        args.canvas_height,
    )
    canvas.Divide(args.canvas_columns, args.canvas_rows, 0.001, 0.001)
    keepalive = []
    shared_color_range = color_range_3d(hist3, hist_name, args)
    max_pads = args.canvas_columns * args.canvas_rows
    for ipad, (fraction, threshold) in enumerate(zip(fractions[:max_pads], thresholds[:max_pads]), start=1):
        keepalive.extend(
            draw_3d_threshold_panel(
                root,
                canvas.cd(ipad),
                hist3,
                hist_name,
                threshold,
                fraction,
                ipad,
                args,
                shared_color_range,
            )
        )
    if not thresholds:
        canvas.cd(1)
        keepalive.append(draw_label(root, f"{quantity_short(hist_name)}: no positive bins", size=0.055))
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def draw_3d_comparison_canvas(root, hists, hist_names, args):
    fraction = min(max(args.comparison_top_fraction, 1.0e-12), 1.0)
    comparison_name = "__".join(safe_name(name) for name in hist_names)
    canvas = root.TCanvas(
        f"canvas_comparison_{comparison_name}_3d",
        "3D voxel comparison",
        args.canvas_width,
        max(1200, args.canvas_height // 2),
    )
    canvas.Divide(len(hists), 1, 0.001, 0.001)
    keepalive = []
    shared_color_range = combined_color_range_3d(hists, hist_names, args)
    for index, (hist, hist_name) in enumerate(zip(hists, hist_names), start=1):
        values = positive_values_3d(hist)
        thresholds = thresholds_for_top_fractions(values, [fraction])
        if not thresholds:
            canvas.cd(index)
            keepalive.append(
                draw_label(root, f"{quantity_short(hist_name)}: no positive bins", size=0.055)
            )
            continue
        keepalive.extend(
            draw_3d_threshold_panel(
                root,
                canvas.cd(index),
                hist,
                hist_name,
                thresholds[0],
                fraction,
                index,
                args,
                shared_color_range,
            )
        )
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def percentile(values, fraction):
    if not values:
        return 0.0
    ordered = sorted(float(value) for value in values)
    position = min(max(float(fraction), 0.0), 1.0) * (len(ordered) - 1)
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def make_xy_log_ratio(root, observed, expected, name, pseudo_count):
    hist = root.TH2F(
        name,
        "",
        observed.GetNbinsX(),
        observed.GetXaxis().GetXmin(),
        observed.GetXaxis().GetXmax(),
        observed.GetNbinsY(),
        observed.GetYaxis().GetXmin(),
        observed.GetYaxis().GetXmax(),
    )
    hist.SetDirectory(0)
    copy_axis_style(observed.GetXaxis(), hist.GetXaxis(), hist_axis_unit(observed))
    copy_axis_style(observed.GetYaxis(), hist.GetYaxis(), hist_axis_unit(observed))
    hist.GetZaxis().SetTitle("ln[(data + 0.5)/(MC + 0.5)]")
    for ix in range(1, observed.GetNbinsX() + 1):
        for iy in range(1, observed.GetNbinsY() + 1):
            observed_sum = 0.0
            expected_sum = 0.0
            for iz in range(1, observed.GetNbinsZ() + 1):
                observed_sum += float(observed.GetBinContent(ix, iy, iz))
                expected_sum += float(expected.GetBinContent(ix, iy, iz))
            hist.SetBinContent(
                ix,
                iy,
                math.log((observed_sum + pseudo_count) / (expected_sum + pseudo_count)),
            )
    return hist


def make_z_ratio_hist(root, observed, expected, name, config):
    hist = root.TH1F(
        name,
        "",
        observed.GetNbinsZ(),
        observed.GetZaxis().GetXmin(),
        observed.GetZaxis().GetXmax(),
    )
    hist.SetDirectory(0)
    copy_axis_style(observed.GetZaxis(), hist.GetXaxis(), hist_axis_unit(observed))
    hist.GetYaxis().SetTitle("ratio")
    unit_to_cm = 100.0 if hist_axis_unit(observed) == "m" else 1.0
    exclude_core = bool(config.get("excluded_cms_core", False))
    core_x_cm = float(config.get("core_x_cm", 900.0))
    core_y_cm = float(config.get("core_y_cm", 1300.0))
    core_z_cm = float(config.get("core_z_cm", 1200.0))
    min_expected = float(config.get("min_expected", 1.0e-6))
    total_numerator = 0.0
    total_denominator = 0.0
    for iz in range(1, observed.GetNbinsZ() + 1):
        numerator = 0.0
        denominator = 0.0
        z_cm = unit_to_cm * observed.GetZaxis().GetBinCenter(iz)
        for ix in range(1, observed.GetNbinsX() + 1):
            x_cm = unit_to_cm * observed.GetXaxis().GetBinCenter(ix)
            for iy in range(1, observed.GetNbinsY() + 1):
                y_cm = unit_to_cm * observed.GetYaxis().GetBinCenter(iy)
                if (
                    exclude_core
                    and abs(x_cm) < core_x_cm
                    and abs(y_cm) < core_y_cm
                    and z_cm < core_z_cm
                ):
                    continue
                if expected.GetBinContent(ix, iy, iz) <= min_expected:
                    continue
                numerator += float(observed.GetBinContent(ix, iy, iz))
                denominator += float(expected.GetBinContent(ix, iy, iz))
        if denominator > 0.0:
            hist.SetBinContent(iz, numerator / denominator)
            hist.SetBinError(iz, math.sqrt(max(numerator, 0.0)) / denominator)
            total_numerator += numerator
            total_denominator += denominator
    shape_norm = total_numerator / total_denominator if total_denominator > 0.0 else 1.0
    if shape_norm > 0.0:
        hist.Scale(1.0 / shape_norm)
    return hist


def make_z_spectrum_hist(root, source, selection_expected, name, config):
    hist = root.TH1F(
        name,
        "",
        source.GetNbinsZ(),
        source.GetZaxis().GetXmin(),
        source.GetZaxis().GetXmax(),
    )
    hist.SetDirectory(0)
    copy_axis_style(source.GetZaxis(), hist.GetXaxis(), hist_axis_unit(source))
    hist.GetYaxis().SetTitle("voxel crossings / z bin")
    unit_to_cm = 100.0 if hist_axis_unit(source) == "m" else 1.0
    exclude_core = bool(config.get("excluded_cms_core", False))
    core_x_cm = float(config.get("core_x_cm", 900.0))
    core_y_cm = float(config.get("core_y_cm", 1300.0))
    core_z_cm = float(config.get("core_z_cm", 1200.0))
    min_expected = float(config.get("min_expected", 1.0e-6))
    for iz in range(1, source.GetNbinsZ() + 1):
        total = 0.0
        z_cm = unit_to_cm * source.GetZaxis().GetBinCenter(iz)
        for ix in range(1, source.GetNbinsX() + 1):
            x_cm = unit_to_cm * source.GetXaxis().GetBinCenter(ix)
            for iy in range(1, source.GetNbinsY() + 1):
                y_cm = unit_to_cm * source.GetYaxis().GetBinCenter(iy)
                if (
                    exclude_core
                    and abs(x_cm) < core_x_cm
                    and abs(y_cm) < core_y_cm
                    and z_cm < core_z_cm
                ):
                    continue
                if selection_expected.GetBinContent(ix, iy, iz) <= min_expected:
                    continue
                total += float(source.GetBinContent(ix, iy, iz))
        hist.SetBinContent(iz, total)
    return hist


def draw_correction_diagnostics_canvas(root, observed, raw_expected, calibrated_expected, config):
    pseudo_count = float(config.get("pseudo_count", 0.5))
    data_spectrum = make_z_spectrum_hist(
        root,
        observed,
        raw_expected,
        "data_depth_spectrum",
        config,
    )
    raw_mc_spectrum = make_z_spectrum_hist(
        root,
        raw_expected,
        raw_expected,
        "raw_mc_depth_spectrum",
        config,
    )
    corrected_mc_spectrum = make_z_spectrum_hist(
        root,
        calibrated_expected,
        raw_expected,
        "corrected_mc_depth_spectrum",
        config,
    )
    raw_ratio = make_z_ratio_hist(
        root, observed, raw_expected, "data_over_raw_mc_z", config
    )
    corrected_ratio = make_z_ratio_hist(
        root,
        observed,
        calibrated_expected,
        "data_over_corrected_mc_z",
        config,
    )
    applied_correction = make_z_ratio_hist(
        root,
        calibrated_expected,
        raw_expected,
        "corrected_over_raw_mc_z",
        config,
    )
    residual_before = make_xy_log_ratio(
        root,
        observed,
        raw_expected,
        "data_mc_log_ratio_xy_before",
        pseudo_count,
    )
    residual_after = make_xy_log_ratio(
        root,
        observed,
        calibrated_expected,
        "data_mc_log_ratio_xy_after",
        pseudo_count,
    )

    residual_values = []
    for hist in (residual_before, residual_after):
        for ix in range(1, hist.GetNbinsX() + 1):
            for iy in range(1, hist.GetNbinsY() + 1):
                residual_values.append(abs(float(hist.GetBinContent(ix, iy))))
    residual_bound = max(0.10, percentile(residual_values, 0.98))
    for hist in (residual_before, residual_after):
        hist.SetMinimum(-residual_bound)
        hist.SetMaximum(residual_bound)
        hist.SetStats(False)
        hist.GetXaxis().SetTitleSize(0.050)
        hist.GetYaxis().SetTitleSize(0.050)
        hist.GetZaxis().SetTitleSize(0.043)
        hist.GetXaxis().SetLabelSize(0.040)
        hist.GetYaxis().SetLabelSize(0.040)
        hist.GetZaxis().SetLabelSize(0.036)
        for ix in range(1, hist.GetNbinsX() + 1):
            for iy in range(1, hist.GetNbinsY() + 1):
                value = float(hist.GetBinContent(ix, iy))
                display_bound = 0.999 * residual_bound
                hist.SetBinContent(ix, iy, min(max(value, -display_bound), display_bound))

    colors = {
        "data": root.kBlack,
        "raw": root.TColor.GetColor("#2563EB"),
        "corrected": root.TColor.GetColor("#15803D"),
        "factor": root.TColor.GetColor("#EA580C"),
    }
    for hist, color, marker, line_style in (
        (data_spectrum, colors["data"], 20, 1),
        (raw_mc_spectrum, colors["raw"], 24, 7),
        (corrected_mc_spectrum, colors["corrected"], 25, 1),
    ):
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetLineWidth(4)
        hist.SetLineStyle(line_style)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(0.75)
        hist.SetStats(False)

    spectrum_values = []
    for hist in (data_spectrum, raw_mc_spectrum, corrected_mc_spectrum):
        spectrum_values.extend(
            float(hist.GetBinContent(ibin))
            for ibin in range(1, hist.GetNbinsX() + 1)
            if hist.GetBinContent(ibin) > 0.0
        )
    spectrum_max = max(spectrum_values) if spectrum_values else 1.0
    spectrum_min = min(spectrum_values) if spectrum_values else 1.0
    data_spectrum.SetMinimum(max(0.5, 0.5 * spectrum_min))
    data_spectrum.SetMaximum(2.5 * spectrum_max)
    data_spectrum.GetXaxis().SetTitle(f"cavern z / up [{hist_axis_unit(observed)}]")
    data_spectrum.GetYaxis().SetTitle("voxel crossings / z bin")
    data_spectrum.GetXaxis().SetTitleSize(0.050)
    data_spectrum.GetYaxis().SetTitleSize(0.050)
    data_spectrum.GetXaxis().SetLabelSize(0.040)
    data_spectrum.GetYaxis().SetLabelSize(0.040)

    for hist, color, marker in (
        (raw_ratio, colors["raw"], 20),
        (corrected_ratio, colors["corrected"], 21),
        (applied_correction, colors["factor"], 24),
    ):
        hist.SetLineColor(color)
        hist.SetMarkerColor(color)
        hist.SetLineWidth(4)
        hist.SetMarkerStyle(marker)
        hist.SetMarkerSize(0.85)
        hist.SetStats(False)

    visible_ratios = []
    for hist in (raw_ratio, corrected_ratio, applied_correction):
        visible_ratios.extend(
            float(hist.GetBinContent(ibin))
            for ibin in range(1, hist.GetNbinsX() + 1)
            if hist.GetBinContent(ibin) > 0.0
        )
    ratio_low = min(0.95, percentile(visible_ratios, 0.02))
    ratio_high = max(1.05, percentile(visible_ratios, 0.98))
    ratio_padding = max(0.08, 0.12 * (ratio_high - ratio_low))
    raw_ratio.SetMinimum(max(0.0, ratio_low - ratio_padding))
    raw_ratio.SetMaximum(ratio_high + ratio_padding)
    raw_ratio.GetXaxis().SetTitle(f"cavern z / up [{hist_axis_unit(observed)}]")
    raw_ratio.GetYaxis().SetTitle("ratio")
    raw_ratio.GetXaxis().SetTitleSize(0.050)
    raw_ratio.GetYaxis().SetTitleSize(0.050)
    raw_ratio.GetXaxis().SetLabelSize(0.040)
    raw_ratio.GetYaxis().SetLabelSize(0.040)

    canvas = root.TCanvas(
        "canvas_data_mc_correction_diagnostics",
        "Data-driven Monte Carlo correction diagnostics",
        3600,
        2400,
    )
    canvas.Divide(2, 2, 0.002, 0.002)
    keepalive = [
        data_spectrum,
        raw_mc_spectrum,
        corrected_mc_spectrum,
        raw_ratio,
        corrected_ratio,
        applied_correction,
        residual_before,
        residual_after,
    ]

    spectrum_pad = canvas.cd(1)
    spectrum_pad.SetLeftMargin(0.13)
    spectrum_pad.SetRightMargin(0.04)
    spectrum_pad.SetBottomMargin(0.13)
    spectrum_pad.SetTopMargin(0.14)
    spectrum_pad.SetLogy(True)
    data_spectrum.Draw("P")
    raw_mc_spectrum.Draw("HIST SAME")
    corrected_mc_spectrum.Draw("HIST SAME")
    data_spectrum.Draw("P SAME")
    spectrum_legend = root.TLegend(0.15, 0.67, 0.72, 0.84)
    spectrum_legend.SetBorderSize(0)
    spectrum_legend.SetFillStyle(0)
    spectrum_legend.SetTextFont(42)
    spectrum_legend.SetTextSize(0.034)
    spectrum_legend.AddEntry(data_spectrum, "data", "p")
    spectrum_legend.AddEntry(raw_mc_spectrum, "MC before correction", "l")
    spectrum_legend.AddEntry(corrected_mc_spectrum, "MC after correction", "l")
    spectrum_legend.Draw()
    keepalive.append(spectrum_legend)
    keepalive.append(draw_label(root, "Depth spectra", size=0.047, x=0.14, y=0.93))
    keepalive.append(
        draw_label(
            root,
            "absolute voxel crossings in the fit region",
            size=0.028,
            x=0.14,
            y=0.885,
            font=42,
        )
    )

    ratio_pad = canvas.cd(2)
    ratio_pad.SetLeftMargin(0.13)
    ratio_pad.SetRightMargin(0.04)
    ratio_pad.SetBottomMargin(0.13)
    ratio_pad.SetTopMargin(0.14)
    raw_ratio.Draw("E1")
    corrected_ratio.Draw("E1 SAME")
    applied_correction.Draw("HIST SAME")
    unity = root.TLine(
        raw_ratio.GetXaxis().GetXmin(),
        1.0,
        raw_ratio.GetXaxis().GetXmax(),
        1.0,
    )
    unity.SetLineColor(root.kGray + 2)
    unity.SetLineStyle(2)
    unity.SetLineWidth(2)
    unity.Draw()
    legend = root.TLegend(0.15, 0.68, 0.91, 0.84)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.034)
    legend.AddEntry(raw_ratio, "data / raw MC (shape)", "lep")
    legend.AddEntry(corrected_ratio, "data / corrected MC (residual)", "lep")
    legend.AddEntry(applied_correction, "applied MC correction (shape)", "l")
    legend.Draw()
    keepalive.extend([unity, legend])
    keepalive.append(draw_label(root, "Depth ratios and correction", size=0.047, x=0.14, y=0.93))

    for ipad, hist, label in (
        (3, residual_before, "Before correction"),
        (4, residual_after, "After correction"),
    ):
        pad = canvas.cd(ipad)
        pad.SetLeftMargin(0.12)
        pad.SetRightMargin(0.17)
        pad.SetBottomMargin(0.13)
        pad.SetTopMargin(0.14)
        hist.Draw("COLZ")
        keepalive.append(draw_label(root, label, size=0.047, x=0.13, y=0.93))
        keepalive.append(
            draw_label(
                root,
                "z-integrated x-y residual; shared 98% color range",
                size=0.028,
                x=0.13,
                y=0.885,
                font=42,
            )
        )

    z_smooth_cm = float(config.get("z_smooth_cm", 0.0))
    xy_order = config.get("xy_poly_order", "?")
    canvas.cd(2)
    keepalive.append(
        draw_label(
            root,
            f"unit-normalized fit region; C(z): {0.01 * z_smooth_cm:g} m smoothing; "
            f"Cxy: order {xy_order}",
            size=0.026,
            x=0.14,
            y=0.875,
            font=42,
        )
    )
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas, {
        "data_depth_spectrum": data_spectrum,
        "raw_mc_depth_spectrum": raw_mc_spectrum,
        "corrected_mc_depth_spectrum": corrected_mc_spectrum,
        "data_over_raw_mc_z": raw_ratio,
        "data_over_corrected_mc_z": corrected_ratio,
        "corrected_over_raw_mc_z": applied_correction,
        "data_mc_log_ratio_xy_before": residual_before,
        "data_mc_log_ratio_xy_after": residual_after,
    }


def slice_axes_for(source_hist, stack_kind):
    coordinate_system = hist_coordinate_system(source_hist)
    if coordinate_system == "cavern":
        if stack_kind == "xy_z":
            return 0, 1, 2
        if stack_kind == "xz_y":
            return 0, 2, 1
    else:
        if stack_kind == "xy_z":
            return 0, 2, 1
        if stack_kind == "xz_y":
            return 0, 1, 2
    raise ValueError(f"unknown slice stack kind: {stack_kind}")


def make_slice_stack_canvas(root, output_dir, hist3, hist_name, stack_kind, values, args):
    axis_a, axis_b, axis_c = slice_axes_for(hist3, stack_kind)
    values = clipped_values(values, axis(hist3, axis_c))
    canvas = root.TCanvas(
        f"canvas_{safe_name(hist_name)}_{stack_kind}_slices",
        f"{hist_name} {stack_kind} slices",
        args.canvas_width,
        args.canvas_height,
    )
    canvas.Divide(args.canvas_columns, args.canvas_rows, 0.001, 0.001)
    keepalive = []
    max_pads = args.canvas_columns * args.canvas_rows
    slice_dir = output_dir.mkdir(f"{stack_kind}_slices")
    for ipad, value in enumerate(values[:max_pads], start=1):
        slice_name = f"{stack_kind}_{value_token(value)}"
        hist2 = make_slice(root, hist3, axis_a, axis_b, axis_c, value, slice_name)
        slice_dir.cd()
        hist2.Write()
        canvas.cd(ipad)
        fixed_bin = axis(hist3, axis_c).FindBin(value)
        fixed_center = axis(hist3, axis_c).GetBinCenter(fixed_bin)
        label = (
            f"{quantity_short(hist_name)}; "
            f"{stack_kind.replace('_', ' ')} = {fixed_center:g} {hist_axis_unit(hist3)}"
        )
        keepalive.extend(
            draw_slice_panel(
                root,
                canvas.cd(ipad),
                hist2,
                hist3,
                hist_name,
                axis_a,
                axis_b,
                axis_c,
                fixed_center,
                label,
                args,
            )
        )
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def copy_hist(root, source_file, hist_name, plot_units="cm", display_rebin=1):
    hist = require_th3(source_file, hist_name)
    display_hist = clone_hist_with_scaled_axes(root, hist, safe_name(hist_name), plot_units)
    display_hist = rebin_display_th3(root, display_hist, hist_name, display_rebin)
    display_hist.SetTitle(root_title_text(quantity_detail(hist_name)))
    return display_hist


def export_canvas(canvas, png_dir):
    if not png_dir:
        return
    os.makedirs(png_dir, exist_ok=True)
    canvas.SaveAs(os.path.join(png_dir, f"{canvas.GetName()}.png"))


def main():
    args = parse_args()

    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)

    source_file = ROOT.TFile.Open(args.input)
    if not source_file or source_file.IsZombie():
        raise RuntimeError(f"Could not open {args.input}")

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create {args.output}")

    copied = {}
    inputs_dir = output.mkdir("inputs")
    inputs_dir.cd()
    for hist_name in split_csv(args.copy_hists):
        hist = copy_hist(ROOT, source_file, hist_name, args.plot_units, args.display_rebin)
        hist.Write(hist_name)
        copied[hist_name] = hist
    output.cd()

    z_values_cm = split_float_csv(args.z_slices_cm)
    y_values_cm = split_float_csv(args.y_slices_cm)
    canvases = []
    requested_hist_names = split_csv(args.hist_names)
    comparison_hist_names = split_csv(args.comparison_hists)

    if not args.no_3d and not args.no_geometry_baseline:
        baseline_name = requested_hist_names[0] if requested_hist_names else comparison_hist_names[0]
        baseline_hist = copied.get(baseline_name) or copy_hist(
            ROOT, source_file, baseline_name, args.plot_units, args.display_rebin
        )
        output.cd()
        canvas = draw_geometry_baseline_canvas(ROOT, baseline_hist, args)
        canvas.Write()
        export_canvas(canvas, args.png_dir)
        canvases.append(canvas)

    if not args.no_3d and comparison_hist_names:
        comparison_hists = [
            copied.get(name) or copy_hist(ROOT, source_file, name, args.plot_units, args.display_rebin)
            for name in comparison_hist_names
        ]
        output.cd()
        canvas = draw_3d_comparison_canvas(ROOT, comparison_hists, comparison_hist_names, args)
        canvas.Write()
        export_canvas(canvas, args.png_dir)
        canvases.append(canvas)

    correction_hist_names = (
        "observed_voxels",
        "expected_voxels_global_scaled",
        "calibrated_expected_voxels",
    )
    if all(source_file.Get(name) for name in correction_hist_names):
        correction_hists = [
            copied.get(name) or copy_hist(ROOT, source_file, name, args.plot_units, args.display_rebin)
            for name in correction_hist_names
        ]
        calibration_config = {}
        calibration_config_object = source_file.Get("smooth_calibration_config_json")
        if calibration_config_object:
            try:
                calibration_config = json.loads(calibration_config_object.GetString().Data())
            except (AttributeError, TypeError, ValueError):
                calibration_config = {}
        output.cd()
        canvas, diagnostic_hists = draw_correction_diagnostics_canvas(
            ROOT,
            correction_hists[0],
            correction_hists[1],
            correction_hists[2],
            calibration_config,
        )
        canvas.Write()
        export_canvas(canvas, args.png_dir)
        canvases.append(canvas)
        diagnostic_dir = output.mkdir("correction_diagnostics")
        diagnostic_dir.cd()
        for name, hist in diagnostic_hists.items():
            hist.Write(name)
        output.cd()

    for hist_name in requested_hist_names:
        hist = copied.get(hist_name) or copy_hist(ROOT, source_file, hist_name, args.plot_units, args.display_rebin)
        hist_dir = output.mkdir(safe_name(hist_name))
        hist_dir.cd()
        hist.Write(hist_name)
        if not args.no_3d and args.draw_3d_no_cut:
            canvas = draw_3d_canvas(ROOT, hist, hist_name, args)
            hist_dir.cd()
            canvas.Write()
            export_canvas(canvas, args.png_dir)
            canvases.append(canvas)
        if not args.no_3d and not args.no_3d_grid:
            canvas = draw_3d_threshold_grid_canvas(ROOT, hist, hist_name, args)
            hist_dir.cd()
            canvas.Write()
            export_canvas(canvas, args.png_dir)
            canvases.append(canvas)
        z_values = [cm_to_hist_units(hist, value) for value in z_values_cm]
        y_values = [cm_to_hist_units(hist, value) for value in y_values_cm]
        for stack_kind, values in (("xy_z", z_values), ("xz_y", y_values)):
            hist_dir.cd()
            canvas = make_slice_stack_canvas(ROOT, hist_dir, hist, hist_name, stack_kind, values, args)
            hist_dir.cd()
            canvas.Write()
            export_canvas(canvas, args.png_dir)
            canvases.append(canvas)
        output.cd()

    metadata = {
        "input": args.input,
        "hist_names": split_csv(args.hist_names),
        "copy_hists": split_csv(args.copy_hists),
        "z_slices_cm": z_values_cm,
        "y_slices_cm": y_values_cm,
        "plot_units": args.plot_units,
        "display_rebin": args.display_rebin,
        "display_rebin_modes": {
            name: display_rebin_mode(name)
            for name in sorted(set(split_csv(args.copy_hists) + split_csv(args.hist_names)))
        },
        "quantity_descriptions": {
            name: quantity_detail(name)
            for name in sorted(set(split_csv(args.copy_hists) + split_csv(args.hist_names)))
        },
        "log_color": not args.linear_color,
        "log_3d_color": args.log_3d_color,
        "cavern_outline": not args.no_cavern_outline,
        "draw_3d_no_cut": args.draw_3d_no_cut,
        "draw_geometry_baseline": not args.no_3d and not args.no_geometry_baseline,
        "draw_3d_threshold_grid": not args.no_3d and not args.no_3d_grid,
        "three_d_top_fractions": split_float_csv(args.three_d_top_fractions),
        "comparison_hists": comparison_hist_names,
        "comparison_top_fraction": args.comparison_top_fraction,
        "shared_3d_color_scales": True,
        "data_mc_correction_diagnostics": True,
        "cavern_model": {
            "hall_length_m": 0.01 * (CAVERN_HALL_BOX_CM["ymax"] - CAVERN_HALL_BOX_CM["ymin"]),
            "hall_width_m": 0.01 * (CAVERN_HALL_BOX_CM["xmax"] - CAVERN_HALL_BOX_CM["xmin"]),
            "hall_height_m": 0.01 * (CAVERN_HALL_BOX_CM["zmax"] - CAVERN_HALL_BOX_CM["zmin"]),
            "floor_z_m": 0.01 * CAVERN_HALL_BOX_CM["zmin"],
            "arch_spring_z_m": 0.01 * CAVERN_HALL_ARCH_SPRING_Z_CM,
            "crown_z_m": 0.01 * CAVERN_HALL_BOX_CM["zmax"],
            "px56_shaft_radius_m": 0.01 * CAVERN_SHAFT_CM["radius"],
            "px56_shaft_center_x_m": 0.01 * CAVERN_SHAFT_CM["x0"],
            "px56_shaft_center_y_m": 0.01 * CAVERN_SHAFT_CM["y0"],
            "note": (
                "Approximate display geometry constrained by public CMS cavern and detector "
                "dimensions; floor is inferred just below the outer CMS radius, and PX56 is "
                "placed on the CMS -z side, which is positive y in the cavern frame."
            ),
        },
        "note": (
            "Display ROOT file for Li-Ma/cavity studies. It keeps the Li-Ma input "
            "TH3s under inputs/ and uses slice stacks rather than thresholded "
            "projection grids."
        ),
    }
    output.WriteObject(ROOT.TObjString(json.dumps(metadata, indent=2, sort_keys=True)), "slice_stack_config_json")
    output.Close()
    source_file.Close()
    print(f"Wrote {args.output} with {len(canvases)} canvases")
    if args.png_dir:
        print(f"Exported PNGs to {args.png_dir}")


if __name__ == "__main__":
    main()
