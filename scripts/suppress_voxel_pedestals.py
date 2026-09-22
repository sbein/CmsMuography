#!/usr/bin/env python3
"""Create staged pedestal-suppressed TH3 voxel maps and 2D views."""

import argparse
import math
import os
from array import array


AXIS_NAMES = ("x", "y", "z")
CMS_COMPONENT_RADII = [
    ("tracker", 120.0, 1),
    ("ecal", 180.0, 1),
    ("hcal", 295.0, 1),
    ("solenoid", 380.0, 1),
    ("muon_mb1", 445.0, 1),
    ("muon_mb2", 535.0, 1),
    ("muon_mb3", 625.0, 1),
    ("muon_mb4", 740.0, 1),
]
CMS_HALF_LENGTH_CM = 1100.0
CMS_Z_PLANES = [-1100.0, -800.0, -600.0, 600.0, 800.0, 1100.0]
CMS_POINT5_LATLON = (46.310819, 6.078311)
ATLAS_POINT1_LATLON = (46.23561, 6.05426)
CARDINAL_BEARINGS = (("N", 0.0), ("E", 90.0), ("S", 180.0), ("W", 270.0))
DEFAULT_GROUND_COEFFICIENTS_CM = (
    10088.8888888889,
    -0.01591579598920878,
    0.03600154133089286,
    -2.318859898172708e-06,
    5.994453725606241e-07,
    1.5621173989662048e-07,
)
CAVERN_HALL_BOX_CM = {
    "xmin": -1350.0,
    "xmax": 1350.0,
    "ymin": -2650.0,
    "ymax": 2650.0,
    "zmin": -800.0,
    "zmax": 1600.0,
}
CAVERN_HALL_ARCH_SPRING_Z_CM = 400.0
CAVERN_SHAFT_CM = {
    "x0": 0.0,
    # PX56 is on the CMS -z side; in the cavern frame y=-CMS z, so this is +y.
    "y0": 1625.0,
    "radius": 1025.0,
    "zmin": CAVERN_HALL_BOX_CM["zmax"],
    "zmax": 11000.0,
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Read a TH3 voxel map, build top-fraction pedestal thresholds, and "
            "write one ROOT directory per suppressed stage."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input ROOT file containing the voxel TH3.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output ROOT file for staged pedestal-suppressed maps.",
    )
    parser.add_argument(
        "--hist-name",
        default="cavern_muon_voxels",
        help="Name/path of the source TH3 in the input file.",
    )
    parser.add_argument(
        "--top-fractions",
        default="0.25,0.10,0.05,0.025,0.01,0.005,0.0025,0.001,0.0001",
        help=(
            "Comma-separated retained top fractions for the default 9 cut sets. "
            "These are quantile-based and need not map to equidistant thresholds."
        ),
    )
    parser.add_argument(
        "--exponential-stages",
        action="store_true",
        help="Use top 10^-N stages instead of --top-fractions.",
    )
    parser.add_argument(
        "--first-exponent",
        type=int,
        default=1,
        help="First top-fraction exponent for --exponential-stages.",
    )
    parser.add_argument(
        "--nstages",
        type=int,
        default=9,
        help="Number of stages for --exponential-stages.",
    )
    parser.add_argument(
        "--include-zero-bins",
        action="store_true",
        help="Compute percentiles over all voxels, including empty bins.",
    )
    parser.add_argument(
        "--mode",
        choices=("keep_above", "subtract"),
        default="keep_above",
        help=(
            "keep_above zeros bins below the pedestal and keeps original counts; "
            "subtract uses max(count - pedestal, 0)."
        ),
    )
    parser.add_argument("--slice-x", type=float, default=0.0, help="X value in cm for the YZ slice.")
    parser.add_argument("--slice-y", type=float, default=0.0, help="Y value in cm for the XZ slice.")
    parser.add_argument("--slice-z", type=float, default=0.0, help="Z value in cm for the XY slice.")
    parser.add_argument(
        "--write-source",
        action="store_true",
        help="Also write a source directory with a copy of the unsuppressed TH3.",
    )
    parser.add_argument(
        "--include-source-in-canvases",
        action="store_true",
        help=(
            "Include the unsuppressed source map as the first pad in comparison canvases. "
            "This is the default when --write-source is used."
        ),
    )
    parser.add_argument(
        "--exclude-source-from-canvases",
        action="store_true",
        help="Do not include the unsuppressed source map in comparison canvases.",
    )
    parser.add_argument(
        "--draw-source-3d-pad",
        action="store_true",
        help=(
            "Draw the unsuppressed source TH3 in the first 3D comparison pad. "
            "By default that pad is left blank to avoid the heavy raw 3D render."
        ),
    )
    parser.add_argument(
        "--skip-stage-dirs",
        action="store_true",
        help="Do not write per-stage directories; only write comparison canvases and metadata.",
    )
    parser.add_argument(
        "--skip-canvases",
        action="store_true",
        help="Do not write comparison canvases.",
    )
    parser.add_argument(
        "--skip-projection-canvases",
        action="store_true",
        help="Do not write xy/xz/yz projection canvases; keep only 3D and fixed-slice canvases.",
    )
    parser.add_argument(
        "--skip-slice-canvases",
        action="store_true",
        help="Do not write fixed xy/xz/yz slice canvases for each threshold; useful when separate slice-stack products exist.",
    )
    parser.add_argument(
        "--skip-source-canvases",
        action="store_true",
        help="Write the unsuppressed source histogram and views, but do not draw source canvases.",
    )
    parser.add_argument(
        "--canvas-columns",
        type=int,
        default=3,
        help="Number of pads along x for threshold comparison canvases.",
    )
    parser.add_argument(
        "--canvas-rows",
        type=int,
        default=3,
        help="Number of pads along y for threshold comparison canvases.",
    )
    parser.add_argument("--canvas-width", type=int, default=3600, help="Canvas width in pixels.")
    parser.add_argument("--canvas-height", type=int, default=3600, help="Canvas height in pixels.")
    parser.add_argument(
        "--plot-units",
        choices=("m", "cm"),
        default="m",
        help="Length units used for drawn canvas axes. Histogram contents are unchanged.",
    )
    parser.add_argument(
        "--no-cms-outline",
        action="store_true",
        help="Do not overlay the simple CMS detector wireframe on canvases.",
    )
    parser.add_argument(
        "--no-cavern-outline",
        action="store_true",
        help="Do not overlay the approximate CMS experimental cavern and shaft wireframe.",
    )
    parser.add_argument(
        "--no-cms-outline-2d",
        action="store_true",
        help="Keep the CMS wireframe on 3D canvases but omit it from 2D projection/slice canvases.",
    )
    parser.add_argument(
        "--no-orientation-overlay",
        action="store_true",
        help="Do not draw the purple north/up compass rose or out-of-page markers.",
    )
    parser.add_argument(
        "--cms-x-bearing-deg",
        type=float,
        default=None,
        help=(
            "Geographic bearing, in degrees clockwise from north, of +CMS x/+cavern x. "
            "By default this is inferred from CMS Point 5 toward ATLAS Point 1, "
            "using +x as the inward LHC-center direction."
        ),
    )
    parser.add_argument(
        "--no-ground-overlay",
        action="store_true",
        help="Do not draw the approximate CMS Point 5 ground surface on canvases.",
    )
    parser.add_argument(
        "--ground-coefficients-cm",
        default=",".join(str(value) for value in DEFAULT_GROUND_COEFFICIENTS_CM),
        help=(
            "Comma-separated coefficients c0,cx,cy,cxx,cxy,cyy for the approximate "
            "surface z_cm=c0+cx*x+cy*y+cxx*x^2+cxy*x*y+cyy*y^2 in cavern cm."
        ),
    )
    parser.add_argument(
        "--ground-grid-step-cm",
        type=float,
        default=500.0,
        help="Grid step in cm for the approximate 3D ground wireframe.",
    )
    parser.add_argument(
        "--png-dir",
        help="Optional directory where comparison and source canvases are also exported as PNG.",
    )
    return parser.parse_args()


def require_th3(root_file, hist_name):
    hist = root_file.Get(hist_name)
    if not hist:
        raise RuntimeError(f"Could not find histogram {hist_name!r} in {root_file.GetName()}")
    if hist.InheritsFrom("TH3") is False:
        raise RuntimeError(f"Object {hist_name!r} is not a TH3")
    hist.SetDirectory(0)
    return hist


def title_without_axis_directives(hist):
    return str(hist.GetTitle()).split(";", 1)[0]


def hist_coordinate_system(hist):
    titles = " ".join(
        str(axis(hist, index).GetTitle()).lower()
        for index in range(3)
        if hasattr(hist, "GetZaxis") or index < 2
    )
    if "cavern" in titles or ("vertical" in titles and "beam" in titles and "cms z" in titles):
        return "cavern"
    return "cms"


def axis_role(coordinate_system, axis_index):
    if coordinate_system == "cavern":
        return ("transverse", "beam", "transverse")[axis_index]
    return ("transverse", "transverse", "beam")[axis_index]


def vertical_axis(coordinate_system):
    return 2 if coordinate_system == "cavern" else 1


def axis_direction_label(coordinate_system, axis_index):
    if coordinate_system == "cavern":
        return ("+x", "+y=-CMS z", "+z up")[axis_index]
    return ("+x", "+y up", "+z beam")[axis_index]


def initial_bearing_degrees(lat1, lon1, lat2, lon2):
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_lon = math.radians(lon2 - lon1)
    y = math.sin(delta_lon) * math.cos(phi2)
    x = math.cos(phi1) * math.sin(phi2) - math.sin(phi1) * math.cos(phi2) * math.cos(delta_lon)
    return (math.degrees(math.atan2(y, x)) + 360.0) % 360.0


def default_cms_x_bearing_degrees():
    # Point 5 and Point 1 are opposite LHC interaction regions, so P5->P1 approximates
    # the inward radial direction used by +CMS x.
    return initial_bearing_degrees(*CMS_POINT5_LATLON, *ATLAS_POINT1_LATLON)


def cms_x_bearing_degrees(args):
    if args.cms_x_bearing_deg is not None:
        return args.cms_x_bearing_deg % 360.0
    return default_cms_x_bearing_degrees()


def horizontal_axis_bearings(coordinate_system, cms_x_bearing_deg):
    # CMS is right-handed: +x points toward the LHC center, +y is vertical/up,
    # and +z is the beam axis.  The Z-up cavern frame used by the muography
    # plots is x=CMS x, y=-CMS z, z=CMS y, so its horizontal map plane is
    # right-handed when viewed from above.
    cms_z_bearing = (cms_x_bearing_deg + 90.0) % 360.0
    if coordinate_system == "cavern":
        return {0: cms_x_bearing_deg, 1: (cms_z_bearing + 180.0) % 360.0}
    return {0: cms_x_bearing_deg, 2: cms_z_bearing}


def geographic_direction_vector(coordinate_system, cms_x_bearing_deg, bearing_deg):
    components = [0.0, 0.0, 0.0]
    for axis_index, axis_bearing in horizontal_axis_bearings(coordinate_system, cms_x_bearing_deg).items():
        components[axis_index] = math.cos(math.radians(bearing_deg - axis_bearing))
    return tuple(components)


def parse_ground_coefficients(text):
    values = [float(piece.strip()) for piece in text.split(",") if piece.strip()]
    if len(values) != 6:
        raise ValueError("--ground-coefficients-cm must have 6 comma-separated numbers")
    return tuple(values)


def ground_z_cavern_cm(args, cavern_x, cavern_y):
    c0, cx, cy, cxx, cxy, cyy = parse_ground_coefficients(args.ground_coefficients_cm)
    return (
        c0
        + cx * cavern_x
        + cy * cavern_y
        + cxx * cavern_x * cavern_x
        + cxy * cavern_x * cavern_y
        + cyy * cavern_y * cavern_y
    )


def cavern_to_hist_point(coordinate_system, cavern_x, cavern_y, cavern_z):
    if coordinate_system == "cavern":
        return (cavern_x, cavern_y, cavern_z)
    return (cavern_x, cavern_z, -cavern_y)


def hist_point_to_cavern(coordinate_system, point):
    if coordinate_system == "cavern":
        return point
    return (point[0], -point[2], point[1])


def cavern_cm_to_hist_units_point(hist, coordinate_system, cavern_x, cavern_y, cavern_z):
    point_cm = cavern_to_hist_point(coordinate_system, cavern_x, cavern_y, cavern_z)
    return tuple(cm_to_hist_units(hist, value) for value in point_cm)


def hist_units_point_to_cavern_cm(hist, coordinate_system, point):
    point_cm = tuple(hist_units_to_cm(hist, value) for value in point)
    return hist_point_to_cavern(coordinate_system, point_cm)


def horizontal_cavern_ranges(hist, coordinate_system):
    unit_to_cm = hist_unit_to_cm(hist)
    ranges = [
        (axis(hist, iaxis).GetXmin(), axis(hist, iaxis).GetXmax())
        for iaxis in range(3)
    ]
    if coordinate_system == "cavern":
        return (
            (ranges[0][0] * unit_to_cm, ranges[0][1] * unit_to_cm),
            (ranges[1][0] * unit_to_cm, ranges[1][1] * unit_to_cm),
        )
    return (
        (ranges[0][0] * unit_to_cm, ranges[0][1] * unit_to_cm),
        tuple(sorted((-ranges[2][0] * unit_to_cm, -ranges[2][1] * unit_to_cm))),
    )


def coordinate_value_to_cavern_horizontal(hist, coordinate_system, axis_index, value, fixed_other=0.0):
    point = [0.0, 0.0, 0.0]
    point[axis_index] = value
    for candidate_axis in horizontal_axis_bearings(coordinate_system, default_cms_x_bearing_degrees()):
        if candidate_axis != axis_index:
            point[candidate_axis] = fixed_other
            break
    cavern_x, cavern_y, _cavern_z = hist_units_point_to_cavern_cm(hist, coordinate_system, tuple(point))
    return cavern_x, cavern_y


def view_axis_order(hist3):
    if hist_coordinate_system(hist3) == "cavern":
        return {
            "xy_projection": (0, 1),
            "xz_projection": (0, 2),
            "yz_projection": (1, 2),
            "xy_slice_z": (0, 1, 2),
            "xz_slice_y": (0, 2, 1),
            "yz_slice_x": (1, 2, 0),
        }
    return {
        "xy_projection": (0, 1),
        "xz_projection": (2, 0),
        "yz_projection": (2, 1),
        "xy_slice_z": (0, 1, 2),
        "xz_slice_y": (2, 0, 1),
        "yz_slice_x": (2, 1, 0),
    }


def axis(hist, axis_index):
    if axis_index == 0:
        return hist.GetXaxis()
    if axis_index == 1:
        return hist.GetYaxis()
    if axis_index == 2:
        return hist.GetZaxis()
    raise ValueError(f"bad axis index: {axis_index}")


def hist_axis_unit(hist):
    titles = " ".join(str(axis(hist, index).GetTitle()).lower() for index in range(3))
    if "[m]" in titles:
        return "m"
    return "cm"


def hist_unit_to_cm(hist):
    return 100.0 if hist_axis_unit(hist) == "m" else 1.0


def cm_to_hist_units(hist, value_cm):
    return value_cm / hist_unit_to_cm(hist)


def hist_units_to_cm(hist, value):
    return value * hist_unit_to_cm(hist)


def axis_title_in_units(title, unit):
    title = title or ""
    if unit == "m":
        title = title.replace("[cm]", "[m]")
        title = title.replace(" cm", " m")
    elif unit == "cm":
        title = title.replace("[m]", "[cm]")
        title = title.replace(" m", " cm")
    return title


def copy_axis_style(source, target, unit):
    target.SetTitle(axis_title_in_units(source.GetTitle(), unit))
    target.SetLabelSize(source.GetLabelSize())
    target.SetTitleSize(source.GetTitleSize())
    target.SetTitleOffset(source.GetTitleOffset())
    target.SetLabelOffset(source.GetLabelOffset())
    target.SetNdivisions(source.GetNdivisions())


def clone_hist_with_scaled_axes(root, hist, name, unit):
    current_unit = hist_axis_unit(hist)
    if current_unit == unit:
        clone = hist.Clone(name)
        clone.SetDirectory(0)
        return clone
    scale = 0.01 if current_unit == "cm" and unit == "m" else 100.0
    if hist.InheritsFrom("TH3"):
        clone = root.TH3F(
            name,
            axis_title_in_units(hist.GetTitle(), unit),
            hist.GetNbinsX(),
            hist.GetXaxis().GetXmin() * scale,
            hist.GetXaxis().GetXmax() * scale,
            hist.GetNbinsY(),
            hist.GetYaxis().GetXmin() * scale,
            hist.GetYaxis().GetXmax() * scale,
            hist.GetNbinsZ(),
            hist.GetZaxis().GetXmin() * scale,
            hist.GetZaxis().GetXmax() * scale,
        )
    elif hist.InheritsFrom("TH2"):
        clone = root.TH2F(
            name,
            axis_title_in_units(hist.GetTitle(), unit),
            hist.GetNbinsX(),
            hist.GetXaxis().GetXmin() * scale,
            hist.GetXaxis().GetXmax() * scale,
            hist.GetNbinsY(),
            hist.GetYaxis().GetXmin() * scale,
            hist.GetYaxis().GetXmax() * scale,
        )
    else:
        clone = hist.Clone(name)
        clone.SetDirectory(0)
        return clone
    clone.SetDirectory(0)
    copy_axis_style(hist.GetXaxis(), clone.GetXaxis(), unit)
    copy_axis_style(hist.GetYaxis(), clone.GetYaxis(), unit)
    if hist.InheritsFrom("TH3"):
        copy_axis_style(hist.GetZaxis(), clone.GetZaxis(), unit)
    for ix in range(0, hist.GetNbinsX() + 2):
        for iy in range(0, hist.GetNbinsY() + 2):
            if hist.InheritsFrom("TH3"):
                for iz in range(0, hist.GetNbinsZ() + 2):
                    source_bin = hist.GetBin(ix, iy, iz)
                    clone.SetBinContent(source_bin, hist.GetBinContent(source_bin))
                    if hist.GetSumw2N() > 0:
                        clone.SetBinError(source_bin, hist.GetBinError(source_bin))
            else:
                source_bin = hist.GetBin(ix, iy)
                clone.SetBinContent(source_bin, hist.GetBinContent(source_bin))
                if hist.GetSumw2N() > 0:
                    clone.SetBinError(source_bin, hist.GetBinError(source_bin))
    clone.SetEntries(hist.GetEntries())
    return clone


def plot_hist(root, hist, name, args):
    return clone_hist_with_scaled_axes(root, hist, name, getattr(args, "plot_units", "cm"))


def nbins(hist, axis_index):
    return axis(hist, axis_index).GetNbins()


def bin_ranges(hist):
    return (
        range(1, hist.GetNbinsX() + 1),
        range(1, hist.GetNbinsY() + 1),
        range(1, hist.GetNbinsZ() + 1),
    )


def collect_scores(hist, include_zero_bins):
    values = []
    xr, yr, zr = bin_ranges(hist)
    for ix in xr:
        for iy in yr:
            for iz in zr:
                value = float(hist.GetBinContent(ix, iy, iz))
                if include_zero_bins or value > 0.0:
                    values.append(value)
    values.sort()
    if not values:
        raise RuntimeError("The selected voxel-score distribution is empty.")
    return values


def stage_name(index, percentile):
    token = f"{percentile:.4g}".replace(".", "p").replace("-", "m")
    return f"cut_{index:02d}_p{token}"


def threshold_for_top_fraction(values, top_fraction):
    target_count = max(1, int(math.ceil(top_fraction * len(values))))
    threshold = values[-target_count]
    return threshold, target_count


def parse_top_fractions(text):
    fractions = []
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        value = float(item)
        if value <= 0.0 or value > 1.0:
            raise ValueError(f"Top fractions must be in (0, 1], got {value}")
        fractions.append(value)
    if not fractions:
        raise ValueError("No top fractions were configured")
    return fractions


def make_stages_from_fractions(values, fractions):
    stages = []
    for index, top_fraction in enumerate(fractions, start=1):
        threshold, target_count = threshold_for_top_fraction(values, top_fraction)
        percentile = 100.0 * (1.0 - top_fraction)
        stages.append(
            {
                "name": stage_name(index, percentile),
                "index": index,
                "exponent": 0,
                "top_fraction": top_fraction,
                "percentile": percentile,
                "threshold": threshold,
                "target_count": target_count,
            }
        )
    return stages


def make_exponential_stages(values, first_exponent, nstages):
    stages = []
    for index, exponent in enumerate(range(first_exponent, first_exponent + nstages), start=1):
        top_fraction = 10.0 ** (-exponent)
        threshold, target_count = threshold_for_top_fraction(values, top_fraction)
        percentile = 100.0 * (1.0 - top_fraction)
        stages.append(
            {
                "name": stage_name(index, percentile),
                "index": index,
                "exponent": exponent,
                "top_fraction": top_fraction,
                "percentile": percentile,
                "threshold": threshold,
                "target_count": target_count,
            }
        )
    return stages


def make_suppressed_hist(source, stage, mode):
    hist = source.Clone("cavern_muon_voxels")
    hist.SetDirectory(0)
    hist.Reset("ICES")
    threshold = stage["threshold"]
    xr, yr, zr = bin_ranges(source)
    kept_bins = 0
    total_content = 0.0
    for ix in xr:
        for iy in yr:
            for iz in zr:
                value = float(source.GetBinContent(ix, iy, iz))
                if mode == "keep_above":
                    new_value = value if value >= threshold and value > 0.0 else 0.0
                else:
                    new_value = max(value - threshold, 0.0)
                if new_value > 0.0:
                    kept_bins += 1
                    total_content += new_value
                    hist.SetBinContent(ix, iy, iz, new_value)
    hist.SetEntries(total_content)
    hist.SetTitle(
        "{title} pedestal={threshold:g} top_fraction={fraction:g} mode={mode}".format(
            title=title_without_axis_directives(source),
            threshold=threshold,
            fraction=stage["top_fraction"],
            mode=mode,
        )
    )
    return hist, kept_bins, total_content


def make_th2(root, name, title, hist3, axis_a, axis_b):
    ax_a = axis(hist3, axis_a)
    ax_b = axis(hist3, axis_b)
    hist2 = root.TH2F(
        name,
        title,
        ax_a.GetNbins(),
        ax_a.GetXmin(),
        ax_a.GetXmax(),
        ax_b.GetNbins(),
        ax_b.GetXmin(),
        ax_b.GetXmax(),
    )
    hist2.GetXaxis().SetTitle(ax_a.GetTitle() or AXIS_NAMES[axis_a])
    hist2.GetYaxis().SetTitle(ax_b.GetTitle() or AXIS_NAMES[axis_b])
    hist2.SetDirectory(0)
    return hist2


def bin_tuple_from_axes(axis_a, bin_a, axis_b, bin_b, axis_c, bin_c):
    bins = [0, 0, 0]
    bins[axis_a] = bin_a
    bins[axis_b] = bin_b
    bins[axis_c] = bin_c
    return bins


def make_projection(root, hist3, axis_a, axis_b, name):
    axis_c = ({0, 1, 2} - {axis_a, axis_b}).pop()
    title = f"{title_without_axis_directives(hist3)} {AXIS_NAMES[axis_a]}{AXIS_NAMES[axis_b]} projection"
    hist2 = make_th2(root, name, title, hist3, axis_a, axis_b)
    for bin_a in range(1, nbins(hist3, axis_a) + 1):
        for bin_b in range(1, nbins(hist3, axis_b) + 1):
            total = 0.0
            for bin_c in range(1, nbins(hist3, axis_c) + 1):
                ix, iy, iz = bin_tuple_from_axes(axis_a, bin_a, axis_b, bin_b, axis_c, bin_c)
                total += float(hist3.GetBinContent(ix, iy, iz))
            if total != 0.0:
                hist2.SetBinContent(bin_a, bin_b, total)
    return hist2


def nearest_bin(axis_obj, value):
    found = axis_obj.FindBin(value)
    if found < 1:
        return 1
    if found > axis_obj.GetNbins():
        return axis_obj.GetNbins()
    return found


def make_slice(root, hist3, axis_a, axis_b, axis_c, fixed_value, name):
    ax_c = axis(hist3, axis_c)
    fixed_bin = nearest_bin(ax_c, fixed_value)
    fixed_center = ax_c.GetBinCenter(fixed_bin)
    unit = hist_axis_unit(hist3)
    title = (
        f"{title_without_axis_directives(hist3)} {AXIS_NAMES[axis_a]}{AXIS_NAMES[axis_b]} "
        f"slice at {AXIS_NAMES[axis_c]}={fixed_center:g} {unit}"
    )
    hist2 = make_th2(root, name, title, hist3, axis_a, axis_b)
    for bin_a in range(1, nbins(hist3, axis_a) + 1):
        for bin_b in range(1, nbins(hist3, axis_b) + 1):
            ix, iy, iz = bin_tuple_from_axes(axis_a, bin_a, axis_b, bin_b, axis_c, fixed_bin)
            value = float(hist3.GetBinContent(ix, iy, iz))
            if value != 0.0:
                hist2.SetBinContent(bin_a, bin_b, value)
    return hist2


def write_views(root, hist3, slice_x, slice_y, slice_z):
    order = view_axis_order(hist3)
    views = [
        make_projection(root, hist3, *order["xy_projection"], "xy_projection"),
        make_projection(root, hist3, *order["xz_projection"], "xz_projection"),
        make_projection(root, hist3, *order["yz_projection"], "yz_projection"),
        make_slice(root, hist3, *order["xy_slice_z"], slice_z, "xy_slice_z"),
        make_slice(root, hist3, *order["xz_slice_y"], slice_y, "xz_slice_y"),
        make_slice(root, hist3, *order["yz_slice_x"], slice_x, "yz_slice_x"),
    ]
    for view in views:
        view.Write()
    return views


def make_views(root, hist3, slice_x, slice_y, slice_z, include_projections=True):
    order = view_axis_order(hist3)
    views = {}
    if include_projections:
        views.update(
            {
                "xy_projection": make_projection(root, hist3, *order["xy_projection"], "xy_projection"),
                "xz_projection": make_projection(root, hist3, *order["xz_projection"], "xz_projection"),
                "yz_projection": make_projection(root, hist3, *order["yz_projection"], "yz_projection"),
            }
        )
    views.update(
        {
            "xy_slice_z": make_slice(root, hist3, *order["xy_slice_z"], slice_z, "xy_slice_z"),
            "xz_slice_y": make_slice(root, hist3, *order["xz_slice_y"], slice_y, "xz_slice_y"),
            "yz_slice_x": make_slice(root, hist3, *order["yz_slice_x"], slice_x, "yz_slice_x"),
        }
    )
    return views


def write_view_dict(views):
    for view in views.values():
        view.Write()


def draw_label(root, text, size=0.055, x=0.12, y=0.92, font=62):
    latex = root.TLatex()
    latex.SetNDC(True)
    latex.SetTextFont(font)
    latex.SetTextSize(size)
    latex.SetTextColor(1)
    latex.DrawLatex(x, y, text)
    return latex


def compass_color(root):
    return getattr(root, "kViolet", 880) + 2


def draw_ndc_text(root, x, y, text, color, size=0.055, align=22):
    latex = root.TLatex()
    latex.SetNDC(True)
    latex.SetTextAlign(align)
    latex.SetTextSize(size)
    latex.SetTextColor(color)
    latex.DrawLatex(x, y, text)
    return latex


def draw_ndc_line(root, x1, y1, x2, y2, color, width=2, style=1):
    line = root.TLine(x1, y1, x2, y2)
    line.SetNDC(True)
    line.SetLineColor(color)
    line.SetLineWidth(width)
    line.SetLineStyle(style)
    line.Draw()
    return line


def draw_ndc_arrow(root, x1, y1, x2, y2, color, width=2):
    arrow = root.TArrow(x1, y1, x2, y2, 0.022, "|>")
    arrow.SetNDC(True)
    arrow.SetLineColor(color)
    arrow.SetFillColor(color)
    arrow.SetLineWidth(width)
    arrow.Draw()
    return arrow


def draw_ndc_marker(root, x, y, style, color, size):
    marker = root.TMarker(x, y, style)
    marker.SetNDC(True)
    marker.SetMarkerStyle(style)
    marker.SetMarkerColor(color)
    marker.SetMarkerSize(size)
    marker.Draw()
    return marker


def plane_orientation_sign(axis_a, axis_b, missing_axis):
    axes = (axis_a, axis_b, missing_axis)
    inversions = 0
    for i in range(3):
        for j in range(i + 1, 3):
            if axes[i] > axes[j]:
                inversions += 1
    return 1 if inversions % 2 == 0 else -1


def draw_out_of_page_marker(root, x, y, outward, color):
    objects = [draw_ndc_marker(root, x, y, 24, color, 1.55)]
    if outward:
        objects.append(draw_ndc_marker(root, x, y, 20, color, 0.55))
    else:
        objects.append(draw_ndc_marker(root, x, y, 5, color, 1.05))
    return objects


def make_polyline_3d(root, points, color, width=2, style=1):
    poly = root.TPolyLine3D(len(points))
    for ipoint, point in enumerate(points):
        poly.SetPoint(ipoint, point[0], point[1], point[2])
    cms_line_style(poly, color, width, style)
    return poly


def axis_unit_vector(axis_index):
    values = [0.0, 0.0, 0.0]
    values[axis_index] = 1.0
    return tuple(values)


def vector_add(*vectors):
    return tuple(sum(vector[index] for vector in vectors) for index in range(3))


def vector_scale(vector, scale):
    return tuple(scale * value for value in vector)


def draw_letter_3d(root, letter, center, right_vector, up_vector, size, color):
    glyphs = {
        "N": [[(-0.5, -0.5), (-0.5, 0.5), (0.5, -0.5), (0.5, 0.5)]],
        "S": [[(0.5, 0.5), (-0.5, 0.5), (-0.5, 0.0), (0.5, 0.0), (0.5, -0.5), (-0.5, -0.5)]],
        "E": [[(0.5, 0.5), (-0.5, 0.5), (-0.5, -0.5), (0.5, -0.5)], [(-0.5, 0.0), (0.3, 0.0)]],
        "W": [[(-0.5, 0.5), (-0.25, -0.5), (0.0, 0.15), (0.25, -0.5), (0.5, 0.5)]],
    }
    objects = []
    for stroke in glyphs[letter]:
        points = [
            vector_add(
                center,
                vector_scale(right_vector, x_value * size),
                vector_scale(up_vector, y_value * size),
            )
            for x_value, y_value in stroke
        ]
        poly = make_polyline_3d(root, points, color, 4, 1)
        poly.Draw("SAME")
        objects.append(poly)
    return objects


def draw_compass_rose_3d(root, hist, cms_x_bearing_deg):
    coordinate_system = hist_coordinate_system(hist)
    color = compass_color(root)
    objects = []

    ranges = [
        (axis(hist, iaxis).GetXmin(), axis(hist, iaxis).GetXmax())
        for iaxis in range(3)
    ]
    spans = [high - low for low, high in ranges]
    up_axis = vertical_axis(coordinate_system)
    axis_bearings = horizontal_axis_bearings(coordinate_system, cms_x_bearing_deg)
    horizontal_axes = sorted(axis_bearings)
    h_axis, d_axis = horizontal_axes

    origin = [0.0, 0.0, 0.0]
    for iaxis, (low, high) in enumerate(ranges):
        if iaxis == up_axis:
            origin[iaxis] = low + 0.08 * spans[iaxis]
        elif iaxis == d_axis:
            origin[iaxis] = high - 0.12 * spans[iaxis]
        else:
            origin[iaxis] = low + 0.33 * spans[iaxis]
    origin = tuple(origin)

    length = 0.15 * min(spans[iaxis] for iaxis in horizontal_axes)
    rose_radius = 0.30 * length
    north_vector = geographic_direction_vector(coordinate_system, cms_x_bearing_deg, 0.0)
    east_vector = geographic_direction_vector(coordinate_system, cms_x_bearing_deg, 90.0)
    up_vector = axis_unit_vector(up_axis)

    rose_points = []
    for ipoint in range(73):
        phi = 2.0 * math.pi * ipoint / 72.0
        rose_points.append(
            vector_add(
                origin,
                vector_scale(north_vector, rose_radius * math.cos(phi)),
                vector_scale(east_vector, rose_radius * math.sin(phi)),
            )
        )
    rose = make_polyline_3d(root, rose_points, color, 2, 1)
    rose.Draw("SAME")
    objects.append(rose)

    for label, bearing in (("N", 0.0), ("E", 90.0)):
        direction = geographic_direction_vector(coordinate_system, cms_x_bearing_deg, bearing)
        spoke = make_polyline_3d(
            root,
            [
                vector_add(origin, vector_scale(direction, -0.45 * length)),
                vector_add(origin, vector_scale(direction, 0.45 * length)),
            ],
            color,
            2,
            2,
        )
        spoke.Draw("SAME")
        objects.append(spoke)

    arrow_tip = vector_add(origin, vector_scale(north_vector, length))
    shaft = make_polyline_3d(root, [origin, arrow_tip], color, 4, 1)
    shaft.Draw("SAME")
    objects.append(shaft)

    head_center = vector_add(origin, vector_scale(north_vector, 0.82 * length))
    for sign in (-1.0, 1.0):
        wing_base = vector_add(
            head_center,
            vector_scale(east_vector, sign * 0.14 * length),
        )
        wing = make_polyline_3d(root, [arrow_tip, wing_base], color, 3, 1)
        wing.Draw("SAME")
        objects.append(wing)

    for label, bearing in CARDINAL_BEARINGS:
        direction = geographic_direction_vector(coordinate_system, cms_x_bearing_deg, bearing)
        label_center = vector_add(
            origin,
            vector_scale(direction, 1.08 * length),
            vector_scale(up_vector, 0.18 * length),
        )
        objects.extend(
            draw_letter_3d(
                root,
                label,
                label_center,
                axis_unit_vector(h_axis),
                up_vector,
                0.24 * length,
                color,
            )
        )

    up_marker = make_polyline_3d(
        root,
        [origin, vector_add(origin, vector_scale(up_vector, 0.35 * length))],
        color,
        1,
        3,
    )
    up_marker.Draw("SAME")
    objects.append(up_marker)
    return objects


def draw_projected_cardinal_2d(root, label, bearing, coordinate_system, axis_a, axis_b, cms_x_bearing_deg, color):
    cx, cy = 0.18, 0.22
    radius = 0.095
    direction = geographic_direction_vector(coordinate_system, cms_x_bearing_deg, bearing)
    vx = direction[axis_a]
    vy = direction[axis_b]
    mag = math.hypot(vx, vy)
    if mag < 0.08:
        return []
    ux = vx / mag
    uy = vy / mag
    x2 = cx + radius * ux
    y2 = cy + radius * uy
    objects = []
    if label == "N":
        objects.append(draw_ndc_arrow(root, cx, cy, x2, y2, color, 2))
    else:
        objects.append(draw_ndc_line(root, cx, cy, x2, y2, color, 1, 2))
    objects.append(draw_ndc_text(root, x2 + 0.020 * ux, y2 + 0.020 * uy, label, color, 0.045, 22))
    return objects


def draw_projected_compass_2d(root, coordinate_system, axis_a, axis_b, cms_x_bearing_deg, color):
    horizontal_axes = set(horizontal_axis_bearings(coordinate_system, cms_x_bearing_deg))
    displayed_axes = {axis_a, axis_b}
    labels = CARDINAL_BEARINGS if displayed_axes == horizontal_axes else (("N", 0.0),)
    objects = []
    for label, bearing in labels:
        objects.extend(
            draw_projected_cardinal_2d(
                root,
                label,
                bearing,
                coordinate_system,
                axis_a,
                axis_b,
                cms_x_bearing_deg,
                color,
            )
        )
    return objects


def draw_orientation_overlay_2d(root, coordinate_system, axis_a, axis_b, cms_x_bearing_deg):
    color = compass_color(root)
    up_axis = vertical_axis(coordinate_system)
    objects = draw_projected_compass_2d(root, coordinate_system, axis_a, axis_b, cms_x_bearing_deg, color)
    cx, cy = 0.18, 0.08

    if axis_b == up_axis:
        objects.append(draw_ndc_arrow(root, cx, cy - 0.005, cx, cy + 0.075, color, 1))
        objects.append(
            draw_ndc_text(root, cx + 0.040, cy + 0.035, axis_direction_label(coordinate_system, up_axis), color, 0.035, 12)
        )
        return objects

    if axis_a == up_axis:
        objects.append(draw_ndc_arrow(root, cx - 0.035, cy, cx + 0.045, cy, color, 1))
        objects.append(
            draw_ndc_text(root, cx + 0.005, cy + 0.045, axis_direction_label(coordinate_system, up_axis), color, 0.035, 22)
        )
        return objects

    sign = plane_orientation_sign(axis_a, axis_b, up_axis)
    outward = sign > 0
    objects.extend(draw_out_of_page_marker(root, cx, cy, outward, color))
    direction = "out" if outward else "in"
    label = f"{axis_direction_label(coordinate_system, up_axis)} {direction}"
    objects.append(draw_ndc_text(root, cx + 0.04, cy - 0.005, label, color, 0.04, 12))
    return objects


def cms_line_style(obj, color, width=1, style=1):
    obj.SetLineColor(color)
    obj.SetLineWidth(width)
    obj.SetLineStyle(style)
    if hasattr(obj, "SetFillStyle"):
        obj.SetFillStyle(0)
    return obj


def cavern_color(root):
    return 1


def cavern_hall_half_width_cm():
    return 0.5 * (CAVERN_HALL_BOX_CM["xmax"] - CAVERN_HALL_BOX_CM["xmin"])


def cavern_hall_half_length_cm():
    return 0.5 * (CAVERN_HALL_BOX_CM["ymax"] - CAVERN_HALL_BOX_CM["ymin"])


def cavern_hall_half_width_at_z_cm(z_value):
    floor = CAVERN_HALL_BOX_CM["zmin"]
    spring = CAVERN_HALL_ARCH_SPRING_Z_CM
    crown = CAVERN_HALL_BOX_CM["zmax"]
    half_width = cavern_hall_half_width_cm()
    if z_value < floor or z_value > crown:
        return None
    if z_value <= spring:
        return half_width
    arch_height = max(crown - spring, 1.0e-9)
    frac = (z_value - spring) / arch_height
    return half_width * math.sqrt(max(1.0 - frac * frac, 0.0))


def cavern_hall_roof_z_at_x_cm(x_value):
    half_width = cavern_hall_half_width_cm()
    if abs(x_value) > half_width:
        return None
    spring = CAVERN_HALL_ARCH_SPRING_Z_CM
    arch_height = CAVERN_HALL_BOX_CM["zmax"] - spring
    frac = x_value / max(half_width, 1.0e-9)
    return spring + arch_height * math.sqrt(max(1.0 - frac * frac, 0.0))


def cavern_hall_cross_section_points(y_value, npoints=49):
    half_width = cavern_hall_half_width_cm()
    floor = CAVERN_HALL_BOX_CM["zmin"]
    spring = CAVERN_HALL_ARCH_SPRING_Z_CM
    arch_height = CAVERN_HALL_BOX_CM["zmax"] - spring
    points = [
        (-half_width, y_value, floor),
        (half_width, y_value, floor),
        (half_width, y_value, spring),
    ]
    for ipoint in range(npoints):
        angle = math.pi * ipoint / (npoints - 1)
        points.append(
            (
                half_width * math.cos(angle),
                y_value,
                spring + arch_height * math.sin(angle),
            )
        )
    points.extend(
        [
            (-half_width, y_value, spring),
            (-half_width, y_value, floor),
        ]
    )
    return points


def cavern_box_hist_ranges(hist, coordinate_system, box):
    corners = []
    for x_value in (box["xmin"], box["xmax"]):
        for y_value in (box["ymin"], box["ymax"]):
            for z_value in (box["zmin"], box["zmax"]):
                corners.append(cavern_cm_to_hist_units_point(hist, coordinate_system, x_value, y_value, z_value))
    return tuple((min(point[axis_index] for point in corners), max(point[axis_index] for point in corners)) for axis_index in range(3))


def draw_rectangle_2d(root, x_low, x_high, y_low, y_high, color, width=2, style=7):
    coords = [
        (x_low, y_low),
        (x_high, y_low),
        (x_high, y_high),
        (x_low, y_high),
        (x_low, y_low),
    ]
    poly = root.TPolyLine(len(coords))
    for ipoint, (x_value, y_value) in enumerate(coords):
        poly.SetPoint(ipoint, x_value, y_value)
    cms_line_style(poly, color, width, style)
    poly.Draw("SAME")
    return poly


def fixed_hist_axis_to_cavern_component(hist, coordinate_system, fixed_axis, fixed_value):
    origin = hist_units_point_to_cavern_cm(hist, coordinate_system, (0.0, 0.0, 0.0))
    unit_point = [0.0, 0.0, 0.0]
    unit_point[fixed_axis] = 1.0
    shifted = hist_units_point_to_cavern_cm(hist, coordinate_system, tuple(unit_point))
    deltas = [abs(shifted[index] - origin[index]) for index in range(3)]
    cavern_component = max(range(3), key=lambda index: deltas[index])

    fixed_point = [0.0, 0.0, 0.0]
    fixed_point[fixed_axis] = fixed_value
    cavern_point = hist_units_point_to_cavern_cm(hist, coordinate_system, tuple(fixed_point))
    return cavern_component, cavern_point[cavern_component]


def draw_cavern_polyline_2d(root, hist, coordinate_system, axis_a, axis_b, cavern_points, color, width=2, style=7):
    poly = root.TPolyLine(len(cavern_points))
    for ipoint, point in enumerate(cavern_points):
        hist_point = cavern_cm_to_hist_units_point(hist, coordinate_system, *point)
        poly.SetPoint(ipoint, hist_point[axis_a], hist_point[axis_b])
    cms_line_style(poly, color, width, style)
    poly.Draw("SAME")
    return poly


def draw_cavern_shaft_circle_2d(root, hist, coordinate_system, axis_a, axis_b, z_value, color):
    points = []
    for ipoint in range(97):
        phi = 2.0 * math.pi * ipoint / 96.0
        points.append(
            (
                CAVERN_SHAFT_CM["x0"] + CAVERN_SHAFT_CM["radius"] * math.cos(phi),
                CAVERN_SHAFT_CM["y0"] + CAVERN_SHAFT_CM["radius"] * math.sin(phi),
                z_value,
            )
        )
    return draw_cavern_polyline_2d(root, hist, coordinate_system, axis_a, axis_b, points, color)


def draw_cavern_shaft_vertical_section_2d(root, hist, coordinate_system, axis_a, axis_b, fixed_component, fixed_value, color):
    center_key = "x0" if fixed_component == 0 else "y0"
    delta = fixed_value - CAVERN_SHAFT_CM[center_key]
    radius = CAVERN_SHAFT_CM["radius"]
    if abs(delta) > radius:
        return None
    other_component = 1 if fixed_component == 0 else 0
    other_center = CAVERN_SHAFT_CM["y0"] if other_component == 1 else CAVERN_SHAFT_CM["x0"]
    half_width = math.sqrt(max(radius * radius - delta * delta, 0.0))
    low = other_center - half_width
    high = other_center + half_width

    def make_point(other_value, z_value):
        coords = [CAVERN_SHAFT_CM["x0"], CAVERN_SHAFT_CM["y0"], z_value]
        coords[fixed_component] = fixed_value
        coords[other_component] = other_value
        return tuple(coords)

    points = [
        make_point(low, CAVERN_SHAFT_CM["zmin"]),
        make_point(high, CAVERN_SHAFT_CM["zmin"]),
        make_point(high, CAVERN_SHAFT_CM["zmax"]),
        make_point(low, CAVERN_SHAFT_CM["zmax"]),
        make_point(low, CAVERN_SHAFT_CM["zmin"]),
    ]
    return draw_cavern_polyline_2d(root, hist, coordinate_system, axis_a, axis_b, points, color)


def draw_cavern_hall_outline_2d(root, hist, coordinate_system, axis_a, axis_b, fixed_axis, fixed_value, color):
    if fixed_axis is None:
        ranges = cavern_box_hist_ranges(hist, coordinate_system, CAVERN_HALL_BOX_CM)
        return [
            draw_rectangle_2d(
                root,
                ranges[axis_a][0],
                ranges[axis_a][1],
                ranges[axis_b][0],
                ranges[axis_b][1],
                color,
                2,
                7,
            )
        ]

    fixed_component, fixed_value_cm = fixed_hist_axis_to_cavern_component(hist, coordinate_system, fixed_axis, fixed_value)
    half_length = cavern_hall_half_length_cm()
    half_width = cavern_hall_half_width_cm()
    floor = CAVERN_HALL_BOX_CM["zmin"]

    if fixed_component == 2:
        half_width_at_z = cavern_hall_half_width_at_z_cm(fixed_value_cm)
        if half_width_at_z is None:
            return []
        points = [
            (-half_width_at_z, -half_length, fixed_value_cm),
            (half_width_at_z, -half_length, fixed_value_cm),
            (half_width_at_z, half_length, fixed_value_cm),
            (-half_width_at_z, half_length, fixed_value_cm),
            (-half_width_at_z, -half_length, fixed_value_cm),
        ]
        return [draw_cavern_polyline_2d(root, hist, coordinate_system, axis_a, axis_b, points, color)]

    if fixed_component == 1:
        if abs(fixed_value_cm) > half_length:
            return []
        return [
            draw_cavern_polyline_2d(
                root,
                hist,
                coordinate_system,
                axis_a,
                axis_b,
                cavern_hall_cross_section_points(fixed_value_cm),
                color,
            )
        ]

    if fixed_component == 0:
        roof = cavern_hall_roof_z_at_x_cm(fixed_value_cm)
        if roof is None:
            return []
        points = [
            (fixed_value_cm, -half_length, floor),
            (fixed_value_cm, half_length, floor),
            (fixed_value_cm, half_length, roof),
            (fixed_value_cm, -half_length, roof),
            (fixed_value_cm, -half_length, floor),
        ]
        return [draw_cavern_polyline_2d(root, hist, coordinate_system, axis_a, axis_b, points, color)]

    return []


def draw_cavern_outline_2d(root, hist, axis_a, axis_b, fixed_axis=None, fixed_value=None):
    coordinate_system = hist_coordinate_system(hist)
    color = cavern_color(root)
    objects = draw_cavern_hall_outline_2d(root, hist, coordinate_system, axis_a, axis_b, fixed_axis, fixed_value, color)

    up_axis = vertical_axis(coordinate_system)
    displayed_axes = {axis_a, axis_b}
    if fixed_axis is None:
        if up_axis not in displayed_axes:
            objects.append(draw_cavern_shaft_circle_2d(root, hist, coordinate_system, axis_a, axis_b, CAVERN_SHAFT_CM["zmin"], color))
        else:
            shaft_box = {
                "xmin": CAVERN_SHAFT_CM["x0"] - CAVERN_SHAFT_CM["radius"],
                "xmax": CAVERN_SHAFT_CM["x0"] + CAVERN_SHAFT_CM["radius"],
                "ymin": CAVERN_SHAFT_CM["y0"] - CAVERN_SHAFT_CM["radius"],
                "ymax": CAVERN_SHAFT_CM["y0"] + CAVERN_SHAFT_CM["radius"],
                "zmin": CAVERN_SHAFT_CM["zmin"],
                "zmax": CAVERN_SHAFT_CM["zmax"],
            }
            shaft_ranges = cavern_box_hist_ranges(hist, coordinate_system, shaft_box)
            objects.append(
                draw_rectangle_2d(
                    root,
                    shaft_ranges[axis_a][0],
                    shaft_ranges[axis_a][1],
                    shaft_ranges[axis_b][0],
                    shaft_ranges[axis_b][1],
                    color,
                    2,
                    7,
                )
            )
        return objects

    fixed_component, fixed_value_cm = fixed_hist_axis_to_cavern_component(hist, coordinate_system, fixed_axis, fixed_value)
    if fixed_component == 2:
        if CAVERN_SHAFT_CM["zmin"] <= fixed_value_cm <= CAVERN_SHAFT_CM["zmax"]:
            objects.append(draw_cavern_shaft_circle_2d(root, hist, coordinate_system, axis_a, axis_b, fixed_value_cm, color))
    elif up_axis in displayed_axes:
        section = draw_cavern_shaft_vertical_section_2d(
            root,
            hist,
            coordinate_system,
            axis_a,
            axis_b,
            fixed_component,
            fixed_value_cm,
            color,
        )
        if section is not None:
            objects.append(section)
    return objects


def draw_cms_outline_2d(root, hist, axis_a, axis_b, fixed_axis=None, fixed_value=None):
    objects = []
    coordinate_system = hist_coordinate_system(hist)
    role_a = axis_role(coordinate_system, axis_a)
    role_b = axis_role(coordinate_system, axis_b)
    role_c = axis_role(coordinate_system, fixed_axis) if fixed_axis is not None else None
    fixed_value_cm = (
        hist_units_to_cm(hist, fixed_value)
        if fixed_axis is not None and fixed_value is not None
        else None
    )

    def u(value_cm):
        return cm_to_hist_units(hist, value_cm)

    if role_a != "beam" and role_b != "beam":
        if role_c == "beam" and abs(fixed_value_cm) > CMS_HALF_LENGTH_CM:
            return objects
        for _name, radius, color in CMS_COMPONENT_RADII:
            ellipse = root.TEllipse(0.0, 0.0, u(radius), u(radius))
            cms_line_style(ellipse, color, 1, 1)
            ellipse.Draw("SAME")
            objects.append(ellipse)
        return objects

    if "beam" not in (role_a, role_b):
        return objects

    beam_axis_on_x = role_a == "beam"
    for _name, radius, color in CMS_COMPONENT_RADII:
        half_transverse = radius
        if role_c != "beam" and fixed_value_cm is not None:
            if abs(fixed_value_cm) > radius:
                continue
            half_transverse = math.sqrt(max(radius * radius - fixed_value_cm * fixed_value_cm, 0.0))
        coords = [
            (-CMS_HALF_LENGTH_CM, -half_transverse),
            (CMS_HALF_LENGTH_CM, -half_transverse),
            (CMS_HALF_LENGTH_CM, half_transverse),
            (-CMS_HALF_LENGTH_CM, half_transverse),
            (-CMS_HALF_LENGTH_CM, -half_transverse),
        ]
        poly = root.TPolyLine(len(coords))
        for ipoint, (beam_value, transverse_value) in enumerate(coords):
            xplot, yplot = (
                (u(beam_value), u(transverse_value))
                if beam_axis_on_x
                else (u(transverse_value), u(beam_value))
            )
            poly.SetPoint(ipoint, xplot, yplot)
        cms_line_style(poly, color, 1, 1)
        poly.Draw("SAME")
        objects.append(poly)

    outer_radius = CMS_COMPONENT_RADII[-1][1]
    half_outer = outer_radius
    if role_c != "beam" and fixed_value_cm is not None:
        if abs(fixed_value_cm) > outer_radius:
            return objects
        half_outer = math.sqrt(max(outer_radius * outer_radius - fixed_value_cm * fixed_value_cm, 0.0))
    for beam_value in CMS_Z_PLANES:
        if beam_axis_on_x:
            line = root.TLine(u(beam_value), u(-half_outer), u(beam_value), u(half_outer))
        else:
            line = root.TLine(u(-half_outer), u(beam_value), u(half_outer), u(beam_value))
        cms_line_style(line, 1, 1, 2)
        line.Draw("SAME")
        objects.append(line)
    return objects


def make_cylinder_ring_3d(root, hist, radius, beam_value, color, npoints=97):
    coordinate_system = hist_coordinate_system(hist)
    poly = root.TPolyLine3D(npoints)
    for ipoint in range(npoints):
        phi = 2.0 * math.pi * ipoint / (npoints - 1)
        transverse_x = radius * math.cos(phi)
        transverse_y = radius * math.sin(phi)
        point = cavern_cm_to_hist_units_point(hist, coordinate_system, transverse_x, beam_value, transverse_y)
        poly.SetPoint(ipoint, point[0], point[1], point[2])
    cms_line_style(poly, color, 1, 1)
    return poly


def make_line_3d(root, x1, y1, z1, x2, y2, z2, color, style=1):
    poly = root.TPolyLine3D(2)
    poly.SetPoint(0, x1, y1, z1)
    poly.SetPoint(1, x2, y2, z2)
    cms_line_style(poly, color, 1, style)
    return poly


def cavern_hist_point(hist, coordinate_system, point):
    return cavern_cm_to_hist_units_point(hist, coordinate_system, point[0], point[1], point[2])


def draw_cavern_line_3d(root, hist, coordinate_system, start, stop, color, width=2, style=7):
    poly = make_polyline_3d(
        root,
        [
            cavern_hist_point(hist, coordinate_system, start),
            cavern_hist_point(hist, coordinate_system, stop),
        ],
        color,
        width,
        style,
    )
    poly.Draw("SAME")
    return poly


def draw_cavern_ring_3d(root, hist, coordinate_system, z_value, color, width=2, style=7):
    points = []
    for ipoint in range(97):
        phi = 2.0 * math.pi * ipoint / 96.0
        points.append(
            cavern_hist_point(
                hist,
                coordinate_system,
                (
                    CAVERN_SHAFT_CM["x0"] + CAVERN_SHAFT_CM["radius"] * math.cos(phi),
                    CAVERN_SHAFT_CM["y0"] + CAVERN_SHAFT_CM["radius"] * math.sin(phi),
                    z_value,
                ),
            )
        )
    poly = make_polyline_3d(root, points, color, width, style)
    poly.Draw("SAME")
    return poly


def draw_cavern_outline_3d(root, hist):
    coordinate_system = hist_coordinate_system(hist)
    color = cavern_color(root)
    objects = []
    half_length = cavern_hall_half_length_cm()
    half_width = cavern_hall_half_width_cm()
    floor = CAVERN_HALL_BOX_CM["zmin"]
    spring = CAVERN_HALL_ARCH_SPRING_Z_CM
    arch_height = CAVERN_HALL_BOX_CM["zmax"] - spring

    for y_value in (-half_length, 0.0, half_length):
        points = [
            cavern_hist_point(hist, coordinate_system, point)
            for point in cavern_hall_cross_section_points(y_value)
        ]
        poly = make_polyline_3d(root, points, color, 2, 7)
        poly.Draw("SAME")
        objects.append(poly)

    for x_value in (-half_width, half_width):
        for z_value in (floor, spring):
            objects.append(
                draw_cavern_line_3d(
                    root,
                    hist,
                    coordinate_system,
                    (x_value, -half_length, z_value),
                    (x_value, half_length, z_value),
                    color,
                )
            )
    for x_value in (-0.5 * half_width, 0.0, 0.5 * half_width):
        objects.append(
            draw_cavern_line_3d(
                root,
                hist,
                coordinate_system,
                (x_value, -half_length, floor),
                (x_value, half_length, floor),
                color,
            )
        )
    for fraction in (0.0, 0.15, 0.30, 0.50, 0.70, 0.85, 1.0):
        angle = math.pi * fraction
        x_value = half_width * math.cos(angle)
        z_value = spring + arch_height * math.sin(angle)
        objects.append(
            draw_cavern_line_3d(
                root,
                hist,
                coordinate_system,
                (x_value, -half_length, z_value),
                (x_value, half_length, z_value),
                color,
            )
        )

    up_axis = vertical_axis(coordinate_system)
    z_low = hist_units_to_cm(hist, axis(hist, up_axis).GetXmin())
    z_high = hist_units_to_cm(hist, axis(hist, up_axis).GetXmax())
    shaft_zmin = max(CAVERN_SHAFT_CM["zmin"], min(z_low, z_high))
    shaft_zmax = min(CAVERN_SHAFT_CM["zmax"], max(z_low, z_high))
    if shaft_zmax > shaft_zmin:
        z_mid = 0.5 * (shaft_zmin + shaft_zmax)
        for z_value in (shaft_zmin, z_mid, shaft_zmax):
            objects.append(draw_cavern_ring_3d(root, hist, coordinate_system, z_value, color))
        for iphi in range(16):
            phi = 2.0 * math.pi * iphi / 16.0
            x_value = CAVERN_SHAFT_CM["x0"] + CAVERN_SHAFT_CM["radius"] * math.cos(phi)
            y_value = CAVERN_SHAFT_CM["y0"] + CAVERN_SHAFT_CM["radius"] * math.sin(phi)
            objects.append(
                draw_cavern_line_3d(
                    root,
                    hist,
                    coordinate_system,
                    (x_value, y_value, shaft_zmin),
                    (x_value, y_value, shaft_zmax),
                    color,
                )
            )
    return objects


def ground_color(root):
    return getattr(root, "kGreen", 416) + 3


def draw_ground_surface_3d(root, hist, coordinate_system, args):
    if args.no_ground_overlay:
        return []
    color = ground_color(root)
    objects = []
    (x_low, x_high), (y_low, y_high) = horizontal_cavern_ranges(hist, coordinate_system)
    step = max(float(args.ground_grid_step_cm), 50.0)

    def values_between(low, high):
        values = []
        nsteps = max(1, int(math.ceil((high - low) / step)))
        for index in range(nsteps + 1):
            values.append(low + (high - low) * index / nsteps)
        return values

    x_values = values_between(x_low, x_high)
    y_values = values_between(y_low, y_high)
    for cavern_y in y_values:
        points = []
        for cavern_x in x_values:
            cavern_z = ground_z_cavern_cm(args, cavern_x, cavern_y)
            points.append(cavern_cm_to_hist_units_point(hist, coordinate_system, cavern_x, cavern_y, cavern_z))
        poly = make_polyline_3d(root, points, color, 2, 2)
        poly.Draw("SAME")
        objects.append(poly)
    for cavern_x in x_values:
        points = []
        for cavern_y in y_values:
            cavern_z = ground_z_cavern_cm(args, cavern_x, cavern_y)
            points.append(cavern_cm_to_hist_units_point(hist, coordinate_system, cavern_x, cavern_y, cavern_z))
        poly = make_polyline_3d(root, points, color, 2, 2)
        poly.Draw("SAME")
        objects.append(poly)
    return objects


def draw_ground_vertical_line_2d(root, hist, coordinate_system, axis_a, axis_b, args):
    up_axis = vertical_axis(coordinate_system)
    if up_axis not in (axis_a, axis_b):
        return []
    color = ground_color(root)
    horizontal_axis = axis_b if axis_a == up_axis else axis_a
    horizontal_range = axis(hist, horizontal_axis)
    npoints = 121
    graph = root.TGraph(npoints)
    for ipoint in range(npoints):
        value = horizontal_range.GetXmin() + (
            horizontal_range.GetXmax() - horizontal_range.GetXmin()
        ) * ipoint / (npoints - 1)
        cavern_x, cavern_y = coordinate_value_to_cavern_horizontal(
            hist, coordinate_system, horizontal_axis, value, 0.0
        )
        hist_point = cavern_cm_to_hist_units_point(
            hist,
            coordinate_system,
            cavern_x,
            cavern_y,
            ground_z_cavern_cm(args, cavern_x, cavern_y),
        )
        if axis_a == up_axis:
            graph.SetPoint(ipoint, hist_point[axis_a], hist_point[axis_b])
        else:
            graph.SetPoint(ipoint, hist_point[axis_a], hist_point[axis_b])
    cms_line_style(graph, color, 3, 2)
    graph.Draw("L SAME")
    label = root.TLatex()
    label.SetTextColor(color)
    label.SetTextSize(0.030)
    label.SetTextFont(42)
    label.DrawLatexNDC(0.60, 0.88, "approx. ground")
    return [graph, label]


def draw_ground_horizontal_contours_2d(root, hist, coordinate_system, axis_a, axis_b, args):
    up_axis = vertical_axis(coordinate_system)
    if up_axis in (axis_a, axis_b):
        return []
    horizontal_axes = set(horizontal_axis_bearings(coordinate_system, cms_x_bearing_degrees(args)))
    if {axis_a, axis_b} != horizontal_axes:
        return []
    color = ground_color(root)
    ax = axis(hist, axis_a)
    ay = axis(hist, axis_b)
    nx = 40
    ny = 40
    draw_ground_horizontal_contours_2d.counter += 1
    contours = root.TH2F(
        f"ground_contours_{hist.GetName()}_{axis_a}_{axis_b}_{draw_ground_horizontal_contours_2d.counter}",
        "",
        nx,
        ax.GetXmin(),
        ax.GetXmax(),
        ny,
        ay.GetXmin(),
        ay.GetXmax(),
    )
    values = []
    for ix in range(1, nx + 1):
        x_value = contours.GetXaxis().GetBinCenter(ix)
        for iy in range(1, ny + 1):
            y_value = contours.GetYaxis().GetBinCenter(iy)
            point = [0.0, 0.0, 0.0]
            point[axis_a] = x_value
            point[axis_b] = y_value
            cavern_x, cavern_y, _cavern_z = hist_units_point_to_cavern_cm(hist, coordinate_system, tuple(point))
            ground_z = ground_z_cavern_cm(args, cavern_x, cavern_y)
            contours.SetBinContent(ix, iy, ground_z)
            values.append(ground_z)
    if not values:
        return []
    low = min(values)
    high = max(values)
    if high <= low:
        levels = array("d", [low])
    else:
        levels = array("d", [low + (high - low) * frac for frac in (0.25, 0.50, 0.75)])
    contours.SetContour(len(levels), levels)
    contours.SetLineColor(color)
    contours.SetLineWidth(2)
    contours.SetLineStyle(2)
    contours.SetStats(False)
    contours.Draw("CONT3 SAME")
    return [contours]


draw_ground_horizontal_contours_2d.counter = 0


def draw_ground_overlay_2d(root, hist, coordinate_system, axis_a, axis_b, args):
    if args.no_ground_overlay:
        return []
    objects = draw_ground_vertical_line_2d(root, hist, coordinate_system, axis_a, axis_b, args)
    if objects:
        return objects
    return draw_ground_horizontal_contours_2d(root, hist, coordinate_system, axis_a, axis_b, args)


def make_beam_line_3d(root, hist, transverse_x, transverse_y, beam_min, beam_max, color, style=1):
    coordinate_system = hist_coordinate_system(hist)
    start = cavern_cm_to_hist_units_point(hist, coordinate_system, transverse_x, beam_min, transverse_y)
    stop = cavern_cm_to_hist_units_point(hist, coordinate_system, transverse_x, beam_max, transverse_y)
    return make_line_3d(root, *start, *stop, color, style)


def draw_cms_outline_3d(root, hist):
    objects = []
    for _name, radius, color in CMS_COMPONENT_RADII:
        for beam_value in (-CMS_HALF_LENGTH_CM, 0.0, CMS_HALF_LENGTH_CM):
            poly = make_cylinder_ring_3d(root, hist, radius, beam_value, color)
            poly.Draw("SAME")
            objects.append(poly)
    outer_radius = CMS_COMPONENT_RADII[-1][1]
    outer_color = CMS_COMPONENT_RADII[-1][2]
    for iphi in range(16):
        phi = 2.0 * math.pi * iphi / 16.0
        x = outer_radius * math.cos(phi)
        transverse_y = outer_radius * math.sin(phi)
        poly = make_beam_line_3d(
            root, hist, x, transverse_y, -CMS_HALF_LENGTH_CM, CMS_HALF_LENGTH_CM, outer_color
        )
        poly.Draw("SAME")
        objects.append(poly)
    for beam_value in CMS_Z_PLANES:
        poly = make_cylinder_ring_3d(root, hist, outer_radius, beam_value, 1)
        poly.SetLineStyle(2)
        poly.Draw("SAME")
        objects.append(poly)
    beam = make_beam_line_3d(root, hist, 0.0, 0.0, -1250.0, 1250.0, 1, 2)
    beam.Draw("SAME")
    objects.append(beam)
    return objects


def draw_hist_for_canvas(root, hist, object_name, draw_cms_outline, draw_orientation, cms_x_bearing_deg, args):
    hist = plot_hist(root, hist, f"{hist.GetName()}_{object_name}_plot", args)
    hist.SetStats(False)
    coordinate_system = hist_coordinate_system(hist)
    objects = [hist]
    if object_name == "cavern_muon_voxels":
        hist.Draw("BOX2Z")
        objects.extend(draw_ground_surface_3d(root, hist, coordinate_system, args))
        if not getattr(args, "no_cavern_outline", False):
            objects.extend(draw_cavern_outline_3d(root, hist))
        if draw_cms_outline:
            objects.extend(draw_cms_outline_3d(root, hist))
        if draw_orientation:
            objects.extend(draw_compass_rose_3d(root, hist, cms_x_bearing_deg))
        return objects
    hist.Draw("COLZ")
    view_axes = view_axis_order(hist)[object_name]
    displayed_axes = view_axes[:2]
    fixed_axis = view_axes[2] if len(view_axes) > 2 else None
    fixed_value = None
    if fixed_axis is not None:
        fixed_value = {"xy_slice_z": args.slice_z, "xz_slice_y": args.slice_y, "yz_slice_x": args.slice_x}[
            object_name
        ]
        fixed_value = cm_to_hist_units(hist, fixed_value)
    if not getattr(args, "no_cavern_outline", False):
        objects.extend(draw_cavern_outline_2d(root, hist, *displayed_axes, fixed_axis, fixed_value))
    if draw_cms_outline and not args.no_cms_outline_2d:
        objects.extend(draw_cms_outline_2d(root, hist, *displayed_axes, fixed_axis, fixed_value))
    objects.extend(draw_ground_overlay_2d(root, hist, coordinate_system, *displayed_axes, args))
    if draw_orientation:
        objects.extend(draw_orientation_overlay_2d(root, coordinate_system, *displayed_axes, cms_x_bearing_deg))
    return objects


def make_comparison_canvas(root, object_name, stage_objects, args):
    canvas = root.TCanvas(
        f"canvas_{object_name}",
        f"{object_name} pedestal comparison",
        args.canvas_width,
        args.canvas_height,
    )
    canvas.Divide(args.canvas_columns, args.canvas_rows, 0.001, 0.001)
    keepalive = []
    max_pads = args.canvas_columns * args.canvas_rows
    cms_bearing = cms_x_bearing_degrees(args)
    for ipad, item in enumerate(stage_objects[:max_pads], start=1):
        pad = canvas.cd(ipad)
        pad.SetRightMargin(0.16 if object_name != "cavern_muon_voxels" else 0.08)
        pad.SetLeftMargin(0.12)
        pad.SetBottomMargin(0.11)
        pad.SetTopMargin(0.08)
        if object_name == "cavern_muon_voxels":
            pad.SetTheta(25)
            pad.SetPhi(35)
        label = item["label"]
        if (
            object_name == "cavern_muon_voxels"
            and item["name"] == "source"
            and not args.draw_source_3d_pad
        ):
            keepalive.append(draw_label(root, f"{label} (not drawn)", size=0.040))
            continue
        hist = item["objects"][object_name]
        drawn = draw_hist_for_canvas(
            root,
            hist,
            object_name,
            not args.no_cms_outline,
            not args.no_orientation_overlay,
            cms_bearing,
            args,
        )
        keepalive.extend(drawn)
        keepalive.append(draw_label(root, label))
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def make_single_object_canvas(root, object_name, item, args):
    canvas = root.TCanvas(
        f"canvas_{item['name']}_{object_name}",
        f"{item['name']} {object_name}",
        args.canvas_width,
        args.canvas_height,
    )
    pad = canvas.cd()
    pad.SetRightMargin(0.16 if object_name != "cavern_muon_voxels" else 0.08)
    pad.SetLeftMargin(0.12)
    pad.SetBottomMargin(0.11)
    pad.SetTopMargin(0.08)
    if object_name == "cavern_muon_voxels":
        pad.SetTheta(25)
        pad.SetPhi(35)
    keepalive = draw_hist_for_canvas(
        root,
        item["objects"][object_name],
        object_name,
        not args.no_cms_outline,
        not args.no_orientation_overlay,
        cms_x_bearing_degrees(args),
        args,
    )
    canvas.cd()
    canvas._codex_keepalive = keepalive
    return canvas


def make_stage_summary(root, stages):
    tree = root.TTree("pedestal_stages", "pedestal suppression stage metadata")
    index = array("i", [0])
    exponent = array("i", [0])
    top_fraction = array("d", [0.0])
    percentile = array("d", [0.0])
    threshold = array("d", [0.0])
    target_count = array("L", [0])
    kept_bins = array("L", [0])
    integral = array("d", [0.0])
    tree.Branch("index", index, "index/I")
    tree.Branch("exponent", exponent, "exponent/I")
    tree.Branch("top_fraction", top_fraction, "top_fraction/D")
    tree.Branch("percentile", percentile, "percentile/D")
    tree.Branch("threshold", threshold, "threshold/D")
    tree.Branch("target_count", target_count, "target_count/l")
    tree.Branch("kept_bins", kept_bins, "kept_bins/l")
    tree.Branch("integral", integral, "integral/D")
    for stage in stages:
        index[0] = stage["index"]
        exponent[0] = stage["exponent"]
        top_fraction[0] = stage["top_fraction"]
        percentile[0] = stage["percentile"]
        threshold[0] = stage["threshold"]
        target_count[0] = stage["target_count"]
        kept_bins[0] = stage.get("kept_bins", 0)
        integral[0] = stage.get("integral", 0.0)
        tree.Fill()
    return tree


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
        raise RuntimeError(f"Could not open input file {args.input}")
    source_hist = require_th3(source_file, args.hist_name)
    source_file.Close()

    values = collect_scores(source_hist, args.include_zero_bins)
    if args.exponential_stages:
        stages = make_exponential_stages(values, args.first_exponent, args.nstages)
    else:
        stages = make_stages_from_fractions(values, parse_top_fractions(args.top_fractions))

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create output file {args.output}")

    stage_objects = []
    source_item = None

    include_source_in_canvases = (
        args.write_source and not args.exclude_source_from_canvases
    ) or args.include_source_in_canvases

    if args.write_source:
        source_display_hist = plot_hist(ROOT, source_hist, args.hist_name, args)
        source_views = make_views(
            ROOT,
            source_display_hist,
            cm_to_hist_units(source_display_hist, args.slice_x),
            cm_to_hist_units(source_display_hist, args.slice_y),
            cm_to_hist_units(source_display_hist, args.slice_z),
            include_projections=not args.skip_projection_canvases,
        )
        source_item = {
            "name": "source",
            "label": "no cut",
            "objects": {"cavern_muon_voxels": source_display_hist, **source_views},
        }
        if include_source_in_canvases:
            stage_objects.append(source_item)
        source_dir = output.mkdir("source")
        source_dir.cd()
        source_display_hist.Write(args.hist_name)
        write_view_dict(source_views)
        output.cd()

    for stage in stages:
        hist, kept_bins, integral = make_suppressed_hist(source_hist, stage, args.mode)
        display_hist = plot_hist(ROOT, hist, args.hist_name, args)
        views = make_views(
            ROOT,
            display_hist,
            cm_to_hist_units(display_hist, args.slice_x),
            cm_to_hist_units(display_hist, args.slice_y),
            cm_to_hist_units(display_hist, args.slice_z),
            include_projections=not args.skip_projection_canvases,
        )
        stage["kept_bins"] = kept_bins
        stage["integral"] = integral
        stage_objects.append(
            {
                "name": stage["name"],
                "label": f"cut {stage['index']} p{stage['percentile']:.4g} thr {stage['threshold']:g}",
                "objects": {"cavern_muon_voxels": display_hist, **views},
            }
        )
        if not args.skip_stage_dirs:
            directory = output.mkdir(stage["name"])
            directory.cd()
            display_hist.Write(args.hist_name)
            write_view_dict(views)
            output.cd()

    summary = make_stage_summary(ROOT, stages)
    output.cd()
    summary.Write()
    if not args.skip_canvases:
        object_names = ["cavern_muon_voxels"]
        if not args.skip_projection_canvases:
            object_names.extend(["xy_projection", "xz_projection", "yz_projection"])
        if not args.skip_slice_canvases:
            object_names.extend(["xy_slice_z", "xz_slice_y", "yz_slice_x"])
        for object_name in object_names:
            if source_item is not None and not args.skip_source_canvases:
                source_dir.cd()
                source_canvas = make_single_object_canvas(ROOT, object_name, source_item, args)
                source_canvas.Write()
                export_canvas(source_canvas, args.png_dir)
                output.cd()
            canvas = make_comparison_canvas(ROOT, object_name, stage_objects, args)
            canvas.Write()
            export_canvas(canvas, args.png_dir)
    output.Close()

    print(f"Wrote {args.output}")
    if args.png_dir:
        print(f"Exported PNGs to {args.png_dir}")
    for stage in stages:
        print(
            "{name}: top_fraction={fraction:g} percentile={percentile:.10g} "
            "threshold={threshold:g} kept_bins={kept_bins} integral={integral:g}".format(
                name=stage["name"],
                fraction=stage["top_fraction"],
                percentile=stage["percentile"],
                threshold=stage["threshold"],
                kept_bins=stage["kept_bins"],
                integral=stage["integral"],
            )
        )


if __name__ == "__main__":
    main()
