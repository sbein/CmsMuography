#!/usr/bin/env python3
"""Build a TH3 voxel line-integral map from CMS cosmic AOD muons/tracks."""

import argparse
import json
import math
import os
from array import array


DEFAULT_VOLUMES = {
    "cms": {
        "xmin": -2000.0,
        "xmax": 2000.0,
        "ymin": -1500.0,
        "ymax": 12000.0,
        "zmin": -3000.0,
        "zmax": 3000.0,
    },
    "cavern": {
        "xmin": -2000.0,
        "xmax": 2000.0,
        "ymin": -3000.0,
        "ymax": 3000.0,
        "zmin": -1500.0,
        "zmax": 12000.0,
    },
}


def resolve_args(args):
    defaults = DEFAULT_VOLUMES[args.coordinate_system]
    for key, value in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
    if args.trace_axis is None:
        args.trace_axis = "z" if args.coordinate_system == "cavern" else "y"
    return args


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Trace quality reconstructed cosmic muons through a CMS-global "
            "cavern volume and fill one TH3 voxel per crossed bin."
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        action="append",
        default=[],
        help="Input ROOT file, /store path, or text file containing inputs. Can be repeated.",
    )
    parser.add_argument("inputs", nargs="*", help="Additional input files or file lists.")
    parser.add_argument(
        "-o",
        "--output",
        default="cavern_voxel_map.root",
        help="Output ROOT file.",
    )
    parser.add_argument("--max-events", type=int, default=-1, help="Maximum events to process.")
    parser.add_argument(
        "--redirector",
        default="root://cms-xrd-global.cern.ch/",
        help="Redirector prepended to plain /store paths.",
    )
    parser.add_argument(
        "--coordinate-system",
        choices=("cavern", "cms"),
        default="cavern",
        help=(
            "Output map coordinate system. 'cms' keeps CMS global axes. "
            "'cavern' writes x=CMS x, y=-CMS z beam, z=CMS y vertical/up "
            "so the horizontal x-y plane is right-handed when viewed from above."
        ),
    )
    parser.add_argument(
        "--product-type",
        choices=("muon", "track"),
        default="muon",
        help="Read std::vector<reco::Muon> or std::vector<reco::Track>.",
    )
    parser.add_argument(
        "--label",
        default=None,
        help=(
            "EDM label as module,instance,process. Defaults to muons,,RECO for "
            "muons and cosmicMuons,,RECO for tracks."
        ),
    )
    parser.add_argument(
        "--muon-track-methods",
        default="outerTrack,standAloneMuon,globalTrack,innerTrack",
        help="Comma-separated reco::Muon track references to try, in order.",
    )
    parser.add_argument("--xmin", type=float, default=None)
    parser.add_argument("--xmax", type=float, default=None)
    parser.add_argument("--ymin", type=float, default=None)
    parser.add_argument("--ymax", type=float, default=None)
    parser.add_argument("--zmin", type=float, default=None)
    parser.add_argument("--zmax", type=float, default=None)
    parser.add_argument(
        "--trace-axis",
        choices=("x", "y", "z"),
        default=None,
        help=(
            "Axis, in the selected output coordinate system, used to orient/extrapolate "
            "candidate trajectories. Defaults to z for cavern coordinates and y for CMS."
        ),
    )
    parser.add_argument(
        "--trace-sign",
        choices=("positive", "negative"),
        default="positive",
        help="Direction along trace-axis used for extrapolation.",
    )
    parser.add_argument(
        "--allow-negative-trace-side",
        action="store_true",
        help="Do not clip the tracing volume to the selected positive/negative half-space.",
    )
    parser.add_argument(
        "--voxel-size",
        type=float,
        default=100.0,
        help="Nominal cubic voxel size in cm; axis bin counts are rounded up.",
    )
    parser.add_argument("--nx", type=int, default=0, help="Override x-axis bin count.")
    parser.add_argument("--ny", type=int, default=0, help="Override y-axis bin count.")
    parser.add_argument("--nz", type=int, default=0, help="Override z-axis bin count.")
    parser.add_argument("--min-p", type=float, default=1.0, help="Minimum track momentum in GeV.")
    parser.add_argument("--min-pt", type=float, default=0.0, help="Minimum track pT in GeV.")
    parser.add_argument(
        "--min-abs-u-trace",
        "--min-abs-uy",
        dest="min_abs_uy",
        type=float,
        default=0.05,
        help="Minimum absolute direction cosine along the selected trace axis.",
    )
    parser.add_argument(
        "--min-valid-hits",
        type=int,
        default=0,
        help="Minimum track valid hits.",
    )
    parser.add_argument(
        "--min-valid-muon-hits",
        type=int,
        default=5,
        help="Minimum valid muon-system hits from the track hitPattern.",
    )
    parser.add_argument(
        "--max-normalized-chi2",
        type=float,
        default=1000.0,
        help="Maximum normalized chi2; negative disables this cut.",
    )
    parser.add_argument(
        "--max-voxels-per-track",
        type=int,
        default=20000,
        help="Safety cap for voxel traversal of one track.",
    )
    parser.add_argument(
        "--weight",
        type=float,
        default=1.0,
        help="Weight added to each crossed voxel.",
    )
    parser.add_argument("--verbose", action="store_true", help="Print per-event progress.")
    parser.add_argument(
        "--progress-every",
        type=int,
        default=0,
        help="Print one progress line every N processed events; 0 disables periodic progress.",
    )
    parser.add_argument(
        "--skip-bad-files",
        action="store_true",
        help="Process inputs one at a time and continue past files that cannot be opened/read.",
    )
    return resolve_args(parser.parse_args())


def normalize_input(path, redirector):
    path = path.strip()
    if not path or path.startswith("#"):
        return None
    if path.startswith("/store/"):
        return redirector + path
    return path


def expand_inputs(items, redirector):
    files = []
    for item in items:
        if not item:
            continue
        item = item[1:] if item.startswith("@") else item
        if os.path.exists(item) and not item.endswith(".root"):
            with open(item, encoding="utf-8") as handle:
                for line in handle:
                    normalized = normalize_input(line, redirector)
                    if normalized:
                        files.append(normalized)
        else:
            normalized = normalize_input(item, redirector)
            if normalized:
                files.append(normalized)
    return files


def setup_fwlite():
    import ROOT

    ROOT.PyConfig.IgnoreCommandLineOptions = True
    for lib in (
        "libFWCoreFWLite",
        "libDataFormatsFWLite",
        "libDataFormatsTrackReco",
        "libDataFormatsMuonReco",
    ):
        ROOT.gSystem.Load(lib)
    ROOT.FWLiteEnabler.enable()

    from DataFormats.FWLite import Events, Handle

    return ROOT, Events, Handle


def parse_label(text):
    parts = [part.strip() for part in text.split(",")]
    if len(parts) == 1:
        return (parts[0], "", "RECO")
    if len(parts) == 2:
        return (parts[0], parts[1], "RECO")
    return (parts[0], parts[1], parts[2])


def label_text(label):
    return ",".join(label)


def get_product(event, label, handle):
    try:
        event.getByLabel(label, handle)
        if hasattr(handle, "isValid") and not handle.isValid():
            return None
        return handle.product()
    except Exception:
        return None


def event_ids(event):
    try:
        aux = event.eventAuxiliary()
    except AttributeError:
        aux = event.object().eventAuxiliary()
    return int(aux.run()), int(aux.luminosityBlock()), int(aux.event())


def safe_call(obj, name, default):
    try:
        return getattr(obj, name)()
    except Exception:
        return default


def finite_or(value, default=-999.0):
    try:
        value = float(value)
        return value if math.isfinite(value) else default
    except Exception:
        return default


def safe_hit_count(track, name):
    try:
        return int(getattr(track.hitPattern(), name)())
    except Exception:
        return -1


def choose_muon_track(muon, methods):
    for method in methods:
        method = method.strip()
        if not method:
            continue
        try:
            ref = getattr(muon, method)()
            if ref.isNonnull() and ref.isAvailable():
                return ref.get(), method
        except Exception:
            continue
    return None, ""


def trace_axis_index(axis_name):
    return {"x": 0, "y": 1, "z": 2}[axis_name]


def transform_coordinates(values, coordinate_system):
    if coordinate_system == "cms":
        return values
    return values[0], -values[2], values[1]


def axis_titles(coordinate_system):
    if coordinate_system == "cms":
        return "x CMS [cm]", "y CMS vertical [cm]", "z CMS beam [cm]"
    return (
        "x cavern = CMS x [cm]",
        "y cavern = -CMS z beam [cm]",
        "z cavern = CMS y vertical [cm]",
    )


def vertical_axis_index(coordinate_system):
    return 1 if coordinate_system == "cms" else 2


def coordinate_transform_note(coordinate_system):
    if coordinate_system == "cms":
        return "output x,y,z are CMS global x,y,z; CMS y is vertical and CMS z is the beam axis"
    return "output x=CMS x, y=-CMS z beam axis, z=CMS y vertical/up; this keeps the Z-up cavern frame right-handed"


def oriented_direction(track, trace_axis, trace_sign, coordinate_system):
    px = finite_or(safe_call(track, "px", -999.0))
    py = finite_or(safe_call(track, "py", -999.0))
    pz = finite_or(safe_call(track, "pz", -999.0))
    norm = math.sqrt(px * px + py * py + pz * pz)
    if norm <= 0.0 or not math.isfinite(norm):
        return None
    ux, uy, uz = px / norm, py / norm, pz / norm
    cms_components = (ux, uy, uz)
    components = transform_coordinates(cms_components, coordinate_system)
    axis_component = components[trace_axis_index(trace_axis)]
    should_flip = axis_component < 0.0 if trace_sign == "positive" else axis_component > 0.0
    if should_flip:
        components = tuple(-value for value in components)
        cms_components = tuple(-value for value in cms_components)
    return components, cms_components


def track_point(track):
    point = (
        finite_or(safe_call(track, "vx", -999.0)),
        finite_or(safe_call(track, "vy", -999.0)),
        finite_or(safe_call(track, "vz", -999.0)),
    )
    if any(value <= -998.0 for value in point):
        return None
    return point


def passes_quality(track, direction, args):
    p = finite_or(safe_call(track, "p", -999.0))
    pt = finite_or(safe_call(track, "pt", -999.0))
    norm_chi2 = finite_or(safe_call(track, "normalizedChi2", -999.0))
    valid_hits = int(safe_call(track, "numberOfValidHits", -1))
    valid_muon_hits = safe_hit_count(track, "numberOfValidMuonHits")

    if p < args.min_p:
        return False, "p", p, pt, norm_chi2, valid_hits, valid_muon_hits
    if pt < args.min_pt:
        return False, "pt", p, pt, norm_chi2, valid_hits, valid_muon_hits
    if abs(direction[trace_axis_index(args.trace_axis)]) < args.min_abs_uy:
        return False, f"u{args.trace_axis}", p, pt, norm_chi2, valid_hits, valid_muon_hits
    if valid_hits < args.min_valid_hits:
        return False, "valid_hits", p, pt, norm_chi2, valid_hits, valid_muon_hits
    if valid_muon_hits < args.min_valid_muon_hits:
        return False, "muon_hits", p, pt, norm_chi2, valid_hits, valid_muon_hits
    if args.max_normalized_chi2 >= 0.0 and norm_chi2 >= 0.0 and norm_chi2 > args.max_normalized_chi2:
        return False, "chi2", p, pt, norm_chi2, valid_hits, valid_muon_hits

    return True, "pass", p, pt, norm_chi2, valid_hits, valid_muon_hits


def ray_box_intersection(point, direction, bounds):
    t_enter = -float("inf")
    t_exit = float("inf")
    mins = (bounds["xmin"], bounds["ymin"], bounds["zmin"])
    maxs = (bounds["xmax"], bounds["ymax"], bounds["zmax"])

    for p0, d0, low, high in zip(point, direction, mins, maxs):
        if abs(d0) < 1.0e-12:
            if p0 < low or p0 > high:
                return None
            continue
        t1 = (low - p0) / d0
        t2 = (high - p0) / d0
        t_low, t_high = (t1, t2) if t1 <= t2 else (t2, t1)
        t_enter = max(t_enter, t_low)
        t_exit = min(t_exit, t_high)
        if t_exit < t_enter:
            return None

    if t_exit < 0.0:
        return None
    return max(t_enter, 0.0), t_exit


def bin_counts(args, bounds):
    if args.voxel_size <= 0.0:
        raise ValueError("--voxel-size must be positive")

    nx = args.nx if args.nx > 0 else int(math.ceil((bounds["xmax"] - bounds["xmin"]) / args.voxel_size))
    ny = args.ny if args.ny > 0 else int(math.ceil((bounds["ymax"] - bounds["ymin"]) / args.voxel_size))
    nz = args.nz if args.nz > 0 else int(math.ceil((bounds["zmax"] - bounds["zmin"]) / args.voxel_size))
    if min(nx, ny, nz) <= 0:
        raise ValueError("Invalid voxel binning or volume bounds")
    return nx, ny, nz


def effective_bounds(args):
    bounds = {
        "xmin": args.xmin,
        "xmax": args.xmax,
        "ymin": args.ymin,
        "ymax": args.ymax,
        "zmin": args.zmin,
        "zmax": args.zmax,
    }
    if not args.allow_negative_trace_side:
        low_key = f"{args.trace_axis}min"
        high_key = f"{args.trace_axis}max"
        if args.trace_sign == "positive":
            bounds[low_key] = max(bounds[low_key], 0.0)
        else:
            bounds[high_key] = min(bounds[high_key], 0.0)
    return bounds


def axis_index(value, low, width, nbins):
    index = int(math.floor((value - low) / width))
    if index < 0:
        return 0
    if index >= nbins:
        return nbins - 1
    return index


def fill_voxels(hist, point, direction, t0, t1, bounds, bins, weight, max_voxels):
    nx, ny, nz = bins
    lows = (bounds["xmin"], bounds["ymin"], bounds["zmin"])
    highs = (bounds["xmax"], bounds["ymax"], bounds["zmax"])
    widths = tuple((highs[i] - lows[i]) / bins[i] for i in range(3))

    eps = 1.0e-7
    t = min(t1, t0 + eps)
    pos = tuple(point[i] + t * direction[i] for i in range(3))
    indices = [
        axis_index(pos[0], lows[0], widths[0], nx),
        axis_index(pos[1], lows[1], widths[1], ny),
        axis_index(pos[2], lows[2], widths[2], nz),
    ]

    t_max = []
    t_delta = []
    steps = []
    for axis in range(3):
        d0 = direction[axis]
        if abs(d0) < 1.0e-12:
            steps.append(0)
            t_max.append(float("inf"))
            t_delta.append(float("inf"))
            continue
        if d0 > 0.0:
            steps.append(1)
            boundary = lows[axis] + (indices[axis] + 1) * widths[axis]
        else:
            steps.append(-1)
            boundary = lows[axis] + indices[axis] * widths[axis]
        t_boundary = (boundary - point[axis]) / d0
        delta = widths[axis] / abs(d0)
        while t_boundary <= t:
            t_boundary += delta
        t_max.append(t_boundary)
        t_delta.append(delta)

    filled = 0
    truncated = False
    while (
        0 <= indices[0] < nx
        and 0 <= indices[1] < ny
        and 0 <= indices[2] < nz
        and t <= t1 + eps
    ):
        global_bin = hist.GetBin(indices[0] + 1, indices[1] + 1, indices[2] + 1)
        hist.AddBinContent(global_bin, weight)
        filled += 1
        if filled >= max_voxels:
            truncated = True
            break

        next_t = min(t_max)
        if next_t > t1:
            break
        tied_axes = [axis for axis in range(3) if abs(t_max[axis] - next_t) < 1.0e-8]
        t = next_t
        for axis in tied_axes:
            indices[axis] += steps[axis]
            t_max[axis] += t_delta[axis]

    return filled, truncated


def make_scalar(tree, name, typecode, leaf):
    value = array(typecode, [0])
    tree.Branch(name, value, f"{name}/{leaf}")
    return value


def make_cutflow(ROOT):
    labels = [
        "events",
        "candidates",
        "no_track",
        "quality_fail",
        "geometry_fail",
        "accepted_tracks",
        "voxel_fills",
        "truncated_tracks",
    ]
    hist = ROOT.TH1D("cutflow", "voxel map cutflow", len(labels), 0, len(labels))
    for ibin, label in enumerate(labels, start=1):
        hist.GetXaxis().SetBinLabel(ibin, label)
    return hist, {label: ibin for ibin, label in enumerate(labels, start=1)}


def add_cut(cutflow, cut_bins, label, amount=1.0):
    ibin = cut_bins[label]
    cutflow.SetBinContent(ibin, cutflow.GetBinContent(ibin) + amount)


def main():
    args = parse_args()
    inputs = expand_inputs(args.input + args.inputs, args.redirector)
    if not inputs:
        raise SystemExit("No input files were provided.")

    if args.label is None:
        label = ("muons", "", "RECO") if args.product_type == "muon" else ("cosmicMuons", "", "RECO")
    else:
        label = parse_label(args.label)
    bounds = effective_bounds(args)
    bins = bin_counts(args, bounds)
    muon_methods = [method.strip() for method in args.muon_track_methods.split(",") if method.strip()]

    ROOT, Events, Handle = setup_fwlite()
    handle_type = "std::vector<reco::Muon>" if args.product_type == "muon" else "std::vector<reco::Track>"
    handle = Handle(handle_type)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")

    nx, ny, nz = bins
    x_title, y_title, z_title = axis_titles(args.coordinate_system)
    voxel_counts = ROOT.TH3F(
        "cavern_muon_voxels",
        (
            f"cosmic-muon voxel crossings along {args.trace_sign} "
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

    tracks = ROOT.TTree("tracks", "accepted tracks used in the voxel map")
    tr_run = make_scalar(tracks, "run", "i", "I")
    tr_lumi = make_scalar(tracks, "lumi", "i", "I")
    tr_event = make_scalar(tracks, "event", "L", "l")
    tr_index = make_scalar(tracks, "index", "i", "I")
    tr_p = make_scalar(tracks, "p", "f", "F")
    tr_pt = make_scalar(tracks, "pt", "f", "F")
    tr_eta = make_scalar(tracks, "eta", "f", "F")
    tr_phi = make_scalar(tracks, "phi", "f", "F")
    tr_vx = make_scalar(tracks, "vx", "f", "F")
    tr_vy = make_scalar(tracks, "vy", "f", "F")
    tr_vz = make_scalar(tracks, "vz", "f", "F")
    tr_cms_vx = make_scalar(tracks, "cms_vx", "f", "F")
    tr_cms_vy = make_scalar(tracks, "cms_vy", "f", "F")
    tr_cms_vz = make_scalar(tracks, "cms_vz", "f", "F")
    tr_ux = make_scalar(tracks, "ux", "f", "F")
    tr_uy = make_scalar(tracks, "uy", "f", "F")
    tr_uz = make_scalar(tracks, "uz", "f", "F")
    tr_cms_ux = make_scalar(tracks, "cms_ux", "f", "F")
    tr_cms_uy = make_scalar(tracks, "cms_uy", "f", "F")
    tr_cms_uz = make_scalar(tracks, "cms_uz", "f", "F")
    tr_abs_uy = make_scalar(tracks, "abs_uy", "f", "F")
    tr_abs_u_vertical = make_scalar(tracks, "abs_u_vertical", "f", "F")
    tr_abs_u_trace = make_scalar(tracks, "abs_u_trace", "f", "F")
    tr_t_enter = make_scalar(tracks, "t_enter", "f", "F")
    tr_t_exit = make_scalar(tracks, "t_exit", "f", "F")
    tr_path_length = make_scalar(tracks, "path_length", "f", "F")
    tr_n_voxels = make_scalar(tracks, "n_voxels", "i", "I")
    tr_valid_hits = make_scalar(tracks, "valid_hits", "i", "I")
    tr_valid_muon_hits = make_scalar(tracks, "valid_muon_hits", "i", "I")
    tr_norm_chi2 = make_scalar(tracks, "norm_chi2", "f", "F")
    tr_track_method = make_scalar(tracks, "track_method", "i", "I")

    track_method_codes = {method: code for code, method in enumerate(muon_methods, start=1)}
    if args.product_type == "track":
        track_method_codes = {"track": 1}

    total_voxel_fills = 0
    accepted_tracks = 0
    processed_events = 0
    skipped_inputs = []

    event_sources = [([input_file], input_file) for input_file in inputs] if args.skip_bad_files else [(inputs, "all inputs")]
    for source_inputs, source_name in event_sources:
        if 0 <= args.max_events <= processed_events:
            break
        try:
            events = Events(source_inputs)
            for event in events:
                if 0 <= args.max_events <= processed_events:
                    break

                run, lumi, event_number = event_ids(event)
                add_cut(cutflow, cut_bins, "events")
                product = get_product(event, label, handle)
                event_candidates = int(product.size()) if product else 0
                event_accepted = 0
                event_voxels = 0

                if product:
                    for index, candidate in enumerate(product):
                        add_cut(cutflow, cut_bins, "candidates")
                        if args.product_type == "muon":
                            track, track_method = choose_muon_track(candidate, muon_methods)
                        else:
                            track, track_method = candidate, "track"
                        if not track:
                            add_cut(cutflow, cut_bins, "no_track")
                            continue

                        directions = oriented_direction(
                            track, args.trace_axis, args.trace_sign, args.coordinate_system
                        )
                        point_cms = track_point(track)
                        if directions is None or point_cms is None:
                            add_cut(cutflow, cut_bins, "quality_fail")
                            continue
                        direction, direction_cms = directions
                        point = transform_coordinates(point_cms, args.coordinate_system)

                        passed, _reason, p, pt, norm_chi2, valid_hits, valid_muon_hits = passes_quality(
                            track, direction, args
                        )
                        if not passed:
                            add_cut(cutflow, cut_bins, "quality_fail")
                            continue

                        segment = ray_box_intersection(point, direction, bounds)
                        if segment is None:
                            add_cut(cutflow, cut_bins, "geometry_fail")
                            continue

                        t_enter, t_exit = segment
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
                            add_cut(cutflow, cut_bins, "geometry_fail")
                            continue

                        accepted_tracks += 1
                        total_voxel_fills += n_voxels
                        event_accepted += 1
                        event_voxels += n_voxels
                        add_cut(cutflow, cut_bins, "accepted_tracks")
                        add_cut(cutflow, cut_bins, "voxel_fills", n_voxels)
                        if truncated:
                            add_cut(cutflow, cut_bins, "truncated_tracks")

                        tr_run[0] = run
                        tr_lumi[0] = lumi
                        tr_event[0] = event_number
                        tr_index[0] = index
                        tr_p[0] = finite_or(p)
                        tr_pt[0] = finite_or(pt)
                        tr_eta[0] = finite_or(safe_call(track, "eta", -999.0))
                        tr_phi[0] = finite_or(safe_call(track, "phi", -999.0))
                        tr_vx[0], tr_vy[0], tr_vz[0] = point
                        tr_cms_vx[0], tr_cms_vy[0], tr_cms_vz[0] = point_cms
                        tr_ux[0], tr_uy[0], tr_uz[0] = direction
                        tr_cms_ux[0], tr_cms_uy[0], tr_cms_uz[0] = direction_cms
                        tr_abs_uy[0] = abs(direction[1])
                        tr_abs_u_vertical[0] = abs(direction[vertical_axis_index(args.coordinate_system)])
                        tr_abs_u_trace[0] = abs(direction[trace_axis_index(args.trace_axis)])
                        tr_t_enter[0] = finite_or(t_enter)
                        tr_t_exit[0] = finite_or(t_exit)
                        tr_path_length[0] = finite_or(t_exit - t_enter)
                        tr_n_voxels[0] = n_voxels
                        tr_valid_hits[0] = valid_hits
                        tr_valid_muon_hits[0] = valid_muon_hits
                        tr_norm_chi2[0] = finite_or(norm_chi2)
                        tr_track_method[0] = track_method_codes.get(track_method, 0)
                        tracks.Fill()

                processed_events += 1
                if args.verbose or (args.progress_every > 0 and processed_events % args.progress_every == 0):
                    print(
                        "event {n}: run={run} lumi={lumi} event={event} "
                        "candidates={candidates} accepted={accepted} voxels={voxels}".format(
                            n=processed_events,
                            run=run,
                            lumi=lumi,
                            event=event_number,
                            candidates=event_candidates,
                            accepted=event_accepted,
                            voxels=event_voxels,
                        )
                    )
        except Exception as exc:
            if not args.skip_bad_files:
                raise
            message = str(exc).splitlines()[-1] if str(exc).splitlines() else repr(exc)
            skipped_inputs.append({"input": source_name, "error": message})
            print(f"Skipping unreadable input: {source_name} :: {message}")
            continue

    config = {
        "inputs": inputs,
        "product_type": args.product_type,
        "label": label_text(label),
        "muon_track_methods": muon_methods,
        "coordinate_system": args.coordinate_system,
        "coordinate_transform": coordinate_transform_note(args.coordinate_system),
        "bounds_cm": bounds,
        "bins": {"nx": nx, "ny": ny, "nz": nz},
        "voxel_size_request_cm": args.voxel_size,
        "trace_axis": args.trace_axis,
        "trace_sign": args.trace_sign,
        "allow_negative_trace_side": args.allow_negative_trace_side,
        "skip_bad_files": args.skip_bad_files,
        "skipped_inputs": skipped_inputs,
        "quality": {
            "min_p": args.min_p,
            "min_pt": args.min_pt,
            "min_abs_uy": args.min_abs_uy,
            "min_valid_hits": args.min_valid_hits,
            "min_valid_muon_hits": args.min_valid_muon_hits,
            "max_normalized_chi2": args.max_normalized_chi2,
        },
        "charge_handling": "positive and negative reconstructed muons are combined in the same histograms",
        "track_method_codes": track_method_codes,
        "note": (
            "Trajectories are transformed to the configured output coordinate system, "
            "then oriented along the configured trace axis/sign before clipping to the "
            "configured volume."
        ),
    }

    voxel_counts.SetEntries(total_voxel_fills)
    output.cd()
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "config_json")
    voxel_counts.Write()
    cutflow.Write()
    tracks.Write()
    output.Close()

    print(
        "Wrote {output} with {events} events, {tracks} accepted tracks, and {voxels} voxel fills.".format(
            output=args.output,
            events=processed_events,
            tracks=accepted_tracks,
            voxels=total_voxel_fills,
        )
    )


if __name__ == "__main__":
    main()
