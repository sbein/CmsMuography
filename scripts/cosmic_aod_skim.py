#!/usr/bin/env python3
"""FWLite/PyROOT baseline skim for CMS cosmic AOD muography studies."""

import argparse
import math
import os
import sys
from array import array


MUON_COLLECTIONS = [
    ("muons", "", "RECO"),
    ("muons1Leg", "", "RECO"),
    ("muonsNoRPC", "", "RECO"),
    ("muonsWitht0Correction", "", "RECO"),
    ("splitMuons", "", "RECO"),
    ("muonsBeamHaloEndCapsOnly", "", "RECO"),
    ("lhcSTAMuons", "", "RECO"),
]

TRACK_COLLECTIONS = [
    ("globalCosmicMuons", "", "RECO"),
    ("globalCosmicMuons1Leg", "", "RECO"),
    ("globalCosmicMuonsNoRPC", "", "RECO"),
    ("globalCosmicMuonsWitht0Correction", "", "RECO"),
    ("globalCosmicSplitMuons", "", "RECO"),
    ("cosmicMuons", "", "RECO"),
    ("cosmicMuons1Leg", "", "RECO"),
    ("cosmicMuonsNoRPC", "", "RECO"),
    ("cosmicMuonsWitht0Correction", "", "RECO"),
    ("cosmicMuonsEndCapsOnly", "", "RECO"),
    ("standAloneMuons", "", "RECO"),
    ("standAloneMuons", "UpdatedAtVtx", "RECO"),
    ("regionalCosmicTracks", "", "RECO"),
    ("cosmictrackfinderP5", "", "RECO"),
    ("ctfWithMaterialTracksP5", "", "RECO"),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create a compact ROOT tree from CMS cosmic AOD muon/track collections."
    )
    parser.add_argument(
        "-i",
        "--input",
        action="append",
        default=[],
        help="Input ROOT file, /store path, or text file containing inputs. Can be repeated.",
    )
    parser.add_argument(
        "inputs",
        nargs="*",
        help="Additional input ROOT files, /store paths, or text files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="cosmic_muography.root",
        help="Output ROOT file.",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=-1,
        help="Maximum events to process; negative means all.",
    )
    parser.add_argument(
        "--redirector",
        default="root://cms-xrd-global.cern.ch/",
        help="Redirector prepended to plain /store paths.",
    )
    parser.add_argument(
        "--track-collections",
        default=",".join(label[0] for label in TRACK_COLLECTIONS),
        help="Comma-separated track collection module labels to keep, or 'all'.",
    )
    parser.add_argument(
        "--muon-collections",
        default=",".join(label[0] for label in MUON_COLLECTIONS),
        help="Comma-separated muon collection module labels to keep, or 'all'.",
    )
    parser.add_argument("--skip-muons", action="store_true", help="Do not read reco::Muon collections.")
    parser.add_argument("--skip-tracks", action="store_true", help="Do not read reco::Track collections.")
    parser.add_argument("--verbose", action="store_true", help="Print per-event collection counts.")
    return parser.parse_args()


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


def selected_collections(defaults, selection):
    wanted = {item.strip() for item in selection.split(",") if item.strip()}
    if not wanted or "all" in wanted:
        return list(defaults)
    return [label for label in defaults if label[0] in wanted]


def clean_name(label):
    module, instance, _process = label
    return module if not instance else f"{module}_{instance}"


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


def safe_hit_count(track, name):
    try:
        return int(getattr(track.hitPattern(), name)())
    except Exception:
        return -1


def finite_or(value, default=-999.0):
    try:
        value = float(value)
        return value if math.isfinite(value) else default
    except Exception:
        return default


def best_muon_track(muon):
    for method in ("globalTrack", "innerTrack", "outerTrack", "standAloneMuon"):
        try:
            ref = getattr(muon, method)()
            if ref.isNonnull() and ref.isAvailable():
                return ref.get()
        except Exception:
            continue
    return None


def make_scalar(tree, name, typecode, leaf):
    value = array(typecode, [0])
    tree.Branch(name, value, f"{name}/{leaf}")
    return value


def main():
    args = parse_args()
    inputs = expand_inputs(args.input + args.inputs, args.redirector)
    if not inputs:
        raise SystemExit("No input files were provided.")

    track_collections = [] if args.skip_tracks else selected_collections(TRACK_COLLECTIONS, args.track_collections)
    muon_collections = [] if args.skip_muons else selected_collections(MUON_COLLECTIONS, args.muon_collections)
    all_collections = [("muon", label) for label in muon_collections] + [
        ("track", label) for label in track_collections
    ]

    ROOT, Events, Handle = setup_fwlite()
    events = Events(inputs)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")

    events_tree = ROOT.TTree("events", "event-level cosmic collection counts")
    tracks_tree = ROOT.TTree("tracks", "per-candidate cosmic muon and track quantities")

    ev_run = make_scalar(events_tree, "run", "i", "I")
    ev_lumi = make_scalar(events_tree, "lumi", "i", "I")
    ev_event = make_scalar(events_tree, "event", "L", "l")
    count_branches = {}
    for kind, label in all_collections:
        branch = f"n_{kind}_{clean_name(label)}"
        count_branches[(kind, label)] = make_scalar(events_tree, branch, "i", "I")

    tr_run = make_scalar(tracks_tree, "run", "i", "I")
    tr_lumi = make_scalar(tracks_tree, "lumi", "i", "I")
    tr_event = make_scalar(tracks_tree, "event", "L", "l")
    tr_collection = make_scalar(tracks_tree, "collection", "i", "I")
    tr_kind = make_scalar(tracks_tree, "kind", "i", "I")
    tr_index = make_scalar(tracks_tree, "index", "i", "I")
    tr_pt = make_scalar(tracks_tree, "pt", "f", "F")
    tr_eta = make_scalar(tracks_tree, "eta", "f", "F")
    tr_phi = make_scalar(tracks_tree, "phi", "f", "F")
    tr_p = make_scalar(tracks_tree, "p", "f", "F")
    tr_px = make_scalar(tracks_tree, "px", "f", "F")
    tr_py = make_scalar(tracks_tree, "py", "f", "F")
    tr_pz = make_scalar(tracks_tree, "pz", "f", "F")
    tr_ux = make_scalar(tracks_tree, "ux", "f", "F")
    tr_uy = make_scalar(tracks_tree, "uy", "f", "F")
    tr_uz = make_scalar(tracks_tree, "uz", "f", "F")
    tr_abs_uy = make_scalar(tracks_tree, "abs_uy", "f", "F")
    tr_theta_vertical = make_scalar(tracks_tree, "theta_from_vertical_abs", "f", "F")
    tr_charge = make_scalar(tracks_tree, "charge", "i", "I")
    tr_norm_chi2 = make_scalar(tracks_tree, "norm_chi2", "f", "F")
    tr_ndof = make_scalar(tracks_tree, "ndof", "f", "F")
    tr_vx = make_scalar(tracks_tree, "vx", "f", "F")
    tr_vy = make_scalar(tracks_tree, "vy", "f", "F")
    tr_vz = make_scalar(tracks_tree, "vz", "f", "F")
    tr_dxy = make_scalar(tracks_tree, "dxy", "f", "F")
    tr_dz = make_scalar(tracks_tree, "dz", "f", "F")
    tr_valid_hits = make_scalar(tracks_tree, "valid_hits", "i", "I")
    tr_lost_hits = make_scalar(tracks_tree, "lost_hits", "i", "I")
    tr_valid_tracker_hits = make_scalar(tracks_tree, "valid_tracker_hits", "i", "I")
    tr_valid_pixel_hits = make_scalar(tracks_tree, "valid_pixel_hits", "i", "I")
    tr_valid_muon_hits = make_scalar(tracks_tree, "valid_muon_hits", "i", "I")
    tr_valid_dt_hits = make_scalar(tracks_tree, "valid_dt_hits", "i", "I")
    tr_valid_csc_hits = make_scalar(tracks_tree, "valid_csc_hits", "i", "I")
    tr_valid_rpc_hits = make_scalar(tracks_tree, "valid_rpc_hits", "i", "I")
    tr_is_global_muon = make_scalar(tracks_tree, "is_global_muon", "i", "I")
    tr_is_tracker_muon = make_scalar(tracks_tree, "is_tracker_muon", "i", "I")
    tr_is_standalone_muon = make_scalar(tracks_tree, "is_standalone_muon", "i", "I")

    track_handle = Handle("std::vector<reco::Track>")
    muon_handle = Handle("std::vector<reco::Muon>")

    labels_text = []
    for idx, (kind, label) in enumerate(all_collections):
        labels_text.append(f"{idx} kind={kind} label={label[0]} instance={label[1]} process={label[2]}")
    output.WriteObject(ROOT.TObjString("\n".join(labels_text)), "collection_labels")

    def fill_candidate(run, lumi, event_number, collection_index, kind, index, candidate, track=None):
        if track is None:
            track = candidate

        p = finite_or(safe_call(candidate, "p", -999.0))
        px = finite_or(safe_call(candidate, "px", -999.0))
        py = finite_or(safe_call(candidate, "py", -999.0))
        pz = finite_or(safe_call(candidate, "pz", -999.0))
        if p > 0:
            ux, uy, uz = px / p, py / p, pz / p
            abs_uy = abs(uy)
            theta_vertical = math.acos(max(-1.0, min(1.0, abs_uy)))
        else:
            ux, uy, uz, abs_uy, theta_vertical = -999.0, -999.0, -999.0, -999.0, -999.0

        tr_run[0] = run
        tr_lumi[0] = lumi
        tr_event[0] = event_number
        tr_collection[0] = collection_index
        tr_kind[0] = 0 if kind == "muon" else 1
        tr_index[0] = index
        tr_pt[0] = finite_or(safe_call(candidate, "pt", -999.0))
        tr_eta[0] = finite_or(safe_call(candidate, "eta", -999.0))
        tr_phi[0] = finite_or(safe_call(candidate, "phi", -999.0))
        tr_p[0] = p
        tr_px[0] = px
        tr_py[0] = py
        tr_pz[0] = pz
        tr_ux[0] = finite_or(ux)
        tr_uy[0] = finite_or(uy)
        tr_uz[0] = finite_or(uz)
        tr_abs_uy[0] = finite_or(abs_uy)
        tr_theta_vertical[0] = finite_or(theta_vertical)
        tr_charge[0] = int(safe_call(candidate, "charge", 0))

        tr_norm_chi2[0] = finite_or(safe_call(track, "normalizedChi2", -999.0)) if track else -999.0
        tr_ndof[0] = finite_or(safe_call(track, "ndof", -999.0)) if track else -999.0
        tr_vx[0] = finite_or(safe_call(track, "vx", -999.0)) if track else -999.0
        tr_vy[0] = finite_or(safe_call(track, "vy", -999.0)) if track else -999.0
        tr_vz[0] = finite_or(safe_call(track, "vz", -999.0)) if track else -999.0
        tr_dxy[0] = finite_or(safe_call(track, "dxy", -999.0)) if track else -999.0
        tr_dz[0] = finite_or(safe_call(track, "dz", -999.0)) if track else -999.0
        tr_valid_hits[0] = int(safe_call(track, "numberOfValidHits", -1)) if track else -1
        tr_lost_hits[0] = int(safe_call(track, "numberOfLostHits", -1)) if track else -1
        tr_valid_tracker_hits[0] = safe_hit_count(track, "numberOfValidTrackerHits") if track else -1
        tr_valid_pixel_hits[0] = safe_hit_count(track, "numberOfValidPixelHits") if track else -1
        tr_valid_muon_hits[0] = safe_hit_count(track, "numberOfValidMuonHits") if track else -1
        tr_valid_dt_hits[0] = safe_hit_count(track, "numberOfValidMuonDTHits") if track else -1
        tr_valid_csc_hits[0] = safe_hit_count(track, "numberOfValidMuonCSCHits") if track else -1
        tr_valid_rpc_hits[0] = safe_hit_count(track, "numberOfValidMuonRPCHits") if track else -1

        tr_is_global_muon[0] = int(safe_call(candidate, "isGlobalMuon", -1)) if kind == "muon" else -1
        tr_is_tracker_muon[0] = int(safe_call(candidate, "isTrackerMuon", -1)) if kind == "muon" else -1
        tr_is_standalone_muon[0] = int(safe_call(candidate, "isStandAloneMuon", -1)) if kind == "muon" else -1
        tracks_tree.Fill()

    processed = 0
    for iev, event in enumerate(events):
        if 0 <= args.max_events <= iev:
            break

        run, lumi, event_number = event_ids(event)
        ev_run[0] = run
        ev_lumi[0] = lumi
        ev_event[0] = event_number
        tr_counts = {}

        for collection_index, (kind, label) in enumerate(all_collections):
            handle = muon_handle if kind == "muon" else track_handle
            product = get_product(event, label, handle)
            count = int(product.size()) if product else 0
            count_branches[(kind, label)][0] = count
            tr_counts[clean_name(label)] = count

            if not product:
                continue
            for index, candidate in enumerate(product):
                track = best_muon_track(candidate) if kind == "muon" else candidate
                fill_candidate(run, lumi, event_number, collection_index, kind, index, candidate, track)

        events_tree.Fill()
        processed += 1
        if args.verbose:
            counts = " ".join(f"{name}={count}" for name, count in sorted(tr_counts.items()))
            print(f"event {processed}: run={run} lumi={lumi} event={event_number} {counts}")

    output.cd()
    events_tree.Write()
    tracks_tree.Write()
    track_entries = int(tracks_tree.GetEntries())
    output.Close()
    print(f"Wrote {args.output} with {processed} events and {track_entries} candidates.")


if __name__ == "__main__":
    main()
