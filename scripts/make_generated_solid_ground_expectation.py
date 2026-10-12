#!/usr/bin/env python3
"""Generate an independent EXPACS/PARMA solid-ground muon null for CMS cavern voxels."""

import argparse
import bisect
import csv
import json
import math
import os
import random
from array import array

from cosmic_voxel_map import fill_voxels, ray_box_intersection
from make_solid_ground_expectation import GroundArgs, ground_intersection_t, hist_bounds_and_bins
from suppress_voxel_pedestals import (
    CMS_COMPONENT_RADII,
    CMS_HALF_LENGTH_CM,
    DEFAULT_GROUND_COEFFICIENTS_CM,
)


MUON_MASS_MEV = 105.6583755
ELECTRON_MASS_MEV = 0.51099895
BETHE_K_MEV_CM2_G = 0.307075


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Throw independent EXPACS/PARMA cosmic muons over a detector impact plane, "
            "require geometric CMS-cylinder crossings with 100% efficiency, propagate "
            "through a solid-ground range model, and fill an expected voxel TH3."
        )
    )
    parser.add_argument("--template-root", required=True, help="ROOT file whose TH3 defines the output binning.")
    parser.add_argument("--template-hist", default="cavern_muon_voxels")
    parser.add_argument("--flux-table", required=True, help="CSV produced from the EXPACS/PARMA flux grid.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--nthrows", type=int, default=2000000)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--chunk-index", type=int, default=0)
    parser.add_argument("--nchunks", type=int, default=1)
    parser.add_argument("--progress-every", type=int, default=50000)
    parser.add_argument(
        "--sampling-mode",
        choices=("forward-plane", "reverse-cms-volume"),
        default="reverse-cms-volume",
        help=(
            "forward-plane samples x/y on a detector plane and rejects non-CMS crossings. "
            "reverse-cms-volume samples points inside the CMS cylinder and weights by 1/chord length."
        ),
    )
    parser.add_argument("--throw-x-half-width-cm", type=float, default=2500.0)
    parser.add_argument("--throw-y-half-width-cm", type=float, default=2500.0)
    parser.add_argument("--cms-radius-cm", type=float, default=CMS_COMPONENT_RADII[-1][1])
    parser.add_argument("--cms-half-length-cm", type=float, default=CMS_HALF_LENGTH_CM)
    parser.add_argument("--min-surface-energy-mev", type=float, default=30000.0)
    parser.add_argument("--min-cosz", type=float, default=0.1, help="Minimum sampled vertical direction cosine.")
    parser.add_argument("--density-g-cm3", type=float, default=2.30)
    parser.add_argument("--material-z-over-a", type=float, default=0.5)
    parser.add_argument("--material-I-eV", type=float, default=136.4)
    parser.add_argument("--radiative-b-per-gcm2", type=float, default=4.0e-6)
    parser.add_argument("--range-grid-points", type=int, default=5000)
    parser.add_argument(
        "--range-smearing-gcm2",
        type=float,
        default=0.0,
        help="If positive, use a logistic survival turn-on with this width in g/cm^2.",
    )
    parser.add_argument(
        "--no-cos-projection",
        action="store_true",
        help="Do not weight thrown surface directions by cos(theta).",
    )
    parser.add_argument(
        "--ground-coefficients-cm",
        default=",".join(str(value) for value in DEFAULT_GROUND_COEFFICIENTS_CM),
        help="Six-coefficient cavern ground model shared with the display scripts.",
    )
    parser.add_argument("--max-voxels-per-track", type=int, default=30000)
    parser.add_argument("--weight-floor", type=float, default=0.0)
    parser.add_argument("--write-tree", action="store_true", help="Store accepted generated rays for debugging.")
    return parser.parse_args()


def split_float_csv(text):
    return [float(item.strip()) for item in text.split(",") if item.strip()]


class FluxSampler:
    def __init__(self, path, min_energy_mev=0.0, min_cosz=0.0, include_cos_projection=True):
        self.cells = []
        cumulative = []
        total = 0.0
        with open(path, encoding="utf-8") as handle:
            reader = csv.DictReader(line for line in handle if not line.startswith("#"))
            for row in reader:
                e_low = max(float(row["energy_low_mev"]), min_energy_mev)
                e_high = float(row["energy_high_mev"])
                if e_high <= e_low:
                    continue
                cos_low = max(float(row["cos_low"]), min_cosz)
                cos_high = float(row["cos_high"])
                if cos_high <= cos_low:
                    continue
                cos_mid = float(row["cos_mid"])
                if cos_mid < cos_low or cos_mid > cos_high:
                    cos_mid = 0.5 * (cos_low + cos_high)
                flux = max(float(row["flux"]), 0.0)
                weight = flux * (e_high - e_low) * (cos_high - cos_low)
                if include_cos_projection:
                    weight *= max(cos_mid, 0.0)
                if weight <= 0.0:
                    continue
                total += weight
                self.cells.append((e_low, e_high, cos_low, cos_high, flux, weight))
                cumulative.append(total)
        if not self.cells:
            raise RuntimeError(f"No usable flux cells found in {path}")
        self.cumulative = cumulative
        self.total_weight = total

    def sample(self, rng):
        value = rng.random() * self.total_weight
        index = bisect.bisect_left(self.cumulative, value)
        e_low, e_high, cos_low, cos_high, _flux, _weight = self.cells[min(index, len(self.cells) - 1)]
        energy = e_low + (e_high - e_low) * rng.random()
        cosz = cos_low + (cos_high - cos_low) * rng.random()
        return energy, min(max(cosz, 0.0), 1.0)


def bethe_bloch_ionization_mev_gcm2(kinetic_mev, z_over_a, mean_excitation_eV):
    kinetic_mev = max(kinetic_mev, 1.0e-6)
    total_energy = kinetic_mev + MUON_MASS_MEV
    gamma = total_energy / MUON_MASS_MEV
    beta2 = max(1.0 - 1.0 / (gamma * gamma), 1.0e-9)
    mass_ratio = ELECTRON_MASS_MEV / MUON_MASS_MEV
    tmax = (
        2.0
        * ELECTRON_MASS_MEV
        * beta2
        * gamma
        * gamma
        / (1.0 + 2.0 * gamma * mass_ratio + mass_ratio * mass_ratio)
    )
    mean_excitation_mev = mean_excitation_eV * 1.0e-6
    argument = max(
        2.0 * ELECTRON_MASS_MEV * beta2 * gamma * gamma * tmax / (mean_excitation_mev * mean_excitation_mev),
        1.0 + 1.0e-9,
    )
    stopping = BETHE_K_MEV_CM2_G * z_over_a / beta2 * (0.5 * math.log(argument) - beta2)
    return max(stopping, 0.1)


class RangeTable:
    def __init__(self, max_energy_mev, points, z_over_a, mean_excitation_eV, radiative_b_per_gcm2):
        points = max(points, 200)
        min_energy = 1.0
        max_energy = max(max_energy_mev, min_energy * 10.0)
        log_min = math.log(min_energy)
        log_max = math.log(max_energy)
        self.energies = [math.exp(log_min + (log_max - log_min) * i / (points - 1)) for i in range(points)]
        self.ranges = [0.0 for _ in self.energies]
        total_range = 0.0
        for index in range(1, points):
            e0 = self.energies[index - 1]
            e1 = self.energies[index]
            emid = math.sqrt(e0 * e1)
            dedx = bethe_bloch_ionization_mev_gcm2(emid, z_over_a, mean_excitation_eV)
            dedx += max(radiative_b_per_gcm2, 0.0) * emid
            total_range += (e1 - e0) / max(dedx, 1.0e-9)
            self.ranges[index] = total_range

    def range_gcm2(self, kinetic_mev):
        if kinetic_mev <= self.energies[0]:
            return 0.0
        if kinetic_mev >= self.energies[-1]:
            return self.ranges[-1]
        index = bisect.bisect_left(self.energies, kinetic_mev)
        e0 = self.energies[index - 1]
        e1 = self.energies[index]
        r0 = self.ranges[index - 1]
        r1 = self.ranges[index]
        fraction = (kinetic_mev - e0) / (e1 - e0)
        return r0 + fraction * (r1 - r0)


def survival_probability(range_gcm2, depth_gcm2, smearing_gcm2):
    if smearing_gcm2 <= 0.0:
        return 1.0 if range_gcm2 >= depth_gcm2 else 0.0
    scaled = (depth_gcm2 - range_gcm2) / smearing_gcm2
    if scaled > 50.0:
        return 0.0
    if scaled < -50.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(scaled))


def line_cms_cylinder_intersection(point, direction, radius_cm, half_length_cm):
    px, py, pz = point
    ux, uy, uz = direction
    intervals = []

    a = ux * ux + uz * uz
    b = 2.0 * (px * ux + pz * uz)
    c = px * px + pz * pz - radius_cm * radius_cm
    if abs(a) < 1.0e-12:
        if c > 0.0:
            return None
        radial_interval = (-float("inf"), float("inf"))
    else:
        discriminant = b * b - 4.0 * a * c
        if discriminant < 0.0:
            return None
        root = math.sqrt(max(discriminant, 0.0))
        radial_interval = ((-b - root) / (2.0 * a), (-b + root) / (2.0 * a))
    intervals.append(radial_interval)

    if abs(uy) < 1.0e-12:
        if abs(py) > half_length_cm:
            return None
        beam_interval = (-float("inf"), float("inf"))
    else:
        t1 = (-half_length_cm - py) / uy
        t2 = (half_length_cm - py) / uy
        beam_interval = (min(t1, t2), max(t1, t2))
    intervals.append(beam_interval)

    low = max(item[0] for item in intervals)
    high = min(item[1] for item in intervals)
    if low > high:
        return None
    return low, high


def sample_point_in_cms_cylinder(rng, radius_cm, half_length_cm):
    radius = radius_cm * math.sqrt(rng.random())
    phi = 2.0 * math.pi * rng.random()
    return (
        radius * math.cos(phi),
        rng.uniform(-half_length_cm, half_length_cm),
        radius * math.sin(phi),
    )


def make_log_hist(root, name, title, bins, low, high):
    edges = array("d")
    log_low = math.log(low)
    log_high = math.log(high)
    for index in range(bins + 1):
        edges.append(math.exp(log_low + (log_high - log_low) * index / bins))
    hist = root.TH1F(name, title, bins, edges)
    hist.SetDirectory(0)
    return hist


def make_debug_tree(root):
    tree = root.TTree("generated_tracks", "accepted generated solid-ground rays")
    values = {}
    branches = [
        ("energy_mev", "f"),
        ("cosz", "f"),
        ("phi", "f"),
        ("x0", "f"),
        ("y0", "f"),
        ("ux", "f"),
        ("uy", "f"),
        ("uz", "f"),
        ("depth_gcm2", "f"),
        ("range_gcm2", "f"),
        ("survival_probability", "f"),
        ("chord_cm", "f"),
        ("ray_weight", "f"),
        ("n_voxels", "i"),
    ]
    for name, typecode in branches:
        values[name] = array(typecode, [0])
        leaf = "F" if typecode == "f" else "I"
        tree.Branch(name, values[name], f"{name}/{leaf}")
    return tree, values


def main():
    args = parse_args()
    if args.nchunks <= 0:
        raise RuntimeError("--nchunks must be positive")
    if not (0 <= args.chunk_index < args.nchunks):
        raise RuntimeError("--chunk-index must be in [0, nchunks)")

    import ROOT

    rng = random.Random(args.seed + 1000003 * args.chunk_index)
    include_cos_projection = not args.no_cos_projection
    flux_sampler = FluxSampler(args.flux_table, args.min_surface_energy_mev, args.min_cosz, include_cos_projection)
    max_energy = max(cell[1] for cell in flux_sampler.cells)
    range_table = RangeTable(
        max_energy,
        args.range_grid_points,
        args.material_z_over_a,
        args.material_I_eV,
        args.radiative_b_per_gcm2,
    )

    template_file = ROOT.TFile.Open(args.template_root)
    if not template_file or template_file.IsZombie():
        raise RuntimeError(f"Could not open {args.template_root}")
    template_hist = template_file.Get(args.template_hist)
    if not template_hist:
        raise RuntimeError(f"Could not find {args.template_hist} in {args.template_root}")
    template_hist.SetDirectory(0)
    template_file.Close()

    bounds, bins = hist_bounds_and_bins(template_hist)
    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    if not output or output.IsZombie():
        raise RuntimeError(f"Could not create {args.output}")

    expected = ROOT.TH3F(
        "generator_expected_voxels",
        "independent EXPACS Bethe-range generated expected voxel crossings;{};{};{}".format(
            template_hist.GetXaxis().GetTitle(),
            template_hist.GetYaxis().GetTitle(),
            template_hist.GetZaxis().GetTitle(),
        ),
        template_hist.GetNbinsX(),
        template_hist.GetXaxis().GetXmin(),
        template_hist.GetXaxis().GetXmax(),
        template_hist.GetNbinsY(),
        template_hist.GetYaxis().GetXmin(),
        template_hist.GetYaxis().GetXmax(),
        template_hist.GetNbinsZ(),
        template_hist.GetZaxis().GetXmin(),
        template_hist.GetZaxis().GetXmax(),
    )
    expected.SetDirectory(0)
    expected.GetXaxis().SetTitle(template_hist.GetXaxis().GetTitle())
    expected.GetYaxis().SetTitle(template_hist.GetYaxis().GetTitle())
    expected.GetZaxis().SetTitle(template_hist.GetZaxis().GetTitle())

    h_energy = make_log_hist(ROOT, "thrown_energy_mev", "thrown surface kinetic energy;E_{#mu} [MeV];throws", 80, args.min_surface_energy_mev, max_energy)
    h_surv_energy = make_log_hist(ROOT, "surviving_energy_mev", "surviving weighted surface kinetic energy;E_{#mu} [MeV];weighted throws", 80, args.min_surface_energy_mev, max_energy)
    h_cos = ROOT.TH1F("thrown_cosz", "thrown #cos#theta_{z};#cos#theta_{z};throws", 100, 0.0, 1.0)
    h_surv_cos = ROOT.TH1F("surviving_cosz", "surviving weighted #cos#theta_{z};#cos#theta_{z};weighted throws", 100, 0.0, 1.0)
    h_depth = ROOT.TH1F("depth_gcm2", "solid-ground depth;X [g/cm^{2}];weighted throws", 120, 0.0, 60000.0)
    h_survprob = ROOT.TH1F("survival_probability", "CSDA survival probability;P_{surv};throws", 100, 0.0, 1.0)
    h_chord = ROOT.TH1F("cms_chord_cm", "CMS chord length;L_{CMS} [cm];weighted throws", 120, 0.0, 3000.0)
    for hist in (h_cos, h_surv_cos, h_depth, h_survprob, h_chord):
        hist.SetDirectory(0)

    debug_tree = None
    debug_values = None
    if args.write_tree:
        debug_tree, debug_values = make_debug_tree(ROOT)

    ground_args = GroundArgs(args.ground_coefficients_cm)
    total_throws = 0
    cms_crossings = 0
    ground_crossings = 0
    weighted_survivors = 0.0
    filled_tracks = 0
    total_voxel_weight = 0.0
    skipped_no_box = 0

    for throw_index in range(args.chunk_index, args.nthrows, args.nchunks):
        total_throws += 1
        energy, cosz = flux_sampler.sample(rng)
        phi = 2.0 * math.pi * rng.random()
        sin_theta = math.sqrt(max(1.0 - cosz * cosz, 0.0))
        direction = (sin_theta * math.cos(phi), sin_theta * math.sin(phi), cosz)
        if args.sampling_mode == "reverse-cms-volume":
            point = sample_point_in_cms_cylinder(rng, args.cms_radius_cm, args.cms_half_length_cm)
        else:
            point = (
                rng.uniform(-args.throw_x_half_width_cm, args.throw_x_half_width_cm),
                rng.uniform(-args.throw_y_half_width_cm, args.throw_y_half_width_cm),
                0.0,
            )
        h_energy.Fill(energy)
        h_cos.Fill(cosz)

        cms_segment = line_cms_cylinder_intersection(
            point,
            direction,
            args.cms_radius_cm,
            args.cms_half_length_cm,
        )
        if cms_segment is None:
            continue
        cms_crossings += 1
        chord_length = max(cms_segment[1] - cms_segment[0], 1.0e-9)
        sampling_weight = 1.0 / chord_length if args.sampling_mode == "reverse-cms-volume" else 1.0
        cms_exit_t = max(cms_segment)
        cms_exit_point = tuple(point[index] + cms_exit_t * direction[index] for index in range(3))
        t_ground = ground_intersection_t(cms_exit_point, direction, ground_args)
        if t_ground is None:
            continue
        ground_crossings += 1
        depth = args.density_g_cm3 * t_ground
        muon_range = range_table.range_gcm2(energy)
        psurv = survival_probability(muon_range, depth, args.range_smearing_gcm2)
        h_survprob.Fill(psurv)
        ray_weight = sampling_weight * psurv
        if ray_weight <= args.weight_floor:
            continue

        segment = ray_box_intersection(point, direction, bounds)
        if segment is None:
            skipped_no_box += 1
            continue
        n_voxels, _truncated = fill_voxels(
            expected,
            point,
            direction,
            segment[0],
            segment[1],
            bounds,
            bins,
            ray_weight,
            args.max_voxels_per_track,
        )
        if n_voxels <= 0:
            continue

        h_surv_energy.Fill(energy, ray_weight)
        h_surv_cos.Fill(cosz, ray_weight)
        h_depth.Fill(depth, ray_weight)
        h_chord.Fill(chord_length, ray_weight)
        weighted_survivors += ray_weight
        filled_tracks += 1
        total_voxel_weight += ray_weight * n_voxels

        if debug_tree is not None:
            debug_values["energy_mev"][0] = energy
            debug_values["cosz"][0] = cosz
            debug_values["phi"][0] = phi
            debug_values["x0"][0] = point[0]
            debug_values["y0"][0] = point[1]
            debug_values["ux"][0] = direction[0]
            debug_values["uy"][0] = direction[1]
            debug_values["uz"][0] = direction[2]
            debug_values["depth_gcm2"][0] = depth
            debug_values["range_gcm2"][0] = muon_range
            debug_values["survival_probability"][0] = psurv
            debug_values["chord_cm"][0] = chord_length
            debug_values["ray_weight"][0] = ray_weight
            debug_values["n_voxels"][0] = n_voxels
            debug_tree.Fill()

        if args.progress_every > 0 and total_throws % args.progress_every == 0:
            print(
                f"throw {total_throws}: cms={cms_crossings} ground={ground_crossings} "
                f"filled={filled_tracks} weighted_survivors={weighted_survivors:g}"
            )

    expected.SetEntries(total_voxel_weight)
    output.cd()
    expected.Write()
    for hist in (h_energy, h_surv_energy, h_cos, h_surv_cos, h_depth, h_survprob, h_chord):
        hist.Write()
    if debug_tree is not None:
        debug_tree.Write()

    config = {
        "template_root": args.template_root,
        "template_hist": args.template_hist,
        "flux_table": args.flux_table,
        "flux_model": {
            "name": "EXPACS/PARMA",
            "url": "https://phits.jaea.go.jp/expacs/",
            "citations": [
                "T. Sato, PLOS ONE 10(12): e0144679 (2015), doi:10.1371/journal.pone.0144679",
                "T. Sato, PLOS ONE 11(8): e0160390 (2016), doi:10.1371/journal.pone.0160390",
            ],
        },
        "nthrows_requested": args.nthrows,
        "chunk_index": args.chunk_index,
        "nchunks": args.nchunks,
        "processed_throws": total_throws,
        "sampling_mode": args.sampling_mode,
        "cms_crossings": cms_crossings,
        "ground_crossings": ground_crossings,
        "filled_tracks": filled_tracks,
        "weighted_survivors": weighted_survivors,
        "voxel_weight_sum": total_voxel_weight,
        "skipped_no_box": skipped_no_box,
        "throw_x_half_width_cm": args.throw_x_half_width_cm,
        "throw_y_half_width_cm": args.throw_y_half_width_cm,
        "cms_radius_cm": args.cms_radius_cm,
        "cms_half_length_cm": args.cms_half_length_cm,
        "min_surface_energy_mev": args.min_surface_energy_mev,
        "min_cosz": args.min_cosz,
        "include_cos_projection": include_cos_projection,
        "density_g_cm3": args.density_g_cm3,
        "material_z_over_a": args.material_z_over_a,
        "material_I_eV": args.material_I_eV,
        "radiative_b_per_gcm2": args.radiative_b_per_gcm2,
        "range_smearing_gcm2": args.range_smearing_gcm2,
        "ground_coefficients_cm": args.ground_coefficients_cm,
        "note": (
            "Independent first-pass generator null. EXPACS/PARMA supplies surface energy "
            "and zenith-angle sampling; azimuth and detector-plane impact parameters are "
            "random in forward-plane mode, while reverse-cms-volume samples points in the "
            "CMS cylinder and applies a 1/chord-length projected-area correction; CMS "
            "crossing is geometric with 100% efficiency; solid-ground survival uses a "
            "Bethe-Bloch plus bE continuous-slowing-down range table."
        ),
    }
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "generator_config_json")
    output.Close()
    print(
        "Wrote {output}: throws={throws} cms={cms} ground={ground} filled={filled} "
        "weighted_survivors={survivors:g} voxel_weight={voxels:g}".format(
            output=args.output,
            throws=total_throws,
            cms=cms_crossings,
            ground=ground_crossings,
            filled=filled_tracks,
            survivors=weighted_survivors,
            voxels=total_voxel_weight,
        )
    )


if __name__ == "__main__":
    main()
