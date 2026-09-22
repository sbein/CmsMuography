#!/usr/bin/env python3
"""Build a smooth calibrated null and cavity-oriented residual volumes."""

import argparse
import json
import math
import os


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Use observed and uniform-earth expected TH3s to fit a deliberately smooth "
            "multiplicative calibration C(z)*Cxy(x,y), then recompute cavity-oriented "
            "significance and transmission residual volumes."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Li-Ma ROOT with observed/expected TH3s.")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--observed-hist", default="observed_voxels")
    parser.add_argument("--expected-hist", default="expected_voxels_scaled")
    parser.add_argument("--method", choices=("lima", "poisson-known"), default="lima")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-expected", type=float, default=1.0e-6)
    parser.add_argument(
        "--calibration-model",
        choices=("zxy-poly", "product1d"),
        default="zxy-poly",
        help=(
            "zxy-poly fits C(z)*Cxy(x,y). product1d uses iterative smooth "
            "Cx(x)*Cy(y)*Cz(z) raking, which is less able to absorb localized "
            "2D/3D cavity patterns."
        ),
    )
    parser.add_argument("--x-smooth-cm", type=float, default=2500.0)
    parser.add_argument("--y-smooth-cm", type=float, default=2500.0)
    parser.add_argument("--z-smooth-cm", type=float, default=1500.0)
    parser.add_argument("--ipf-iterations", type=int, default=8)
    parser.add_argument("--xy-poly-order", choices=("0", "1", "2"), default="2")
    parser.add_argument("--huber-k", type=float, default=1.5)
    parser.add_argument("--robust-iterations", type=int, default=8)
    parser.add_argument("--correction-min", type=float, default=0.25)
    parser.add_argument("--correction-max", type=float, default=4.0)
    parser.add_argument("--pseudo-count", type=float, default=0.5)
    parser.add_argument(
        "--exclude-cms-core",
        action="store_true",
        default=True,
        help="Exclude the approximate CMS detector envelope from calibration fitting.",
    )
    parser.add_argument("--include-cms-core", dest="exclude_cms_core", action="store_false")
    parser.add_argument("--core-x-cm", type=float, default=900.0)
    parser.add_argument("--core-y-cm", type=float, default=1300.0)
    parser.add_argument("--core-z-cm", type=float, default=1200.0)
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
    hist.SetDirectory(0)
    return hist


def poisson_known_z(observed, expected):
    if expected <= 0.0:
        return 0.0
    if observed <= 0.0:
        q = 2.0 * expected
    else:
        q = 2.0 * (observed * math.log(observed / expected) - (observed - expected))
    q = max(q, 0.0)
    return (1.0 if observed >= expected else -1.0) * math.sqrt(q)


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
    return (1.0 if observed >= expected else -1.0) * math.sqrt(q)


def hist_to_array(hist):
    import numpy as np

    array = np.zeros((hist.GetNbinsX(), hist.GetNbinsY(), hist.GetNbinsZ()), dtype="float64")
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            for iz in range(1, hist.GetNbinsZ() + 1):
                array[ix - 1, iy - 1, iz - 1] = hist.GetBinContent(ix, iy, iz)
    return array


def fill_hist_from_array(hist, values):
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            for iz in range(1, hist.GetNbinsZ() + 1):
                value = float(values[ix - 1, iy - 1, iz - 1])
                if value != 0.0:
                    hist.SetBinContent(ix, iy, iz, value)


def axis_centers(axis_obj):
    import numpy as np

    return np.array(
        [axis_obj.GetBinCenter(index) for index in range(1, axis_obj.GetNbins() + 1)],
        dtype="float64",
    )


def gaussian_smooth_1d(values, sigma_bins):
    import numpy as np

    if sigma_bins <= 0.0:
        return values.copy()
    half_width = max(1, int(math.ceil(4.0 * sigma_bins)))
    offsets = np.arange(-half_width, half_width + 1, dtype="float64")
    kernel = np.exp(-0.5 * (offsets / sigma_bins) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(values, (half_width, half_width), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def fit_z_correction(observed, expected, fit_mask, z_smooth_bins):
    import numpy as np

    nz = observed.shape[2]
    ratio = np.ones(nz, dtype="float64")
    for iz in range(nz):
        mask = fit_mask[:, :, iz]
        exp_sum = expected[:, :, iz][mask].sum()
        obs_sum = observed[:, :, iz][mask].sum()
        ratio[iz] = obs_sum / exp_sum if exp_sum > 0.0 else 1.0
    ratio = np.clip(ratio, 1.0e-6, 1.0e6)
    return np.exp(gaussian_smooth_1d(np.log(ratio), z_smooth_bins))


def polynomial_design(x, y, order):
    import numpy as np

    columns = [np.ones_like(x)]
    if order >= 1:
        columns.extend([x, y])
    if order >= 2:
        columns.extend([x * x, y * y, x * y])
    return np.vstack(columns).T


def robust_weighted_fit(design, target, weights, huber_k, iterations):
    import numpy as np

    coeff = np.zeros(design.shape[1], dtype="float64")
    fit_weights = weights.copy()
    for _iteration in range(max(1, iterations)):
        sqrt_w = np.sqrt(np.clip(fit_weights, 0.0, None))
        lhs = design * sqrt_w[:, None]
        rhs = target * sqrt_w
        coeff = np.linalg.lstsq(lhs, rhs, rcond=None)[0]
        residual = target - design.dot(coeff)
        mad = np.median(np.abs(residual - np.median(residual)))
        sigma = max(1.4826 * mad, 1.0e-6)
        robust = np.minimum(1.0, huber_k * sigma / np.maximum(np.abs(residual), 1.0e-12))
        fit_weights = weights * robust
    residual = target - design.dot(coeff)
    return coeff, residual


def fit_xy_correction(observed, expected_z, fit_mask, x_centers, y_centers, order, huber_k, iterations):
    import numpy as np

    obs_xy = (observed * fit_mask).sum(axis=2)
    exp_xy = (expected_z * fit_mask).sum(axis=2)
    valid = exp_xy > 0.0
    ratio = np.ones_like(exp_xy)
    ratio[valid] = obs_xy[valid] / exp_xy[valid]
    ratio = np.clip(ratio, 1.0e-6, 1.0e6)

    xx, yy = np.meshgrid(x_centers / max(abs(x_centers).max(), 1.0), y_centers / max(abs(y_centers).max(), 1.0), indexing="ij")
    design = polynomial_design(xx[valid], yy[valid], order)
    target = np.log(ratio[valid])
    weights = np.sqrt(np.clip(exp_xy[valid], 1.0, None))
    coeff, residual = robust_weighted_fit(design, target, weights, huber_k, iterations)

    all_design = polynomial_design(xx.ravel(), yy.ravel(), order)
    log_corr = all_design.dot(coeff).reshape(exp_xy.shape)
    return np.exp(log_corr), coeff, residual


def fit_axis_factor(observed, expected_corrected, fit_mask, axis_index, sigma_bins):
    import numpy as np

    moved_obs = np.moveaxis(observed, axis_index, 0)
    moved_exp = np.moveaxis(expected_corrected, axis_index, 0)
    moved_mask = np.moveaxis(fit_mask, axis_index, 0)
    ratio = np.ones(moved_obs.shape[0], dtype="float64")
    for ibin in range(moved_obs.shape[0]):
        mask = moved_mask[ibin]
        exp_sum = moved_exp[ibin][mask].sum()
        obs_sum = moved_obs[ibin][mask].sum()
        ratio[ibin] = obs_sum / exp_sum if exp_sum > 0.0 else 1.0
    ratio = np.clip(ratio, 1.0e-6, 1.0e6)
    return np.exp(gaussian_smooth_1d(np.log(ratio), sigma_bins))


def broadcast_axis_factor(factor, axis_index):
    shape = [1, 1, 1]
    shape[axis_index] = factor.shape[0]
    return factor.reshape(shape)


def fit_product1d_correction(observed, expected, fit_mask, smooth_bins, iterations, correction_min, correction_max):
    import numpy as np

    correction = np.ones_like(expected, dtype="float64")
    axis_factors = [None, None, None]
    for _iteration in range(max(1, iterations)):
        for axis_index in (2, 1, 0):
            factor = fit_axis_factor(
                observed,
                expected * correction,
                fit_mask,
                axis_index,
                smooth_bins[axis_index],
            )
            axis_factors[axis_index] = factor
            correction *= broadcast_axis_factor(factor, axis_index)
            correction = np.clip(correction, correction_min, correction_max)
    return correction, axis_factors


def make_fit_mask(observed, expected, x_centers, y_centers, z_centers, args):
    import numpy as np

    mask = expected > args.min_expected
    if args.exclude_cms_core:
        xx, yy, zz = np.meshgrid(x_centers, y_centers, z_centers, indexing="ij")
        core = (
            (np.abs(xx) < args.core_x_cm)
            & (np.abs(yy) < args.core_y_cm)
            & (zz < args.core_z_cm)
        )
        mask = mask & (~core)
    return mask


def main():
    args = parse_args()

    import numpy as np
    import ROOT

    source = ROOT.TFile.Open(args.input)
    if not source or source.IsZombie():
        raise RuntimeError(f"Could not open {args.input}")
    observed_hist = require_th3(source, args.observed_hist)
    expected_hist = require_th3(source, args.expected_hist)
    source.Close()

    observed = hist_to_array(observed_hist)
    expected = hist_to_array(expected_hist)
    x_centers = axis_centers(observed_hist.GetXaxis())
    y_centers = axis_centers(observed_hist.GetYaxis())
    z_centers = axis_centers(observed_hist.GetZaxis())
    x_smooth_bins = args.x_smooth_cm / observed_hist.GetXaxis().GetBinWidth(1)
    y_smooth_bins = args.y_smooth_cm / observed_hist.GetYaxis().GetBinWidth(1)
    z_smooth_bins = args.z_smooth_cm / observed_hist.GetZaxis().GetBinWidth(1)

    fit_mask = make_fit_mask(observed, expected, x_centers, y_centers, z_centers, args)
    axis_factor_ranges = {}
    if args.calibration_model == "product1d":
        correction, axis_factors = fit_product1d_correction(
            observed,
            expected,
            fit_mask,
            (x_smooth_bins, y_smooth_bins, z_smooth_bins),
            args.ipf_iterations,
            args.correction_min,
            args.correction_max,
        )
        xy_coeff = np.array([], dtype="float64")
        xy_residual = np.array([], dtype="float64")
        for name, factor in zip(("x", "y", "z"), axis_factors):
            if factor is not None:
                axis_factor_ranges[name] = [float(factor.min()), float(factor.max())]
    else:
        z_corr = fit_z_correction(observed, expected, fit_mask, z_smooth_bins)
        expected_z = expected * z_corr[None, None, :]
        xy_corr, xy_coeff, xy_residual = fit_xy_correction(
            observed,
            expected_z,
            fit_mask,
            x_centers,
            y_centers,
            int(args.xy_poly_order),
            args.huber_k,
            args.robust_iterations,
        )
        correction = z_corr[None, None, :] * xy_corr[:, :, None]
        correction = np.clip(correction, args.correction_min, args.correction_max)
    calibrated_expected = expected * correction
    norm = observed.sum() / calibrated_expected.sum()
    calibrated_expected *= norm
    correction *= norm

    signed = np.zeros_like(observed)
    excess = np.zeros_like(observed)
    ratio = np.zeros_like(observed)
    log_transmission = np.zeros_like(observed)
    valid = calibrated_expected > args.min_expected
    z_min = 999.0
    z_max = -999.0
    for ix in range(observed.shape[0]):
        for iy in range(observed.shape[1]):
            for iz in range(observed.shape[2]):
                if not valid[ix, iy, iz]:
                    continue
                obs = float(observed[ix, iy, iz])
                exp = float(calibrated_expected[ix, iy, iz])
                if args.method == "lima":
                    z_value = lima_z(obs, exp, args.alpha)
                else:
                    z_value = poisson_known_z(obs, exp)
                signed[ix, iy, iz] = z_value
                excess[ix, iy, iz] = max(z_value, 0.0)
                ratio[ix, iy, iz] = obs / exp if exp > 0.0 else 0.0
                log_transmission[ix, iy, iz] = math.log(
                    (obs + args.pseudo_count) / (exp + args.pseudo_count)
                )
                z_min = min(z_min, z_value)
                z_max = max(z_max, z_value)

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output = ROOT.TFile.Open(args.output, "RECREATE")
    output.cd()
    observed_hist.Write("observed_voxels")
    expected_hist.Write("expected_voxels_global_scaled")

    products = {
        "calibrated_expected_voxels": (
            calibrated_expected,
            "smooth-calibrated expected voxel crossings",
        ),
        "smooth_calibration_factor": (
            correction,
            "smooth multiplicative calibration C(z)*Cxy(x,y)",
        ),
        "signed_significance_calibrated": (
            signed,
            "signed observed-minus-smooth-null significance",
        ),
        "cavity_excess_significance_calibrated": (
            excess,
            "cavity-candidate excess significance after smooth null calibration",
        ),
        "observed_over_expected_calibrated": (
            ratio,
            "observed / smooth-calibrated expected",
        ),
        "log_transmission_excess": (
            log_transmission,
            "log((observed+pseudo)/(smooth-calibrated expected+pseudo)); positive is cavity-like",
        ),
    }
    for name, (values, title) in products.items():
        hist = clone_empty(
            ROOT,
            observed_hist,
            name,
            f"{title};{observed_hist.GetXaxis().GetTitle()};{observed_hist.GetYaxis().GetTitle()};{observed_hist.GetZaxis().GetTitle()}",
        )
        fill_hist_from_array(hist, values)
        hist.SetEntries(float(np.count_nonzero(values)))
        hist.Write()

    config = {
        "input": args.input,
        "observed_hist": args.observed_hist,
        "expected_hist": args.expected_hist,
        "method": args.method,
        "alpha": args.alpha,
        "min_expected": args.min_expected,
        "calibration_model": args.calibration_model,
        "x_smooth_cm": args.x_smooth_cm,
        "y_smooth_cm": args.y_smooth_cm,
        "z_smooth_cm": args.z_smooth_cm,
        "ipf_iterations": args.ipf_iterations,
        "xy_poly_order": int(args.xy_poly_order),
        "xy_coefficients": [float(value) for value in xy_coeff],
        "xy_residual_median": float(np.median(xy_residual)) if xy_residual.size else 0.0,
        "xy_residual_mad": float(np.median(np.abs(xy_residual - np.median(xy_residual)))) if xy_residual.size else 0.0,
        "axis_factor_ranges": axis_factor_ranges,
        "global_post_calibration_norm": float(norm),
        "correction_range": [float(correction.min()), float(correction.max())],
        "observed_integral": float(observed.sum()),
        "expected_integral_before": float(expected.sum()),
        "expected_integral_after": float(calibrated_expected.sum()),
        "signed_z_range": [float(z_min), float(z_max)],
        "fit_voxels": int(np.count_nonzero(fit_mask)),
        "excluded_cms_core": bool(args.exclude_cms_core),
        "note": (
            "Prototype smooth sideband calibration. The correction is intentionally low "
            "frequency, C(z)*Cxy(x,y), to absorb bulk survival/normalization biases "
            "without fitting away localized hall/shaft-scale residuals."
        ),
    }
    output.WriteObject(ROOT.TObjString(json.dumps(config, indent=2, sort_keys=True)), "smooth_calibration_config_json")
    output.Close()
    print(
        "Wrote {output}: correction=[{cmin:g},{cmax:g}] z_range=[{zmin:g},{zmax:g}]".format(
            output=args.output,
            cmin=config["correction_range"][0],
            cmax=config["correction_range"][1],
            zmin=config["signed_z_range"][0],
            zmax=config["signed_z_range"][1],
        )
    )


if __name__ == "__main__":
    main()
