#!/usr/bin/env python3
"""Mirror old Z-up cavern ROOT histograms into the right-handed cavern frame.

Older local products used x=CMS x, y=+CMS z, z=CMS y.  The preferred cavern
frame is now x=CMS x, y=-CMS z, z=CMS y, which keeps z vertical/up and makes
top-down x-y displays map-handed.  This script mirrors the TH3 y axis so old
histogram-only products can be compared with newly produced outputs.
"""

import argparse
import json
import os

import ROOT


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--input", required=True, help="Input ROOT file in the old cavern frame.")
    parser.add_argument("-o", "--output", required=True, help="Output ROOT file in the right-handed cavern frame.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite the output file if it already exists.",
    )
    return parser.parse_args()


def update_axis_title(axis):
    title = axis.GetTitle()
    title = title.replace("y cavern = CMS z beam", "y cavern = -CMS z beam")
    title = title.replace("output x=CMS x, y=CMS z beam axis", "output x=CMS x, y=-CMS z beam axis")
    axis.SetTitle(title)


def copy_axis_style(source, target):
    target.SetTitle(source.GetTitle())
    target.SetLabelSize(source.GetLabelSize())
    target.SetTitleSize(source.GetTitleSize())
    target.SetTitleOffset(source.GetTitleOffset())
    target.SetLabelOffset(source.GetLabelOffset())
    target.SetNdivisions(source.GetNdivisions())


def clone_empty_th3(source):
    xaxis = source.GetXaxis()
    yaxis = source.GetYaxis()
    zaxis = source.GetZaxis()
    hist = source.__class__(
        source.GetName(),
        source.GetTitle(),
        xaxis.GetNbins(),
        xaxis.GetXmin(),
        xaxis.GetXmax(),
        yaxis.GetNbins(),
        yaxis.GetXmin(),
        yaxis.GetXmax(),
        zaxis.GetNbins(),
        zaxis.GetXmin(),
        zaxis.GetXmax(),
    )
    copy_axis_style(xaxis, hist.GetXaxis())
    copy_axis_style(yaxis, hist.GetYaxis())
    copy_axis_style(zaxis, hist.GetZaxis())
    update_axis_title(hist.GetYaxis())
    hist.GetZaxis().SetTitle(source.GetZaxis().GetTitle())
    hist.SetDirectory(0)
    if source.GetSumw2N() > 0:
        hist.Sumw2()
    return hist


def mirror_th3_y(source):
    hist = clone_empty_th3(source)
    nx = source.GetNbinsX()
    ny = source.GetNbinsY()
    nz = source.GetNbinsZ()
    for ix in range(0, nx + 2):
        for iy in range(0, ny + 2):
            if iy == 0 or iy == ny + 1:
                new_iy = iy
            else:
                new_iy = ny + 1 - iy
            for iz in range(0, nz + 2):
                old_bin = source.GetBin(ix, iy, iz)
                new_bin = hist.GetBin(ix, new_iy, iz)
                hist.SetBinContent(new_bin, source.GetBinContent(old_bin))
                if source.GetSumw2N() > 0:
                    hist.SetBinError(new_bin, source.GetBinError(old_bin))
    return hist


def write_reoriented_object(source_dir, output_dir, key):
    key_name = key.GetName()
    obj = key.ReadObj()
    output_dir.cd()
    if obj.InheritsFrom("TDirectory"):
        subdir = output_dir.mkdir(key_name)
        copy_directory(obj, subdir)
        return
    if obj.InheritsFrom("TH3"):
        mirrored = mirror_th3_y(obj)
        mirrored.Write(key_name)
        return
    obj.Write(key_name)


def copy_directory(source_dir, output_dir):
    for key in source_dir.GetListOfKeys():
        write_reoriented_object(source_dir, output_dir, key)


def main():
    args = parse_args()
    if os.path.exists(args.output) and not args.overwrite:
        raise SystemExit(f"Output already exists: {args.output}")
    ROOT.gROOT.SetBatch(True)
    infile = ROOT.TFile.Open(args.input)
    if not infile or infile.IsZombie():
        raise SystemExit(f"Could not open input ROOT file: {args.input}")
    mode = "RECREATE" if args.overwrite else "CREATE"
    outfile = ROOT.TFile.Open(args.output, mode)
    if not outfile or outfile.IsZombie():
        raise SystemExit(f"Could not create output ROOT file: {args.output}")
    copy_directory(infile, outfile)
    note = {
        "source_file": os.path.abspath(args.input),
        "old_convention": "x=CMS x, y=+CMS z, z=CMS y",
        "new_convention": "x=CMS x, y=-CMS z, z=CMS y",
        "operation": "TH3 y-axis bin mirror: new_y_bin = nbins_y + 1 - old_y_bin",
    }
    outfile.WriteObject(ROOT.TObjString(json.dumps(note, indent=2, sort_keys=True)), "cavern_reorientation_note_json")
    outfile.Close()
    infile.Close()
    print(args.output)


if __name__ == "__main__":
    main()
