CMS Cavern Muography presentation bundle

Google Slides deck:
https://docs.google.com/presentation/d/1bzIsWHNY5xSYaxeBXLMEuMU-wxgaOkStQTgCUWG4_u8/edit

Each slide_NN directory contains:
- content.txt: the slide title, bullets, captions, and links
- slide.png: a complete 1600 x 899 render of that slide
- visual_*.png: the original scientific raster used on the slide, when applicable

The standalone directory contains a reusable simulation-chain table in SVG
and PNG formats. It separates the implemented fast null from the future
Geant4 + CMSSW DIGI/RECO extension.

It also contains a presentation-format CMS geometry data sheet summarizing
the detector wireframe, cavern coordinates, voxel grid, hall, shaft, ground
surface, and the model's interpretive limitations.

Interactive ROOT file on lxplus:
/afs/cern.ch/work/s/sbein/Muography/CmsCavern/CosmicMuography/test_output/presentation_highstat_2M_balanced/cms_muography_highstat_canvases.root

Open it with:
root -l test_output/presentation_highstat_2M_balanced/cms_muography_highstat_canvases.root

Then open the ROOT browser:
new TBrowser;

The file contains the copied input TH3 histograms and 21 TCanvas objects, including
the geometry baseline, observed/raw-null/corrected-null comparison, the 1D/2D
data-driven correction diagnostic, threshold grids, and XY/XZ slice stacks. The 3D
threshold canvases use shared color ranges within each quantity, the 3D comparison
uses one shared count scale, and the before/after correction maps use one shared
residual scale.
