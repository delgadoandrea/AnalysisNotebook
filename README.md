# AnalysisNotebook

This repository assumes you already have Heppy Ntuples for signal / MC.

1. The first step is tro create weights for data/MC.
- For BTagging, run the Notebook under DeepCSVReweight/notebooks/MCEff.ipynb
- For the PU weight distribution for data/MC, this plot was obtained running the Heppy NTuplizer twice. The distribution was created automatically.

2. For the trigger efficiency estimation and weights:
- The efficiency estimation scripts can be found at Scripts/TriggerEstimate.C
- The weights were derived by fitting a distribution to the efficiency plots and adding an up/down variation. The ROOT file can be found under ROOT Files/TriggerFits4b.root

3. Once you have these weight files, you can go to Scripts/Uncertainties.C to get the mjj shapes with variations up/down of the uncertainties considered.
