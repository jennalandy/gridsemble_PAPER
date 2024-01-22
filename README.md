# Gridsemble: Selective Ensembling for False Discovery Rates

This repository contains the code to replicate all results reported in the paper *Gridsemble: Selective Ensembling for False Discovery Rates*. We use [quarto](https://quarto.org/) documents which can be edited and run within RStudio, Jupyter Lab, or Visual Studio Code.

## Simulation Studies

## Experimental Application

Our experimental application relies on the Platinum Spike dataset$^1$.

- `PAPER_platinum_data.qmd`: load and pre-process Platinum Spike data. This needs to be run before either of the analysis documents.
- `PAPER_platinum_run_subsets.qmd`: analyses on subsets of Platinum Spike data with $\pi_0 \in [0.6, 0.95]$.
- `PAPER_platinum_run_all_data.qmd`: analyses on the full Platinum Spike dataset, used to create Supplementary Figure 3 and Supplementary Table 2.
- `PAPER_metrics_helpers.R`: functions to calculate metrics and helper functions.

## References

[1] Q. Zhu, J.C. Miecznikowski, and M.S. Halfon. Preferred analysis methods for affymetrix genechips. II. an expanded, balanced, wholly-defined spike-in dataset. *BMC Bioinformatics*, 11:285, 2010. doi:https://doi.org/10.1186/1471-2105-11-285.
