# Gridsemble: Selective Ensembling for False Discovery Rates

This repository contains the code to replicate all results reported in the paper *Gridsemble: Selective Ensembling for False Discovery Rates*. See details in the Simulation Study and Experimental Application sections of our paper.

### Simulation Studies

Our simulation studies are in R scripts. When each script is run, it will log progress and results in a new sub-directory. Scripts assume you are in the `simulation_studies` directory. We use [`sink`](https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/sink) to log progress; if you terminate a test early you will need to run `sink()` for output to show up in the console again. Note that these scripts take many hours to run.

**Simulation studies presented in Figure 1**

These tests compare `gridsemble`, partial implementations, and benchmarks on Symmetric, Asymmetric, and Curated Ovarian Data-Based simulation studies.

- `symmetric_test.R`
- `asymmetric_test.R`
- `cod_based_test.R`

**Simulation studies presented in Figure 3**

These tests compare `gridsemble` and `ensemble` with varying number of synthetic datasets and model size in the Symmetric and Asymmetric simulation studies.

- `inc_n_symmetric.R`
- `inc_n_asymmetric.R`

**Simulation studies presented in Supplementary Figure 2**

These tests compare `gridsemble`, partial implementations, and benchmarks on Symmetric, Asymmetric, and COD-Based simulation studies using random search in place of grid search. 

These scripts can only be run after their grid search counterparts.

- `symmetric_test_random.R`
- `asymmetric_test_random.R`
- `cod_based_test_random.R`
- `get_random_grid.R`: functions to construct grids for a random search.

**Scripts sourced by the above**

- `test.R`: defines wrapper functions to run simulation study given a data generating function.
- `evaluate.R`: functions to compute metrics given fdr estimates and ground truth.
- `simulate.R`: functions to simulate each type of data.
- `utils.R`: other utility functions.

### Experimental Application

Our experimental application relies on the Platinum Spike dataset<sup>1</sup>. We use [quarto](https://quarto.org/) documents which can be edited and run with RStudio, Jupyter Lab, or Visual Studio Code.

**Notebooks**

- [`PAPER_platinum_data`](https://github.com/jennalandy/gridsemble_PAPER/blob/main/experimental/PAPER_platinum_data.pdf): download and pre-process Platinum Spike data. This needs to be run before either of the analysis documents.
- [`PAPER_platinum_run_subsets`](https://github.com/jennalandy/gridsemble_PAPER/blob/main/experimental/PAPER_platinum_run_subsets.pdf): analyses on subsets of Platinum Spike data with $\pi_0 \in [0.6, 0.95]$, used to create Figure 2.
- [`PAPER_platinum_run_all_data`](https://github.com/jennalandy/gridsemble_PAPER/blob/main/experimental/PAPER_platinum_run_all_data.pdf): analyses on the full Platinum Spike dataset, used to create Supplementary Figure 3 and Supplementary Table 2.

**Other**

- `PAPER_metrics_helpers.R`: functions to calculate metrics and helper functions.

### References

[1] Q. Zhu, J.C. Miecznikowski, and M.S. Halfon. Preferred analysis methods for affymetrix genechips. II. an expanded, balanced, wholly-defined spike-in dataset. *BMC Bioinformatics*, 11:285, 2010. doi:https://doi.org/10.1186/1471-2105-11-285.
