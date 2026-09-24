# ESMC Pathway

ESMC Pathway is a multi-regional, multi-year extension of **EnergyScope Multi-Cell (ESMC)**.
It optimizes the energy transition pathway of a multi-regional energy system from 2015 to 2050
with a linear programming model written in AMPL.

Main features:
- Multi-regional whole-energy system with exchanges between regions.
- Transition pathway over 2015-2050 (phases of 4-6 years), with investment, decommissioning and change-rate limits.
- Rolling horizon (myopic) optimization; a single window covering all years gives perfect foresight.
- Typical days selected by k-medoids clustering, with an optional intra-day segmentation (tsam).
- Post-processing to CSV files and interactive plots (Plotly).

The current case study is Bolivia with six regions: the national grid (SIN) split into
Highlands, Lowlands and Valleys, and three isolated systems (Norte, Tarija, Santa Cruz).

## Requirements

- Python >= 3.9
- [AMPL](https://ampl.com/) with the CPLEX solver (license required)
- Python packages: numpy, pandas, amplpy, plotly, kaleido
- Optional: tsam (only for intra-day segmentation, `--n-seg` < 24)

Install with conda:

```bash
conda env create -f environment.yml
conda activate esmc_pathway
```

or with pip:

```bash
pip install -e .
```

## How to run

From the repository root:

```bash
python run.py                          # full run (optimization + post-processing)
python run.py --window-size 15 --overlap 5
python run.py --skip-optimize          # only post-processing, from _Results.pkl
python run.py --n-seg 12               # intra-day segmentation (12 segments per typical day)
python run.py --help                   # all options
```

Main settings (years, phases, number of typical days, window size, case name) are at the top of `run.py`.
Results are written to `case_studies/<regions>/PATHWAY_<case>/`.

To regenerate only the yearly .dat files: `python -m esmc.preprocessing.generate_year_dat_files`

## Repository structure

```
run.py                             Main script (rolling horizon)
esmc/
  __init__.py                      Logging setup (console + one log file per run in logs/)
  common.py                        Global variables (repository root, region codes, colors, names)
  utils/
    esmcpt.py                      EsmcPathway class (model setup, solving, results, plots)
    esmc_preprocessor.py           Optimization windows management
    esmc_collector.py              Collects results across windows
    region.py                      Region data
    opti_probl.py                  AMPL interface
    df_utils.py                    DataFrame helpers
  preprocessing/
    generate_year_dat_files.py     Builds the yearly .dat files from the CSV data
    segmentation_tsam.py           Intra-day segmentation of typical days
    temporal_aggregation.py        Typical days selection (k-medoids)
    dat_print.py                   .dat file writing
    kmedoid_clustering/            AMPL model for the typical days selection
  postprocessing/                  amplpy/pandas helpers
  energy_model/                    AMPL model (.mod) and headers
Data/                              Input data per year and region
case_studies/                      Outputs (generated)
logs/                              Run logs (generated)
```

## Logs

- `logs/<date>_<time>_esmc.log`: full log of each run (Python side).
- `case_studies/.../PATHWAY_<case>/log.txt`: AMPL/CPLEX solver log.
- `case_studies/<regions>/00_td_dat/log_<N>.txt`: AMPL log of the typical days selection.

## Authors

- Pablo Jimenez Zabalaga - ESMC Pathway (pathway extension, rolling horizon, intra-day segmentation, Bolivian case study)

ESMC Pathway is built on EnergyScope Multi-Cell and EnergyScope. Previous versions and authors:
- Stefano Moret, EPFL (Switzerland)
- Gauthier Limpens, UCLouvain (Belgium)
- Paolo Thiran, UCLouvain (Belgium) - EnergyScope Multi-Cell code
- Aurélia Hernandez, UCLouvain (Belgium)
- Noé Cornet, UCLouvain (Belgium)
- Pauline Eloy, UCLouvain (Belgium)

Previous releases:
- EnergyScope (v1, v2): https://github.com/energyscope/EnergyScope
- First Multi-Cell release: https://github.com/pathiran22/EnergyScope/tree/Hernandez_Thiran_Multi_cell_2020
- EnergyScope MC for Western Europe: https://github.com/16NoCo/EnergyScope/tree/Multi_cell_West-Eu_2021

Please report bugs through GitHub issues.

## How to cite

Please cite the references below when using this code (see also `NOTICE`):

[1] G. Limpens, S. Moret, H. Jeanmart, F. Maréchal (2019). EnergyScope TD: a novel open-source model for regional energy systems and its application to the case of Switzerland. https://doi.org/10.1016/j.apenergy.2019.113729

[2] V. Codina Gironès, S. Moret, F. Maréchal, D. Favrat (2015). Strategic energy planning for large-scale energy systems: A modelling framework to aid decision-making. Energy, 90(PA1), 173–186. https://doi.org/10.1016/j.energy.2015.06.008

[3] S. Moret, M. Bierlaire, F. Maréchal (2016). Strategic Energy Planning under Uncertainty: a Mixed-Integer Linear Programming Modeling Framework for Large-Scale Energy Systems. https://doi.org/10.1016/B978-0-444-63428-3.50321-0

[4] Limpens, G. (2021). Generating energy transition pathways: application to Belgium.

[5] Hernandez, A., Thiran, P., Jeanmart, H., & Limpens, G. (2020). EnergyScope Multi-Cell: A novel open-source model for multi-regional energy systems and application to a 3-cell, low-carbon energy system [UCLouvain]. http://hdl.handle.net/2078.1/thesis:25229

[6] Cornet, N., Eloy, P., Jeanmart, H., & Limpens, G. (2021). Energy Exchanges between Countries for a Future Low-Carbon Western Europe By merging cells in EnergyScope MC to handle wider regions.

[7] Thiran, P., Hernandez, A., Limpens, G., Prina, M. G., Jeanmart, H., & Contino, F. (2021). Flexibility options in a multi-regional whole-energy system: the role of energy carriers in the Italian energy transition. Proceedings of ECOS 2021, 1–12.

## License

Licensed under the Apache License, Version 2.0. See `LICENSE` and `NOTICE`.
