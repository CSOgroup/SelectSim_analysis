# SelectSim_analysis

## Description
- The github repos consist of scripts to perform the analysis related to SelectSim paper.
- The contents of the repository is mentioned below.
    ```
    ├── analysis
    │   ├── catalogue
    │   │   └── notebook
    │   │       ├── catalouge_consistent.ipynb
    │   │       ├── dfci_merge_results.ipynb
    │   │       ├── dfci_runs.ipynb
    │   │       ├── merge_all_results.ipynb
    │   │       ├── msk_merge_results.ipynb
    │   │       ├── msk_runs.ipynb
    │   │       ├── tcga_merge_results.ipynb
    │   │       └── tcga_runs.ipynb
    │   ├── normal_mutation
    │   │   ├── notebook
    │   │   │   └── data_processing.ipynb
    │   │   └── run_selectSim.r
    │   ├── primary_met
    │   │   └── notebook
    │   │       ├── data_processing.ipynb
    │   │       ├── msk_p_m_enrichment.ipynb
    │   │       ├── msk_p_m_gene_frequency.ipynb
    │   │       └── msk_p_m_runs.ipynb
    │   ├── primary_met_elastic_net
    │   │   ├── EN_functions.r
    │   │   └── run_EN.r
    │   └── roubustness_benchmarking
    │       ├── coselens.r
    │       ├── notebook
    │       │   ├── cooccur.ipynb
    │       │   ├── Discover.ipynb
    │       │   ├── Discover_SelectSim.ipynb
    │       │   ├── GamTOC.ipynb
    │       │   ├── MCC.ipynb
    │       │   ├── mutational_load.ipynb
    │       │   ├── sampling_recovery.ipynb
    │       │   ├── select.ipynb
    │       │   ├── SelectSim.ipynb
    │       │   ├── tcga_robustness_analysis.ipynb
    │       │   └── WeSME_CO.ipynb
    │       └── time_analysis_berewire_simulation.r
    ├── data_processing
    │   └── notebook
    │       ├── genie_data_processing.ipynb
    │       ├── msk_dfci_gam_create.ipynb
    │       ├── msk_dfci_sample_filtering_annotation.ipynb
    │       └── tcga_gam_create.ipynb
    └── README.md
    ```
- The analysis folder consists of scripts (jupyter R-notebooks and .R script) files related to different analysis (catalogue, normal mutation, primary vs metastasis etc.)
- The data processing folder consists of scripts (jupyter R-notebooks) to process the data to generate the `SelectSim` algorithm input run objects.

## Data Availablity

- Data needed for the scripts are avialabe from [zenodo](https://zenodo.org/records/13752870) and [zenodo](https://zenodo.org/records/13769786) on request.

## SelectSim Algorithm

- SelectSim algorithm is avilable from [github](https://github.com/CSOgroup/SelectSim).

## Who do I talk to?

- For any bugs or feature support in using SelectSim, please use the [issue tracker](https://github.com/CSOgroup/SelectSim_analysis/issues).
- For any other question related to SelectSim, please contact Prof Giovanni Ciriello (<giovanni.ciriello@unil.ch>) or Arvind Iyer (<arvind.iyer@unil.ch>).
