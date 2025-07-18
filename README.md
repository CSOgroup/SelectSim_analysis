# SelectSim_analysis

## Description
- The GitHub repository consists of scripts to perform the analysis related to SelectSim paper.
- The contents of the repository is mentioned below.
 ```
  .
├── analysis
│   ├── catalogue
│   │   └── notebook
│   │       ├── catalouge_consistent.ipynb
│   │       ├── dfci_merge_results.ipynb
│   │       ├── dfci_runs.ipynb
│   │       ├── merge_all_results.ipynb
│   │       ├── msk_merge_results.ipynb
│   │       ├── msk_runs.ipynb
│   │       ├── tcga_merge_results.ipynb
│   │       └── tcga_runs.ipynb
│   ├── ground_truth_benchmarking
│   │   ├── 1_design_ground_truth_modelB.R
│   │   └── 2_run_gt_all_methods.R
│   ├── hotspot_analysis
│   │   ├── data_processing.ipynb
│   │   ├── hotspot_pan_can_analysis.ipynb
│   │   └── SelectSim_runs.ipynb
│   ├── normal_mutation
│   │   ├── 1_preprocess_Abby2023.R
│   │   ├── 1_preprocess_ESOPHAGUS.R
│   │   ├── 1_preprocess_Fowler.R
│   │   ├── 1_preprocess_Lawson.R
│   │   ├── 1_preprocess_Martincorena2015.R
│   │   ├── 1_preprocess_Martincorena2018.R
│   │   ├── 1_preprocess_SKIN.R
│   │   ├── 2_run_ESOPHAGUS.R
│   │   ├── 2_run_Lawson.R
│   │   ├── 2_run_SKIN.R
│   │   ├── notebook
│   │   │   └── data_processing.ipynb
│   │   └── run_selectSim.r
│   ├── oncogenic_non_oncogenic
│   │   └── mix_analysis.ipynb
│   ├── oncokb_random_analysis
│   │   └── selectsim_oncokb_random_analysis.ipynb
│   ├── primary_met
│   │   └── notebook
│   │       ├── data_processing.ipynb
│   │       ├── msk_p_m_enrichment.ipynb
│   │       ├── msk_p_m_gene_frequency.ipynb
│   │       └── msk_p_m_runs.ipynb
│   ├── primary_met_elastic_net
│   │   ├── 1_preprocessTreatment.R
│   │   ├── 2_run_EN.R
│   │   ├── 3_merge_results.R
│   │   └── EN_functions.R
│   └── roubustness_benchmarking
│       ├── coselens.r
│       ├── notebook
│       │   ├── cooccur.ipynb
│       │   ├── Discover.ipynb
│       │   ├── Discover_SelectSim.ipynb
│       │   ├── GamTOC.ipynb
│       │   ├── MCC.ipynb
│       │   ├── mutational_load.ipynb
│       │   ├── sampling_recovery.ipynb
│       │   ├── select.ipynb
│       │   ├── SelectSim.ipynb
│       │   ├── tcga_robustness_analysis.ipynb
│       │   └── WeSME_CO.ipynb
│       └── time_analysis_berewire_simulation.r
├── data_processing
│   └── notebook
│       ├── genie_data_processing.ipynb
│       ├── msk_dfci_gam_create.ipynb
│       ├── msk_dfci_sample_filtering_annotation.ipynb
│       └── tcga_gam_create.ipynb
└── README.md

16 directories, 53 files
```
- The analysis folder consists of scripts (Jupyter R-notebooks and.R script) files related to different analysis (catalogue, normal mutation, primary vs metastasis etc.)
- The data processing folder consists of scripts (Jupyter R-notebooks) to process the data to generate the `SelectSim` algorithm input run objects.

## Data Availability

- Data needed for the scripts are available from [zenodo](https://zenodo.org/records/13752870) and [zenodo](https://zenodo.org/records/13769786) on request.

## SelectSim Algorithm

- SelectSim algorithm is available from [github](https://github.com/CSOgroup/SelectSim).

## Who do I talk to?

- For any bugs or question/data request (missing), please use the [issue tracker](https://github.com/CSOgroup/SelectSim_analysis/issues) to raise your query.
- For any other question related to SelectSim, please contact Prof Giovanni Ciriello (<giovanni.ciriello@unil.ch>) or Arvind Iyer (<arvind.iyer@unil.ch>), Miljan Petrovic (miljan.petrovic@unil.ch).
