Metaproteomics Workflow
==============================

Summary
-------
The metaproteomics workflow/pipeline is an end-to-end data processing workflow for protein identification and characterization using MS/MS data. Briefly, mass spectrometry instrument generated data files(.RAW) are converted to mzML, an open data format, using MSConvert. Peptide identification is achieved using MSGF+ and the associated metagenomic information in the FASTA (protein sequences) file format. Intensity information for identified species is extracted using MASIC and combined with protein information.

Workflow Diagram
------------------

.. image:: workflow_diagram.png

Workflow Dependencies
---------------------

Third party software
~~~~~~~~~~~~~~~~~~~~
|                            |                                          |
|----------------------------|------------------------------------------|
| MSGFPlus                   | v20190628                                |
| Mzid-To-Tsv-Converter      | v1.3.3                                   |
| PeptideHitResultsProcessor | v1.5.7130                                |
| pwiz-bin-windows           | x86_64-vc141-release-3_0_20149_b73158966 |
| MASIC                      | v3.0.7235                                |
| sqlite-netFx-full-source   | 1.0.111.0                                |
| Conda                      | (3-clause BSD)                           |
|                            |                                          |

Workflow Availability
---------------------

The workflow is available in GitHub:
https://github.com/microbiomedata/metaPro

The container is available at Docker Hub (microbiomedata/mg-annotation):
https://hub.docker.com/r/microbiomedata/mepro

Inputs
~~~~~~~~

- .raw, FASTA, parameter files

Outputs
~~~~~~~~

1. Processing multiple datasets.
.
├── Data/
├── FDR_table.csv
├── Plots/
├── dataset_job_map.csv
├── peak_area_crosstab_by_dataset_id.csv
├── protein_peptide_map.csv
├── specID_table.csv
└── spectra_count_crosstab_by_dataset_id.csv

2. Processing single FICUS dataset.

| DatasetName | PeptideSequence | FirstHitProtein | SpectralCount | sum(MasicAbundance) | GeneCount | FullGeneList | FirstHitDescription | DescriptionList | min(Qvalue) |
|-------------|-----------------|-----------------|---------------|---------------------|-----------|--------------|---------------------|-----------------|-------------|


Requirements for Execution
--------------------------

- Docker or other Container Runtime

Version History
---------------

- 1.0.0

Point of contact
----------------

Package maintainer: Anubhav <anubhav@pnnl.gov>
