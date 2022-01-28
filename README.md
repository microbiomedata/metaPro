# Meta-proteomics workflow.
<hr>

![workflow](docs/workflow_diagram.png)  

<hr>

[![Battelle Memorial Institute license](https://img.shields.io/badge/license-Battelle%20Memorial%20Institute-green)](https://github.com/microbiomedata/metaPro/blob/master/license.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5748584.svg)](https://doi.org/10.5281/zenodo.5748584)
![GitHub issues](https://img.shields.io/github/issues-raw/microbiomedata/metaPro)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/microbiomedata/metaPro)
![Lines of code](https://img.shields.io/tokei/lines/github/microbiomedata/metaPro)
<hr>

## About
Meta-proteomics workflow is an end-to-end data processing and analyzing pipeline for studying proteomes i.e studying [protein identification and characterization using MS/MS data](https://www.sciencedirect.com/science/article/pii/S1874391911002053).

We identify the active organisms/species in a metagenome corresponding to a wet-lab sample obtained from JGI after gene sequencing. Then the researchers at PNNL culture these samples and make it appropriate to study it as a protein sample. This protein sample may have a single protein or a complex mixture of proteins. Later, this sample is passed through a [mass spectrometry](https://nationalmaglab.org/user-facilities/icr/techniques/tandem-ms) instrument to obtain a proprietary data format .RAW file. This file contains MS/MS spectrum i.e mass analysis(mass-to-charge (*m/z*) ratios) for each peptide sequences identified in the sample. 
### How to run the workflow:

- Python codebase:
    
  - `Virtual machine`:
    1. Make your input datasets ready- [as described here](src/prepare_input/prepare_input.md) 
       - Make your input `storage/` folder visible to workflow. [You need to provide path in docker-compose.yml](docker-compose.yml)
         
         Note: 
           - I left `./storage/` already configured assuming you kept inputs in the project directory itself.
           - Typically, a Study(such as `stegen`) has more than 1 datasets(RAW files-MSMS spectra) and multiple fastas to search against.
            This is information is must and a sample is provide [here](storage/mappings/EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping_2021-01-25.xlsx)
             
    2. Configure workflow as per need. Typically, we run in following ways:
       1. Fully-Tryptic with No modifications (recommended for large datasets such as Prosser Soil.)
       2. Fully-Tryptic with modifications    
       3. partially-tryptic with Modification( such as MetOx).
       4. partially-tryptic No Modification.
         
          Notes:
    
           - User need to tweek [configuration file](metaPro-config.env). To reproduce results achieved for FICUS dataset studies(`Hess`, `Stegen`, `Blanchard`)
           - we provided [parameter files](storage/parameters) and a [pre-configured env file](conf.env) that could be use to run the workflow.
    4. Must have installed [docker](https://docs.docker.com/get-docker/) and [docker-compose](https://docs.docker.com/compose/install/) on your system.
       
    5. To run workflow, From [project directory](./):
           1. `make build_unified` to start services.
              Notes: - (to take containers down and remove volumes: `docker-compose down -v`)
           2. `make run_workflow` 
              It will create a `storage/results` folder and create all the necessary files.

    
  - HPC([Tahoma- high performance computing cluster at EMSL](https://www.emsl.pnnl.gov/MSC/UserGuide/tahoma/tahoma_overview.html)):
      1. Pull image from dockerhub:
         ```
         singularity pull docker://microbiomedata/mepro:2.0.0
         singularity sif list mepro_2.0.0.sif
         singularity run docker://microbiomedata/mepro:2.0.0
         ```
      2. 
<hr>

- WDL support codebase:

    1. prepare you input.json
        `make prepare-your-input`
       Note: User need to generate the input.json file based on the 
             - mapping (dataset(raw) to annotations(.faa & .gff ))
             - actual files respective file locations. 
             For you help, a [script](wdl/scripts/prepare_input.sh) has been provided.
    2. run the WDL:
       Need an 
       - execution engine(tested with [cromwell-66](https://github.com/broadinstitute/cromwell/releases)) to run WDL 
       - along with Java runtime(tested with `openjdk 12.0.1`)
             1. if docker support 
                1. `make run_wdl`
             2. if shifter support to run on cori:
                1. `make run_wdl_on_cori`


[More about workflow...](docs/index.rst)

[Documentation](https://nmdc-proteomics-workflow.readthedocs.io/en/latest/)