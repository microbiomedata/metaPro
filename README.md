# Meta-proteomics workflow.
<hr>

![workflow](docs/workflow_diagram.png)  

<hr>

[![Battelle Memorial Institute license](https://img.shields.io/badge/license-Battelle%20Memorial%20Institute-green)](https://github.com/microbiomedata/metaPro/blob/master/license.txt)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5748584.svg)](https://doi.org/10.5281/zenodo.5748584)
![GitHub issues](https://img.shields.io/github/issues-raw/microbiomedata/metaPro)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/microbiomedata/metaPro)
<hr>

## About
<!-- Meta-proteomics workflow is an end-to-end data processing and analyzing pipeline for studying proteomes i.e studying [protein identification and characterization using MS/MS data](https://www.sciencedirect.com/science/article/pii/S1874391911002053). -->

The NMDC metaproteomics workflow is an end-to-end data processing pipeline for studying the proteomes of complex communities using data dependent LC-MS/MS data. The workflow matches the sequences from a metagenome annotation to fragmentation spectra (MS2) in the mass spectrometry data to identify peptide sequences, which are fragments of proteins contained in a processed proteomic sample. The metagenome annotation can either come from genomic sequencing that has been processed using the NMDC metagenomics pipeline, or for the more recent development that removes the paired NMDC metagenome requirement, by utilizing a neural network model to identify peptide sequences de novo from the fragmentation spectra (MS2). This sequence information is then used to identify model organisms with genome sequences archived on UniProt. The matched organisms genomes are then concatenated to create a pseudo-metagenome for in-depth peptide identification. Relative abundance measurements are derived by taking the area under the LC elution curve for the peptide from the intact spectra (MS1).
Metaproteomics analyses identify the active organisms/species in a microbial community through the identification and relative quantification of protein expression.

### How to run the workflow:

<!-- - Python codebase:
    
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
<hr> -->

- WDL support codebase:

    <!-- 1. prepare you input.json
        `make prepare-your-input`
       Note: User need to generate the input.json file based on the 
             - mapping (dataset(raw) to annotations(.faa & .gff ))
             - actual files respective file locations. 
             For you help, a [script](wdl/scripts/prepare_input.sh) has been provided. -->

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
