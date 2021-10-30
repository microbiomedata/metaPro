### Meta-proteomics data processing workflow for Mass Spectrometry identifications.

Run the workflow:
- Using docker:
    1. Make your input datasets ready- [as described here](src/prepare_input/prepare_input.md) 
       - Make you input `storage/` folder visible to workflow. [You need to provide path in docker-compose.yml](docker-compose.yml)
         Note: 
           - I left `./storage/` already configured assuming you kept inputs in the project directory itself.
           - Typically, a Study(such as `stegen`) has more than 1 datasets(RAW files-MSMS spectra) and multiple fastas to search against.
             This is information is must and a sample is provide [here](storage/mappings/EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping_2021-01-25.xlsx)
             
    2. Configure workflow as per need. Typically, we run in following ways:
        1. Fully-Tryptic with No modifications (recommended for large datasets such as Prosser Soil.)
        2. Fully-Tryptic with modifications    
        3. partially-tryptic with Modification( such as MetOx).
        4. partially-tryptic No Modification.
         Notes: - User need to tweek [configuration file](conf.env). To reproduce results achieved for FICUS dataset studies(`Hess`, `Stegen`, `Blanchard`)
                - we provided [parameter files](storage/parameters) and a [pre-configured env file](conf.env) that could be use to run the workflow.
    3. Must have installed [docker](https://docs.docker.com/get-docker/) and [docker-compose](https://docs.docker.com/compose/install/) on your system.
       
    4. To run workflow, From [project directory](./):
        1. `make build_unified` to start services.
           Notes: - (to take containers down and remove volumes: `docker-compose down -v`)
        2. `make run_workflow` 
           It will create a `storage/results` folder and create all the necessary files.

- Using WDL
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