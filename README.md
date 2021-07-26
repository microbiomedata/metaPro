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
        1. Run 
           `docker-compose up --build` to start services.
           Notes: - (to take containers down and remove volumes: `docker-compose down -v`)
        2. Run 
           `bash run_task.sh` 
           It will create a `storage/results` folder and create all the necessary files.

- Using Cromwell-WDL:


[More about workflow...](docs/index.rst)