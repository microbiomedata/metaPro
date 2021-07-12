# A Meta-proteomics data processing pipeline(Workflow) for Mass Spectrometry identifications.

Run the workflow:
- Using docker:
    1. Structure your input datasets following as [described here](src/prepare_input/prepare_input.md) 
       and add your `storage` in [docker-compose.yml file](docker-compose.yml) to work as docker volume.
         Notes: - Typically, a Study(such as `Stegen`) has more than 1 datasets(RAW files-MSMS specra) and multiple fastas to search against.
                  You can [mount your network-share drives using this](utility/mount_network_drive.sh)(optional!)
    2. Configure workflow as per need, How would you like to run the workflow:
        1. Fully-Tryptic with No modifications (recommened for large datasets such as Prosser Soil.)
        2. partially-tryptic with Modification( such as MetOx).
         Notes: - User need to tweek [configuration file](conf.env). To reproduce results achieved for FICUS dataset studies(`Hess`, `Stegen`, `Blanchard`) we have provided parameter files that could be use to run the workflow in above configurations.
                - TODO, add a example usecase for each configuration.
    3. Must have installed docker and docker-compose.
       - from project directory, run `docker-compose up --build`
         Notes: (to take containers down and remove volumes: `docker-compose down -v`)
    

- Using Cromwell-WDL:
