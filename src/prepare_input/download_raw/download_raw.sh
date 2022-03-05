#!/bin/sh
export DATABASE_USERNAME=dmsreader DATABASE_PASSWORD=dms4fun DATABASE_SERVER=gigasax DATABASE_NAME=DMS5

storage="/app/storage"

projectName="spruce"
datasets="598849, 598852, 598515, 598853, 598850, 598509, 598517, 598847, 598511, 598506, 598846, 598516, 598513, 598520, 598848, 598851, 598512, 598508, 598510"

# Based on dataset_ids | DMS or NMDC-FICUS
python3.8 /metaProNew/src/prepare_input/download_raw/via_DMS/MetProWorkflowApp.py \
    --Mode developer \
    --InputType 2 \
    --Storage $storage \
    --ProjectName $projectName \
    --Input "$datasets"