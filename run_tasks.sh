# download study raw input files.
docker exec -it downloadDMSDataContainer bash /metaProNew/src/prepare_input/download_raw/download_raw.sh
# scp fastas from NERSC-cori

# mapping.xls to mapping.json
docker exec -it analysisJobContainer python3.8 /metaProNew/src/prepare_input/emsl_to_jgi.py
# run job analysis
docker exec -it analysisJobContainer python3.8 /metaProNew/src/analysis_jobs/run_analysis_job.py

# run report generation.
docker exec -it postProcessingContainer python ./metaPro/src/post_processing/run_fa.py

# generate NMDC-schema complaint metadata
docker exec -it postProcessingContainer python ./metaPro/src/metadata_collection/gen_meta_data.py

# validate against NMDC-schema.
docker exec -it postProcessingContainer python ./metaPro/src/metadata_collection/validate.py
