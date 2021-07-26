
docker exec -it analysisJobContainer python3.8 ./metaPro/src/prepare_input/emsl_to_jgi.py

docker exec -it analysisJobContainer python3.8 ./metaPro/src/analysis_jobs/run_analysis_job.py

docker exec -it postProcessingContainer python ./metaPro/src/post_processing/run_fa.py

docker exec -it postProcessingContainer python ./metaPro/src/metadata_collection/gen_meta_data.py

#docker exec -it downloadDMSDataContainer python ./metaPro/src/post_processing/run_fa.py
