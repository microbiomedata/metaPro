docker exec -it metapro_analysis-job_1 python3.8 ./metaPro/src/prepare_input/emsl_to_jgi.py
# echo "========"
docker exec -it metapro_analysis-job_1 python3.8 ./metaPro/src/analysis_jobs/run_analysis_job.py
# echo "========"
docker exec -it postProcessingContainer python ./metaPro/src/post_processing/run_fa.py
#echo "========"
docker exec -it postProcessingContainer python ./metaPro/src/metadata_collection/gen_meta_data.py
#echo "========"
docker exec -it downloadDMSDataContainer python ./metaPro/src/post_processing/run_fa.py
