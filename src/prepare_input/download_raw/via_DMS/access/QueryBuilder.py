from src.data_access.via_DMS.DMSDatabase import DMSDatabase
from src.data_access.via_DMS.secure import Config
from src.data_access.via_DMS.Query import Query

import pandas as pd
import os
import logging
import sys

from utility.utils import logger


class QueryBuilder:
    """1. Load MS-SQl Queries.
    2. Execute them.
    3. Save to job_query_info.csv to disk
     Job |	Dataset | Experiment  | OrganismDBName |ProteinCollectionList  |ParameterFileName
    4. Save to start_file.csv to disk and keeps in analysis_jobs object.
     Dataset_ID  | MSGFPlusJob  | Data Folder Link  | NewestMasicJob  |	Results Folder Path |
    """

    def __init__(self, user_input=None, storage=None, project_name=None):
        """

        :param user_input: parsed list of job_nums/dataset_ids/datapackage_id From Module:Input.py
        :param storage: User-defined data storage location.
        :param project_name: User-defined Study name.
        """
        self.db = DMSDatabase(Config)
        self.user_input = user_input
        self.analysis_jobs = None
        self.job_info = None
        self.parent_data_folder = storage
        self.project_name = project_name
        # TODO: NMDC: 18 Create a Summary_Stats.txt to distinguish between proteomics & Meta-proteomics study!
        self.Summary_Stats = None

    # def create_job_info_query(self):

    def save_to_disk(self, data, data_path, msgf_job_list: list):
        """

        :param data: analysis_jobs object
        :param data_path: User-defined storage followed by
                          data/dpkgs/{}/ or
                          data/set_of_Dataset_IDs/{}/ or
                          data/set_of_Jobs/
        :param msgf_job_list: list of DMS_MSGF+ jobs.
        :return:
        """
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        data.to_csv(data_path + "start_file.csv")
        logger.info(
            msg="@{}: start_file.csv/analysis_jobs_obj shape: {} size:{}".format(
                data_path, data.shape, sys.getsizeof(data)
            )
        )

        query = Query.JOB_INFO.format(",".join(str(job) for job in msgf_job_list))
        result_set = self.db.run_query(query).fetchall()
        df = pd.DataFrame(result_set)
        self.job_info = df
        df.to_csv(data_path + "job_query_info.csv")
        logger.info(msg="@{}: job_query_info.csv shape: {}".format(data_path, df.shape))

    def start_with_datapackage_id(self, id: int):
        """Find me all the jobs for each dataset in this datapackage!
        1 Given a datapackage_id.
        1.0 Find out the Dataset_ID , MSGFPlusJob.
        1.0.0 Using MSGFPlusJob, findout  "Data Folder Link".
        1.0.1 Using Dataset_ID,  findout  NewestMasicJob.
        1.0.1.0 Using NewestMasicJob findout "Results Folder Path".

        2 Merge results to create "analysis_jobs".
        :param id: datapackage_id from Module:Input.py
        :return:
        """
        query = Query.DATASET_MSFG.format(id, id)
        result_set = self.db.run_query(query).fetchall()
        MSGFPlusJobs = pd.DataFrame(result_set)["MSGFPlusJob"]
        Dataset_ID = pd.DataFrame(result_set)["Dataset_ID"]

        query = Query.MSGF_loc.format(
            ",".join(str(job) for job in MSGFPlusJobs.to_list())
        )
        result_set = self.db.run_query(query).fetchall()
        MSGFPlusJob_loc = pd.DataFrame(result_set)["MSGFplus_loc"]

        query = Query.DATASET_MASIC.format(
            ",".join(str(job) for job in Dataset_ID.to_list())
        )
        result_set = self.db.run_query(query).fetchall()
        df_dataset_newest_MASIC = pd.DataFrame(result_set)
        NewestMasicJob = pd.DataFrame(result_set)["NewestMasicJob"]

        query = Query.MASIC_loc.format(
            ",".join(str(job) for job in NewestMasicJob.to_list())
        )
        result_set = self.db.run_query(query).fetchall()
        df_MASIC_loc = pd.DataFrame(result_set)

        df = pd.concat([Dataset_ID, MSGFPlusJobs, MSGFPlusJob_loc], axis=1)

        first = df.merge(
            df_dataset_newest_MASIC, how="left", on=["Dataset_ID", "Dataset_ID"]
        )
        second = first.merge(
            df_MASIC_loc, how="left", on=["NewestMasicJob", "NewestMasicJob"]
        )

        self.analysis_jobs = second
        self.parent_data_folder = (
            self.parent_data_folder + "/" + "data/dpkgs/{}/".format(id)
        )
        self.save_to_disk(second, self.parent_data_folder, MSGFPlusJobs.to_list())

    def start_with_dataset_ids(self, id_list: list):
        """find me all the jobs for each dataset!

        1 Given set of dataset-IDs
        1.0 findout MSGFPlusJob, "Results Folder Path"
        1.1 Using Dataset_ID,  findout  NewestMasicJob
        1.1.0 Using NewestMasicJob findout "Results Folder Path"

        2 Merge results to create "analysis_jobs".
        :param id_list: set of dataset-IDs from Module:Input.py
        :return:
        """
        query = Query.DATASET.format(",".join(str(job) for job in id_list))
        result_set = self.db.run_query(query).fetchall()
        df = pd.DataFrame(result_set)
        MSGFPlusJobs = pd.DataFrame(result_set)["MSGFPlusJob"]

        query = Query.DATASET_MASIC.format(",".join(str(job) for job in id_list))
        result_set = self.db.run_query(query).fetchall()
        df_dataset_newest_MASIC = pd.DataFrame(result_set)
        NewestMasicJob = pd.DataFrame(result_set)["NewestMasicJob"]

        query = Query.MASIC_loc.format(
            ",".join(str(job) for job in NewestMasicJob.to_list())
        )
        result_set = self.db.run_query(query).fetchall()
        df_MASIC_loc = pd.DataFrame(result_set)

        first = df.merge(
            df_dataset_newest_MASIC, how="left", on=["Dataset_ID", "Dataset_ID"]
        )
        second = first.merge(
            df_MASIC_loc, how="left", on=["NewestMasicJob", "NewestMasicJob"]
        )

        self.analysis_jobs = second
        self.parent_data_folder = (
            self.parent_data_folder
            + "/"
            + "data/set_of_Dataset_IDs/{}/".format(self.project_name)
        )

        self.save_to_disk(second, self.parent_data_folder, MSGFPlusJobs.to_list())

    def start_with_job_nums(self, id_list: list):
        """
        1 Given set of MSGFJobs
        1.0 Find the Dataset_ID, & "Results Folder Path"
        1.0.0 Using Dataset_ID, findout MASIC
        1.0.0.0  Using MASIC, findout "Results Folder Path"

        2 Merge results to create "analysis_jobs".
        :param id_list:  set of JobNums from Module:Input.py
        :return:
        """
        query = Query.MSGF.format(",".join(str(job) for job in id_list))
        result_set = self.db.run_query(query).fetchall()
        msgf_loc_dataset = pd.DataFrame(result_set)
        MSGFPlusJobs = pd.DataFrame(result_set)["MSGFPlusJob"]
        Dataset_ID = pd.DataFrame(result_set)["Dataset_ID"]

        query = Query.DATASET_MASIC.format(
            ",".join(str(job) for job in Dataset_ID.to_list())
        )
        result_set = self.db.run_query(query).fetchall()
        df_dataset_newest_MASIC = pd.DataFrame(result_set)
        NewestMasicJob = pd.DataFrame(result_set)["NewestMasicJob"]

        query = Query.MASIC_loc.format(
            ",".join(str(job) for job in NewestMasicJob.to_list())
        )
        result_set = self.db.run_query(query).fetchall()
        df_MASIC_loc = pd.DataFrame(result_set)

        first = msgf_loc_dataset.merge(
            df_dataset_newest_MASIC, how="left", on=["Dataset_ID", "Dataset_ID"]
        )
        second = first.merge(
            df_MASIC_loc, how="left", on=["NewestMasicJob", "NewestMasicJob"]
        )

        self.analysis_jobs = second
        self.parent_data_folder = self.parent_data_folder + "/" + "data/set_of_Jobs/"
        self.save_to_disk(second, self.parent_data_folder, MSGFPlusJobs.to_list())

    def execute(self):
        """Based on the parsed input from Module:Input.py
        It determines which UserType has to be executed and starts it's execution."""
        try:
            # TODO: Remove query redundancy!
            dpkg_id = self.user_input.datapackage_id
            if dpkg_id:
                self.start_with_datapackage_id(dpkg_id)
            dataset_ids = self.user_input.dataset_ids
            if dataset_ids:
                self.start_with_dataset_ids(dataset_ids)
            msgf_jobs = self.user_input.job_nums
            if msgf_jobs:
                # find me all the datasets and their respective jobs!
                self.start_with_job_nums(msgf_jobs)

        except Exception as error:
            raise error
            logger.error(msg="Query Builder Failed: User must provide correct input!")
