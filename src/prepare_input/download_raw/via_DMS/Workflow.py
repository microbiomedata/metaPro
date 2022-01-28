from src.data_access.via_DMS.Input import Input
from src.data_access.via_DMS.QueryBuilder import QueryBuilder
from src.data_access.via_DMS.FileOperations import FileOperations
from src.merge.DatasetsMerger import DatasetsMerger
from src.analysis.FicusAnalysis import FicusAnalysis
from src.analysis.InternalAnalysis import InternalAnalysis

import os
import fnmatch
from utility.utils import logger


class Workflow:
    """Automate the Meta-proteomics workflow"""

    def __init__(
        self,
        mode=None,
        InputType=None,
        path_to_data=None,
        project_name=None,
        UserInput=None,
        CombineDatasets=None,
        SelectAnalysis=None,
        PipeLineMode=None,
    ):
        self.Mode = mode
        self.InputType = InputType
        self.Storage = path_to_data
        self.Project = project_name
        self.UserInput = UserInput
        self.CombineDatasets = CombineDatasets
        self.SelectAnalysis = SelectAnalysis
        self.PipeLineMode = PipeLineMode

    def run_Analysis(self, on_file, analysis_type):
        if analysis_type == "internal":
            print("Run Internal analysis on {}".format(on_file))
            # internalAnalysis(on_file)
        elif analysis_type == "ficus":
            print("Run Ficus analysis on {}".format(on_file))
            fa = FicusAnalysis(on_file, self.Storage)
            fa.create_DatasetOutputTable()
        else:  # "both"
            print("Run Internal & Ficus analysis on {}".format(on_file))
            # internalAnalysis(on_file)
            # ficusAnalysis(on_file)

    def start_downStreamAnalysis(self, result_path):

        if self.CombineDatasets:
            # generate report on single file.
            self.run_Analysis(result_path + "resultants_df.tsv", self.SelectAnalysis)
        else:
            # generate report on multiple files.
            # stop =0
            for path, subdirs, files in os.walk(result_path):
                for file in files:
                    if fnmatch.fnmatch(file, "MSGFjobs_MASIC_resultant.tsv"):

                        # if stop!=1:
                        self.run_Analysis(os.path.join(path, file), self.SelectAnalysis)
                        # else:
                        #     break
                        # stop += 1

    def start_merging(self, folder):
        merge = DatasetsMerger(folder, self.CombineDatasets)
        result_path = merge.merge_all_jobs_in_UserInput()
        return result_path

    def download_data(self, user_obj):
        """
        If PNNL
        1. build & execute SQL queries to dowload data from DMS.
        2. Data: 1. PNNL Fasta & parameter files
                 2. PNNL .raw files
                 3. PNNL, DMS processed MASIC and MSGF+ jobs results.
        If NMDC
        :param user_obj: parsed list of job_nums/dataset_ids/datapackage_id
        :return: path to each data.
        """
        myQuery = QueryBuilder(user_obj, self.Storage, self.Project)
        myQuery.execute()
        analysis_jobs, parent_data_folder, job_info = (
            myQuery.analysis_jobs,
            myQuery.parent_data_folder,
            myQuery.job_info,
        )
        # TODO: Make it run for PNNl/NMDC
        # 1. Class Internal
        #     Download DMS processed MASIC, MSGF+jobs

        # 2. Class External
        #     Download DMS .raw, parameter files and , FASTA from NERSC!
        file_obj = FileOperations(analysis_jobs, parent_data_folder, job_info)
        file_obj.get_files()
        return parent_data_folder

    def start_workflow(self):
        """
        Runs the workflow in 3 stages
        1. Download relevant datasets from specified source.
        2. Aggregation of analysis tools{MSGF+, MASIC} results: to extract useful data from datasets.
        3. Generation of experimental report.

        :return:
        """
        # TODO: Use User-mode: to suppress file creations & Developer-mode: to Generate files!
        ## prepare user's input
        user_obj = Input()
        if self.InputType is None:
            # user_obj.user_input() # nomore manual execution
            pass
        else:
            user_obj.other_input(self.InputType, self.UserInput)
            ## 1.
            logger.info("1. Start Download relevant datasets from specified source.")
            data_parent_folder = self.download_data(user_obj)
            logger.info(msg="Input Data located at  @:{}".format(data_parent_folder))
            # ## 2.
            # logger.info(msg="2. Start Aggregation of analysis tools{MSGF+, MASIC} results")
            # result_path= self.start_merging(data_parent_folder)
            # logger.info(msg="Merged jobs located at @:{}".format(result_path))
            # ## 3.
            # logger.info(msg="3. Generation of experimental report.")
            # self.start_downStreamAnalysis(result_path)
            # logger.info(msg="Generated reports      @:{}".format(result_path))

        logger.info(msg="```````````````")
        logger.info(msg="Finished running Meta-proteomics pipeline!")
