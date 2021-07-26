from src.analysis_jobs.merge_jobs.MSGFplusMerger import MSGFplusMerger
from src.analysis_jobs.merge_jobs.MASICmerger import MASICmerger
import os
import pandas as pd
from utility.utils import stats,logger
import fnmatch

class DatasetsMerger(MSGFplusMerger):
    ''' MANDATORY_INPUT : [../data/set_of_Dataset_IDs/<STUDY>/*]
        Run for each dataset in [../data/set_of_Dataset_IDs/<STUDY>/*]
                  and generate [../results/set_of_Dataset_IDs/<STUDY>/*] keeping file hierarchy as of data folder.
             1. Run for UserInput:
                 a datapackage or
                 a set of datasets or
                 a set of MSGFJobNums
             2. create a crossTab object
    '''
    def __init__(self, results= None, combineDatasets=None):
        self.resultants = []
        self.save_jobs_at= results
        self.resultants_df= None
        self.crossTab = None
        self.combineDatasets= combineDatasets
    @stats
    def merge_all_jobs_in_UserInput(self):
        '''
        1. Run for each dataset.
        2. Merge all MSGFjobs_MASIC_resultant objects.
        :return:
        '''
        msfg_obj= MSGFplusMerger(self.save_jobs_at)
        msfg_obj.consolidate_syn_files()
        masic = MASICmerger(self.save_jobs_at)
        masic.merge_msgfplus_msaic(msfg_obj.MSGFjobs_Merged)
        return masic.MSGFjobs_MASIC_resultant
