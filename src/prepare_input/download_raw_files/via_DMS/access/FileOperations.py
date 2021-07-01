import os
import fnmatch
import requests
import logging
from bs4 import BeautifulSoup
from urllib.parse import urlparse

from utility.utils import stats,logger

class FileOperations:
    ''' MANDATORY_INPUT : "analysis_jobs object" a.k.a "start_file.csv" at
                          [../data/set_of_Dataset_IDs/<STUDY>/*] created by DMS which has information from where to
                          - download DMS processed: MSGFplus, MASIC jobs.
                            "start_file.csv" :: [Dataset_ID,MSGFPlusJob,MSGFplus_loc,NewestMasicJob,MASIC_loc]
                          - Download DMS processed: .raw files:
                            "start_file.csv" :: "MSGFplus_loc", one level up in file hierarchy:
                          - Download DMS fasta files from either
                            url1="http://gigasax/DMS_Organism_Files/Microbial_Communities/FASTA/" or
                            url2="http://gigasax/DMS_FASTA_File_Archive/dynamic/forward/"
                         - Downlaod DMS parameter files from
                            http://gigasax/DMS_Parameter_Files/MSGFPlus/
                            http://gigasax/DMS_Parameter_Files/MASIC/

        Read each data set in [../data/set_of_Dataset_IDs/<STUDY>/*]
        and download
                    1. PNNL Fasta & parameter files
                    2. PNNL .raw files
                    3. PNNL, DMS processed MASIC and MSGF+ jobs results.
    '''

    def __init__(self, analysis_jobs= None, parent_folder= None, job_info=None):
        '''

        :param analysis_jobs: set in Module: QueryBuilder.py
        :param parent_folder: set in Module: QueryBuilder.py
                              User-defined storage followed by
                              data/dpkgs/{}/ or
                              data/set_of_Dataset_IDs/{}/ or
                              data/set_of_Jobs/
        :param job_info: set in Module: QueryBuilder.py
        '''
        self.Input = analysis_jobs
        self.parent_folder = parent_folder
        self.job_info= job_info
        self.url = None
        self.started_from= parent_folder
        self.file_pattern_types = [
                                   "*syn.txt",
                                   "*SeqToProteinMap.txt",
                                   "*ResultToSeqMap.txt",
                                   "*_SICstats.txt",
                                   "*.raw"]
        self.row=None

    def create_dir(self,folder):
        if not os.path.exists(folder):
            os.makedirs(folder)
        os.chdir(folder)

    def write_to_disk(self, url: str):
        '''
        :param url: Job's file path on DMS.
        :return:
        '''
        if not os.path.isfile(url.split('/')[-1]):
            try:
                os.system('wget %s' % url)
                # logging.info("Files transferred successfully!")
            except Exception as e:
                logger.info("FAILED to download file!")

    def check_url(self, url):
        response = requests.get(url)
        try:
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            print("Error: " + str(e)+ ", if .fasta, searching on another location!")
            return False
        return True

    def download_over_http(self):
        '''
        Given a url, copy files from DMS to disk!
        :return:
        '''
        if self.check_url(self.url):
            response = requests.get(self.url)
            soup = BeautifulSoup(response.text, 'html.parser')
            filenames= [ link.get('href') for link in soup.find_all('a')]
            for file in filenames:
                for p in self.file_pattern_types:
                    if fnmatch.fnmatch(file, p):
                        parsed_uri = urlparse(self.url)
                        domain_name = '{uri.scheme}://{uri.netloc}/'.format(uri=parsed_uri)
                        file_url = domain_name + file
                        self.write_to_disk(file_url)

    def parse_fileserverpath_to_web_url(self, file_server_path):
        ''' Converts Windows FileSever path to webURL.
        :param file_server_path: windows server file path.
        :return:
        '''
        folders = file_server_path.split('\\')
        self.url = 'http://' + folders[2] + '.pnl.gov/' + '/'.join(folders[3:])

    def download_msgf_jobs(self,df):
        '''
        Downloads DMS_MSGF+ plus.
        :param df: row corresponding to a dataset
        :return:
        '''
        self.create_dir(self.parent_folder + '/' + 'DMS_MSGFjobs' + '/' + str(df['MSGFPlusJob']) )
        path_or_url = str(df['MSGFplus_loc'])
        if not path_or_url.startswith("http"):
            self.parse_fileserverpath_to_web_url(path_or_url)
            self.download_over_http()
        else:
            self.url = path_or_url
            self.download_over_http()
        return path_or_url

    def download_masic_jobs(self,df):
        '''
        Downloads DMS_MASIC plus.
        :param df: row corresponding to a dataset
        :return:
        '''
        self.create_dir(self.parent_folder + '/' + 'DMS_MASICjob' + '/' + str(df['NewestMasicJob']) )
        # print("donwload_masic", os.getcwd())
        path_or_url = str(df['MASIC_loc'])
        if not path_or_url.startswith("http"):
            self.parse_fileserverpath_to_web_url(path_or_url)
            self.download_over_http()
        else:
            self.url= path_or_url
            self.download_over_http()

    def download_raw_files(self ,df , path_or_url):
        '''
        Downloads RAW plus.
        :param df: row corresponding to a dataset
        :param path_or_url:
        :return:
        '''
        os.chdir(self.parent_folder)
        # print("donwload_raw", os.getcwd())
        # print("!!!", path_or_url)
        if not path_or_url.startswith("http"): # works with datasets | jobs
            split_path = path_or_url.split("\\")
            path_or_url = '\\'.join(split_path[:-1])
            self.parse_fileserverpath_to_web_url(path_or_url)
            self.download_over_http()
        else:
            # print(">>", path_or_url)
            split_path = path_or_url.split("/")
            path_or_url = '/'.join(split_path[:-2])
            self.url = path_or_url
            # print("<<", self.url)
            self.download_over_http()

    def download_fasta_param_files(self):
        '''
        Parses job_info object and Downlaod
        1. parameter files
        2. fasta files.
        from DMS.
        :return:
        '''
        self.create_dir(self.started_from + '/' + 'DMS_fasta_param' )
        fasta_file = list(set(self.job_info["OrganismDBName"]))
        param_file = list(set(self.job_info["ParameterFileName"]))
        for file in fasta_file:
            url = "http://gigasax/DMS_Organism_Files/Microbial_Communities/FASTA/" + file
            if self.check_url(url):
                # print("MC downlaod", url)
                self.write_to_disk(url)
            url = "http://gigasax/DMS_FASTA_File_Archive/dynamic/forward/" + file
            if self.check_url(url):
                # print("Forward downlaod", url)
                self.write_to_disk(url)
            else:
                print("Can't find FASTA!")

        for file in param_file:
            url= "http://gigasax/DMS_Parameter_Files/MSGFPlus/" + file
            if self.check_url(url):
                self.write_to_disk(url)
            else:
                print("Failed to grab params file")

    def use_df(self, df):
        '''
        Called for each dataset in the dataFrame!
        :param df: analysis_jobs object.
        :return:
        '''
        self.row=df
        self.parent_folder = self.started_from + '/' + str(self.row['Dataset_ID'])
        self.create_dir(self.parent_folder)
        path_or_url = self.download_msgf_jobs(self.row)
        self.download_masic_jobs(self.row)
        self.download_raw_files(self.row, path_or_url)


    @stats
    def get_files(self):
        '''
        :return:
        '''
        # TODO: Check if .raw files exisits in subdirs!: if not os.path.exists(self.started_from):
        # FIXME: Enabling above stmt, doesn't allow to download datasets when they aren't on the disk!
        os.chdir(self.started_from)
        logger.info(msg="get_files()  @:".format(os.getcwd()))
        self.download_fasta_param_files()
        self.Input.apply( lambda x: self.use_df(x), axis=1)
        logger.info(msg="````````")
        logger.info(msg="Finished downloading data at loc:{}".format(self.started_from))
        logger.info(msg="````````")
        # logger.info(msg="Data already exist at loc:{}".format(self.started_from))


    # def download_over_ftp(self):
    #     '''
    #     TODO : Directly from Proto-X windows file-server
    #     Eg. folder_path :"\\proto-6\QExactHF03\2015_2\MinT_Kans_No_Gly_pool_19_Qexactive_22May15_Arwen_14-12-03\SIC201505251246_Auto1197920"
    #     :return:
    #     '''
    #
    # def download_using_DMS_api(self):
    #     '''
    #     TODO: need to explore!
    #     Source: https://prismwiki.pnl.gov/wiki/DMS_Data_Export#Advanced_DMS_Data_Export
    #     url= https://dms2.pnl.gov/data/ax/tsv/aux_info_categories/aux_info_def/501
    #     :return:
    #     '''
    # def handle_workflow_failure(self):
    #     # TODO:  NMDC-10, add logic to check which files were downloaded before pipeline failure. It should only download the ones which weren't successful and start from there.
    #     #        Using MD5 checksum.
    #     pass