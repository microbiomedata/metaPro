import os
import hashlib
import fnmatch
import json
import pandas as pd
from datetime import datetime, timezone
from pymongo import MongoClient
import csv
from utility.utils import timeit
'''
TODO: planning to split this class calling instances of this class from 
      different components of pipeline.(only Definations will live here!)
      Then we won't require below Global variable.
'''


STORAGE="/Volumes/MSSHARE/Anubhav/storage"
STUDY= "stegen"
RESULT_LOC= os.path.join(STORAGE,"results/set_of_Dataset_IDs", STUDY)

xlsx_file = os.path.join(STORAGE, "data/mappings", "stegen_sample.xlsx")
ROOT_LOC = "/Users/anub229/PycharmProjects/nmdc-proteomics-workflow/benchmark/metadata"
# FASTA_LOC= os.path.join(STORAGE,"fasta", STUDY)
# FASTA_LOC= os.path.join(STORAGE,"fastas/Stegen_IMG_annotations/1781_100336/")
FASTA_LOC= os.path.join(STORAGE,"fastas/test_new/")
CONTAMINANT_FILE= os.path.join(STORAGE,"fastas/Tryp_Pig_Bov.fasta")

PIPELINE_TYPE= "nmdc:MetaProteomicAnalysis"
EXECUTION_RESOURCE = "EMSL"
TYPE= "nmdc:DataObject"
GIT_URL = "https://github.com/microbiomedata/metaPro/releases/tag/1.0.0"

class GenMetadata:
    '''
    Generate metadata for the pipeline!
    '''
    def __init__(self):

        self.dataset_id=None,
        self.genome_directory=None
        self.activity={}
        self.uri = 'mongodb://localhost:27017/'
        self.activity_coll= None
        self.data_obj_coll= None
        self.new=[]
        self.data_object_file=[]
        self.quant_bucket=[]
        self.iso_date_with_microseconds=None

    def write_to_json_file(self, filenames):
        '''

        :param filenames:
        :return:
        '''
        with open(filenames[0], 'w') as fptr1 ,\
            open(filenames[1], 'w') as fptr2:
            fptr1.write('[\n')
            fptr2.write('[\n')
            activity_count= self.activity_coll.count()
            data_obj_count= self.activity_coll.count()
            for index, (doc1, doc2) in enumerate(zip(self.activity_coll.find(),self.data_obj_coll.find()) ):
                doc1["id"] = doc1.pop("_id")
                doc2["id"] = doc2.pop("_id")
                fptr1.write(json.dumps(doc1, indent=4))
                fptr2.write(json.dumps(doc2, indent=4))
                if index != activity_count-1: # skip for last document.
                    fptr1.write(",\n")
                    fptr2.write(',\n')
            fptr1.write('\n]')
            fptr2.write('\n]')

        print("Finished parsing {}=={}".format(activity_count,data_obj_count))

    def make_connection(self, db_name, coll_names):
        '''
        1. Make connection to mongodb database
        2. Make cursors available.

        :param db_name: database name
        :param coll_names: list of collection names: _activity then _data_obj
        :return:
        '''
        client = MongoClient(self.uri)
        # makes db.coll if not exists.
        self.activity_coll = client[db_name][coll_names[0]]
        self.data_obj_coll = client[db_name][coll_names[1]]

    def get_md5(self, file):
        '''
        generates MD5 checksum of a file.

        :param file: filename with absolute path.
        :return: chechsum or empty string.
        '''
        if not file in (None, ""):
            md5 = hashlib.md5(open(file, 'rb').read()).hexdigest()
            return md5
        else:
            return ""

    def grab_fasta_file(self, fasta_file_name):
        '''
        Search for desired fasta file

        :return: absolute path to a fasta.
        '''

        for path,subdirs,files in os.walk(FASTA_LOC):
            for file in files:
                # looking for exact match!
                if fnmatch.fnmatch(file, '.'.join([fasta_file_name,'faa'])):
                    return os.path.join(path, file)
                else:
                    print("{}.faa doesn't exist at {}".format(self.genome_directory, os.path.join(path)))

    def gettime(self):
        '''
        Start and end time of running the pipeline.
        '''
        datetime_now = datetime.now(tz=timezone.utc)
        self.iso_date_with_microseconds = datetime_now.strftime('%Y-%m-%dT%H:%M:%S.%fZ%z')
        return self.iso_date_with_microseconds

    def gen_id(self, nersc_seq_id ):
        '''
        - Generate unique ID foreach dataset in the _activity.json
        '''
        txt = "{}\n{}\n{}\n{}\n".format( str(self.dataset_id),
                                     str(self.genome_directory),
                                     str(nersc_seq_id),
                                     self.iso_date_with_microseconds
                                    )
        return 'nmdc:{}'.format(hashlib.md5(txt.encode('utf-8')).hexdigest())

    def prepare_quant_activity_object(self, file):

        with open(file, encoding='utf-8-sig') as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            header = next(tsvreader)
            PeptideSequence = header[1]
            BestProtein = header[2]
            # all_proteins=header[9] "FullGeneList" aka "all_proteins"
            min_QValue = header[11]
            SpectralCount = header[12]
            sum_MASICAbundance = header[13]

            for line in tsvreader:
                quant_dict = {}
                quant_dict[PeptideSequence] = line[1]
                quant_dict[BestProtein] = 'nmdc:' + line[2].replace(" ", "")
                quant_dict["all_proteins"] = ['nmdc:' + protein.replace(" ", "") for protein in line[9].split(',')]
                quant_dict[min_QValue] = line[11]
                quant_dict[SpectralCount] = line[12]
                quant_dict[sum_MASICAbundance] = line[13]
                self.quant_bucket.append(quant_dict)


    def prepare_file_data_object_(self, file_path, file_name , description):
        '''
        - Makes entry in _emsl_analysis_data_objects.json
        - Contains information about Pipeline's analysis results E.g.
            1. resultant.tsv
            2. data_out_table.tsv
        '''
        checksum= self.get_md5(file_path)
        if checksum:
            file_id = 'nmdc:'+ checksum
            print("{} : {}".format(checksum, file_name))
            data_object={}
            data_object['id']= file_id
            data_object['name'] = self.dataset_id+'_'+self.genome_directory+'_'+file_name
            data_object['description'] = description
            data_object['file_size_bytes'] =os.stat(file_path).st_size
            data_object['type'] =TYPE
            data_object['url']= "https://nmdcdemo.emsl.pnnl.gov/proteomics/results/"+ self.dataset_id +'_'+ self.genome_directory+'_'+file_name
            self.data_object_file.append(data_object)
            return file_id
        else:
            print("Found HASH empty for {}".format(file_name))

    def create_has_output(self):
        '''
        Files:
            MSGFjobs_MASIC_resultant.tsv
            Peptide_Report.tsv
            Protein_Report.tsv
            QC_Metrics.tsv
        Quantification:
            peptide
            #TODO: protein
        '''
        has_output=[]
        # add MSGFjobs_MASIC_resultant
        # TODO: Remove temporary Sams_results dependency
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, "MSGFjobs_MASIC_resultant.tsv"),
                                                            'MSGFjobs_MASIC_resultant.tsv',
                                                            "Aggregation of analysis tools{MSGFplus, MASIC} results"
                                                            ))
        # add Peptide_Report
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, "Sams_results/500088_1781_100336_Peptide_Report.tsv"),
                                                            'Peptide_Report.tsv',
                                                            "Aggregated peptides sequences from MSGF+ search results filtered to ~5% FDR"
                                                            ))
        # add Protein_Report
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, "Sams_results/500088_1781_100336_Protein_Report.tsv"),
                                                            'Protein_Report.tsv',
                                                            "Aggregated protein lists from MSGF+ search results filtered to ~5% FDR"
                                                            ))
        # add QC_Metrics
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, "Sams_results/500088_1781_100336_QC_metrics.tsv"),
                                                            'QC_Metrics.tsv',
                                                            "Overall statistics from MSGF+ search results filtered to ~5% FDR"
                                                            ))

        self.activity["has_output"] = has_output
        pass

    def create_has_input(self):
        '''
        "has_input":
              - created for .RAW file pointer to Bill's JSON
              - fasta_checksum
              - contaminant_checksum
        '''

        emsl_to_jgi = {}
        has_input=[]
        # add .RAW
        ptr_to_raw_in_bills_json=  "emsl:output_"+self.dataset_id
        has_input.append(ptr_to_raw_in_bills_json)
        # add JGI fasta
        fasta = self.grab_fasta_file("Ga0482236_proteins")
        emsl_to_jgi[self.dataset_id]= fasta
        self.new.append(emsl_to_jgi)
        emsl_to_jgi.clear()
        fasta_checksum= self.get_md5(fasta)
        if fasta_checksum:
            print("{} : {}".format(fasta_checksum, fasta))
            has_input.append('nmdc:' + fasta_checksum)
            self.activity["has_input"] = has_input
        else:
            print("Found HASH empty for {}".format(fasta))
        # add EMSL contaminants
        contaminant_checksum= self.get_md5(CONTAMINANT_FILE)
        if contaminant_checksum:
            print("{} : {}".format(contaminant_checksum, CONTAMINANT_FILE))
            has_input.append('nmdc:' + contaminant_checksum)
            self.activity["has_input"] = has_input
        else:
            print("Found HASH empty for {}".format(fasta))

        # add Quantifications

        self.prepare_quant_activity_object( os.path.join(RESULT_LOC,self.dataset_id,  "Sams_results/500088_1781_100336_Peptide_Report.tsv"))
        self.activity["has_peptide_quantifications"] = self.quant_bucket
        pass

    def prepare_activity(self, nersc_seq_id):
        '''
        - Makes entry in _MetaProteomicAnalysis_activity.json
        - Foreach dataset, a pointer:
            from : _MetaProteomicAnalysis_activity.json.has_output.[*_file_id]
            to   : _emsl_analysis_data_objects.json."id"

        '''
        self.gettime()
        self.activity["id"] = self.gen_id(nersc_seq_id)
        self.activity["name"]= ":".join(["Metaproteome",self.dataset_id, self.genome_directory])
        self.activity["was_informed_by"]= ":".join(["emsl", self.dataset_id])
        self.activity["started_at_time"]= self.iso_date_with_microseconds
        self.gettime()
        self.activity["ended_at_time"]= self.iso_date_with_microseconds
        self.activity["type"]=PIPELINE_TYPE
        self.activity["execution_resource"]=EXECUTION_RESOURCE
        self.activity["git_url"]=GIT_URL
        self.create_has_input()
        self.create_has_output()


    def on_each_row(self, row):
        '''
        Runs foreach dataset and make entry in a collection.

        '''
        self.dataset_id = str(row['Dataset ID'])
        nersc_seq_id = str(row['sequencing_project_extid'])
        self.genome_directory = str(row['genome directory'])
        print(">>Prepare activity for datasetID:{} genome_directory:{}".format(self.dataset_id,self.genome_directory))
        if not self.genome_directory in ( "", "missing"):
            # skip empty or missing genome_directory
            self.prepare_activity(nersc_seq_id)
            print(self.activity)
            print('*'*50)
            print(self.data_object_file)
            self.activity_coll.insert_one(self.activity)
            self.data_obj_coll.insert_one(self.data_object)
            self.activity.clear()
            self.data_object.clear()
        else:
            print("genome_directory:{} can't be empty/missing!".format(self.genome_directory))

    @timeit
    def start(self):
        '''
        Beging parsing EMSL-JGI mapper file.

        '''
        df_xlsx = pd.read_excel(xlsx_file)
        print("Parsing file {} :: {}".format(xlsx_file, df_xlsx.shape))
        df_xlsx.apply(lambda row: self.on_each_row(row), axis=1)

if __name__ == "__main__":

    meta_file= GenMetadata()
    db_name= "mp_metadata"
    coll_names= ["{}_MetaProteomicAnalysis_activity".format(STUDY),
                 "{}_emsl_analysis_data_objects".format(STUDY)]
    #setup the cursors
    meta_file.make_connection(db_name, coll_names)

    meta_file.start()

    # database to json dump
    activity= os.path.join(ROOT_LOC , STUDY+'_MetaProteomicAnalysis_activity.json')
    data_obj= os.path.join(ROOT_LOC , STUDY+'_emsl_analysis_data_objects.json')
    meta_file.write_to_json_file([activity, data_obj])
    print(meta_file.new)

    # Pipeline uses this file to process datasets.
    # emsl_to_jgi_file= os.path.join(DATA_LOC, "emsl_to_jgi.json")
    # if not os.path.exists(emsl_to_jgi_file):
    #     with open(emsl_to_jgi_file , 'w' ) as fptr:
    #         fptr.write(json.dumps(meta_file.new, indent=2))

    print("Hit the bottom!")