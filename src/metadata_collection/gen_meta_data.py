import os
import hashlib
import fnmatch
import json
import pandas as pd
from datetime import datetime
from pymongo import MongoClient
from bson import json_util
import csv

import json
from bson import ObjectId


'''
TODO: planning to split this class calling instances of this class from 
      different components of pipeline.(only Definations will live here!)
      Then we won't require below Global variable.
'''


STORAGE="/Volumes/storage/"
STUDY= "stegen"
RESULT_LOC= os.path.join(STORAGE,"results", STUDY)
RESULT_PICK_LOC =os.path.join('/Volumes/Anubhav/storage/',"results/set_of_Dataset_IDs/", STUDY)
xlsx_file = os.path.join('/Volumes/Anubhav/storage/data/mappings', "EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping_2021-01-25.xlsx")
ROOT_LOC = "./benchmark/metadata"
FASTA_LOC= '/Volumes/Anubhav/storage/fastas/stegen_copy/Stegen_IMG_annotations'
CONTAMINANT_FILE= os.path.join('/Volumes/Anubhav/storage/',"fastas/Tryp_Pig_Bov.fasta")

PIPELINE_TYPE= "nmdc:MetaProteomicAnalysis"
EXECUTION_RESOURCE = "EMSL"
TYPE= "nmdc:DataObject"
GIT_URL = "https://github.com/microbiomedata/metaPro/releases/tag/1.0.0"

class GenMetadata:
    '''
    Generate metadata for the pipeline!
    '''
    def __init__(self):

        self.annotation_file_name=None
        self.dataset_id=None,
        self.genome_directory=None
        self.activity={}
        self.data_object = {}
        self.uri = 'mongodb://localhost:27017/'
        self.activity_coll= None
        self.data_obj_coll= None
        self.new=[]
        self.quant_bucket=[]


    def write_to_json_file(self, filenames):
        '''

        :param filenames:
        :return:
        '''
        activity_count = self.activity_coll.count()
        data_obj_count = self.activity_coll.count()
        print(f"_MetaProteomicAnalysis_activity count: {activity_count} \n _emsl_analysis_data_objects count: {data_obj_count}")

        with open(filenames[0], 'w') as fptr1 , open(filenames[1], 'w') as fptr2:

            fptr1.write('[\n')
            for doc1 in self.activity_coll.find():
                doc1["id"] = doc1.pop("_id")
                fptr1.write( json.dumps(doc1, default=str, indent=4) )
                fptr1.write(",\n")
            fptr1.write('\n]')

            fptr2.write('[\n')
            for doc2 in self.data_obj_coll.find():
                doc2["id"] = doc2.pop("_id")
                fptr2.write( json.dumps(doc2, default=str, indent=4))
                fptr2.write(',\n')
            fptr2.write('\n]')



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

    def grab_fasta_file(self):
        '''
        Search for desired fasta file

        :return: absolute path to a fasta.
        '''
        fasta_file=''
        for path,subdirs,files in os.walk( os.path.join(FASTA_LOC, self.genome_directory) ):
            for file in files:
                # looking for exact match!
                if fnmatch.fnmatch(file, self.annotation_file_name ):
                    fasta_file= os.path.join(path, file)
                    return fasta_file
            if not fasta_file:
                print("{}.faa doesn't exist at {}".format(self.genome_directory, os.path.join(path)))


    def gettime(self):
        '''
        Start and end time of running the pipeline.
        '''
        return datetime.today().strftime('%Y-%m-%d')

    def gen_id(self, nersc_seq_id ):
        '''
        - Generate unique ID foreach dataset in the _activity.json
        '''
        txt = "{}\n{}\n{}\n".format( str(self.dataset_id),
                                     str(self.genome_directory),
                                     str(nersc_seq_id))
        return 'nmdc:{}'.format(hashlib.md5(txt.encode('utf-8')).hexdigest())

    def prepare_quant_activity_object(self, file):

        with open(file, encoding='utf-8-sig') as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            header = next(tsvreader)
            #reading the header names
            PeptideSequence = header[2]
            BestProtein = header[3]
            # all_proteins=header[10] "FullGeneList" aka "all_proteins"
            # min_QValue = header[12]
            # SpectralCount = header[13]
            # sum_MASICAbundance = header[14]

            for line in tsvreader:
                quant_dict = {}
                quant_dict[PeptideSequence] = line[2]
                quant_dict[BestProtein] = 'nmdc:' + line[3].replace(" ", "")
                quant_dict["all_proteins"] = ['nmdc:' + protein.replace(" ", "") for protein in line[10].split(',')]
                quant_dict['min_QValue'] = line[12]
                quant_dict['SpectralCount'] = line[13]
                quant_dict['sum_MASICAbundance'] = line[14]
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

            self.data_object['id']= file_id
            self.data_object['name'] = file_name
            self.data_object['description'] = description
            self.data_object['file_size_bytes'] =os.stat(file_path).st_size
            self.data_object['type'] =TYPE
            self.data_object['url']= "https://nmdcdemo.emsl.pnnl.gov/proteomics/"+file_name
            # self.data_object_file.append(self.data_object)
            self.data_obj_coll.insert_one(self.data_object)
            self.data_object.clear()
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
        resultant_filename= f"{self.dataset_id}_{self.genome_directory}_MSGFjobs_MASIC_resultant.tsv"
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_PICK_LOC, self.dataset_id, resultant_filename),
                                                            resultant_filename,
                                                            "Aggregation of analysis tools{MSGFplus, MASIC} results"
                                                            ))
        # add Peptide_Report
        peptide_report_filename=f"{self.dataset_id}_{self.genome_directory}_Peptide_Report.tsv"
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, peptide_report_filename),
                                                            peptide_report_filename,
                                                            "Aggregated peptides sequences from MSGF+ search results filtered to ~5% FDR"
                                                            ))
        # add Protein_Report
        protein_report_filename =f"{self.dataset_id}_{self.genome_directory}_Protein_Report.tsv"
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, protein_report_filename),
                                                            protein_report_filename,
                                                            "Aggregated protein lists from MSGF+ search results filtered to ~5% FDR"
                                                            ))
        # add QC_Metrics
        qc_metric_filename =f"{self.dataset_id}_{self.genome_directory}_QC_Metrics.tsv"
        has_output.append(self.prepare_file_data_object_(   os.path.join(RESULT_LOC, self.dataset_id, qc_metric_filename),
                                                            qc_metric_filename,
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
        fasta = self.grab_fasta_file()
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

        self.prepare_quant_activity_object( os.path.join(RESULT_LOC,self.dataset_id,  f"{self.dataset_id}_{self.genome_directory}_Peptide_Report.tsv"))
        self.activity["has_peptide_quantifications"] = self.quant_bucket
        pass

    def prepare_activity(self, nersc_seq_id):
        '''
        - Makes entry in _MetaProteomicAnalysis_activity.json
        - Foreach dataset, a pointer:
            from : _MetaProteomicAnalysis_activity.json.has_output.[*_file_id]
            to   : _emsl_analysis_data_objects.json."id"

        '''
        self.activity["id"] = self.gen_id(nersc_seq_id)
        self.activity["name"]= ":".join(["Metaproteome",self.dataset_id, self.genome_directory])
        self.activity["was_informed_by"]= ":".join(["emsl", self.dataset_id])
        self.activity["started_at_time"]= self.gettime()
        self.activity["ended_at_time"]= self.gettime()
        self.activity["type"]=PIPELINE_TYPE
        self.activity["execution_resource"]=EXECUTION_RESOURCE
        self.activity["git_url"]=GIT_URL

        self.create_has_input()
        self.create_has_output()

        # add to database!
        self.activity_coll.insert_one(self.activity)
        self.activity.clear()
        pass


    def on_each_row(self, row):
        '''
        Runs foreach dataset and make entry in a collection.
        '''

        self.dataset_id = str(row['Dataset ID'])
        nersc_seq_id = str(row['sequencing_project_extid'])
        self.genome_directory = str(row['genome directory'])
        self.annotation_file_name= str(row['annotation file name'])

        print(">>Prepare activity for datasetID:{} genome_directory:{}".format(self.dataset_id,self.genome_directory))

        if not self.genome_directory in ( "", "missing"):
            # skip empty or missing genome_directory
            self.prepare_activity(nersc_seq_id)
            # print(self.activity)
            print('*'*50)
            # print(self.data_object_file)
        else:
            print("genome_directory:{} can't be empty/missing!".format(self.genome_directory))

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
    coll_names= [f"{STUDY}_MetaProteomicAnalysis_activity",
                 f"{STUDY}_emsl_analysis_data_objects"]

    # setup the cursors
    meta_file.make_connection(db_name, coll_names)

    # 1. Make collection and populate them.
    # meta_file.start()

    # 2. dump collections in json files
    if not os.path.exists(ROOT_LOC):
        os.makedirs(ROOT_LOC)
    activity = os.path.join(ROOT_LOC, STUDY + '_MetaProteomicAnalysis_activity.json')
    data_obj = os.path.join(ROOT_LOC, STUDY + '_emsl_analysis_data_objects.json')
    meta_file.write_to_json_file([activity, data_obj])

