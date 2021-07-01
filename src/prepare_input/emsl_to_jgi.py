import os
import json
import pandas as pd
# from utility.utils import stats

"""
Program reads in mapper file and creates "emsl_to_jgi_file".
Pipeline uses this file to process datasets.

Choose:
    "STUDY" name such as "stegen" / "hess" / "blanchard"
    and 
    " MAPPER_FILE file: @Sam
        "EMSL48099_JGI1393_Hess_DatasetToMetagenomeMapping.xlsx"
        "EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping.xlsx"
        "EMSL49483_JGI503125_Blanchard_DatasetToMetagenomeMapping.xlsx"
    and
    "STORAGE" location such as a share drive:
        "/mnt/anub229_msshare/anubhav/storage/"
        "/Volumes/MSSHARE/Anubhav/storage/"
"""

def prepare_activity( dataset_id, genome_directory, emsl_to_jgi):

    # faa_file_loc = os.path.join(FASTA_LOC, STUDY, genome_directory, "annotation")
    faa_file_loc = os.path.join(FASTA_LOC, "stegen01152021")
    faa_file_loc = faa_file_loc + '/' + '.'.join([genome_directory, 'faa'])
    if dataset_id not in emsl_to_jgi:
        emsl_to_jgi[dataset_id] = [faa_file_loc]
    else:
        emsl_to_jgi[dataset_id].append(faa_file_loc)

def on_each_row( row, emsl_to_jgi):
    dataset_id = str(row['Dataset ID'])
    genome_directory = str(row['genome directory'])
    # print(">>Prepare activity for datasetID:{} genome_directory:{}".format(dataset_id, genome_directory))
    if not genome_directory in ("", "missing"):
        # skip empty or missing genome_directory
        prepare_activity(dataset_id, genome_directory, emsl_to_jgi)
    else:
        print("genome_directory:{} can't be empty/missing!".format(genome_directory))

def create_mapper(xlsx_file):
    '''
    Beging parsing EMSL-JGI mapper file.
    :return:
    '''
    emsl_to_jgi={}

    df_xlsx = pd.read_excel(xlsx_file)
    print("Parsing file {} :: {}".format(xlsx_file, df_xlsx.shape))
    df_xlsx.apply(lambda row: on_each_row(row, emsl_to_jgi), axis=1)

    print(emsl_to_jgi)
    emsl_to_jgi_file = os.path.join(DATA_LOC, "emsl_to_jgi_Stegen01212021.json")
    if not os.path.exists(emsl_to_jgi_file):
        with open(emsl_to_jgi_file, 'w') as fptr:
            fptr.write(json.dumps(emsl_to_jgi))
    print("Look at {}".format(emsl_to_jgi_file))
    return emsl_to_jgi_file

STORAGE="/Volumes/MSSHARE/Anubhav/storage/"
STUDY= "stegen"
DATA_LOC= os.path.join(STORAGE,"data/set_of_Dataset_IDs", STUDY)
FASTA_LOC= os.path.join(STORAGE,"fastas")

if __name__ == "__main__":

    xlsx_file = os.path.join(STORAGE, "data/mappings", "EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping.xlsx")
    create_mapper(xlsx_file)

