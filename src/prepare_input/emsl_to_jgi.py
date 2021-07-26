import os
import json
import pandas as pd
import fnmatch
from utility.utils import logger

"""
- Program reads in mapper file provided by user and
- create emsl_to_jgi.json having file location information of user's input data.

Note: workflow uses this file to process datasets.

Choose:
    "STUDY" name such as "stegen" / "hess" / "blanchard"
    and 
    " MAPPER_FILE file: @Sam
        "EMSL48099_JGI1393_Hess_DatasetToMetagenomeMapping.xlsx"
        "EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping.xlsx"
        "EMSL49483_JGI503125_Blanchard_DatasetToMetagenomeMapping.xlsx"
"""

def search_file_loc( dataset_id, dataset_name, genome_directory, nersc_seq_id, emsl_to_jgi):
    '''

    :param dataset_id:
    :param dataset_name:
    :param genome_directory:
    :param emsl_to_jgi:
    :return:
    '''

    data_loc=  os.path.join('storage', 'data')
    fasta_loc = os.path.join('storage', 'fastas')

    raw_file_loc=''
    # search for raw_file.
    for path, subdirs, files in os.walk(data_loc):
        for name in files:
            if fnmatch.fnmatch(name, f'{dataset_name}.raw'):
                # set path
                raw_file_loc = os.path.join(path, name)
    if not raw_file_loc:
        logger.info('rawfile not present', extra={'dataset_id':f'{dataset_id}',
                                                  'dataset_name':f'{dataset_name}'})

    faa_file_loc=''
    gff_file_loc=''
    # search for fasta file.
    for root, dirnames, filenames in os.walk(fasta_loc):
        if genome_directory in os.path.basename(root):
            for filename in filenames:
                if fnmatch.fnmatch(filename, '*_functional_annotation.gff'):
                    # set path
                    gff_file_loc = os.path.join(root, filename)
                if fnmatch.fnmatch(filename, '*_proteins.faa'):
                    # set path
                    faa_file_loc = os.path.join(root, filename)
    if not gff_file_loc:
        logger.info('gff_file not present', extra={'genome_directory':f'{genome_directory}'})
    if not faa_file_loc:
        logger.info('faa_file not present', extra={'genome_directory': f'{genome_directory}'})
    # add to dict.
    if raw_file_loc and faa_file_loc and gff_file_loc:
        if dataset_id not in emsl_to_jgi:
            emsl_to_jgi[dataset_id] = {
                "raw_file_loc": raw_file_loc,
                "dataset_name": dataset_name,
                'genome_directory': { f"{genome_directory}": { "nersc_seq_id": nersc_seq_id,
                                                               "faa_file_loc": faa_file_loc,
                                                               "gff_file_loc": gff_file_loc
                                                               }
                                      }
            }
        else:
            emsl_to_jgi[dataset_id]['genome_directory'][genome_directory]= { "nersc_seq_id": nersc_seq_id,
                                                               "faa_file_loc": faa_file_loc,
                                                               "gff_file_loc": gff_file_loc
                                                               }

def on_each_row( row, emsl_to_jgi):
    '''

    :param row:
    :param emsl_to_jgi:
    :return:
    '''
    dataset_id = str(row['Dataset ID'])
    dataset_name = str(row['Dataset Name'])
    genome_directory = str(row['genome directory'])
    nersc_seq_id = str(row['sequencing_project_extid'])

    if not genome_directory in ("", "missing"):
        # skip empty or missing genome_directory
        search_file_loc(dataset_id, dataset_name, genome_directory, nersc_seq_id, emsl_to_jgi)
    else:
        logger.info('Workflow will not run for missing faa files from JGI', extra={'dataset_id': f'{dataset_id}',
                                                                                   'genome_directory':f'{genome_directory}'})

def create_mapper(xlsx_file):
    '''
    Beging parsing EMSL-JGI mapper file.
    :return:
    '''
    emsl_to_jgi={}
    df_xlsx = pd.read_excel(xlsx_file)
    df_xlsx.apply(lambda row: on_each_row(row, emsl_to_jgi), axis=1)
    return emsl_to_jgi

def write_to_json(file, object_to_write):
    '''

    :param file:
    :param object_to_write:
    :return:
    '''
    with open( file , 'w') as fptr:
        fptr.write(json.dumps(object_to_write, default=str, indent=4))
    logger.info('Workflow driver file', extra={'emsl_to_jgi_file': file})

if __name__ == "__main__":

    os.environ['STUDY']= 'stegen'
    os.environ['MAPPING_FILENAME']= 'EMSL48473_JGI1781_Stegen_DatasetToMetagenomeMapping_2021-01-25.xlsx'
    #----

    mapper_file=os.path.join('storage','mappings', os.environ.get('MAPPING_FILENAME'))
    if os.path.isfile(mapper_file):
        # parse mapping file.
        mapper = create_mapper(mapper_file)
        # add contaminant file loc
        contaminant_file_loc = os.path.join('storage', 'parameters', os.environ.get('CONTAMINANT_FILENAME'))
        mapper['contaminant_file_loc']= contaminant_file_loc
        # add study info.
        mapper['STUDY'] = os.environ['STUDY']
        # dump it.
        results_loc= os.path.join('storage', 'results', os.environ.get('STUDY'))
        if not os.path.exists(results_loc):
            os.makedirs(results_loc)
        write_to_json( os.path.join(results_loc,'emsl_to_jgi.json'), mapper )
    else:
        logger.info('Mapping file not found.', extra={"MAPPING_FILENAME": os.environ.get('MAPPING_FILENAME')})


