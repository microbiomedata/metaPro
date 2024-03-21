import os
import json
import csv
import sys
import hashlib

from urllib.parse import urljoin
from datetime import datetime
from linkml_runtime.dumpers import json_dumper
from nmdc_schema.nmdc import Database, MetaproteomicsAnalysisActivity, DataObject, ProteinQuantification, PeptideQuantification, FileTypeEnum
from nmdc_schema.nmdc_data import get_nmdc_jsonschema_string
from nmdc_id_source import NmdcIdSource, NmdcFakeIdSource
from typing import List


WORKFLOW_METADATA_VERSION = 1

class GenMetadata:
    """
    Generate metadata for the pipeline!
    """
    def __init__(self, analysis_type: str, execution_resource: str, git_url: str, results_url: str,
                id_source: NmdcIdSource,
                activity_id: str,
                contaminate_file: str,
                masic_param_file: str,
                msgfplus_param_file: str,
                masic_param_id: str,
                msgf_param_id: str,
                contam_id: str,
                version: str):

        self.execution_resource = execution_resource
        self.git_url = git_url
        self.type = analysis_type
        self.results_url = results_url
        self.id_source = id_source
        self.activity_id = activity_id
        self.contaminant_file = contaminate_file
        self.masic_param_file = masic_param_file
        self.msgfplus_param_file = msgfplus_param_file
        self.masic_param_id = masic_param_id
        self.msgf_param_id = msgf_param_id
        self.contam_id = contam_id
        self.version = version

        self.dataset_id = None
        self.genome_directory = None
        self.annotation_file_name = None
        self.resultant_file = None
        self.fasta_file = None
        self.peptide_report = None
        self.protein_report = None
        self.qc_metric_report = None
        self.started_at_time = None
        self.ended_at_time = None
        self.fasta_id = None
        self.gff_id = None

    def get_md5(self, file):
        """
        generates MD5 checksum of a file.

        :param file: filename with absolute path.
        :return: chechsum or empty string.
        """
        if not file in (None, ""):
            md5 = hashlib.md5(open(file, "rb").read()).hexdigest()
            return md5
        else:
            return ""

    def get_ProteinQuantification(self):
        protein_arr = []

        with open(self.protein_report, encoding="utf-8-sig") as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            header = next(tsvreader)

            for line in tsvreader:
                protein_quantification = ProteinQuantification(
                    peptide_sequence_count = line[1],
                    best_protein = line[2].replace(" ", ""),
                    all_proteins = [protein.replace(" ", "") for protein in line[9].split(",")],
                    protein_spectral_count = line[12],
                    protein_sum_masic_abundance = line[13],
                )

                protein_arr.append(protein_quantification)

    def get_PeptideQuantification(self):
        peptide_quantifications_arr = []

        with open(self.peptide_report, encoding="utf-8-sig") as tsvfile:
            tsvreader = csv.reader(tsvfile, delimiter="\t")
            header = next(tsvreader)
            for line in tsvreader:
                peptide_quantification = PeptideQuantification(
                peptide_sequence = line[1],
                best_protein = line[2].replace(" ", ""),
                all_proteins = [protein.replace(" ", "") for protein in line[9].split(",")],
                min_q_value = line[11],
                peptide_spectral_count = line[12],
                peptide_sum_masic_abundance = int(float(line[13]))
                )

                peptide_quantifications_arr.append(peptide_quantification)
        
        return peptide_quantifications_arr

    def get_file_data_object(self, file_path: str, file_name: str, description: str, activity_id: str, type: str) -> DataObject:
        """
        Populates and returns a DataObject object for a metaproteomics pipeline results file.
        """

        file_id = self.id_source.get_ids("nmdc:DataObject", 1)[0]
        
        data_object = DataObject(
            id=file_id,
            name=file_name,
            description=description,
            file_size_bytes=os.stat(file_path).st_size,
            md5_checksum=hashlib.md5(open(file_path,'rb').read()).hexdigest(),
            was_generated_by=None,
            data_object_type=type,
            url=urljoin(self.results_url, file_name),
        )

        return data_object

    def get_data_objects(self) -> List[DataObject]:
        data_objects = []
        data_objects.append(
            self.get_file_data_object(
                self.resultant_file,
                os.path.basename(self.resultant_file),
                "Aggregation of analysis tools {MSGFplus, MASIC} results", 
                self.activity_id,
                "Unfiltered Metaproteomics Results"
            )
        )

        data_objects.append(
            self.get_file_data_object(
                self.peptide_report,
                os.path.basename(self.peptide_report),
                "Aggregated peptides sequences from MSGF+ search results filtered to ~5% FDR", 
                self.activity_id,
                "Peptide Report"
            )
        )

        data_objects.append(
            self.get_file_data_object(
                self.protein_report,
                os.path.basename(self.protein_report),
                "Aggregated protein lists from MSGF+ search results filtered to ~5% FDR", 
                self.activity_id,
                "Protein Report"
            )
        )

        data_objects.append(
            self.get_file_data_object(
                self.qc_metric_report,
                os.path.basename(self.qc_metric_report),
                "Overall statistics from MSGF+ search results filtered to ~5% FDR", 
                self.activity_id,
                "Metaproteomics Workflow Statistics"
            )
        )

        return data_objects

    def get_has_input(self) -> List[str]:
        '''
            Returns the DataObject IDs for pipeline input files
        '''

        has_input = [
            self.dataset_id,
            self.gff_id,
            self.fasta_id,
            self.msgf_param_id,
            self.masic_param_id,
            self.contam_id
        ]
        
        return has_input

    def get_metaproteomics_analysis_activity(self, data_objects: list[DataObject]) -> MetaproteomicsAnalysisActivity:
        has_input_arr = self.get_has_input()
        has_output_arr = [data_object.id for data_object in data_objects]

        mp_analysis_activity_obj = MetaproteomicsAnalysisActivity(
            id=self.activity_id,
            execution_resource=self.execution_resource,
            git_url=self.git_url,
            version=self.version,
            name=f"Metaproteomics Analysis Activity for {self.activity_id}",
            was_informed_by=self.genome_directory,
            type=self.type,
            has_output=has_output_arr,
            has_input=has_input_arr,
            started_at_time=self.started_at_time,
            ended_at_time=self.ended_at_time
            )

        mp_analysis_activity_obj.has_peptide_quantifications = self.get_PeptideQuantification()

        return mp_analysis_activity_obj
    
    def set_keys(
        self,
        dataset_id,
        genome_directory,
        annotation_file_name,
        resultant_file,
        fasta_file,
        peptide_report,
        protein_report,
        qc_metric_report,
        started_at_time,
        ended_at_time,
        gff_id,
        fasta_id,
    ):
        self.dataset_id = dataset_id
        self.genome_directory = genome_directory
        self.annotation_file_name = annotation_file_name
        self.resultant_file = resultant_file
        self.fasta_file = fasta_file
        self.peptide_report = peptide_report
        self.protein_report = protein_report
        self.qc_metric_report = qc_metric_report
        self.started_at_time = started_at_time
        self.ended_at_time = ended_at_time
        self.gff_id = gff_id
        self.fasta_id = fasta_id


    @staticmethod
    def create(execution_resource: str, git_url: str, results_url: str,
                  id_source: NmdcIdSource,
                  contaminate_file: str,
                  masic_param_file: str,
                  msgfplus_param_file: str,
                  masic_param_id: str,
                  msgf_param_id: str,
                  contam_id: str,
                  version: str) -> 'GenMetadata':
        '''
        Create a GenMetadata object with newly minted activity ID with OOP in mind.
        '''
        analysis_type = "nmdc:MetaproteomicsAnalysisActivity"
        activity_id = id_source.get_ids(analysis_type, 1)[0]
        activity_id = activity_id + "." + str(WORKFLOW_METADATA_VERSION)

        return GenMetadata(analysis_type, execution_resource, git_url, results_url, id_source, activity_id,
                           contaminate_file,
                           masic_param_file,
                           msgfplus_param_file,
                           masic_param_id,
                           msgf_param_id,
                           contam_id,
                           version)


def underscore_to_colon(curie: str):
    return curie.replace('_', ':')


if __name__ == "__main__":
    mapper_file = sys.argv[1]
    study = sys.argv[2]
    execution_resource = sys.argv[3]
    git_url = sys.argv[4]
    results_url = sys.argv[5]
    contaminate_file = sys.argv[6]
    masic_param_file = sys.argv[7]
    msgfplus_param_file = sys.argv[8]
    masic_param_id = sys.argv[9]
    msgf_param_id = sys.argv[10]
    contam_id = sys.argv[11]
    version = sys.argv[12]

    if results_url[-1] != '/':
        results_url = results_url + '/'

    data_objects_arr: List[DataObject] = []
    metaproteomics_analysis_activity_arr: List[MetaproteomicsAnalysisActivity] = []

    # Minting source
    id_source = NmdcIdSource(os.environ['CLIENT_ID'], os.environ['CLIENT_SECRET'])

    # 1. Make collection and populate them.
    with open(mapper_file, "r+") as json_file:
        mapper_json = json.load(json_file)

        for mapping in mapper_json:

            meta_file = GenMetadata.create(
                execution_resource=execution_resource,
                git_url=git_url,
                results_url=results_url,
                id_source=id_source,
                contaminate_file=contaminate_file,
                masic_param_file=masic_param_file,
                msgfplus_param_file=msgfplus_param_file,
                masic_param_id=underscore_to_colon(masic_param_id),
                msgf_param_id=underscore_to_colon(msgf_param_id),
                contam_id=underscore_to_colon(contam_id),
                version=version)

            meta_file.set_keys(
                underscore_to_colon(mapping["dataset_id"]),
                underscore_to_colon(mapping["genome_directory"]),
                os.path.basename(mapping["txt_faa_file"]),
                mapping["resultant_file"],
                mapping["faa_file"],
                mapping["peptide_report_file"],
                mapping["protein_report_file"],
                mapping["qc_metric_report_file"],
                mapping["start_time"],
                mapping["end_time"],
                underscore_to_colon(mapping["gff_id"]),
                underscore_to_colon(mapping["fasta_id"])
            )

            data_objects = meta_file.get_data_objects()
            metaproteomics_analysis_activity = meta_file.get_metaproteomics_analysis_activity(data_objects)

            metaproteomics_analysis_activity_arr.append(metaproteomics_analysis_activity)
            data_objects_arr.extend(data_objects)
            

        # 2. dump collections in json files
        schema = get_nmdc_jsonschema_string()
        activity_file = f"{study}_MetaProteomicAnalysis_activity.json"
        data_obj_file = f"{study}_analysis_data_objects.json"

        activity_db = Database()
        activity_db.metaproteomics_analysis_activity_set = metaproteomics_analysis_activity_arr

        data_obj_db = Database()
        data_obj_db.data_object_set = data_objects_arr

        # json_dumper.dump(activity_db, contexts=schema, inject_type=False, to_file=activity_file)
        # json_dumper.dump(data_obj_db, contexts=schema, inject_type=False, to_file=data_obj_file)        
        json_dumper.dump(activity_db, inject_type=False, to_file=activity_file)
        json_dumper.dump(data_obj_db, inject_type=False, to_file=data_obj_file)