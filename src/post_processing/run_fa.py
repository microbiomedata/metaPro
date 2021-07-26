from src.post_processing.ficus_analysis import DataOutputtable

import json
import os
import fnmatch


def register_in_emsl_to_jgi(dataset_id, genome_directory, key, value, emsl_to_jgi_copy):
    locations = emsl_to_jgi_copy[dataset_id]['genome_directory'][genome_directory]
    if key not in locations:
        locations[key] = value
    else:
        print(f'{key} already present in {genome_directory}')

if __name__ == '__main__':

    result_loc = os.path.join('storage', 'results', os.environ.get('STUDY'))
    mapper_file= os.path.join(result_loc, 'emsl_to_jgi.json' )

    with open(mapper_file, 'r+') as json_file:
        emsl_to_jgi = json.load(json_file)
        emsl_to_jgi_copy = emsl_to_jgi

        # run for each dataset
        for dataset_id, values in emsl_to_jgi.items():
            if dataset_id not in ['contaminant_file_loc', 'analysis_activity_file_loc', 'data_objects_file_loc',
                                  'STUDY']:
                # dataset search against a fasta file
                for genome_directory, locations in values['genome_directory'].items():
                    resultant_file = locations['resultant_file_loc']
                    gff_file = locations['gff_file_loc']
                    fasta_txt_file = locations['txt_faa_file_loc']

                    data_obj = DataOutputtable(gff_file, resultant_file, fasta_txt_file, os.environ.get('QVALUE_THRESHOLD'), dataset_id, genome_directory)
                    peptide_report,protein_report,qc_metrics_report= data_obj.gen_reports()

                    save_job_results = os.path.join(result_loc, dataset_id, genome_directory)
                    save_at = os.path.join(save_job_results, 'reports')
                    if not os.path.exists(save_at):
                        os.makedirs(save_at)

                    pep_rpt_file= os.path.join(save_at, f"{dataset_id}_{genome_directory}_Peptide_Report.tsv" )
                    register_in_emsl_to_jgi(dataset_id, genome_directory, 'peptide_report_loc', pep_rpt_file,
                                            emsl_to_jgi_copy)
                    pro_rpt_file =os.path.join(save_at, f"{dataset_id}_{genome_directory}_Protein_Report.tsv" )
                    register_in_emsl_to_jgi(dataset_id, genome_directory, 'protein_report_loc', pro_rpt_file,
                                            emsl_to_jgi_copy)
                    qc_m_rpt_file =os.path.join(save_at, f"{dataset_id}_{genome_directory}_QC_metrics.tsv" )
                    register_in_emsl_to_jgi(dataset_id, genome_directory, 'qc_metric_report_loc', qc_m_rpt_file,
                                            emsl_to_jgi_copy)

                    # write dfs to file.
                    peptide_report.to_csv( pep_rpt_file, sep="\t")
                    protein_report.to_csv(pro_rpt_file, sep="\t")
                    qc_metrics_report.to_csv( qc_m_rpt_file, sep="\t")

        json_file.seek(0)  # move back to BOF.
        json_file.truncate()
        json_file.write(json.dumps(emsl_to_jgi_copy, default=str, indent=4))
    pass