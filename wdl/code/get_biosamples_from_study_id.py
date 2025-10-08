from nmdc_api.schema_api import *
import click
import json
from dataclasses import dataclass
import itertools
from pathlib import Path


def dash_to_underscore(s):
    return s.replace('-', '_')

@click.command()
@click.option('--study_id', prompt='Your study id', help='The study id to get biosample IDs from.')
@click.option('--output_dir', prompt='Output directory', help='The output directory to write the biosample IDs file to.')
def main(study_id, output_dir):

    set = NMDCCollection.biosample_set
    filters = [
        filter_field_value_equals("associated_studies", study_id),
        filter_field_value_matches_all("analysis_type", ["metaproteomics", "metagenomics"])
    ]

    records = get_records(set, filters)

    if records is None or len(records) == 0:
        print(f"No records found for study id '{study_id}'")
        return

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    with open(output_path / f"{dash_to_underscore(study_id)}_biosample_records.json", 'w') as file:
        json.dump(records, file, indent=4)

    with open(output_path / f"{dash_to_underscore(study_id)}_biosample_ids.csv", 'w') as file:
        file.write('id\n')
        for record in records:
            file.write(f"{record['id']}\n")

if __name__ == "__main__":
    main()