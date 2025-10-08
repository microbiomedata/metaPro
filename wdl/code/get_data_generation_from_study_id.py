from nmdc_api.schema_api import *
import click
import json
from collections import defaultdict
from dataclasses import dataclass

from pathlib import Path


def dash_to_underscore(s):
    return s.replace('-', '_')

@click.command()
@click.option('--study_id', prompt='Your study id', help='The study id to get biosample IDs from.')
@click.option('--output_dir', prompt='Output directory', help='The output directory to write the biosample IDs file to.')
def main(study_id, output_dir):

    set = NMDCCollection.data_generation_set
    filters = [
        filter_field_value_equals("associated_studies", study_id)
    ]

    records = get_records(set, filters)

    if records is None or len(records) == 0:
        print(f"No records found in '{set}' for study id '{study_id}'")
        return
    
    grouped_records = defaultdict(list)
    for record in records:
        grouped_records[record["has_input"][0]].append(record)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    with open(output_path / f"{dash_to_underscore(study_id)}_data_generation_records.json", 'w') as file:
        json.dump(grouped_records, file, indent=4)

    with open(output_path / f"{dash_to_underscore(study_id)}_data_generation_ids.csv", 'w') as file:
        file.write('id\n')
        for id in grouped_records.keys():
            file.write(f"{id}\n")

if __name__ == "__main__":
    main()