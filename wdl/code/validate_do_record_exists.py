from linkml_runtime.loaders import json_loader
from nmdc_schema.schema_nmdc import Database
from dataclasses import dataclass, asdict, fields
from typing import List
from pathlib import Path
import click
import csv


@dataclass
class FileRecord:
    name: str = ""
    is_found: bool = False


csv_out: List[FileRecord] = []


@click.command()
@click.option('--metadata', prompt='Your metadata file', help='The metadata file to validate.')
@click.option('--compare-dir', required=False, help='File directory to compare filenames against')
@click.option('--compare-file', required=False, help='A text containing a list of files to compare filenames against')
def main(metadata, compare_dir, compare_file):
    db = json_loader.load(metadata, Database)
    do_filenames = [do.name for do in db.data_object_set]
    filenames = []
    if compare_file:
        with open(compare_file, 'r') as f:
            for line in f:
                file = line.strip()
                filenames.append(file)
    elif compare_dir:
        filenames = [p.name for p in Path(compare_dir).iterdir() if p.is_file()]

    for file in filenames:
        if file in do_filenames:
            file_record = FileRecord(name=file, is_found=True)
        else:
            file_record = FileRecord(name=file, is_found=False)
        csv_out.append(file_record)

    with open(f"do_file_exists_report.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[field.name for field in fields(FileRecord)])
        writer.writeheader()
        for row in csv_out:
            writer.writerow(asdict(row))
        

if __name__ == "__main__":
    main()