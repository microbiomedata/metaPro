from linkml_runtime.loaders import json_loader
from nmdc_schema.nmdc import Database
from dataclasses import dataclass, asdict, fields
from typing import List
from pathlib import Path
import click
import requests
import csv


@dataclass
class FileRecord:
    id: str = ""
    name: str = ""
    is_found: bool = False


csv_out: List[FileRecord] = []


@click.command()
@click.option('--metadata', prompt='Your metadata file', help='The metadata file to validate.')
@click.option('--compare-dir', required=False, help='Option file directory to compare filenames against')
def main(metadata, compare_dir):
    db = json_loader.load(metadata, Database)
    
    for do in db.data_object_set:
        file_record = FileRecord(id=do.id, name=do.name)
        csv_out.append(file_record)
        try:
            print(f"checking for {do.name} at {do.url}")
            response = requests.head(do.url)

            if (response.status_code == 200):
                file_size_bytes = int(response.headers.get("Content-Length", 0))
                print(f"Found {do.name} at {do.url} with size {file_size_bytes} bytes (expected {do.file_size_bytes} bytes)")
                if do.file_size_bytes == file_size_bytes:
                    file_record.is_found = True
        except requests.exceptions.RequestException as e:
            print(f"Error: {e}")
            print(f"Data object URL: {do.url}")

    with open(f"do_url_report.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[field.name for field in fields(FileRecord)])
        writer.writeheader()
        for row in csv_out:
            writer.writerow(asdict(row))
        

if __name__ == "__main__":
    main()