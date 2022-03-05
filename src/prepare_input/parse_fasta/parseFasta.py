# # algorithm:
# 1. read the .fasta file
# 2. create a sha1 hash and keep it in workflow.json.
# 2. take the "description" and "protein_seq" and model in a dict object!
# {
#     "protein_seq" : {
#         "matches" : 3,
#         "protein_names" : ["protein_name1", "protein_name3", "protein_name4"]
#     }
# }
#     2.1 while parsing set the flag self.yes_to_new_fasta=1 if found any entry of "matches">1
# 3. output new fasta file
#     if any of the self.yes_to_new_fasta=1
#     print(Input fasta sha1, and new sha1)
#     else no need!
#     print(Input fasta sha1 == output fasta-sha1)
#
# 4. Add to workflow.json file,
#     {
#         "fasta": {
#                 "found_at": "file_loc_on_server",
#                 "name": "filename.fasta",
#                 "ip_checksum": "garegq4gq34g5",
#                 "out_checksum": "garegq4gq34g5",
#                 "duplicates_found": "count the number of duplicates in file!"
#         }
#     }
# https://stackoverflow.com/questions/66462611/remove-duplicated-sequences-in-fasta-with-python
#
# https://biopython.org/wiki/Download
# https://biopython.org/wiki/Sequence_Cleaner
from Bio import SeqIO
import time

start = time.time()

seen = []
records = []

for record in SeqIO.parse("INPUT-FILE", "fasta"):
    if str(record.seq) not in seen:
        seen.append(str(record.seq))
        records.append(record)


# writing to a fasta file
SeqIO.write(records, "OUTPUT-FILE", "fasta")
end = time.time()

print(f"Run time is {(end- start)/60}")
