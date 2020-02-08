import subprocess
import json
import os

with open("metadata.json") as mdfile:
    metadata = json.load(mdfile)

# Get the affiliations orders
affiliations = [c["affiliation"] for c in metadata["authors"]]
affiliations = [i for sl in affiliations for i in sl]

affiliations_data = {}
counter = 1
for affiliation in affiliations:
    if not (affiliation in affiliations_data.values()):
        affiliations_data[counter] = affiliation
        counter += 1

for author in metadata["authors"]:
    auth_affil = [k for k,v in affiliations_data.items() if v in author["affiliation"]]
    author["affiliation_index"] = auth_affil

metadata["affiliation_name"] = affiliations_data

# Fields to remove

with open("metadata.json", "w") as outfile:
    json.dump(metadata, outfile, sort_keys=True, indent=4, separators=(',', ': '))
    outfile.write('\n')
