# The purpose of this file is to fetch homologous IREDs predicted to stereoselectively reduce
# 1-phenyl-dihydroisoquinoline (1-DHIQ) to 1-(S)-phenyl-tetrahydroisoquinoline (1S-THIQ)

import requests
import json
from Bio.Blast.Applications import NcbiblastpCommandline

# IREDs reported in the literature to reduce 1-DHIQ to 1S-THIQ
IREDs = [
    "WP_013019548.1",
    "MF540819",
    "MF540870",
    "MF540873",
    ]


# Define BLAST parameters
blast_db = "nr"
num_aligns = 200
pident_cutoff = 50.0


    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        return response.text
    else:
        print(accID)
        print("FATAL: Bad eFetch request "+ str(response.status_code))
        return None


    # Input protein sequence. Output cached blast results
def blast(acc: str, num_aligns=num_aligns):

    print("Starting BLAST")

    seq = accID2sequence(acc)

    if seq != None:
            # Must have BLAST+ executables in PATH to run this
        blast_cline = NcbiblastpCommandline(db=blast_db, outfmt="6 sseqid pident qcovs", \
            num_alignments=num_aligns, remote=True)

        results, err = blast_cline(stdin=seq)

        results = results.split("\n")[:-1]
        
        homologs = [{"accession": r.split("|")[1], \
                "identity": r.split("|")[2].split("\t")[1], \
                "coverage": r.split("|")[2].split("\t")[2].strip(), \
                "sequence": accID2sequence(r.split("|")[1])} \
                for r in results ]

            # filter out homologs that don't meet the percent identity cutoff
        homologs = [h for h in homologs if float(h["identity"]) >= pident_cutoff]

        #alignment = json.dumps(homologs, indent=4)
        return homologs

    else:
        return None


if __name__ == "__main__":

    # ired_data = []

    # for ired in IREDs:
    #     alignment = blast(ired)
    #     ired_data.append({"ired_acc": ired, "alignment": alignment})

    # with open("ired_homologs.json", "w+") as f:
    #     f.write(json.dumps(ired_data))

    print(accID2sequence("".join([i for i in accID2sequence("BAE95580.1").split("\n")][1:])))


    #print(blast("WP_013019548.1"))