# The purpose of this file is to fetch homologous IREDs predicted to stereoselectively reduce
# 1-phenyl-dihydroisoquinoline (1-DHIQ) to 1-(S)-phenyl-tetrahydroisoquinoline (1S-THIQ)

import requests
import json
import time
from Bio.Blast.Applications import NcbiblastpCommandline

# IREDs reported in the literature to reduce 1-DHIQ to 1S-THIQ
IREDs = [
    "ATW48861.1",
    "GCD40644.1 ",
    "PFG99428.1 ",
    "PFG98216.1",
    "PFG94490.1",  
    "EMF00717.1",
    "EME99395.1",
    "EMF02685.1",
    "EMF01349.1",
    "WP_073486480.1",
    "AUD39487.1",
    "AUD39493.1",
    "AUD39502.1",
    "AUD39504.1",
    "AUD39505.1",
    "AUD39517.1",

    "WP_013019548.1",
    "AUD39529.1",
    "AUD39580.1",
    "AUD39583.1",
    ]


# Define BLAST parameters
blast_db = "nr"
num_aligns = 500
max_pident_cutoff = 70.0
min_pident_cutoff = 50.0
folder_name = "ired_70_50"


    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    # need to wait for 1 second between each request so I don't get a 429 HTTP error
    time.sleep(1)
    response = requests.get(URL)
    if response.ok:
        print(accID+" seq fetched.")
        return response.text
    else:
        if response.status_code == 429:
            print("429 for "+str(accID)+". Retrying ..."+ str(response.status_code))
            time.sleep(5)
            response = requests.get(URL)
            if response.ok:
                print(accID+" seq fetched.")
                return response.text
            elif response.status_code == 429:
                print("Got ANOTHER 429! What's going on?")


    # Input protein sequence. Output cached blast results
def blast(acc: str, num_aligns=num_aligns):

    print("Starting BLAST for "+str(acc))

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
                # "sequence": accID2sequence(r.split("|")[1])} \
                }
                for r in results ]

            # filter out homologs that don't meet the percent identity cutoff
        homologs = [h for h in homologs if float(h["identity"]) <= max_pident_cutoff]
        homologs = [h for h in homologs if float(h["identity"]) >= min_pident_cutoff]

        homologs_with_seq = []

        # limit to max 50 homologs per protein
        c = 0
        for i in homologs:
            if c < 50:
                i["sequence"] = accID2sequence(i["accession"])
                homologs_with_seq.append(i)
            c += 1

        homologs = homologs_with_seq

        with open(folder_name+"/"+str(acc), "w") as f:
            f.write(json.dumps(homologs))
        print("cached seq results for "+str(acc))

    else:
        print("no data returned for "+str(acc))


# def create_fasta_db(acc_list):

#     fasta = ""
#     for i in acc_list:
#         seq = accID2sequence(i)
#         fasta += seq + "\n"



if __name__ == "__main__":


    # for i in IREDs:
    #     blast(i)

    import os



    # Create fasta for clustering
    # seqs = ""
    # for file_name in os.listdir("ired_70_50"):
    #     f = os.path.join("ired_70_50", file_name)

    #     with open(f, "r") as data:
    #         ireds = json.load(data)

    #     for ired in ireds:
    #         full_seq = "".join(ireds[0]["sequence"].split("\n")[1:-2])
    #         if ired["sequence"] != None and len(full_seq) < 320:
    #             seqs += ired["sequence"]


    # # add in documented ireds
    # for i in IREDs:
    #     seq = accID2sequence(i)
    #     seqs += seq

    # # save all ireds in fasta file
    # with open("all_ireds.fasta", "w+") as out:
    #     out.write(seqs)
    #     print("added all ireds")


###############################################


    # Create dictionary for sorting R vs S
    seqs = []
    for file_name in os.listdir("ired_70_50"):
        f = os.path.join("ired_70_50", file_name)

        with open(f, "r") as data:
            ireds = json.load(data)

        for ired in ireds:
            full_seq = "".join(ireds[0]["sequence"].split("\n")[1:-2])
            if ired["sequence"] != None and len(full_seq) < 320:
                entry = ired
                entry["full_seq"] = full_seq
                seqs.append(entry)


    # add in documented ireds
    for i in IREDs:
        ired_seq = accID2sequence(i)
        full_seq = "".join(ired_seq.split("\n")[1:-2])
        entry = {
            "accession": i,
            "coverage": 100,
            "identity": 100,
            "sequence": ired_seq,
            "full_seq": full_seq,
        }
        seqs.append(entry)

    # save all ireds in fasta file
    with open("all_ireds.json", "w+") as out:
        out.write(json.dumps(seqs))
        print("wrote all ireds to json file")