from dnachisel import *
import pandas as pd
import requests
import time

def fetch_protein_seq(accession: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accession+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        seq = "".join(i for i in response.text.split("\n")[1:])
        return seq
    else:
        print("Bad eFetch request "+ str(response.status_code))
        print(accession)
        return None




def codon_opt(protein_seq:str):

    # Create a random DNA seq given the protein seq. Append a stop codon.
    protein_dna_seq = reverse_translate(protein_seq+"*")

    # DEFINE THE OPTIMIZATION PROBLEM
    problem = DnaOptimizationProblem(
        sequence=protein_dna_seq,
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceGCContent(mini=0.4, maxi=0.6, window=50),
            EnforceTranslation(location=(0, len(protein_dna_seq)))
        ],
        objectives=[CodonOptimize(species='e_coli', location=(0, len(protein_dna_seq)))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

    final_sequence = problem.sequence  # string
    #print(final_sequence)
    return final_sequence



def generate_CDSs(file_name):
    # file_name = "MDRR_48.xlsx"

    with open(file_name, "rb+") as f:
        df = pd.read_excel(f)

    refseqs = df.loc[:,"Accession"].values

    for i in range(0, len(refseqs)):
        
        if pd.isna(df.loc[i,"CDS"]):

                # Added a 1 second delay between eFetch requests to prevent HTTP 502 errors
            time.sleep(1)
            acc = df.loc[i, "Accession"]
            prot_seq = fetch_protein_seq(acc)
            if prot_seq != None:
                CDS = codon_opt(prot_seq)
            else:
                CDS = None
            print("CDS generated for "+str(acc))
            df.loc[i, "CDS"] = CDS

    df.to_excel(file_name)




if __name__ == "__main__":
    
    seq = '''
        MARKTKQQALETRQHILDVALRLFSQQGVSATSLAEIANAAGVTRGAIYWHFKNKSDLFSEIWELSESNIGELEIEYQAKFPDDPLSVLREILVHILEATVTEERRRLLMEIIFHKCEFVGEMVVVQQAQRSLCLESYDRIEQTLKHCINAKMLPENLLTRRAAILMRSFISGLMENWLFAPQSFDLKKEARAYVTILLEMYQLCPTLRASTVNGSP    '''
    seq = seq.replace(" ","")
    seq = seq.replace('\n', '')
    seq = seq.upper()

    CDS = codon_opt(seq)

    print(CDS)

