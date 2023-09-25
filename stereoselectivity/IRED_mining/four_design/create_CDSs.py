from dnachisel import *
import requests


def codon_opt(protein_seq:str):

    # Create a random DNA seq given the protein seq. Append a stop codon.
    protein_dna_seq = reverse_translate(protein_seq+"*")

    # DEFINE THE OPTIMIZATION PROBLEM
    problem = DnaOptimizationProblem(
        sequence=protein_dna_seq,
        constraints=[
            AvoidPattern("BsaI_site"),
            EnforceGCContent(mini=0.35, maxi=0.65, window=50),
            EnforceTranslation(location=(0, len(protein_dna_seq)))
        ],
        objectives=[CodonOptimize(species='e_coli', location=(0, len(protein_dna_seq)))]
    )

    # SOLVE THE CONSTRAINTS, OPTIMIZE WITH RESPECT TO THE OBJECTIVE

    problem.resolve_constraints()
    problem.optimize()

    # GET THE FINAL SEQUENCE (AS STRING OR ANNOTATED BIOPYTHON RECORDS)

    final_sequence = problem.sequence  # string
    return final_sequence



    # Input protein accession ID, output sequence in fasta format
def accID2sequence(accID: str):
    URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi/?db=protein&id="+accID+"&rettype=fasta"
    response = requests.get(URL)
    if response.ok:
        print(accID+" seq fetched.")
        full_seq = "".join(response.text.split("\n")[1:-2])
        return full_seq
    else:
        if response.status_code == 429:
            print("429 for "+str(accID))


S = [
    "WP_240033132.1",
    "AUD39583.1",
    "WP_232847731.1",
    "AUD39529.1",
    "WP_277212655.1",
    "WP_253870842.1",
    "WP_123993384.1",
    "EMF01349.1",
    "WP_013019548.1",
    "WP_123102485.1",
    "NUT93668.1",
    "WP_083466898.1",
    "WP_132480478.1",
    "WP_009080771.1"
]

R = [
    "OMI34305.1",
    "AUD39517.1",
    "WP_230421476.1",
    "GCD40644.1",
    "ATL31837.1",
    "ATW48861.1",
    "EME99395.1",
    "WP_112256886.1",
    "WP_116022473.1",
    "PRX45020.1",
    "WP_220260546.1",
    "WP_240598438.1",
    "WP_037254083.1",
    "WP_111548309.1",
    "WP_030629722.1",
]

U = [
    "GHE75528.1",
    "WP_172381658.1",
    "WP_266533875.1",
    "WP_073486480.1",
    "WP_189878949.1",
    "WP_073484372.1",
    "WP_242906400.1",
    "WP_184965005.1",
    "WP_084468547.1",
    "MCP2166807.1",
    "WP_187234818.1",
    "WP_253669224.1",
    "PFG99428.1",
    "WP_029384055.1",
    "WP_229347914.1",
    "WP_191257941.1",
]


with open("S.txt", "w+") as f:
    out = ""
    for i in S:
        seq = accID2sequence(i)
        CDS = codon_opt(seq)
        out += i + "\n" + CDS + "\n\n"

    f.write(out)


with open("R.txt", "w+") as f:
    out = ""
    for i in R:
        seq = accID2sequence(i)
        CDS = codon_opt(seq)
        out += i + "\n" + CDS + "\n\n"

    f.write(out)


with open("U.txt", "w+") as f:
    out = ""
    for i in U:
        seq = accID2sequence(i)
        CDS = codon_opt(seq)
        out += i + "\n" + CDS + "\n\n"

    f.write(out)