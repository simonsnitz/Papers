# The purpose of this script is to filter IRED sequences by amino acid motifs.
# For example, one motif is found in proteins that bind to NADPH

import json

# with open("ired_homologs.json", "r") as f:
#     ireds = json.load(f)

# Protein set as reference for numbering
    #: RefSeq - BAE95580.1
    #: Uniprot - Q1EQE0
    #: Species - Streptomyces kanamyceticus
ref = "MPDNPSTKGRMMRNQQAEHTPVTVIGLGLMGQALAGAFLGAGHPTTVWNRTAAKAEPLVARGAKSAGSVAEAVAASPLVVVCVSDYDAVHALLDPLDGTALQGRTLVNLTSGTSAQARERAAWADGRGADYLDGAILAGPAAIGTADAVVLLSGPRSAFDPHASALGGLGAGTTYLGADHGLASLYDAAGLVMMWSILNGFLQGAALLGTAGVDATTFAPFITQGIGTVADWLPGYARQIDDGAYPADDAAIDTHLATMEHLIHESEFLGVNAELPRFIKALADRAVADGHGGSGYPALIEQFRTHSGK"


# # get the first sequence among the homologs for testing
# test = ireds[0]["alignment"][0]["sequence"]
# test = "".join([i for i in test.split("\n")][1:])

# ired_seqs = []

# print(ireds[3])
# #print(ireds[0]["alignment"][0])
# for i in ireds[0]["alignment"]:
#     try:
#         seq = i["sequence"]
#         seq = "".join([k for k in seq.split("\n")][1:])
#         ired_seqs.append(seq)
#     except:
#         pass




from io import StringIO
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO



def blast_align(sequence1, sequence2):

    # Create two sequence files
    seq1 = SeqRecord(Seq(sequence1),
                    id="seq1")
    seq2 = SeqRecord(Seq(sequence2),
                    id="seq2")
    SeqIO.write(seq1, "seq1.fasta", "fasta")
    SeqIO.write(seq2, "seq2.fasta", "fasta")

    # Run BLAST and parse the output as XML
    output = NcbiblastpCommandline(query="seq1.fasta", subject="seq2.fasta", outfmt=5, word_size=4)()[0]
    
    #print(output)
    blast_result_record = NCBIXML.read(StringIO(output))



    # Print some information on the result
    data = []
    for alignment in blast_result_record.alignments:
        for hsp in alignment.hsps:

            identity = round((hsp.identities/hsp.align_length)*100, 2)
            start_pos = hsp.query_start
            ref_seq = hsp.query
            ired_seq = hsp.sbjct

            coverage = round((hsp.align_length/blast_result_record.query_length)*100, 2)

            entry = {"length": alignment.length,
                    "e value": hsp.expect,
                    "identity": identity,
                    "coverage": coverage}
            data.append(entry)


    aligned_seq = {}
    for i in range(0, len(ref_seq)):
        ref_a = ref_seq[i]
        ired_a = ired_seq[i]
        if ref_a != "-":
            aligned_seq[str(start_pos)] = ired_a
            start_pos += 1
    
    #print(aligned_seq)

    key_residues = {
        "50": aligned_seq["50"],
        "108": aligned_seq["108"],
        "136": aligned_seq["136"],
        "137": aligned_seq["137"],
        "139": aligned_seq["139"],
        "187": aligned_seq["187"],
        "191": aligned_seq["191"],
        "194": aligned_seq["194"],
        "228": aligned_seq["228"],
        "255": aligned_seq["255"],
        "258": aligned_seq["258"]
    }

    return key_residues


    # print(aligned_seq["187"])   
    #     # Should be Y for S-specificity, and D for R-specificity
    # print(aligned_seq["191"])
    #     # Should NOT be N (indicator that it's actually a hydrolase)
    # print(aligned_seq["50"])
    #     # Should be R, a conserved residue
    # print(aligned_seq["108"])
    #     # Should be N, usually conserved
    # print(aligned_seq["139"])
    #     # Should be P, usually conserved
    # print(aligned_seq["194"])
    #     # Should be F, which is conserved for S-specific IREDs

    '''
    Important residues based on structure:
    191 (included), 194 (included), 187 (included)
    136, 137, 258, 255, 228
    '''

#print(ired_seqs)
#for i in ired_seqs:
#    try:
#        blast_align(ref, i)
#    except: 
#        print("Failed")
#blast_align(ref, test)


if __name__ == "__main__":

    with open("all_ireds.json", "rb") as f:
        data = json.load(f)

    new_data = []

    num_ireds = str(len(data))
    c = 0
    f = 0
    for i in data:
        test = i["full_seq"]
        try:
            key_residues = blast_align(ref, test)
            i["key_residues"] = key_residues
            new_data.append(i)
            print("added ired: "+str(c)+" of "+num_ireds)
        except: 
            print("Failed")
            f += 1
        c += 1

    print("number of failed entries: "+str(f))

    with open("labeled_ireds.json", "w") as f:
        f.write(json.dumps(new_data))
