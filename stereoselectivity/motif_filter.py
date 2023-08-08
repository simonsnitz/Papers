# The purpose of this script is to filter IRED sequences by amino acid motifs.
# For example, one motif is found in proteins that bind to NADPH

import json

with open("ired_homologs.json", "r") as f:
    ireds = json.load(f)

# Protein set as reference for numbering
    #: RefSeq - BAE95580.1
    #: Uniprot - Q1EQE0
    #: Species - Streptomyces kanamyceticus
ref = "MPDNPSTKGRMMRNQQAEHTPVTVIGLGLMGQALAGAFLGAGHPTTVWNRTAAKAEPLVARGAKSAGSVAEAVAASPLVVVCVSDYDAVHALLDPLDGTALQGRTLVNLTSGTSAQARERAAWADGRGADYLDGAILAGPAAIGTADAVVLLSGPRSAFDPHASALGGLGAGTTYLGADHGLASLYDAAGLVMMWSILNGFLQGAALLGTAGVDATTFAPFITQGIGTVADWLPGYARQIDDGAYPADDAAIDTHLATMEHLIHESEFLGVNAELPRFIKALADRAVADGHGGSGYPALIEQFRTHSGK"

# get the first sequence among the homologs for testing
test = ireds[0]["alignment"][0]["sequence"]
test = "".join([i for i in test.split("\n")][1:])


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
            #print(hsp.query)
            #print(hsp.sbjct)
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
    print(aligned_seq["187"])
    print(aligned_seq["191"])
    print(aligned_seq["108"])
    print(aligned_seq["139"])
    print(aligned_seq["194"])

blast_align(ref, test)