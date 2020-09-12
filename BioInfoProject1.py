import os

def get_gene(start, end):
    file = open('SARS_COV2.fasta', 'r')
    file.readline()
    SARS_COV2_genome = file.read()
    file.close()
    SARS_COV2_genome = SARS_COV2_genome.replace('\n','')
    start = start - 1
    gene = SARS_COV2_genome[start:end]
    return gene

# param: a single strand 5" to 3" dna
# return: a single strand 5" to 3" dna that is reverse complement to the input strand
def reverse_complement(dna5_3):  # input is a strand from 5" to 3"
    rc_dna5_3 = dna5_3.replace('A', 'X')  # replace A by X
    rc_dna5_3 = rc_dna5_3.replace('T', 'A')  # replace T by A
    rc_dna5_3 = rc_dna5_3.replace('X', 'T')  # replace X by A
    rc_dna5_3 = rc_dna5_3.replace('C', 'X')  # replace C by X
    rc_dna5_3 = rc_dna5_3.replace('G', 'C')
    rc_dna5_3 = rc_dna5_3.replace('X', 'G')
    rc_dna5_3 = rc_dna5_3[::-1]  # reverse the complementary strand to have a strand from 5" to 3"
    return rc_dna5_3

# Read nsp3 gene from 2720:8554
nsp3_gene = get_gene(2720, 8554)
cDNA_S = reverse_complement(nsp3_gene)
rccDNA_S = nsp3_gene

# It is the double strand DNA {cDNA_S, rccDNA_S} that can be amplified.
DNA_S = (cDNA_S, rccDNA_S)


print(DNA_S)