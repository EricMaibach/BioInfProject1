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

def get_complement(strand):
    strandr = strand.replace('A', 'X')  # replace A by X
    strandr = strandr.replace('T', 'A')  # replace T by A
    strandr = strandr.replace('X', 'T')  # replace X by A
    strandr = strandr.replace('C', 'X')  # replace C by X
    strandr = strandr.replace('G', 'C')
    strandr = strandr.replace('X', 'G')
    return strandr

# param: a single strand 5" to 3" dna
# return: a single strand 5" to 3" dna that is reverse complement to the input strand
def reverse_complement(dna5_3):  # input is a strand from 5" to 3"
    rc_dna5_3 = get_complement(dna5_3)
    rc_dna5_3 = rc_dna5_3[::-1]  # reverse the complementary strand to have a strand from 5" to 3"
    
    return rc_dna5_3

def create_dna(stranda, strandb):
    return list(zip(stranda, strandb))

def denaturation(dna):
    sa = [item[0] for item in dna]
    sb = [item[1] for item in dna]
    return sa, sb

def get_primers():
    return 'TGACGCGCACTACAGTCAAT', 'GCCCCGATTTCAGCTATGGT'

def anneal(strand, primer):
    strandp = [None]*len(strand)
    pstart = ''.join(strand).find(primer)
    pend = pstart + len(primer)
    strandp[pstart:pend] = string_to_list(get_complement(primer))
    return create_dna(strand, strandp)

def string_to_list(str):
    return [char for char in str]


    

# Read nsp3 gene from 2720:8554
nsp3_gene = get_gene(2720, 8554)
cDNA_S = reverse_complement(nsp3_gene)
rccDNA_S = nsp3_gene

DNA_S = create_dna(cDNA_S, rccDNA_S)
# d = 200
# r = 50

stranda, strandb = denaturation(DNA_S)

primera, primerb = get_primers()

DNA_N1 = anneal(stranda, primera)
DNA_N2 = anneal(strandb, primerb)