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

def denaturation(dnalist):
    dnastrands = []
    for dna in dnalist:
        dnastrands.append(dna[0])
        dnastrands.append(dna[1])
    return dnastrands

def get_primers():
    return ('GCCCCGATTTCAGCTATGGT', 'TGACGCGCACTACAGTCAAT')

def anneal(strands, primers):
    dnal = []
    for strand in strands:
        for primer in primers:
            primerrc = reverse_complement(primer)
            pstart = strand.find(primerrc)
            if pstart != -1:
                originalstrand = strand[pstart::]
                dnal.append((originalstrand, primer))

    return dnal

# Read nsp3 gene from 2720:8554
nsp3_gene = get_gene(2720, 8554)
cDNA_S = reverse_complement(nsp3_gene)
rccDNA_S = nsp3_gene

dnalist = [(''.join(cDNA_S), ''.join(rccDNA_S))]

primers = get_primers()

strands = denaturation(dnalist)

dnalist = anneal(strands, primers)

for dna in dnalist:
    print()
    print(dna[0])
    print()
    print(dna[1])