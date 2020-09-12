import os

def GetGene(start, end):
    file = open('SARS_COV2.fasta', 'r')
    file.readline()
    SARS_COV2_genome = file.read()
    file.close()
    SARS_COV2_genome = SARS_COV2_genome.replace('\n','')
    start = start - 1
    gene = SARS_COV2_genome[start:end]
    return gene

# Read nsp3 gene from 2720:8554
nsp3_gene = GetGene(2720, 8554)
print(nsp3_gene)