import os
import random

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

def denaturation(dnalist):
    dnastrands = []
    for dna in dnalist:
        dnastrands.append(dna[0])
        dnastrands.append(dna[1])
    return dnastrands

# primers are from 5' to 3'
def get_primers():
    return ('GACGCGCAGGGAATGGATAA', 'TGTGTGGCCAACCTCTTCTG')

def anneal(strands, primers):
    dnal = []
    for strand in strands:
        for primer in primers:
            primerrc = get_complement(primer)
            pstart = strand.find(primerrc)
            if pstart != -1:
                originalstrand = strand[pstart::]
                dnal.append((originalstrand, primer))
                break

    return dnal

def elongation(dnalist):
    newdnalist = []
    for dna in dnalist:
        originalstrand = dna[0]
        primerstrand = dna[1]
        primerstrand = primerstrand + get_complement(originalstrand[len(primerstrand)::])
        r = int(round((random.random() * r_range)-(r_range/2), 0))
        cutoff = d + r
        primerstrand = primerstrand[:cutoff:]
        originalstrand = originalstrand[:cutoff:]
        newdnalist.append((originalstrand, primerstrand[::-1]))
    return newdnalist
        

d = 200
r_range = 30


# Read nsp3 gene from 2720:8554
nsp3_gene = get_gene(2720, 8554)

cDNA_S = get_complement(nsp3_gene)
rccDNA_S = nsp3_gene[::-1]

dnalist = [(''.join(cDNA_S), ''.join(rccDNA_S))]

primers = get_primers()

count = 0

while (count < 20):
    count = count + 1
    strands = denaturation(dnalist)
    dnalist = anneal(strands, primers)
    dnalist = elongation(dnalist)

size = 0
averageGC = 0
for dna in dnalist:
    size = size + len(dna[0])
    gcContent = dna[0].count('G')
    gcContent = gcContent + dna[0].count('C')
    averageGC = averageGC + (gcContent / len(dna[0]))

print("dnalist length")
print(len(dnalist))

print()

print("Max size of DNA Fragments")
print(len(max(dnalist)[0]))

print()

print("Min size of DNA Fragments")
print(len(min(dnalist)[0]))

print()

print("Avg size of DNA Fragments")
print(size / len(dnalist))

print()

print("Avg GC content")
print(averageGC / len(dnalist))