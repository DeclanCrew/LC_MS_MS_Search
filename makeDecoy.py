from Bio import SeqIO
import copy

proteinFile = open("Uniprot-TGondii-Reference.fasta", "rb")
proteins  = SeqIO.parse(proteinFile, "fasta")

def decoy (protein):
    dec = copy.deepcopy(protein)
    dec.id = str("DECOY_"+protein.id)
    dec.name = str("DECOY_"+protein.name)
    dec.description = str("DECOY_"+protein.description)
    dec.seq = i.seq[::-1]
    return dec

output = []
for i in proteins:
    output.append(i)
    output.append(decoy(i))

SeqIO.write(output, "target-decoy-gondii.fasta", "fasta")
