class protein ():
    pass

def fastaParser(pFile, seq_format):
    if seq_format != "fasta":
        print "Cannot read non-fasta protein records without biopython."
        exit()
        yield -1
    newProt = protein()
    sequence = []
    for line in pFile.readlines():
        if line[0] == ";":
            continue
        if line[0] == ">":
            if hasattr(newProt, "name"):
                newProt.seq = "".join(sequence)
                yield newProt
                newProt = protein()
                sequence = []
            newProt.name = line.split(" ")[0][1:]
        else:
            sequence.append(line.strip("\n "))
    newProt.seq = "".join(sequence)
    yield newProt


            
            
