import sys, getopt
import re
import itertools
import csv
from Bio import SeqIO

def inputOutput (argv):
    options = {"maxMissedCleave":1}
    try:
        opts, args = getopt.getopt(argv, "hi:o:m:pc")
    except getopt.GetoptError:
        print "PeptideDB.py -i <inputfasta> -o <outputPep> -m <AAmassref>"
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print "PeptideDB.py -i <inputfasta> -o <outputPep> -m <AAmassref> -p <PostTranslationalMods> -c <maxMissedCleaves>"
            sys.exit()
        elif opt == "-i":
            options["infile"] = arg
        elif opt == "-o":
            options["outfile"] = arg
        elif opt == "-m":
            options["massRef"] = arg
        elif opt == "-c":
            options["maxMissedCleave"] = arg
    return options

#Reads amino acid masses from tab-delimited file to dictionary, for use as reference
def readAAMasses (iFile):
    iData = open(iFile, "rb").readlines()
    oData = {}
    for line in iData:
        entry = line.split("\t")
        oData[str(entry[0])] = float(str(entry[1]).rstrip("\n\r"))
    return oData

#Returns mass of AA string using dictionary of AA mass values, with optional adjustment
def returnPeptideMass (peptide, pKey, pTMs=[], adjustment=0):
    totalMass = 0
    for i in peptide:
        totalMass += pKey[i]
    for i in pTMs:
        totalMass += i
    totalMass += adjustment
    return totalMass

def returnPeptideDict (peptide, pTMs, protein, pKey):
    output = {"peptide":peptide,"pTMs":pTMs,"protein":protein}
    output["Mass"]= returnPeptideMass(peptide, pKey, [i[1] for i in pTMs], float(18.0153))
    bIons =[]
    yIons =[]
    for i in xrange(len(peptide)):
        bIons.append((peptide[:(i+1)], [j[1] for j in pTMs if j[0] < (i+1)]))
        yIons.append((peptide[i:], [j[1] for j in pTMs if j[0] >= i]))
    output["bIons"]= [returnPeptideMass(bIon[0], pKey, bIon[1], float(1.00727)) for bIon in bIons]
    output["yIons"]= [returnPeptideMass(yIon[0], pKey, yIon[1], float(19.02257)) for yIon in yIons]
    output["maxScore"]= (len(peptide)*7)
    return output

#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, MMC):
    mid = []
    peptides = (re.sub(r"(?<=[KR])(?=[^P])",'\n', str(protein))).split()
    for x in xrange(MMC+1):
        mid += ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
    output = []
    for i in mid:
        [output.append(j) for j in postTransMods(i)]
    return output

def postTransMods (peptide, regex=0, nRegex=0, cRegex=0):
    output = []
    if peptide[0] == "n":
        peptide = peptide[1:]
    if peptide[len(peptide)-1] == "c":
        peptide = peptide[:(len(peptide)-1)]
    mods = [(m.start(), int(16)) for m in re.finditer("M",peptide)]
    modPerm =[]
    for i in xrange(len(mods)):
        modPerm += itertools.combinations(mods, i+1)
    for m in modPerm:
        output.append((peptide, m))
    output.append((peptide, ()))
    return output

def writeToCSV (ofile, oArray):
    keys = set()
    for k in oArray[0]:
        keys.add(str(k))
    out = open(ofile,"wb")
    dict_writer = csv.DictWriter(out, list(keys), dialect="excel-tab")
    dict_writer.writer.writerow(list(keys))
    dict_writer.writerows(oArray)
    return

options = inputOutput(sys.argv[1:])
print options
proteinFile = open(options["infile"], "rb")
proteins = SeqIO.parse(proteinFile, "fasta")
AAMassKey = readAAMasses(options["massRef"])
peptideDB = []
print ("[+]Generating peptide database")
for iProtein in proteins:
    print ("Reading %s." % (iProtein.name))
    mProtein = str("n"+iProtein.seq+"c")
    peptides = iSilSpectra(mProtein, options["maxMissedCleave"])
    for i in peptides:
        peptideDB.append(returnPeptideDict(i[0], i[1], iProtein.name, AAMassKey))

writeToCSV(options["outfile"], peptideDB)
