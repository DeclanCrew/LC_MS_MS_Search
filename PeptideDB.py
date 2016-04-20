import sys, getopt
import re
import itertools
import csv
from Bio import SeqIO

#Reads amino acid masses from tab-delimited file to dictionary, for use as reference
def readAAMasses (iFile):
    iData = open(iFile, "rb").readlines()
    oData = {}
    for line in iData:
        entry = line.split("\t")
        oData[str(entry[0])] = float(str(entry[1]).rstrip("\n\r"))
    return oData

#Returns mass of AA string using dictionary of AA mass values, with optional adjustment
def returnPeptideMass (peptide, pKey, adjustment=0):
    totalMass = adjustment
    for i in peptide:
        totalMass += pKey[i]
    return totalMass

#Bins ion mass predictions to the closest 0.5, to enable quicker searching later
def coarsen (invalue, tol=0.5):
    return int(round(invalue/tol, 0)*tol*10)

#Returns list of ion masses
def returnIons (peptide, pKey, adjustment=0):
    output = []
    mass = adjustment
    for i in peptide:
        mass += pKey[i]
        output.append(coarsen(mass))
    return output

def postTranslationallyMod (peptide, ptms):
    return

#Returns peptideDictionary
def returnPeptideDict (peptide, proteins, pKey):
    output = {"peptide":peptide, "proteins":proteins}
    output["mass"]= returnPeptideMass(peptide, pKey, float(18.009467553))
    output["bIons"]= returnIons(peptide, pKey, float(1.00727))
    output["yIons"]= returnIons(peptide[::-1], pKey, float(19.02257))
    return output

#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, MMC, decoy=False):
    output = []
    if decoy == True:
        mProtein = str("n%sc" %(protein.seq[::-1]))
        peptides = (re.sub(r"(?<=[KR])(?=[^P])",'\n', str(mProtein))).split()
        for x in xrange(MMC+1):
            mid = ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
            for i in mid:
                for j in postTransMods(i):
                    output.append((j, str("DECOY_" + protein.name)))
    mProtein = str("n%sc" %(protein.seq))
    peptides = (re.sub(r"(?<=[KR])(?=[^P])",'\n', str(mProtein))).split()
    for x in xrange(MMC+1):
        mid = ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
        for i in mid:
            for j in postTransMods(i):
                output.append((j, protein.name))
    return output

def postTransMods (peptide, regex=0, nRegex=0, cRegex=0):
    output = []
    if peptide[0] == "n":
        peptide = peptide[1:]
    if peptide[len(peptide)-1] == "c":
        peptide = peptide[:(len(peptide)-1)]
   # mods = [(m.start(), int(16)) for m in re.finditer("M",peptide)]
    #modPerm =[]
    #for i in xrange(2):
    #    modPerm += itertools.combinations(mods, i+1)
    #for m in modPerm:
    #    output.append((peptide, m))]
    if len(peptide) > 3:
        output.append(peptide)
    return output

def returnPeptides (protFile, massRef, mmc, useDecoy):
    proteinFile = open(protFile, "rb")
    proteins  = SeqIO.parse(proteinFile, "fasta")
    AAMassKey = readAAMasses(massRef)
    peptideDB = {}
    print ("[+]Digesting proteins in Silico")
    for iProtein in proteins:
        if "X" in iProtein.seq:
            continue
        #print ("Reading %s." % (iProtein.name))
        peptides = iSilSpectra(iProtein, mmc, useDecoy)
        for i in peptides:
            if "B" in i[0] or "Z" in i[0]:
                continue
            if i[0] in peptideDB:
                #print "Appending to match: "+i[1]
                peptideDB[i[0]].append(i[1])
            else:
                #print "Found new peptide: "+i[0]
                peptideDB[i[0]] = [i[1]]
    print "[+]Generating peptide spectra"
    peptideList = [returnPeptideDict(key, peptideDB[key], AAMassKey) for key in peptideDB]
    print "[+]Sorting peptides"
    outHash = {}
    for i in sorted(peptideList, key=lambda entry: entry["mass"]):
        if int(i["mass"]) in outHash:
            outHash[int(i["mass"])].append(i)
        else:
            outHash[int(i["mass"])] = [i]
    return outHash
                    

