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

#Returns list of ion masses
def returnIons (peptide, pKey, adjustment=0):
    output = []
    mass = adjustment
    for i in peptide:
        mass += pKey[i]
        output.append(mass)
    return output

#Returns peptideDictionary
def returnPeptideDict (peptide, protein, pKey):
    output = {"peptide":peptide, "proteins":[]}
    output["proteins"].append(protein)
    output["mass"]= returnPeptideMass(peptide, pKey, float(18.009467553))
    output["bIons"]= returnIons(peptide, pKey, float(1.00727))
    output["yIons"]= returnIons(peptide[::-1], pKey, float(19.02257))
    return output

#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, MMC, pKey, decoy=False):
    output = []
    if decoy == True:
        mProtein = str("n%sc" %(protein.seq[::-1]))
        peptides = (re.sub(r"(?<=[KR])(?=[^P])",'\n', str(mProtein))).split()
        for x in xrange(MMC+1):
            mid = ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
            for i in mid:
                for j in postTransMods(i):
                    output.append(returnPeptideDict(j, str("DECOY_" + protein.name), pKey))
    mProtein = str("n%sc" %(protein.seq))
    peptides = (re.sub(r"(?<=[KR])(?=[^P])",'\n', str(mProtein))).split()
    for x in xrange(MMC+1):
        mid = ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
        for i in mid:
            for j in postTransMods(i):
                output.append(returnPeptideDict(j, protein.name, pKey))
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
    peptideDB = []
    print ("[+]Generating peptide database")
    for iProtein in proteins:
        #print ("Reading %s." % (iProtein.name))
        peptides = iSilSpectra(iProtein, mmc, AAMassKey, useDecoy)
        for i in peptides:
            peptideDB.append(i)
    print "[+]Sorting peptides"
    return sorted(peptideDB, key=lambda entry: entry["mass"])

