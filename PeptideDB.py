import re
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
    totalMass = 0
    for i in peptide:
        totalMass += pKey[i]
    totalMass += adjustment
    return totalMass

def returnPeptideDict (peptide, pTMs, protein, pKey):
    output = {"peptide":peptide,"pTMs":pTMs,"protein":protein,"pKey":pKey}
    output["Mass"]= returnPeptideMass(output["peptide"], output["pKey"])
    
    bIons = [self.Peptide[:(i+1)] for i in range(len(self.Peptide)-1)]
    output["bIons"]= [returnPeptideMass(bIon, output["pKey"], 1) for bIon in bIons]
    
    yIons = [self.Peptide[(i+1):] for i in range(len(self.Peptide)-1)]
    output["yIons"]= [returnPeptideMass(yIon, output["pKey"], 18) for yIon in yIons]
    
    output["maxScore"]= ((len(bIons)+len(yIons))*7)+1
    return output

#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, MMC):
    output = []
    peptides = (re.sub(r"(?<=[KR])(?=[^P])",'\n', str(protein))).split()
    for x in xrange(MMC+1):
        output += ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
    return output

def returnProteinDict(protein, peptides):
    outData = {"Protein":str(protein.name)}
    outData["MaxTotalScore"] = (sum([i.maxScore() for i in peptides])+1)
    return outData

def postTransMods (peptide, regex=0, nRegex=0, cRegex=0):
    output = []
    if peptide[0] = "n":
        peptide = peptide[1:]
    if peptide[len(peptide)] = "c":
        peptide = peptide[:(len(peptide)-1)]
    
##
##proteinFile = open("Uniprot-TGondii-Reference.fasta", "rb")
##proteins = SeqIO.parse(proteinFile, "fasta")
proteins = ["PEREPRTITIE"] 
AAMassKey = readAAMasses("AAMassRef.txt")
peptideDB = []
proteinDB = []
maxMissedCleave = 1
print ("[+]Generating peptide database")
for iProtein in proteins:
    mProtein = str("n"+iProtein+"c")
    print mProtein
    peptides = iSilSpectra(mProtein, maxMissedCleave)
    print peptides
    #peptideDB += [i for i in peptides]
    #proteinDB.append(returnProteinDict(iProtein, peptides))
