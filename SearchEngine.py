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

#N.B need to add adjustment to peptide Mass calculation, to account for asymmetry in cleaving
def returnBIonPeaks (peptide, pKey):
    bIons = [peptide[:(i+1)] for i in range(len(peptide)-1)]
    return [returnPeptideMass(bIon, pKey) for bIon in bIons]

def returnYIonPeaks (peptide, pKey):
    yIons = [peptide[(i+1):] for i in range(len(peptide)-1)]
    return [returnPeptideMass(yIon, pKey) for yIon in yIons]

class Peptide:
    "Predicted peptide object."
    def __init__(self, peptide, protein, pKey):
        self.Peptide = peptide
        self.Protein = protein
        self.pKey = pKey
        self.Mass = returnPeptideMass(peptide, self.pKey)
        self.Score = 0
    def BIons(self):
        return returnBIonPeaks(self.Peptide, self.pKey)
    def YIons(self):
        return returnYIonPeaks(self.Peptide, self.pKey)

#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, regex, pKey):
    peptides = re.sub(r'(?<=[RK])(?=[^P])','\n', str(protein.seq))
    for peptide in peptides.split():
        yield Peptide(str(peptide), str(protein.name), pKey)

#returns true mass of peptide from peak, treating proton mass as 1 for now
def returnTrueMass (massToCharge, Charge):
    return float((massToCharge*Charge)-Charge)

#returns sublist of prediction dataset that has M1 (peptide mass) peaks in the right places
def matchMasses (iPeak, searchIter, tolerance):
    for searchTerm in searchIter:
        if abs(returnTrueMass(iPeak[0], iPeak[1])- float(searchTerm.Mass)) < tolerance:
            yield searchTerm

#returns scores for matches in data at the M2 level
def scoreM2Peaks (iPeakList, m1Match, tolerance):
    score = int(0)
    for iPeak in iPeakList:
        for bIon in m1Match.BIons():
            if abs(returnTrueMass(iPeak[0], iPeak[1])- float(bIon)) < tolerance:
                score +=2
        for yIon in m1Match.YIons():
            if abs(returnTrueMass(iPeak[0], iPeak[1])- float(yIon)) < tolerance:
                score +=5
    return score                

#Sample prediction
samplePeaks = [{"M1Peak":[float(244),2],"M2Peaks":[[float(39),3],[float(68),3],[float(331),1],[float(187.5),2],[float(287),1],[float(53.3),3]]}]
proteinFile = open("Uniprot-TGondii-Reference.fasta", "rb")
proteins = SeqIO.parse(proteinFile, "fasta")
AAMassKey = readAAMasses("AAMassRef.txt")
peakCount = 0
for peak in samplePeaks:
    peakCount += 1
    print ("[+]Searching peak #%s" % (peakCount))
    highestScore = 0
    highestMatches = []
    for iProtein in proteins:
        peptideMatches = matchMasses(peak["M1Peak"], iSilSpectra(iProtein, "(?<=[KR])(?=[^P])", AAMassKey), 0.5)
        for i in peptideMatches:
            i.Score = scoreM2Peaks(peak["M2Peaks"], i, 1)
            if i.Score > highestScore:
                highestScore = i.Score
                highestMatches = []
                highestMatches.append(i)
            elif i.Score == highestScore:
                highestMatches.append(i)
    for highestMatch in highestMatches:
        print ("Matched peptide %s on protein %s with score %s." % (highestMatch.Peptide, highestMatch.Protein, highestMatch.Score))

        
