import re
from Bio import SeqIO
from multiprocessing import Pool
from functools import partial

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

class Peptide:
    "Predicted peptide object."
    def __init__(self, peptide, protein, pKey):
        self.Peptide = peptide
        self.Protein = protein
        self.pKey = pKey
        self.Mass = returnPeptideMass(peptide, self.pKey)
        self.Score = 0
        self.Match = []
    def BIons(self):
        try:
            return self.bIons
        except:
            bIons = [self.Peptide[:(i+1)] for i in range(len(self.Peptide)-1)]
            self.bIons =[returnPeptideMass(bIon, self.pKey, 1) for bIon in bIons]
            return self.bIons
    def YIons(self):
        try:
            return self.yIons
        except:
            yIons = [self.Peptide[(i+1):] for i in range(len(self.Peptide)-1)]
            self.yIons = [returnPeptideMass(yIon, self.pKey, 18) for yIon in yIons]
            return self.yIons
    def addMatch(self, peakName, score):
        if score > self.Score:
            self.Score = score
            self.Match = []
            self.Match.append(peakName)
        elif score == self.Score:
            self.Match.append(peakName)
        return
    def maxScore(self):
        return ((len(self.BIons())+len(self.YIons()))*7)+1
    def relativeScore(self):
        return float(self.Score/self.maxScore())
    def csvOutput(self):
        return str("%s\t%s\t%s\t%s\n" % (self.Peptide,self.Match,self.Score,(self.maxScore()-1)))
    
def IterAddMatch(inIter,peakDict):
    output = []
    for i in inIter:
        i.addMatch(peakDict["name"],scoreM2Peaks(peakDict["m2Peaks"], i, 1))
        output.append(i)
    return output
        
        
#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, regex, pKey):
    output = []
    peptides = re.sub(r'(?<=[RK])(?=[^P])','\n', str(protein.seq))
    for peptide in peptides.split():
         output.append(Peptide(str(peptide), str(protein.name), pKey))
    return output

#returns scores for matches in data at the M2 level
def scoreM2Peaks (iPeakList, m1Match, tolerance):
    score = int(0)
    for iPeak in iPeakList:
        for bIon in m1Match.BIons():
            if abs(iPeak - float(bIon)) < tolerance:
                score +=2
        for yIon in m1Match.YIons():
            if abs(iPeak - float(yIon)) < tolerance:
                score +=5
    return score

#returns sublist of peptide dataset that has correct mass
def MatchMasses (searchIter, mass, tolerance):
    output = []
    for searchTerm in searchIter:
        if abs(mass - float(searchTerm.Mass)) < tolerance:
            output.append(searchTerm)
    return output

#returns Peak dictionaries for entries in MGFfile
def MGFReader (mgfFile):
    outDict = {"m2Peaks":[]}
    for line in mgfFile:
        if "TITLE" in line:
            outDict["name"] = line.split(" ")[0][6:]
        if "PEPMASS" in line:
            outDict["mass"] = float(line.split(" ")[0][8:])
        if "CHARGE" in line:
            outDict["charge"] = int(line[7])
        try:
            outDict["m2Peaks"].append(float(line.split(" ")[0]))
        except:
            next     
        if "END" in line:
            outDict["trueMass"] = float((outDict["mass"]*outDict["charge"])-outDict["charge"])
            yield outDict
            outDict = {"m2Peaks":[]}

def Chunk (iArray, length):
    for i in xrange(0, len(iArray), length):
        yield iArray[i:i+length]

def returnProteinDict(protein, peptides):
    outData = {"Protein":str(protein.name)}
    outData["MaxTotalScore"] = (sum([i.maxScore() for i in peptides])+1)
    return outData

def returnResult(matchDB, protein):
    protein["Peptides"] =[]
    for i in matchDB:
        if i.Protein == protein["Protein"]:
            protein["Peptides"].append(i)
    protein["Peptides"].sort(key=lambda pep: pep.Score, reverse=True)
    protein["TotalScore"] = sum([i.Score for i in protein["Peptides"]]) 
    protein["RelativeScore"] = float(protein["TotalScore"]/protein["MaxTotalScore"])
    return protein
    
        
#Sample prediction

procNumber = 6
p = Pool(procNumber)

MGFFile = open("combined.mgf","rb")
proteinFile = open("Uniprot-TGondii-Reference.fasta", "rb")
proteins = SeqIO.parse(proteinFile, "fasta")
AAMassKey = readAAMasses("AAMassRef.txt")
peptideDB = []
proteinSet = set()
peakNumberCounter = 0
proteinDB = []
print ("[+]Generating peptide database")
for iProtein in proteins:
    peptides = iSilSpectra(iProtein, "(?<=[KR])(?=[^P])", AAMassKey)
    peptideDB += [i for i in peptides]
    proteinDB.append(returnProteinDict(iProtein, peptides))

matchDB = []
splitPeptideDB = list(Chunk(peptideDB, int(len(peptideDB)/procNumber)))
for peak in MGFReader(MGFFile):
    peakNumberCounter += 1
    if peakNumberCounter < 10:
        print ("[+]Searching %s" % (peak["name"]))
        partMatchMasses = partial(MatchMasses, mass=peak["trueMass"], tolerance=0.5)
        massMatches = p.map(partMatchMasses,splitPeptideDB)
        partAddMatch = partial(IterAddMatch, peakDict=peak)
        peptideMatches = p.map(partAddMatch, massMatches)
        for i in peptideMatches :
            for j in i:
                matchDB.append(j)

    else:
        break
for i in matchDB:
    proteinSet.add(i.Protein)

results = []
for p in proteinDB:
    if p["Protein"] in proteinSet:
        results.append(returnResult(matchDB, p))
    

print ("[+]%s results found." % (len(results)))
results.sort(key=lambda result: result["RelativeScore"],reverse=True)
outfile = open("Results.csv","wb")
outfile.write("RESULTS")
for result in results:
    outfile.write("%s\t%s\t%s\t%s\n" % (result["Protein"],result["TotalScore"],result["MaxTotalScore"],result["RelativeScore"]))
    for peptide in result["Peptides"]:
        outfile.write(str(peptide.csvOutput()))
    outfile.write("====\n")
                   

        

