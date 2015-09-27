import re

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

#Generates inSilico dictionaries of peptides with spectra predictions for input protein
def iSilSpectra (protein, regex, pKey):
    peptides = re.sub(r'(?<=[RK])(?=[^P])','\n', protein)
    for peptide in peptides.split():
        mid = {}
        mid["Protein"] = str(protein)
        mid["Peptide"] = str(peptide)
        mid["Mass"]= returnPeptideMass(peptide, pKey)
        mid["BIons"]= returnBIonPeaks (peptide, pKey)
        mid["YIons"]= returnYIonPeaks (peptide, pKey)
        yield mid

#returns true mass of peptide from peak, treating proton mass as 1 for now
def returnTrueMass (massToCharge, Charge):
    return float((massToCharge*Charge)-Charge)

#returns sublist of prediction dataset that has M1 (peptide mass) peaks in the right places
def matchMasses (iPeak, searchIter, tolerance):
    output= []
    for searchTerm in searchIter:
        if abs(returnTrueMass(iPeak[0], iPeak[1])- float(searchTerm["Mass"])) < tolerance:
            output.append(searchTerm)
    return output

#returns scores for matches in data at the M2 level
def scoreM2Peaks (iPeakList, m1Match, tolerance):
    score = int(0)
    for iPeak in iPeakList:
        for bIon in m1Match["BIons"]:
            if abs(returnTrueMass(iPeak[0], iPeak[1])- float(bIon)) < tolerance:
                score +=2
        for yIon in m1Match["YIons"]:
            if abs(returnTrueMass(iPeak[0], iPeak[1])- float(yIon)) < tolerance:
                score +=5
    return score                

#Sample prediction
samplePeaks = [{"M1Peak":[float(173),2],"M2Peaks":[[float(35),3],[float(79.5),2],[float(190),1],[float(123),2]]},
               {"M1Peak":[float(244),2],"M2Peaks":[[float(39),3],[float(68),3],[float(331),1],[float(187.5),2],[float(287),1],[float(53.3),3]]}]
proteins = ["PREDICT","TSRISER","STRAWLIGHT","CITED","CARPENTER"]
AAMassKey = readAAMasses("AAMassRef.txt")
for peak in samplePeaks:
    for iProtein in proteins:
        peptideMatches = matchMasses(peak["M1Peak"], iSilSpectra(iProtein, "(?<=[KR])(?=[^P])", AAMassKey), 2)
        for i in peptideMatches:
            i["Score"] = scoreM2Peaks(peak["M2Peaks"], i, 2)
            print i

        
