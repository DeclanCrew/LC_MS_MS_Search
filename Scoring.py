import sys, getopt
from multiprocessing import Pool
from functools import partial

def inputOutput (argv):
    options = {"maxMissedCleave":1}
    try:
        opts, args = getopt.getopt(argv, "hi:o:m:pc")
    except getopt.GetoptError:
        print "PeptideDB.py -i <inputpep> -o <output> -m <AAmassref>"
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

def csvToDictRead (icsv):
    ifile = open(icsv, "rb")
    idata = csv.DictReader(ifile, dialect="excel-tab")
    return idata

def IterAddMatch(inIter,peakDict):
    output = []
    for i in inIter:
        i.addMatch(peakDict["name"],scoreM2Peaks(peakDict["m2Peaks"], i, 1))
        output.append(i)
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
            outDict["m2Peaks"].append((float(line.split(" ")[0]),float(line.split(" ")[1])))
        except:
            next     
        if "END" in line:
            outDict["trueMass"] = float((outDict["mass"]*outDict["charge"])-outDict["charge"])
            yield outDict
            outDict = {"m2Peaks":[]}

def Chunk (iArray, length):
    for i in xrange(0, len(iArray), length):
        yield iArray[i:i+length]

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

peakNumberCounter = 0
matchDB = []
proteinSet = set()
splitPeptideDB = list(Chunk(peptideDB, int(len(peptideDB)/procNumber)))
for peak in MGFReader(MGFFile):
    peakNumberCounter += 1
    if peakNumberCounter < 3:
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
                   

        

