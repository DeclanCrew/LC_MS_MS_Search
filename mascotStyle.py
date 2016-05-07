'''
MascotStyle, contains the main loop and the core scoring components of the
search engine
'''
import MGFParse
import PeptideDB
import configs
import csv
from itertools import groupby

def matchMasses(searchIter, mass, tolerance):
    '''Returns the subset of the searchDB that has a mass within tolerance'''
    output = []
    for searchTerm in searchIter[int(mass)]:
        if float(searchTerm["mass"]) < (mass - tolerance):
            continue
        elif float(searchTerm["mass"]) < (mass + tolerance):
            searchTerm["delta"] = mass - float(searchTerm["mass"])
            searchTerm["bCount"] = 0
            searchTerm["bSum"] = 0
            searchTerm["yCount"] = 0
            searchTerm["ySum"] = 0
            output.append(searchTerm)
        else:
            return output

def grabIons(matchList, conf):
    '''Takes peptides and returns a sorted list of their predicted ions'''
    output = []
    for i in enumerate(matchList):
        consts = conf["other_constants"]
        masses = i[1]["orderedMasses"]
        bIons = PeptideDB.coarsen(PeptideDB.returnIons(masses, consts["B+"]))
        yIons = PeptideDB.coarsen(PeptideDB.returnIons(masses[::-1], consts["Y+"]))
        for j in bIons:
            output.append([int(j), int((i[0]*2)+1)])
        for j in yIons:
            output.append([int(j), int(i[0])*2])
    return sorted(output, key=lambda entry: entry[0])

def countMatches(matchList, spectra, conf):
    '''Updates peptides with counts and intensity of matches against spectra'''
    if matchList == None:
        return []
    ions = grabIons(matchList, conf)
    matchIndex = 0
    peptides = matchList
    for i in spectra["m2Peaks"]:
        if matchIndex == len(ions):
            break
        while ions[matchIndex][0] < i[0] and matchIndex < (len(ions)-1):
            matchIndex += 1
        if matchIndex == len(ions):
            break
        if ions[matchIndex][0] > i[0]:
            continue
        if ions[matchIndex][1] % 2:
            peptides[ions[matchIndex][1]/2]["yCount"] += 1
            peptides[ions[matchIndex][1]/2]["ySum"] += int(i[1])
        else:
            peptides[(ions[matchIndex][1]+1)/2]["bCount"] += 1
            peptides[(ions[matchIndex][1]+1)/2]["bSum"] += int(i[1])
        matchIndex += 1
    output = []
    for entry in peptides:
        entry["spec"] = spectra["name"]
        output.append(entry)
    return output

def scoreMethod(entry, bias=2.5):
    '''Simple mascot style scoring system, weights yIons stronger than b'''
    return ((bias*entry["yCount"]*entry["ySum"])+(entry["bCount"]*entry["bSum"]))

def returnMaxima(data, value):
    '''Returns the highest scoring peptide for each spectrum'''
    sortedData = sorted(data, key=lambda entry: entry[value], reverse=True)
    if len(sortedData) > 1:
        penultimate = sortedData[1]
    else:
        penultimate = sortedData[0]
    maximum = sortedData[0]
    maximum["diff"] = maximum[value] - penultimate[value]
    return maximum

def generateWriter(dataFileName, keys):
    dataFile = open(dataFileName, "wb")
    writerObject = csv.DictWriter(dataFile, keys, dialect="excel-tab",
                                  extrasaction="ignore")
    writerObject.writer.writerow(keys)
    while True:
        yield writerObject.writerow

configurations = configs.readConfigs("mascotStyle.cfg")
fileConfs = configurations["data"]
peptideSet = PeptideDB.returnPeptides(configurations)
print "[+]"+str(len(peptideSet))+" mass entries generated."
spectra_file = fileConfs["spectra_file"]
MGFFile = open(spectra_file, "rb")
spectraGen = MGFParse.MGFReader(MGFFile, configurations)
full_results = bool(fileConfs["write_full_scores"])

if full_results == True:
    fullKeys = ["peptide", "spec", "proteins", "delta",
                "bCount", "bSum", "yCount", "ySum"]
    full_writer = generateWriter(fileConfs["full_scores_file"], fullKeys)

topKeys = ["peptide", "spec", "proteins", "score", "diff", "delta",
           "bCount", "bSum", "yCount", "ySum"]
top_writer = generateWriter(fileConfs["top_scores_file"], topKeys)
topScores = []

for spectrum in spectraGen:
    print ("[+]Searching %s. \r" % (spectrum["name"])),
    m1Matches = matchMasses(peptideSet, spectrum["trueMass"], 0.001)
    counter = countMatches(m1Matches, spectrum, configurations)
    if len(counter) < 2:
        continue
    for result in counter:
        if full_results == True:
            next(full_writer)(result)
        result["score"] = scoreMethod(result)
    topScore = returnMaxima(counter, "score")
    next(top_writer)(topScore)
    topScores.append(topScore)
print ""
print "[+]Generating protein level scores."

clearPepSet = set()
proteinScores = {}
for i in sorted(topScores, key=lambda entry: entry["diff"], reverse=True):
    if i["peptide"] in clearPepSet:
        continue
    clearPepSet.add(i["peptide"])
    for j in i["proteins"]:
        if j in proteinScores:
            proteinScores[j] += i["diff"]
        else:
            proteinScores[j] = i["diff"]
protFile = open(fileConfs["prot_scores_file"], "wb")
for protein in sorted(proteinScores.items(), key=lambda prot: prot[1], reverse=True):
    protFile.write(str(protein[0]+"\t"+str(protein[1])+"\n"))
