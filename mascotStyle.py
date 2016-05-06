'''
MascotStyle, contains the main loop and the core scoring components of the
search engine
'''
import MGFParse
import PeptideDB
import configs
import csv

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
    if len(sortedData) < 2:
        yield sortedData[-1]
    penultimate = sortedData[1]
    maximum = sortedData[0]
    maximum["diff"] = maximum[value] - penultimate[value]
    yield maximum

configurations = configs.readConfigs("mascotStyle.cfg")
peptideSet = PeptideDB.returnPeptides("target-decoy-gondii.fasta",
                                    configurations)
print "[+]"+str(len(peptideSet))+" mass entries generated."
MGFFile = open("combined.mgf", "rb")
spectraGen = MGFParse.MGFReader(MGFFile, configurations)
##full_results = True

##if full_results == True:
fullFile = open("scoreData.csv", "wb")
fullKeys = ["peptide", "spec", "proteins", "delta",
            "bCount", "bSum", "yCount", "ySum"]
full_writer = csv.DictWriter(fullFile, fullKeys, dialect="excel-tab",
                             extrasaction="ignore")
full_writer.writer.writerow(fullKeys)

##maxFile = open("topScores.csv", "wb")
##maxKeys = ["peptide", "spec", "proteins", "score", "diff", "delta", "bCount", "bSum", "yCount", "ySum"]

for spectrum in spectraGen:
    print ("[-]Searching %s \r" % (spectrum["name"])),
    m1Matches = matchMasses(peptideSet, spectrum["trueMass"], 0.001)
    counter = countMatches(m1Matches, spectrum, configurations)
    for result in counter:
##        if full_results == True:
        full_writer.writerow(result)
##        result["score"] = scoreMethod(result)
##    maxScore = returnMaxima(counter, "score")
print ""
