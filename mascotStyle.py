import MGFParse
import PeptideDB
import csv

#returns sublist of peptide dataset that has correct mass
def matchMasses (searchIter, mass, tolerance, roadMap):
    output = []
    startPoint = 0
    for sign in roadMap:
        if float(sign[0]) < (mass - tolerance):
            startPoint = sign[1]
        else:
            break
    for searchTerm in searchIter[startPoint:]:
        if float(searchTerm["mass"]) < (mass - tolerance):
            continue
        elif float(searchTerm["mass"]) < (mass + tolerance):
            searchTerm["delta"] = mass - float(searchTerm["mass"])
            output.append(searchTerm)
        else:
            return output

def makeRoadMap (peptides, granularity):
    high = int(peptides[-1]["mass"])+1
    low = int(peptides[0]["mass"])
    roadMap = []
    index = 0
    mapPoints = xrange(low,high,granularity)
    mapIndex = 0
    for peptide in peptides:
        if peptide["mass"] > mapPoints[mapIndex]:
            roadMap.append((mapPoints[mapIndex], (index - 1)))
            mapIndex += 1
        index += 1
    return roadMap
        
    for i in xrange(5):
        guess = (high+low)/2
        if peptides[guess]["mass"] > query:
            high = guess
        elif peptides[guess]["mass"] < query:
            low = guess
        else:
            low = guess
            break
    return (low-2)

def simplifyPeptideDB (peptideData):
    output = []
    while len(peptideData) > 0:
        startPep = peptideData.pop(0)
        while peptideData[0]["mass"] == startPep["mass"]:
            nextPep = peptideData.pop(0)
            if nextPep["peptide"] == startPep["peptide"]:
                startPep["proteins"].append(nextPep["proteins"][0])
            else:
                peptideData.append(nextPep)
        output.append(startPep)
        print startPep["mass"]
    return output
            
def countMatches (prediction, spectra, tolerance):
    output = {"peptide":prediction["peptide"], "spec":spectra["name"] ,"proteins":prediction["proteins"],"delta":prediction["delta"], "bCount":0, "yCount":0, "bSum":0.0, "ySum":0.0}
    yIndex = 0
    bIndex = 0
    for i in spectra["m2Peaks"]:
        if prediction["yIons"][yIndex] < (i[0]-tolerance) and yIndex < len(prediction["yIons"])-1:
            yIndex +=1
        elif prediction["yIons"][yIndex] < (i[0]+tolerance) and yIndex < len(prediction["yIons"])-1:
            yIndex +=1
            output["yCount"] += 1
            output["ySum"] += i[1]
        if prediction["bIons"][bIndex] < (i[0]-tolerance) and bIndex < len(prediction["bIons"])-1:
            bIndex +=1
        elif prediction["bIons"][bIndex] < (i[0]+tolerance) and bIndex < len(prediction["bIons"])-1:
            bIndex +=1
            output["bCount"] += 1
            output["bSum"] += i[1]
    return output

def scoreCount (counter, yWeight, bWeight):
    score = float((counter["yCount"]*yWeight)+(counter["bCount"]*bWeight))/len(counter["peptide"]*7)
    return {"peptide": counter["peptide"], "proteins": counter["proteins"],"score": score*len(counter["peptide"])}

peptides = PeptideDB.returnPeptides("Uniprot-TGondii-Reference.fasta", "AAMassRef.txt", 2, True)
print "[+]"+str(len(peptides))+" peptides generated."
#peptides = sorted(simplifyPeptideDB(peptides), key=lambda entry: entry["mass"])
#print "[+]"+str(len(peptides))+" remaining after simplification."
minPeptideMass = peptides[0]["mass"]
MGFFile = open("combined.mgf", "rb")
spectraGen = MGFParse.MGFReader(MGFFile)
roadMap = makeRoadMap(peptides, 20)
#score = []

outFile = open("scoreData.csv", "wb")
keys = ["peptide", "spec", "proteins","delta", "bCount", "bSum", "yCount", "ySum"]
dict_writer = csv.DictWriter(outFile, keys, dialect="excel-tab")
dict_writer.writer.writerow(keys)

for i in spectraGen:
    print "[-]Searching " + i["name"]
    m1Matches = matchMasses(peptides, i["trueMass"], 0.0001, roadMap)
    for j in m1Matches:
        counter = countMatches(j, i, 0.5)
        dict_writer.writerow(counter)
        #scoreCounter = scoreCount(countMatches(j, i, 0.5), 5, 2)
        #score.append(scoreCounter)
        #print scoreCounter
        
