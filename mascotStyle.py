import MGFParse
import PeptideDB
import csv

#returns sublist of peptide dataset that has correct mass
def matchMasses (searchIter, mass, tolerance):
    output = []
    peptides = {}
    for searchTerm in searchIter[int(mass)]:
        if float(searchTerm["mass"]) < (mass - tolerance):
            continue
        elif float(searchTerm["mass"]) < (mass + tolerance):
            searchTerm["delta"] = mass - float(searchTerm["mass"])
            searchTerm["bCount"] = 0
            searchTerm["yCount"] = 0
            searchTerm["bSum"] = 0
            searchTerm["ySum"] = 0
            peptides[searchTerm["peptide"]] = searchTerm
            for i in searchTerm["yIons"]:
                output.append((i, "y", searchTerm["peptide"]))
            for j in searchTerm["bIons"]:
                output.append((j, "b", searchTerm["peptide"]))
        else:
            return [sorted(output, key=lambda entry: entry[0]), peptides]

def countMatches (matchList, spectra):
    if matchList == None:
        return []
    matchIndex = 0
    ions = matchList[0]
    peptides = matchList[1]
    for i in spectra["m2Peaks"]:
        if matchIndex == len(ions):
            break
        while ions[matchIndex][0] < i and matchIndex < (len(ions)-1):
            matchIndex += 1
        if matchIndex == len(ions):
            break
        if ions[matchIndex][0] > i:
            continue
        peptides[ions[matchIndex][2]][str(ions[matchIndex][1]+"Count")] += 1
        peptides[ions[matchIndex][2]][str(ions[matchIndex][1]+"Sum")] += int(spectra["m2Peaks"][i])
        matchIndex += 1
    output = []
    for entry in peptides:
        peptides[entry]["spec"] = spectra["name"]
        output.append(peptides[entry])
    return output

def scoreCount (counter, yWeight, bWeight):
    score = float((counter["yCount"]*yWeight)+(counter["bCount"]*bWeight))/len(counter["peptide"]*7)
    return {"peptide": counter["peptide"], "proteins": counter["proteins"],"score": score*len(counter["peptide"])}

peptides = PeptideDB.returnPeptides("target-decoy-gondii.fasta", "AAMassRef.txt", 2, False)
print "[+]"+str(len(peptides))+" mass entries generated."
MGFFile = open("combined.mgf", "rb")
spectraGen = MGFParse.MGFReader(MGFFile)

outFile = open("scoreData.csv", "wb")
keys = ["peptide", "spec", "proteins","delta", "bCount", "bSum", "yCount", "ySum"]
dict_writer = csv.DictWriter(outFile, keys, dialect="excel-tab", extrasaction="ignore")
dict_writer.writer.writerow(keys)

for i in spectraGen:
    print "[-]Searching " + i["name"]
    m1Matches = matchMasses(peptides, i["trueMass"], 0.0001)
    counter = countMatches(m1Matches, i)
    for result in counter:
        dict_writer.writerow(result)
        #scoreCounter = scoreCount(countMatches(j, i, 0.5), 5, 2)
        #score.append(scoreCounter)
        #print scoreCounter
        
