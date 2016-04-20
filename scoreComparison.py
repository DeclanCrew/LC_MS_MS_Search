import csv
import itertools
from math import factorial, log10
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np
from PIL import Image
import MGFParse


def returnMaxima (generator, value):
    groups = itertools.groupby(generator, lambda entry: entry["spec"])
    for spec, data in groups:
        sortedData = sorted(data, key=lambda entry: entry[value], reverse=True)
        if len(sortedData) < 2:
            continue
        penultimate = sortedData[1]
        maximum = sortedData[0]
        maximum["diff"] = maximum[value] - penultimate[value]
        yield maximum

##MGFFile = open("combined.mgf", "rb")
##spectraGen = MGFParse.MGFReader(MGFFile)
##spectraNames = {}
##counter = 1
##for i in spectraGen:
##    spectraNames[counter] = i["name"]
##    counter += 1
##print counter
tandemData = csv.DictReader(open("newTandemData.csv","rb"), dialect="excel-tab")
filteredTandemData = []

for i in tandemData:
##    if int(i["spectrum"]) > counter:
##        continue
##    i["spec"] = spectraNames[int(i["spectrum"])]
    i["spec"] = int(i["spectrum"])
    i["score"] = i["hyperscore"]
    filteredTandemData.append(i)
##    if i["delta"] == None:
##        continue
##    if abs(float(i["delta"]))< 0.0001:
##        i["proteins"] = [i["proteins"]]
##        filteredTandemData.append(i)
##sortedTandemData = sorted(tandemData, key=lambda j : j["hyperscore"], reverse=True)

def scoreMethod (entry):
    return ((5*entry["yCount"]*entry["ySum"])+(2*entry["bCount"]*entry["bSum"]))
    #return (factorial(entry["yCount"]+entry["bCount"])*(entry["ySum"]+entry["bSum"]))

def scoreProcess (record):
    output = {"bCount":int(record["bCount"]),"spec":record["spec"],
              "bSum":float(record["bSum"]),"yCount":int(record["yCount"]),
              "ySum":float(record["ySum"]),"delta":float(record["delta"])}
    output["proteins"] = [i.strip(r" '").lstrip(r" '") for i in record["proteins"].strip("[]").split(",")]
    output["length"] = len(record["peptide"])
    output["score"] = scoreMethod(output)
    output["peptide"] = record["peptide"]
    return output

def stage2Score (record):
    scoring_coeff = [-2.02849873e-05, -4.28916096e-06, 1.24617099e-02]
    record["2score"] = (scoring_coeff[0]*record["ySum"]*record["yCount"]+
                        scoring_coeff[1]*record["bSum"]*record["bCount"]+
                        scoring_coeff[2]*record["diff"])
    return record

    
mascotStyleData = csv.DictReader(open("scoreData.csv","rb"), dialect="excel-tab")
scoredMascotData = [i for i in returnMaxima(map(scoreProcess, mascotStyleData), "score")]
simpDataFile = open("MascScoreData.csv", "wb")
keys = ["peptide", "spec", "proteins", "delta", "score", "yCount", "ySum", "bCount", "bSum", "diff"]
dict_writer = csv.DictWriter(simpDataFile, keys, dialect="excel-tab", extrasaction="ignore")
dict_writer.writer.writerow(keys)
for i in scoredMascotData:
    if (i["yCount"]+i["bCount"]) > 3:
        dict_writer.writerow(i)
        for j in filteredTandemData:
            if j["peptide"] == i["peptide"]:
                dict_writer.writerow(j)


#for i in sorted(list(scoredMascotData + filteredTandemData), key= lambda entry: int(entry["spec"][8:].split("-")[0])):
#    dict_writer.writerow(i)
##
##scoredMascotData = map(stage2Score, scoredMascotData)
##sortedMascotData = sorted(scoredMascotData, key=lambda j: j["2score"], reverse=True)
##
##def returnHits (sortedData, maxFDR):
##    outSet = set()
##    decoyCount = 0
##    count = 0
##    for i in sortedData:
##        for j in i["proteins"]:
##            count += 1
##            if j[0] == "D":
##                #decoyCount += 100.0
##                if decoyCount/count > maxFDR:
##                    print i
##                    prompt = raw_input("Continue? Y/N \n")
##                    if prompt[0] == "N":
##                        return outSet
##                    else:
##                        decoyCount -= 100.0
##            #print j
##            outSet.add(j)
##    return outSet
##
##mascotStyleSet = returnHits(sortedMascotData, 5)
##tandemSet = returnHits(sortedTandemData, 5)
##with sns.axes_style("darkgrid"):
##    plt.figure(figsize=(4,4))
##    venn2([mascotStyleSet, tandemSet], set_labels = ("mascotStyle", "X!Tandem"))
##    plt.title("Venn diagram of mascotStyle hits versus X!Tandem peptide hits.")
##    plt.show()
##
##print ("Number of mascotStyle Hits: "+str(len(mascotStyleSet)))
##print ("Number of tandem Hits: "+str(len(tandemSet)))
##print ("Number of common Hits: "+str(len(tandemSet & mascotStyleSet)))
##
##outfile = open("tandemunique.txt", "wb")
##for i in (tandemSet - mascotStyleSet):
##    outfile.write(i+"\n")
##    
##with sns.axes_style("darkgrid"):
##    mascotTargetScores = [log10(i["score"]+1) for i in sortedMascotData if i["proteins"][0][0] == "t"]
##    mascotDecoyScores = [log10(i["score"]+1) for i in sortedMascotData if i["proteins"][0][0] == "D"]
##    targetScores = [i["hyperscore"] for i in sortedTandemData]
##    plt.figure()
##    sns.distplot(mascotTargetScores, bins=20)
##    sns.distplot(mascotDecoyScores, bins=20)
##    plt.legend()
##    plt.show()
##
##maxValues = [0,0]
##pixels = []
##for i in sortedMascotData:
##    pixel = {"x": int(i["bSum"]), "y": int(i["ySum"]), "proteins": i["proteins"]}
##    pixels.append(pixel)
##    if pixel["x"] > maxValues[0]:
##        maxValues[0] = pixel["x"]+1
##    if pixel["y"] > maxValues[1]:
##        maxValues[1] = pixel["y"]+1
##
##size = (maxValues[0], maxValues[1])
##print str(size)
####outputImage = Image.new("L", size, 0)
####for j in pixels:
####    outputImage.putpixel((j["x"],j["y"]), 255)
####outputImage.save("heatMap2.png")
##
##outputImage = Image.new("RGB", size)
##for j in pixels:
##    if j["proteins"][0][0] == "D":
##        outputImage.putpixel((j["x"],j["y"]), (255,0,0))
##    else:
##        outputImage.putpixel((j["x"],j["y"]), (0,0,255))
##outputImage.save("scoremap.png")


