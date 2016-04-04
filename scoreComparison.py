import csv
import itertools
from math import factorial
from matplotlib import pyplot as plt
from matplotlib_venn import venn2
import seaborn as sns
import numpy as np
from PIL import Image


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

tandemData = csv.DictReader(open("tandemData.csv","rb"), dialect="excel-tab")
filteredTandemData = []
for i in tandemData:
    if i["delta"] == None:
        continue
    if abs(float(i["delta"]))< 0.0005:
        i["proteins"] = [i["proteins"]]
        filteredTandemData.append(i)
sortedTandemData = sorted(filteredTandemData, key=lambda j : j["hyperscore"], reverse=True)

def scoreMethod (entry):
    return ((5*entry["yCount"]*entry["ySum"])+(2*entry["bCount"]*entry["bSum"]))
    #return (factorial(entry["yCount"]+entry["bCount"])*(entry["ySum"]+entry["bSum"]))

def scoreProcess (record):
    output = {"bCount":int(record["bCount"]),"spec":record["spec"],"bSum":float(record["bSum"]),"yCount":int(record["yCount"]),"ySum":float(record["ySum"]),"delta":float(record["delta"])}
    output["proteins"] = [i.strip(r" '").lstrip(r" '") for i in record["proteins"].strip("[]").split(",")]
    output["length"] = len(record["peptide"])
    output["score"] = scoreMethod(output)
    output["peptide"] = record["peptide"]
    return output

mascotStyleData = csv.DictReader(open("scoreData.csv","rb"), dialect="excel-tab")
scoredMascotData = returnMaxima([scoreProcess(i) for i in mascotStyleData], "score")
sortedMascotData = sorted(scoredMascotData, key=lambda j: j["score"], reverse=True)

def returnHits (sortedData, maxFDR):
    outSet = set()
    decoyCount = 0
    count = 0
    for i in sortedData:
        if "diff" in i and i["diff"] < i["score"]/10:
            continue
        for j in i["proteins"]:
            count += 1
            if j[0] == "D":
                decoyCount += 100.0
                if decoyCount/count > maxFDR:
                    print i
                    prompt = raw_input("Continue? Y/N \n")
                    if prompt[0] == "N":
                        return outSet
                    else:
                        decoyCount -= 100.0
            print j
            outSet.add(j)
    return outSet

mascotStyleSet = returnHits(sortedMascotData, 5)
tandemSet = returnHits(sortedTandemData, 5)
with sns.axes_style("darkgrid"):
    plt.figure(figsize=(4,4))
    venn2([mascotStyleSet, tandemSet], set_labels = ("mascotStyle", "X!Tandem"))
    plt.title("Venn diagram of mascotStyle hits versus X!Tandem peptide hits.")
    plt.show()

print ("Number of mascotStyle Hits: "+str(len(mascotStyleSet)))
print ("Number of tandem Hits: "+str(len(tandemSet)))
print ("Number of common Hits: "+str(len(tandemSet & mascotStyleSet)))

outfile = open("tandemunique.txt", "wb")
for i in (tandemSet - mascotStyleSet):
    outfile.write(i+"\n")
    
with sns.axes_style("darkgrid"):
    mascotTargetScores = [(i["diff"]) for i in sortedMascotData if i["proteins"][0][0] == "t"]
    mascotDecoyScores = [(i["diff"]) for i in sortedMascotData if i["proteins"][0][0] == "D"]
    targetScores = [i["hyperscore"] for i in sortedTandemData]
    plt.figure()
    sns.distplot(mascotTargetScores, bins=20)
    sns.distplot(mascotDecoyScores, bins=20)
    plt.legend()
    plt.show()

maxValues = [0,0]
pixels = []
for i in sortedMascotData:
    pixel = {"x": int(i["bSum"]), "y": int(i["ySum"]), "proteins": i["proteins"]}
    pixels.append(pixel)
    if pixel["x"] > maxValues[0]:
        maxValues[0] = pixel["x"]+1
    if pixel["y"] > maxValues[1]:
        maxValues[1] = pixel["y"]+1

size = (maxValues[0], maxValues[1])
print str(size)
##outputImage = Image.new("L", size, 0)
##for j in pixels:
##    outputImage.putpixel((j["x"],j["y"]), 255)
##outputImage.save("heatMap2.png")

outputImage = Image.new("RGB", size)
for j in pixels:
    if j["proteins"][0] == "D":
        outputImage.putpixel((j["x"],j["y"]), (255,0,0))
    else:
        outputImage.putpixel((j["x"],j["y"]), (0,0,255))
outputImage.save("scoremap.png")


