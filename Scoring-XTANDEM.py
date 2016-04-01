import sys, getopt
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import ast
from MGFParse import MGFReader
from scipy import stats

def inputOutput (argv):
    options = {"maxMissedCleave":1}
    print argv
    try:
        opts, args = getopt.getopt(argv, "hi:o:m:p:c:")
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
    odata = []
    for i in idata:
        odata.append(i)
    return odata
             
#returns hyper scores for matches in data at the M2 level
def returnHyperScore (iPeakList, m1Match, tolerance):
    score = 0
    bCount = 0
    yCount = 0
    for iPeak in iPeakList:
        for bIon in m1Match["bIons"]:
            if abs(iPeak[0] - float(bIon)) < tolerance:
                ##print iPeak[1]
              ##  print m1Match["mass"]
                score +=(iPeak[1])
                bCount +=1
        for yIon in m1Match["yIons"]:
            if abs(iPeak[0] - float(yIon)) < tolerance:
                score +=(iPeak[1])
                yCount +=1
    hyperScore = (float(score)*(bCount)*(yCount))
    return hyperScore

#returns sublist of peptide dataset that has correct mass
def matchMasses (searchIter, mass, tolerance):
    output = []
    for searchTerm in searchIter:
        if abs(mass - float(searchTerm["mass"])) < tolerance:
            output.append(searchTerm)
    return output

def returnResult(matchDB, protein):
    protein["Peptides"] =[]
    for i in matchDB:
        if i.Protein == protein["Protein"]:
            protein["Peptides"].append(i)
    protein["Peptides"].sort(key=lambda pep: pep.Score, reverse=True)
    protein["TotalScore"] = sum([i.Score for i in protein["Peptides"]]) 
    protein["RelativeScore"] = float(protein["TotalScore"]/protein["MaxTotalScore"])
    return protein

def returnHistogram(Input_array):
    histogram ={}
    histogram["n"], histogram["bins"] = np.histogram(Input_array, int(max(Input_array)/30), density=False)
    return histogram

def returnSurvival(histogram):
    zeroes = 0
    output = []
    for i in xrange(len(histogram["n"])):
        if histogram["n"][i] == 0:
            zeroes += 1
        elif zeroes == 1:
            histogram["n"][i]= 1
        elif zeroes > 1:
            histogram["n"][i] =0

    maxIndex = 0
    maxValue = 0
    zeroIndex = 0
    for i in xrange(len(histogram["n"])):
        if histogram["n"][i] > maxValue:
            maxValue = histogram["n"][i]
            maxIndex = i            
    for i in [(len(histogram["n"])-(j+1)) for j in xrange(len(histogram["n"]))]:
        if histogram["n"][i] != 0:
            zeroIndex = i
            break
    histogram["bins"] = histogram["bins"][maxIndex:(zeroIndex+2)]
    histogram["n"] = histogram["n"][maxIndex:(zeroIndex+1)]

    totalValue = 0
    for i in xrange(len(histogram["n"])):
        totalValue += histogram["n"][-(i+1)]
        histogram["n"][-(i+1)] = totalValue

    histogram["n"] = [convolScore(i) for i in histogram["n"]]
    return histogram

def convolScore(iScore):
    if (iScore <= 0):
        return 1
    else:
        return np.math.log10(iScore)
        
def returnEvalue(histogram, maxScore):

  #  print str( histogram["bins"][maxIndex:zeroIndex]) +str( histogram["n"][maxIndex:zeroIndex])
    regressLine = stats.linregress(histogram["bins"][1:(len(histogram["bins"]))], histogram["n"])
    #print regressLine
    return ((regressLine[0]*(maxScore)) + regressLine[1])
    
def make_Histogram (histogram, Set_name):
    fig, ax = plt.subplots()
    n = histogram["n"]
    bins = histogram["bins"]
    left = np.array(bins[:-1])
    right = np.array(bins[1:])
    bottom = np.zeros(len(left))
    top = bottom + n
    XY = np.array([[left,left,right,right], [bottom,top,top,bottom]]).T
    barpath = path.Path.make_compound_path_from_polys(XY)
    patch = patches.PathPatch(barpath, facecolor='blue', edgecolor='gray', alpha=0.8)
    ax.add_patch(patch)
    ax.set_xlim(left.min(), right.max())
    ax.set_ylim(bottom.min(), top.max())
    ax.set_title(Set_name)
    ax.set_xlabel("HyperScore")
    ax.set_ylabel("Frequency")
    plt.savefig("Images/" + Set_name + ".png")
    return

def writeToCSV (ofile, oArray):
    keys = set()
    for k in oArray[0]:
        keys.add(str(k))
    out = open(ofile,"wb")
    dict_writer = csv.DictWriter(out, list(keys), dialect="excel-tab")
    dict_writer.writer.writerow(list(keys))
    dict_writer.writerows(oArray)
    return

#Sample prediction
MGFFile = open("combined.mgf","rb")

proteinSet = set()
print "Reading peptideDB"
peptideDB = csvToDictRead("FullReference.pep")
peakListFile = open("HighestXTandemResults.list","rb")
peakList = [i.strip("\n") for i in peakListFile.readlines()]
peakList = ["Spectrum%s" % (i) for i in peakList]
print "Reading MGF file"
peaksRead = [i for i in MGFReader(MGFFile)]
print str("MGF read %s entries found." % (len(peaksRead)))
output = []
counter = 0

for peak in peaksRead:
##    print ("[+]Searching %s" % (peak["name"]))
##    print len(peak["m2Peaks"])
    counter += 1
    if counter <= 1000:
        outDict = {}
        m1Matches = matchMasses(peptideDB, peak["trueMass"], float(0.1))
        hyperScoreList = []
        maxScore = 0
        maxEntry = dict()
        for i in m1Matches:
           ## try:
            i["bIons"]= ast.literal_eval(i["bIons"])
            i["yIons"]= ast.literal_eval(i["yIons"])
            proteinSet.add(i["protein"])
            hS = returnHyperScore(peak["m2Peaks"],i, 0.5)
            hyperScoreList.append(hS)
            if hS > maxScore:
                maxScore = hS
                maxEntry = i
         ##   except:
           ##     continue
        if len(hyperScoreList) > 3:
            try:
                rawHistogram = returnHistogram(hyperScoreList)
                #make_Histogram(rawHistogram, str("ScoreHistogramDec" + peak["name"]))
                survival = returnSurvival(rawHistogram)
                #make_Histogram(survival, str("SurvivalFuncDec" + peak["name"]))
                if len(survival["bins"] > 3):
                    outDict["proteinName"] = maxEntry["protein"]
                    outDict["peptide"] = maxEntry["peptide"]
                    outDict["evalue"] = returnEvalue(survival, maxScore)
                    print str("%s\t%s"%(outDict["proteinName"],outDict["evalue"]))
                    output.append(outDict)
            except:
                print "Failed to Generate."
                continue


outFile = open("InitialResults.res","wb")
writeToCSV(outFile, output)
