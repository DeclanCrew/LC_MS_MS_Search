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
                score +=iPeak[1]
                bCount +=1
        for yIon in m1Match["yIons"]:
            if abs(iPeak[0] - float(yIon)) < tolerance:
                score +=iPeak[1]
                yCount +=1
    hyperScore = (score)*(np.math.factorial(bCount))*(np.math.factorial(yCount))
    return hyperScore

#returns sublist of peptide dataset that has correct mass
def matchMasses (searchIter, mass, tolerance):
    output = []
    for searchTerm in searchIter:
        if abs(mass - float(searchTerm["Mass"])) < tolerance:
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

def returnEvalue (Input_array):
    maxScore = max(Input_array)
    n, bins = np.histogram(Input_array, len(Input_array))
    n = [np.log((i+1)) for i in n]
    bins = [i for i in bins]
    while len(n) < len(bins):
        n.append(0)
    regressLine = stats.linregress(bins, n)
    return np.exp((regressLine[0]*maxScore) + regressLine[1])
    
def make_Histogram (Input_array, Set_name):
    fig, ax = plt.subplots()
    n, bins = np.histogram(Input_array, len(Input_array))
    n = [np.log(i) for i in n]
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

#Sample prediction
MGFFile = open("combined.mgf","rb")

proteinSet = set()
print "Reading peptideDB"
peptideDB = csvToDictRead("Ref-forward.pep")
peakListFile = open("HighestXTandemResults.list","rb")
peakList = [i.strip("\n") for i in peakListFile.readlines()]
peakList = ["Spectrum%s" % (i) for i in peakList]

for peak in MGFReader(MGFFile):
    print ("[+]Searching %s" % (peak["name"]))
    m1Matches = matchMasses(peptideDB, peak["trueMass"], float(0.1))
    hyperScoreList = []
    maxScore = 0
    maxEntry = dict()
    for i in m1Matches:
        try:
            i["bIons"]= ast.literal_eval(i["bIons"])
            i["yIons"]= ast.literal_eval(i["yIons"])
            proteinSet.add(i["protein"])
            hS = returnHyperScore(peak["m2Peaks"],i, 0.5)
            hyperScoreList.append(hS)
            if hS > maxScore:
                maxScore = hS
                maxEntry = i
        except:
            continue
    if len(hyperScoreList) > 3:
        try:
            print maxEntry["peptide"]
            print str(returnEvalue(hyperScoreList))
        except:
            continue

