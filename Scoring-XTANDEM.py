import sys, getopt
from multiprocessing import Pool
from functools import partial
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as path
import ast

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
        i["bIons"]= ast.literal_eval(i["bIons"])
        i["yIons"]= ast.literal_eval(i["yIons"])
        odata.append(i)
    return odata
             
#returns hyper scores for matches in data at the M2 level
def returnHyperScore (iPeakList, norm, m1Match, tolerance):
    score = 0
    bCount = 0
    yCount = 0
    for iPeak in iPeakList:
        for bIon in m1Match["bIons"]:
            if abs(iPeak[0] - float(bIon)) < tolerance:
                m1Match["bIons"].remove(bIon)
                score +=iPeak[1]
                bCount +=1
        for yIon in m1Match["yIons"]:
            if abs(iPeak[0] - float(yIon)) < tolerance:
                m1Match["yIons"].remove(yIon)
                score +=iPeak[1]
                yCount +=1
    hyperScore = (score/norm)*(np.math.factorial(bCount))*(np.math.factorial(yCount))
    return hyperScore

#returns sublist of peptide dataset that has correct mass
def MatchMasses (searchIter, mass, tolerance):
    output = []
    for searchTerm in searchIter:
        if abs(mass - float(searchTerm["Mass"])) < tolerance:
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
            outDict["m2Peaks"].append(list((float(line.split(" ")[0]),float(line.split(" ")[1]))))
        except:
            next     
        if "END" in line:
            outDict["trueMass"] = float((outDict["mass"]*outDict["charge"])-outDict["charge"])
            m2norm = np.linalg.norm([i[0] for i in outDict["m2Peaks"]])
            outDict["norm"] = m2norm
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
    
def Make_Histogram (Input_array, Set_name):
    fig, ax = plt.subplots()
    n, bins = np.histogram(Input_array, 100)
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
print peakList

for peak in MGFReader(MGFFile):
    if peak["name"] in peakList:
        print ("[+]Searching %s" % (peak["name"]))
        m1Matches = MatchMasses(peptideDB, peak["trueMass"], float(0.1))
        hyperScoreList = []
        maxScore = 0
        maxEntry = dict()
        for i in m1Matches:
            proteinSet.add(i["protein"])
            hS = returnHyperScore(peak["m2Peaks"],peak["norm"], i, 0.5)
            hyperScoreList.append(hS)
            if hS > maxScore:
                maxScore = hS
                maxEntry = i
        print maxEntry
        Make_Histogram(hyperScoreList, peak["name"])
    else:
        break
   

