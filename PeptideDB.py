import re
import itertools
from Bio import SeqIO
import copy
import operator
import numpy as np

#Bins ion mass predictions to the closest 0.5, to enable quicker searching later
def coarsen (inArray, tol=0.5):
    return (np.around(inArray/tol))*tol*10

#Returns list of ion masses
def returnIons (peptideValues, adjustment=0):
    return ((np.cumsum(peptideValues))+adjustment)

def returnPostTranslationMods (peptideDict, mods):
    output = []
    for i in mods:
        index = peptideDict["peptide"].find(i[0])
        while( index >= 0 ):
            index += len(i[0])
            output.append((index, i[1]))
            index = peptideDict["peptide"].find(i[0], index)
    return output

def genModdedPeptide (peptide, mods):
    clean = cleanPeptide(peptide)
    output["mass"]= returnPeptideMass(clean, pKey, float(18.009467553))
    output["bIons"]= returnIons(clean, pKey, float(1.00727))
    output["yIons"]= returnIons(clean[::-1], pKey, float(19.02257))
    outPep = ""
    beginning = 0
    for i in mods:
        substring = peptideDict["peptide"][beginning:i[0]]
        outPep = outPep+substring
        insert = str("[%s]" % (int(round(i[1], 0))))
        outPep = outPep+insert
        beginning = i[0]
        for j in xrange(i[0]):
            output["bIons"][j] += i[1]
            output["yIons"][-(1+j)] += i[1]
        output["mass"] += i[1]
    substring = peptideDict["peptide"][beginning:]
    outPep = outPep+substring
    output["bIons"]= coarsen(returnIons(clean, pKey, float(1.00727)))
    output["yIons"]= coarsen(returnIons(clean[::-1], pKey, float(19.02257)))
    return outPep

def refine (peptide, conf):
    modifications = [(x, conf["variable_ptms"][x]) for x in conf["variable_ptms"]]
    return

#Returns peptideDictionary
def returnPeptideDict (peptide, proteins, conf):
    output = {"peptide":peptide, "proteins":proteins}
    output["orderedMasses"] = np.array([conf["AAMassRef"][acid] for acid in cleanPeptide(peptide)])
    output["mass"]= np.sum(output["orderedMasses"]) + conf["other_constants"]["Mass+"]
    return output

#Generates inSilico peptide objects for input protein SeqRecord
def iSilSpectra (protein, regEx, conf):
    output = []
    mProtein = str("n%sc" %(protein.seq))
    peptides = regEx.sub('\n', str(mProtein)).split()
    for x in xrange(conf["search_options"]["maximum_missed_cleavages"]+1):
        mid = ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
        for i in mid:
            if len(i) >= conf["search_options"]["min_peptide_length"]:
                output.append((i, protein.name))
    return output

def cleanPeptide (peptide):
    if peptide[0] == "n":
        peptide = peptide[1:]
    if peptide[len(peptide)-1] == "c":
        peptide = peptide[:(len(peptide)-1)]
    return peptide

def proteinPreprocessing (proteins, conf, maxX=1, maxB=4):
    output = []
    useDecoy = bool(conf["search_options"]["include_decoy"])
    append = output.append
    for protein in proteins:
        if "X" in protein.seq:
            continue
        elif "B" in protein.seq or "Z" in protein.seq:
            continue
        else:
            append(protein)
        if useDecoy == True:
            decoy = copy.deepcopy(protein)
            decoy.seq = protein.seq[::-1]
            decoy.name = str("DECOY_%s" % (protein.name))
            append(decoy)
    return output
            
def returnPeptides (protFile, conf):
    proteinFile = open(protFile, "rb")
    proteins  = SeqIO.parse(proteinFile, "fasta")
    peptideDB = {}
    print ("[+]Digesting proteins in Silico")
    matchingRegex = re.compile(r"(?<=[KR])(?=[^P])")
    for iProtein in proteinPreprocessing(proteins, conf):
        peptides = iSilSpectra(iProtein, matchingRegex, conf)
        for i in peptides:
            try:
                peptideDB[i[0]].append(i[1])
            except KeyError:
                peptideDB[i[0]] = [i[1]]
    print "[+]Generating peptide spectra"
    peptideList = [returnPeptideDict(key, peptideDB[key], conf) for key in peptideDB]
    print "[+]Sorting peptides"
    outHash = {}
    for i in sorted(peptideList, key=lambda entry: entry["mass"]):
        if int(i["mass"]) in outHash:
            outHash[int(i["mass"])].append(i)
        else:
            outHash[int(i["mass"])] = [i]
    return outHash
