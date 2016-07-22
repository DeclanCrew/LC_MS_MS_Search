'''
Contains a number of functions for generating virtual spectra from a protein
file, makes use of biopython for protein file parsing.
'''

import re
import itertools
import copy
import numpy as np

def coarsen(inArray, tol=0.5):
    '''Bins ion mass predictions to integers within the closest 0.5.'''
    return (np.around(inArray/tol))*tol*10

def returnIons(peptideValues, adjustment=0):
    '''Returns numpy array of ion masses from AA mass array.'''
    return (np.cumsum(peptideValues))+adjustment

def returnPostTranslationMods(peptideDict, mods):
    '''Takes peptide and finds all modifications that can be applied to it.'''
    output = {}
    for i in mods:
        index = peptideDict["peptide"].find(i[0])
        while index >= 0:
            index += len(i[0])
            output[index] = i[1]
            index = peptideDict["peptide"].find(i[0], index)
    return output

def genModdedPeptide(peptideDict, conf, mods):
    '''Modifies peptide entry to incorporate modifications.'''
    output = copy.deepcopy(peptideDict)
    adjust = []
    for i in xrange(len(output["orderedMasses"])):
        if i in mods:
            adjust.append(mods[i])
        else:
            adjust.append(0)
    output["orderedMasses"] += np.array(adjust)
    pepName = list(output["peptide"])
    for mod in mods:
        pepName.insert(mod, str("[%i]" % (mods[mod])))
    output["peptide"] = "".join(pepName)
    output["mass"] = np.sum(output["orderedMasses"]) + conf["other_constants"]["Mass+"]
    return output

def refine(peptideDict, conf):
    '''Generates modified peptide entries for a given peptide,
       for second pass search.'''
    modifications = conf["variable_ptms"].items()
    validM = returnPostTranslationMods(peptideDict, modifications)
    if len(validM) > 0:
        combos = []
        for length in xrange(conf["search_options"]["ptm_number"]):
            combos.append([{j: validM[j] for j in i}
                                 for i in itertools.combinations(validM, length)])
        for subset in combos:
            for combo in subset:
                yield genModdedPeptide(peptideDict, conf, combo)
    else:
        yield peptideDict

def returnPeptideDict(peptide, proteins, conf):
    '''Returns dictionary of peptide characteristics'''
    output = {"peptide":peptide, "proteins":proteins}
    output["orderedMasses"] = np.array([conf["AAMassRef"][acid]
                                        for acid in cleanPeptide(peptide)])
    output["mass"] = np.sum(output["orderedMasses"]) + conf["other_constants"]["Mass+"]
    return output

def iSilSpectra(protein, regEx, conf):
    '''Splits a protein entry into a series of peptide strings for processing.'''
    output = []
    mProtein = str("n%sc" %(protein.seq))
    peptides = regEx.sub('\n', str(mProtein)).split()
    for x in xrange(conf["search_options"]["maximum_missed_cleavages"]+1):
        mid = ["".join(peptides[i:i+(x+1)]) for i in xrange(len(peptides)-x)]
        for i in mid:
            if len(i) >= conf["search_options"]["min_peptide_length"]:
                output.append((i, protein.name))
    return output

def cleanPeptide(peptide):
    '''Removes sentinel values from a peptide string before mass analysis.'''
    if peptide[0] == "n":
        peptide = peptide[1:]
    if peptide[len(peptide)-1] == "c":
        peptide = peptide[:(len(peptide)-1)]
    return peptide

def proteinPreprocessing(proteins, conf):
    '''Removes proteins with unsupported characters, also generates decoy proteins.'''
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
        if useDecoy:
            decoy = copy.deepcopy(protein)
            decoy.seq = protein.seq[::-1]
            decoy.name = str("DECOY_%s" % (protein.name))
            append(decoy)
    return output

def peptideDatabaseParse():
    try:
        from Bio import SeqIO
        return SeqIO.parse
    except ImportError:
        import fastaParse
        return fastaParse.fastaParser
        
def returnPeptides(conf):
    '''
    Generates dictionary of unique peptide entries from a given reference sequence
    dataset, returns dictionary of mass-sorted peptides, with each key holding all
    peptides with the same dalton mass. Also implements protein processing from confs.
    '''
    protFile = conf["data"]["reference_sequences"]
    proteinFile = open(protFile, "rb")
    protparser = peptideDatabaseParse()
    proteins = protparser(proteinFile, conf["data"]["sequence_format"])
    peptideDB = {}
    print "[+]Digesting proteins in Silico."
    matchingRegex = re.compile(r"(?<=[KR])(?=[^P])")
    for iProtein in proteinPreprocessing(proteins, conf):
        peptides = iSilSpectra(iProtein, matchingRegex, conf)
        for i in peptides:
            try:
                peptideDB[i[0]].append(i[1])
            except KeyError:
                peptideDB[i[0]] = [i[1]]
    print "[+]Generating peptide spectra."
    peptideList = [returnPeptideDict(key, peptideDB[key], conf) for key in peptideDB]
    print "[+]Sorting peptides."
    outHash = {}
    for i in sorted(peptideList, key=lambda entry: entry["mass"]):
        if int(i["mass"]) in outHash:
            outHash[int(i["mass"])].append(i)
        else:
            outHash[int(i["mass"])] = [i]
    return outHash
