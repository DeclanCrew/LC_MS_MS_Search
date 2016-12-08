'''
MascotStyle, contains the main loop and the core scoring components of the
search engine
'''
import MGFParse
import PeptideDB
import configs
import scorer
import csv

def match_masses(search_iter, mass, tolerance):
    '''Returns the subset of the searchDB that has a mass within tolerance'''
    output = []
    for search_term in search_iter[int(mass)]:
        if float(search_term["mass"]) < (mass - tolerance):
            continue
        elif float(search_term["mass"]) < (mass + tolerance):
            search_term["delta"] = mass - float(search_term["mass"])
            search_term["bCount"] = 0
            search_term["bSum"] = 0
            search_term["yCount"] = 0
            search_term["ySum"] = 0
            output.append(search_term)
        else:
            return output

def grab_ions(match_list, conf):
    '''Takes peptides and returns a sorted list of their predicted ions'''
    output = []
    for match in enumerate(match_list):
        consts = conf["other_constants"]
        masses = match[1]["orderedMasses"]
        bIons = PeptideDB.coarsen(PeptideDB.returnIons(masses, consts["B+"]))
        yIons = PeptideDB.coarsen(PeptideDB.returnIons(masses[::-1], consts["Y+"]))
        for bIon in bIons:
            output.append([int(bIon), int((match[0]*2)+1)])
        for yIon in yIons:
            output.append([int(yIon), int(match[0])*2])
    return sorted(output, key=lambda entry: entry[0])

def count_matches(match_list, spectra, conf):
    '''Updates peptides with counts and intensity of matches against spectra'''
    if match_list == None:
        return []
    ions = grab_ions(match_list, conf)
    match_index = 0
    peptides = match_list
    for peak in spectra["m2Peaks"]:
        if match_index == len(ions):
            break
        while ions[match_index][0] < peak[0] and match_index < (len(ions)-1):
            match_index += 1
        if match_index == len(ions):
            break
        if ions[match_index][0] > peak[0]:
            continue
        while ions[match_index][0] == peak[0]:
            if ions[match_index][1] % 2:
                peptides[ions[match_index][1]/2]["yCount"] += 1
                peptides[ions[match_index][1]/2]["ySum"] += int(peak[1])
            else:
                peptides[(ions[match_index][1]+1)/2]["bCount"] += 1
                peptides[(ions[match_index][1]+1)/2]["bSum"] += int(peak[1])
            match_index += 1
            if match_index == len(ions):
                break
    output = []
    for entry in peptides:
        entry["spec"] = spectra["name"]
        output.append(entry)
    return output

configurations = configs.readConfigs("mascotStyle.cfg")
file_confs = configurations["data"]

peptide_set = PeptideDB.returnPeptides(configurations)
print "[+]"+str(len(peptide_set))+" mass entries generated."

spectra_file = file_confs["spectra_file"]
MGFFile = open(spectra_file, "rb")
spectra_gen = MGFParse.MGFReader(MGFFile, configurations)

score_mod = scorer.Score(configurations)

top_scores = []

for spectrum in spectra_gen:
    print ("[+]Searching %s. \r" % (spectrum["name"])),
    tol = configurations["search_options"]["initial_tolerance"]*spectrum["trueMass"]
    m1Matches = match_masses(peptide_set, spectrum["trueMass"], tol)
    counter = count_matches(m1Matches, spectrum, configurations)
    score = score_mod.score_series(counter)
    if score:
        top_scores.append(score)

print ""
print "[+]Generating protein level scores."

clear_pep_set = set()
protein_scores = {}
for i in sorted(top_scores, key=lambda entry: entry["mainscore"], reverse=True):
    if i["peptide"] in clear_pep_set:
        continue
    clear_pep_set.add(i["peptide"])
    for j in i["proteins"]:
        if j in protein_scores:
            protein_scores[j] += i["mainscore"]
        else:
            protein_scores[j] = i["mainscore"]
prot_file = open(file_confs["prot_scores_file"], "wb")
for protein in sorted(protein_scores.items(), key=lambda prot: prot[1], reverse=True):
    prot_file.write(str(protein[0]+"\t"+str(protein[1])+"\n"))
