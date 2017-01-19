'''
MascotStyle, contains the main loop and the core matching component of the
search engine
'''
import MGFParse
import PeptideDB
import configs
import scorer

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
    for index, match in enumerate(match_list):
        if "b" not in match:
            match["b"] = PeptideDB.coarsen(PeptideDB.returnIons(
                match["orderedMasses"], conf["other_constants"]["B+"]))
        if "y" not in match:
            match["y"] = PeptideDB.coarsen(PeptideDB.returnIons(
                match["orderedMasses"][::-1], conf["other_constants"]["Y+"]))
        for bIon in match["b"]:
            output.append([int(bIon), index, "b"])
        for yIon in match["y"]:
            output.append([int(yIon), index, "y"])
    return output

def count_matches(peptides, spectra, conf):
    '''Updates peptides with counts and intensity of matches against spectra'''
    if not peptides:
        return []
    ions = grab_ions(peptides, conf)
    for peak, index, iontype in ions:
        if peak in spectra["m2Peaks"]:
            peptides[index][iontype+"Count"] += 1
            peptides[index][iontype+"Sum"] += spectra["m2Peaks"][peak]
    for peptide in peptides:
        peptide["spec"] = spectra["name"]
    return peptides

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
