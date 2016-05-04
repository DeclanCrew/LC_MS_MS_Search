import ConfigParser

#Reads amino acid masses from tab-delimited file to dictionary, for use as reference
def readAAMasses (iFile, index):
    iData = open(iFile, "rb").readlines()
    oData = {}
    for line in iData:
        entry = line.split("\t")
        oData[str(entry[0])] = float(str(entry[index]).rstrip("\n\r"))
    return oData

def readConfigs (configFile):
    config = ConfigParser.ConfigParser()
    config.optionxform = str
    config.read(configFile)
    output = {}
    units= {"Daltons": 1, "ppm": 0.000001}
    isotope_columns= {"Monoisotopic": 1, "Mixed": 2}
    for i in config.sections():
        output[i] = {}
        for j in config.items(i):
            if j[0][:3] == "int":
                output[i][j[0][4:]]= int(j[1])
            elif i[-4:] == "ptms":
                output[i][j[0]]= float(j[1])
            elif j[0][-9:] == "tolerance":
                output[i][j[0]]= float(j[1])
            else:
                output[i][j[0]]= j[1]
    output["AAMassRef"]= readAAMasses (output["amino_acid_options"]["aa_mass_file"],
                                       isotope_columns[output["amino_acid_options"]["isotope_type"]])    
    for i in output["search_options"]:
        if i[-6:] == "_units":
            output["search_options"][i[:-6]] *= units[output["search_options"][i]]
    for i in output["constitutive_ptms"]:
        output["AAMassRef"][i.upper()] += output["constitutive_ptms"][i]
    for i in output["other_constants"]:
        output["other_constants"][i] = float(output["other_constants"][i])
    return output
