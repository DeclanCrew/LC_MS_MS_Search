from collections import OrderedDict

#Normalises to up to 100 range peaks in input
def normalise (iM2array):
    output = []
    iMax = max(iM2array, key=lambda entry: entry[1])[1]
    append = output.append
    for i in iM2array:
        append((i[0],((i[1]/iMax)*100)))        
    return output

def coarsen (invalue, tol=0.5):
    return int(round(invalue/tol, 0)*tol*10)

def simplifyIons (iM2array, threshold):
    output = OrderedDict()
    iMax = 0
    for i in iM2array:
        tollBin = coarsen(i[0])
        if tollBin in output:
            output[tollBin] += int(i[1])
        else:
            output[tollBin] = int(i[1])
    return {key: value for key, value in output.iteritems() if value >= threshold}

#Returns Peak dictionaries for entries in MGFfile
def MGFReader (mgfFile, conf):
    outDict = {"m2Peaks":[]}
    protonMass = float(conf["other_constants"]["H+"])
    noiseThreshold = conf["spectrum_options"]["noiseThreshold"]
    maxCharge = conf["spectrum_options"]["maxCharge"]
    minMass = conf["spectrum_options"]["minMass"]
    minLength = conf["spectrum_options"]["minLength"]
    while True:
        line = mgfFile.readline()
        if not line:
            break
        delimited = line.split(" ")
        if "BEGIN" == line[:5]:
            continue
        elif "TITLE" == line[:5]:
            outDict["name"] = delimited[0][6:]
        elif "PEPMA" == line[:5]:
            outDict["mass"] = float(delimited[0][8:])
        elif "CHARG" == line[:5]:
            outDict["charge"] = int(line[7])
        elif "RTINS" == line[:5]:
            outDict["name"] += "-"+line[12:-2]
        elif "END" == line[:3]:
            outDict["trueMass"] = float((outDict["mass"]*outDict["charge"])-(outDict["charge"]*protonMass))
            outDict["m2Peaks"] = simplifyIons(normalise(outDict["m2Peaks"]), noiseThreshold)
            if outDict["charge"] > maxCharge or outDict["trueMass"] < minMass or len(outDict["m2Peaks"]) < minLength:
                outDict = {"m2Peaks":[]}
                continue
            yield outDict
            outDict = {"m2Peaks":[]}
        elif(len(delimited) == 2):
            outDict["m2Peaks"].append(list((float(delimited[0]),float(delimited[1]))))

##
#MGFFile = open("combined.mgf", "rb")
#parser = MGFReader(MGFFile)
#for entry in parser:
#    print entry
##    print findDisplacements(entry["m2Peaks"])
#    raw_input("Continue? /n")
