
#Doesn't actually remove isotopes, just removes peaks within 1 dalton of each other
def isotopeRemove (iM2array):
    output = []
    for i in xrange(len(iM2array)):
        mid =[iM2array[i]]
        for j in xrange(i, len(iM2array)):
            if abs(iM2array[i][0]-iM2array[j][0]) < 1:
                mid.append(iM2array[j])
            else:
                output.append(max(mid,key=lambda x: x[1]))
                break
    return output

#Normalises to 0,1 range peaks in input
def normalise (iM2array):
    output = []
    iMax = 0
    iMin = max
    for i in xrange(len(iM2array)):
        if iM2array[i][1] > iMax:
            iMax = iM2array[i][1]
        if iM2array[i][1] < iMin:
            iMin = iM2array[i][1]
    for i in xrange(len(iM2array)):
        output.append((iM2array[i][0],((iM2array[i][1]-iMin)/(iMax-iMin))))        
    return output

def removeNoise (iM2array, threshold):
    output = []
    for i in xrange(len(iM2array)):
        if iM2array[i][1] > threshold:
            output.append(iM2array[i])
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
            outDict["m2Peaks"] = removeNoise(normalise(isotopeRemove(outDict["m2Peaks"])), 0.05)
            yield outDict
            outDict = {"m2Peaks":[]}
