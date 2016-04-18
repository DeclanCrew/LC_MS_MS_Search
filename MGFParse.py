#Normalises to up to 100 range peaks in input
def normalise (iM2array):
    output = []
    iMax = max(iM2array, key=lambda entry: entry[1])[1]
    for i in xrange(len(iM2array)):
        output.append((iM2array[i][0],((iM2array[i][1]/iMax)*100)))        
    return output

def removeNoise (iM2array, threshold):
    iM2array = sorted(iM2array, key=lambda entry: entry[1])[-threshold:]
    return sorted(iM2array, key=lambda entry: entry[0])

def findDisplacements (peaks):
    dispSet = set([28,18,17])
    outPeaks = {"unSorted":[],"ionPeaks":[],"dispPeaks":[]}
    frame = []
    for i in xrange(len(peaks)):
        size = int(round(peaks[i][0]))
        frame.append(size)
        print frame
        while (size - frame[0]) > 28:
            outPeaks["unSorted"].append(peaks[i-len(frame)+1])
            frame.pop(0)
        if (size - frame[0]) < 17:
            continue
        for j in xrange(len(frame)):
            if (size - frame[j]) in dispSet:
                intensityRatio = peaks[i][1]/peaks[i-len(frame)+j+1][1]
                if intensityRatio > 1.0:
                    outPeaks["ionPeaks"].append(peaks[i])
                    outPeaks["dispPeaks"].append(peaks[i-len(frame)+j+1])
        print outPeaks
    return outPeaks

def drawPeaks (peaks, colour):
    return

#returns Peak dictionaries for entries in MGFfile
def MGFReader (mgfFile):
    outDict = {"m2Peaks":[]}
    protonMass = 1.007276466
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
            outDict["m2Peaks"] = removeNoise(normalise(outDict["m2Peaks"]), 50)
            if outDict["charge"] > 4 or outDict["trueMass"] < 500 or len(outDict["m2Peaks"]) < 15:
                outDict = {"m2Peaks":[]}
                continue
            yield outDict
            outDict = {"m2Peaks":[]}
        elif(len(delimited) == 2):
            outDict["m2Peaks"].append(list((float(delimited[0]),float(delimited[1]))))

##
##MGFFile = open("combined.mgf", "rb")
##parser = MGFReader(MGFFile)
##for entry in parser:
##    print entry
##    print findDisplacements(entry["m2Peaks"])
##    raw_input("Continue? /n")
