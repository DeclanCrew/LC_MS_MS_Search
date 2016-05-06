'''
MGFParser, contains functions for extracting and processing spectral
data from Mascot Generic Format files.
'''
import numpy as np

def normalise(iM2array):
    '''Normalises input peak array to maximum intensity of 100'''
    iMax = np.amax(iM2array, axis=0)[1]
    iM2array[:,1] *= 100/iMax
    return iM2array

def simplifyIons(iM2array, threshold):
    '''Coalesces peaks to half dalton ranges, removes peaks with
       Intensities below threshhold'''
    output = {}
    for i in iM2array:
        tollBin = int(i[0])
        if tollBin in output:
            output[tollBin] += int(i[1])
        else:
            output[tollBin] = int(i[1])
    filtered = [[key, value] for key, value in output.iteritems()
                if value >= threshold]
    return sorted(filtered, key=lambda t: t[0])


def MGFReader(mgfFile, conf):
    '''Generates peak dictionaries from entries in MGFfile if in thresholds'''
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
            outDict["trueMass"] = float((outDict["mass"]
                                         - protonMass)
                                         * outDict["charge"])
            outDict["m2Peaks"] = simplifyIons(normalise(np.array(
                outDict["m2Peaks"])), noiseThreshold)
            if outDict["charge"] > maxCharge or outDict["trueMass"] < minMass \
               or len(outDict["m2Peaks"]) < minLength:
                outDict = {"m2Peaks": []}
                continue
            yield outDict
            outDict = {"m2Peaks":[]}
        elif len(delimited) == 2:
            peakBin = int(round(float(delimited[0])/0.5, 0)*5)
            if len(outDict["m2Peaks"]) > 0 and peakBin == outDict["m2Peaks"][-1][0]:
                outDict["m2Peaks"][-1][1] += float(delimited[1])
            else:
                outDict["m2Peaks"].append([peakBin, float(delimited[1])])
