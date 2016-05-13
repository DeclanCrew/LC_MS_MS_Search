'''
MGFParser, contains functions for extracting and processing spectral
data from Mascot Generic Format files.
'''
import numpy as np

def normalise(iM2array):
    '''Normalises input peak array to maximum intensity of 100'''
    iMax = np.amax(iM2array, axis=0)[1]
    iM2array[:, 1] *= 100/iMax
    return iM2array

def simplify_ions(iM2array, threshold):
    '''Coalesces peaks to half dalton ranges, removes peaks with
       Intensities below threshhold'''
    output = {}
    for i in iM2array:
        toll_bin = int(i[0])
        if toll_bin in output:
            output[toll_bin] += int(i[1])
        else:
            output[toll_bin] = int(i[1])
    filtered = [[key, value] for key, value in output.iteritems()
                if value >= threshold]
    return sorted(filtered, key=lambda t: t[0])


def MGFReader(mgf_file, conf):
    '''Generates peak dictionaries from entries in mgf_file if in thresholds'''
    out_dict = {"m2Peaks":[]}
    proton_mass = float(conf["other_constants"]["H+"])
    noise_threshold = conf["spectrum_options"]["noiseThreshold"]
    max_charge = conf["spectrum_options"]["maxCharge"]
    min_mass = conf["spectrum_options"]["minMass"]
    min_length = conf["spectrum_options"]["minLength"]
    while True:
        line = mgf_file.readline()
        if not line:
            break
        delimited = line.split(" ")
        if "BEGIN" == line[:5]:
            continue
        elif "TITLE" == line[:5]:
            out_dict["name"] = delimited[0][6:]
        elif "PEPMA" == line[:5]:
            out_dict["mass"] = float(delimited[0][8:])
        elif "CHARG" == line[:5]:
            out_dict["charge"] = int(line[7])
        elif "RTINS" == line[:5]:
            out_dict["name"] += "-"+line[12:-2]
        elif "END" == line[:3]:
            out_dict["trueMass"] = float((out_dict["mass"]
                                         - proton_mass)
                                         * out_dict["charge"])
            out_dict["m2Peaks"] = simplify_ions(normalise(np.array(
                out_dict["m2Peaks"])), noise_threshold)
            if out_dict["charge"] > max_charge or out_dict["trueMass"] < min_mass \
               or len(out_dict["m2Peaks"]) < min_length:
                out_dict = {"m2Peaks": []}
                continue
            yield out_dict
            out_dict = {"m2Peaks":[]}
        elif len(delimited) == 2:
            peak_bin = int(round(float(delimited[0])/0.5, 0)*5)
            if len(out_dict["m2Peaks"]) > 0 and peak_bin == out_dict["m2Peaks"][-1][0]:
                out_dict["m2Peaks"][-1][1] += float(delimited[1])
            else:
                out_dict["m2Peaks"].append([peak_bin, float(delimited[1])])
