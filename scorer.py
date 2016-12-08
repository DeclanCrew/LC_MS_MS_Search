import csv
import numpy as np
from scipy import stats

def generate_writer(data_file_name, keys):
    '''Generates results writer, creating a file with header and then
       returning a writer function for each row that needs writing.'''
    data_file = open(data_file_name, "wb")
    writer_object = csv.DictWriter(data_file, keys, dialect="excel-tab",
                                  extrasaction="ignore")
    writer_object.writer.writerow(keys)
    while True:
        yield writer_object.writerow

class Score:
    def __init__(self, confs):
        self.full = bool(confs["data"]["write_full_scores"])
        if self.full == True:
            full_keys = ["peptide", "spec", "proteins", "delta",
                         "bCount", "bSum", "yCount", "ySum"]
            self.full_writer = generate_writer(confs["data"]["full_scores_file"], full_keys)

        top_keys = ["peptide", "spec", "proteins", "score", "stage2", "delta",
                    "bCount", "bSum", "yCount", "decoy", "mainscore", "ySum"]
        self.top_writer = generate_writer(confs["data"]["top_scores_file"], top_keys)
        self.mascot_bias = 1.9
        self.comp_bias = 20

   
    def decoy_query (self, top):
        for protein in top["proteins"]:
            if "DECOY" not in protein:
                return False
        return True 
        
    def score_method(self, entry):
        '''Simple mascot style scoring system, weights yIons stronger than b'''
        return ((self.mascot_bias*entry["ySum"])+entry["bSum"])*entry["yCount"]*entry["bCount"]

    def secondary_scoring(self, score_series):
        '''XTandem style stage 2 scoring, requires numpy for linear regression'''
        if len(score_series) <= 2:
            return 0
        y, x = np.histogram(score_series[:-1])
        y = np.cumsum(y[::-1])[::-1]
        slope, intercept, r_value, p_value, std_err = stats.linregress(x[1:], y)
        return (score_series[-1]*slope + intercept)*-1

    def composite_score(self, score_record):
        '''Consolidates two scoring methods'''
        return score_record["score"]+(self.comp_bias*score_record["stage2"])

    def score_series(self, series):
        top_score = 0
        top_result = None
        for result in series:
            if self.full:
                next(self.full_writer)(result)
            result["score"] = self.score_method(result)
            if result["score"] > top_score:
                top_score = result["score"]
                top_result = result
        if not top_result:
            return None
        primary_scores = np.log10(np.array(sorted([i["score"] for i  in series if i["score"] > 0])))
        top_result["stage2"] = self.secondary_scoring(primary_scores)
        
        top_result["decoy"] = self.decoy_query(top_result)
        top_result["mainscore"] = self.composite_score(top_result)
        next(self.top_writer)(top_result)
        return top_result



