from pytorch_tabnet.tab_model import TabNetClassifier
from pytorch_tabnet.augmentations import ClassificationSMOTE
from sklearn.preprocessing import LabelEncoder

import pandas as pd
import numpy as np
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", type=str, default=None)
parser.add_argument("-o", "--output", dest="output", type=str,default=None)
args = parser.parse_args()


def data_pro(Data_in, featuse):
    target = "class"
    unused_feat = ["chrom", "begin", "end", "ref", "alt", "id", "gene", 
                   "type", "length", "istv", "aacomb", "Set"]

    unused_feat = unused_feat + [target]
    features = [ col for col in Data_in.columns if col not in unused_feat ] 

    X_test = Data_in[features].values
    y_test = Data_in[target].values

    return([X_test, y_test])


Data_in = pd.read_csv(args.input, sep = "\t", header = 0)

OHcol = ['consequence', 'domain', 'dst2spltype', 'oaa', 'naa', 
            'polyphencat', 'siftcat', 'segway', 'nuccomb']
for col in OHcol:
    Data_in[col] = Data_in[col].astype("str")
    catfile = f"./FeatureOH/{col}_EnCode.tsv"
    catdict = pd.read_csv(catfile, sep = "\t", header = 0)
    catdict["CodeNumb"] = catdict["CodeNumb"].astype(int)
    catdict = dict(zip(catdict["ClassName"], catdict["CodeNumb"]))
    Data_in[col] = Data_in[col].apply(lambda x:catdict.get(x))

X_test, y_test = data_pro(Data_in, "FeaALL")

Model_in = f"../../trained_model/FIND_model.zip"
clf = TabNetClassifier()
clf.load_model(Model_in)
pred_proba = clf.predict_proba(X_test)

pd.DataFrame(pred_proba).to_csv(args.output, index = False, sep = "\t")