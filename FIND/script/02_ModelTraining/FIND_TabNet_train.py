from pytorch_tabnet.tab_model import TabNetClassifier
from pytorch_tabnet.pretraining import TabNetPretrainer
from pytorch_tabnet.augmentations import ClassificationSMOTE

import torch

from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
from sklearn import metrics

import pandas as pd
import numpy as np

from pathlib import Path
from matplotlib import pyplot as plt
import os
import gzip

import argparse

np.random.seed(0)



parser = argparse.ArgumentParser()
parser.add_argument("--ModelP", help="Training data file path")
parser.add_argument("--ModelF", help="Training data file name")
parser.add_argument("--ModelO", help="Model putput path")

args = parser.parse_args()


Monm = args.ModelF.split("_")[0]
## Model Out
Model_out = os.path.join(args.ModelO, Monm + "_model")
## Model Data IN
file_in = os.path.join(args.ModelP, args.ModelF)

train = pd.read_csv(file_in, header = 0, sep = "\t")
train = train.sample(frac=1).reset_index(drop=True)
n_total = len(train)

target = "class"

if "Set" not in train.columns:
    train["Set"] = np.random.choice(["train", "valid"], p =[.9, .1], size=(train.shape[0],))

train_indices = train[train.Set=="train"].index
valid_indices = train[train.Set=="valid"].index

categorical_columns = []
categorical_dims =  {}
unused_feat = ["chrom", "begin", "end", "ref", "alt", "id", "gene", 
               "type", "length", "istv", "aacomb", "Set"]

OHcol = ['consequence', 'domain', 'dst2spltype', 'oaa', 'naa', 
         'polyphencat', 'siftcat', 'segway', 'nuccomb']

for col in OHcol:
    train[col] = train[col].astype("str")
    print(col, train[col].nunique())

    catfile = f"./FeatureOH/{col}_EnCode.tsv"
    catdict = pd.read_csv(catfile, sep = "\t", header = 0)
    catdict["CodeNumb"] = catdict["CodeNumb"].astype(int)
    catdict = dict(zip(catdict["ClassName"], catdict["CodeNumb"]))

    train[col] = train[col].apply(lambda x:catdict.get(x))
    
    categorical_columns.append(col)
    categorical_dims[col] = len(catdict)


unused_feat = unused_feat + [target]

features = [ col for col in train.columns if col not in unused_feat] 
cat_idxs = [ i for i, f in enumerate(features) if f in categorical_columns]
cat_dims = [ categorical_dims[f] for i, f in enumerate(features) if f in categorical_columns]

if os.getenv("CI", False):
# Take only a subsample to run CI
    X_train = train[features].values[train_indices][:1000,:]
    y_train = train[target].values[train_indices][:1000]
else:
    X_train = train[features].values[train_indices]
    y_train = train[target].values[train_indices]

X_valid = train[features].values[valid_indices]
y_valid = train[target].values[valid_indices]

## TabNet Model
clf = TabNetClassifier(
    n_d=64, n_a=64, n_steps=5,
    gamma=1.5, n_independent=2, n_shared=2,
    cat_idxs=cat_idxs,
    cat_dims=cat_dims,
    cat_emb_dim=1,
    lambda_sparse=1e-4, momentum=0.3, clip_value=2.,
    optimizer_fn=torch.optim.Adam,
    optimizer_params=dict(lr=2e-2),
    scheduler_params = {"gamma": 0.95,
                     "step_size": 20},
    scheduler_fn=torch.optim.lr_scheduler.StepLR, epsilon=1e-15
)

max_epochs = 400 if not os.getenv("CI", False) else 2

clf.fit(
    X_train=X_train, y_train=y_train,
    eval_set=[(X_train, y_train), (X_valid, y_valid)],
    eval_name=['train', 'valid'],
    max_epochs=max_epochs, patience=32,
    batch_size=1024 * 4, 
    virtual_batch_size=256
)

# To get final results you may need to use a mapping for classes 
# as you are allowed to use targets like ["yes", "no", "maybe", "I don't know"]
preds_mapper = { idx : class_name for idx, class_name in enumerate(clf.classes_)}
print(f"BEST VALID SCORE FOR {Monm} : {clf.best_cost}")

saved_filename = clf.save_model(Model_out)
