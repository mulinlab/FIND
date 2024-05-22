import pandas as pd
import numpy as np
import configparser
import math
import sys
import random
import gzip

from argparse import ArgumentParser
from columnInfo import get_columns


parser = ArgumentParser(description="%prog name")
parser.add_argument("-i", "--input", dest="input", type=str,
                    default=None,
                    help="File location of variant data in tsv format. If not specified use stdin")
parser.add_argument("-o", "--output", dest="output", type=str,
                    default=None,
                    help="Location were the generated file is stored. If not specified use stdout")
parser.add_argument("-c", "--config", dest="config", type=str,
                    default='./impute_GRCh37_v1.6.cfg',
                    help="Config file that specifies used tracks")
parser.add_argument("-b", "--cat2bool", dest="cat2bool", action='store_true',
                    help="Specify whether categories are split into multiple boolean classifier")
parser.add_argument("--noheader", dest="noheader", default=False,
                    action="store_true",
                    help="Do not print header line")
args = parser.parse_args()


# define input and output sources
if args.input is None:
    stdin = sys.stdin
else:
    stdin = open(args.input, 'r')

if args.output is None:
    stdout = sys.stdout
else:
    stdout = open(args.output, 'w')


ImpFile = "./ScoreMetrix_ALL_change.tsv"
ImpData = pd.read_csv(ImpFile, sep="\t", header=0, index_col=0)

# MissZero (0 or mean)
def score_change(score_list, MissValue = 0):
    for i in range(len(score_list)):
        try:
            score_list[i] = float(score_list[i])
        except:
            score_list[i] = MissValue

    return(score_list)


def score_change_TD(score_value, MissValue = 0):
    try:
        score_value = float(score_value)
    except:
        score_value = MissValue

    return(score_value)


def Enhencer_S(enhancerTSS_list):
    if "NA" in enhancerTSS_list:
        enhancerTSS_score = math.log2(158)
    else:
        enhancerTSS_score = [int(x) for x in enhancerTSS_list]
        enhancerTSS_score = math.log2(max(enhancerTSS_score))
        
    return(enhancerTSS_score)


def TDGenome_Pro(TDLine_in):
    TDLine_out = TDLine_in[0:12]

    TDLine_out.append(round(math.log2(score_change_TD(TDLine_in[12], MissValue = 158)), 6))
    TDLine_out.append(score_change_TD(TDLine_in[13], MissValue = -0.04717))
    TDLine_out.append(score_change_TD(TDLine_in[14], MissValue = 2.09427))
    TDLine_out.append(score_change_TD(TDLine_in[15], MissValue = 0))

    return(TDLine_out)


def ColValue_split(Value_list):
    for i in range(len(Value_list)):
        try:
            Value_list[i] = Value_list[i].split(",")[0]
            if Value_list[i] == ".":
                Value_list[i] = "NA"
        except:
            Value_list[i] = "NA"
            
    return(Value_list)


def regBase_Pro(regLine_in, colnms_in, ImpData):
    for ColNum in range(len(colnms_in)):
        colnm = colnms_in[ColNum]
        # MissValue = ImpData.loc[colnm][0]
        MissValue = ImpData.loc[colnm, "Mean"]
        regLine_in[ColNum] = max(score_change(regLine_in[ColNum].split(","), MissValue))
    
    return(regLine_in)


def EpiData_Pro(EpiLine_in):
    EpiLine_in = list(map(lambda x: float(x.replace("NA", "0")), EpiLine_in))
    EpiLine_in = list(map(lambda x: round(math.log10(x+1), 6), EpiLine_in))
    
    return(EpiLine_in)


def CADD_Pro(CADDL_in, colnms_in):
    columnNames = list(map(str.lower, colnms_in))

    # associate fields with column names
    CADDL_in = ColValue_split(CADDL_in)
    fieldsDict = dict(zip(columnNames, CADDL_in))
    outFields = []
    indicatorFields = []

    for trackName, status in config['Tracks']:
        if 'colname' in trackData[trackName].keys():
            trackName = trackData[trackName]['colname']
        track = trackData[trackName]

        if status == 'Ignore':
            continue
        if status != 'True':
            continue

        if track['type'] == 'combined':
            if args.cat2bool:
                i = trackData[track['base']]['id']
                baseArray = np.array(outFields[i:i+len(trackData[track['base']]['categories'])])
            else:
                baseValue = outFields[trackData[track['base']]['id']]
                baseArray = (np.array(trackData[track['base']]['categories']) == baseValue).astype(int)
            values = []
            for child in track['child']:
                values.extend(baseArray * outFields[trackData[child]['id']])
            outFields.extend([value for value in values])
        else:  
            try:
                if 'derive' in track.keys():
                    value = track['derive'](fieldsDict)
                else:
                    value = fieldsDict[trackName]
                if track['type'] in [float, int]:
                    value = track['type'](value)
                if 'transformation' in track.keys(): # transform is slightly redundant to derive
                    value = track['transformation'](value)
                if track['type'] is list:
                    assert(value in track['categories'])
                if 'indicator' in track.keys():
                    indicatorFields.append('0')
            except:
                value = track['na_value']
                if 'indicator' in track.keys():
                    indicatorFields.append('1')

            if args.cat2bool and track['type'] is list:
                values = (np.array(track['categories']) == value).astype(int)
                outFields.extend(values)
            else:
                outFields.append(value)

    # minimize zeros and stringify
    outFields = ['0' if f == 0 else str(f) for f in outFields]
    outFields.extend(indicatorFields)
    
    return(outFields)


### TRACK PREPARATION ###
newColumnNames, _, trackData, config = \
    get_columns(args.config, None, args.cat2bool, False)


### DATA PROCESSING ###
for line in stdin:
    line = line.strip("\n")
    if line.startswith('#') or line.startswith('CHROM'):
        columnNames = line.strip('#').split('\t')
        CADD_colnm = columnNames[0:117]
        Epi_colnm = columnNames[117:169]
        reg_colnm = columnNames[169:291]
        TDLine_colnm = columnNames[291:307]
        newColumnNames = newColumnNames + Epi_colnm + reg_colnm + TDLine_colnm
        if not args.noheader:
            stdout.write('\t'.join(newColumnNames) + '\n')

        continue
        
    # associate fields with column names
    line = line.split('\t')
    
    CADDLine = line[0:117]
    EpiLine = line[117:169]
    regLine = line[169:291]
    TDLine_in = line[291:307]

    try:
        CADD_outl = CADD_Pro(CADDLine, CADD_colnm)
        
    except:
        continue

    else:
        Epi_outl = EpiData_Pro(EpiLine)
        reg_outl = regBase_Pro(regLine, reg_colnm, ImpData)
        TD_outl = TDGenome_Pro(TDLine_in)

        out_line = CADD_outl + Epi_outl + reg_outl + TD_outl
        out_line = [str(x) for x in out_line]
        stdout.write('\t'.join(out_line) + '\n')

