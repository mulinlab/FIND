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


CA = ["DNase-seq", "ATAC-seq"]
HM = ["H2AK5ac", "H2AK9ac", "H2BK120ac", "H2BK12ac", "H2BK15ac", "H2BK20ac", 
      "H2BK5ac", "H3F3A", "H3K14ac", "H3K18ac", "H3K23ac", "H3K23me2", "H3K27ac", 
      "H3K27me3", "H3K36me3", "H3K4ac", "H3K4me1", "H3K4me2", "H3K4me3", "H3K56ac", 
      "H3K79me1", "H3K79me2", "H3K9ac", "H3K9me1", "H3K9me2","H3K9me3", "H3T11ph", 
      "H4K12ac", "H4K20me1", "H4K5ac", "H4K8ac", "H4K91ac"]
TF = ["CTCF", "POLR2A", "RAD21", "SMC3", "EP300", "H2AFZ"]

RM = ["DNase", "H3K4me3", "H3K27ac", "H3K4me1", "H3K36me3", "H3K9me3", 
      "H3K27me3", "H3K9ac", "H3K4me2", "H2AFZ", "H3K79me2", "H4K20me1"]


CA_h = list(map(lambda x: "EpiMap_" + x, CA))
HM_h = list(map(lambda x: "EpiMap_" + x, HM))
TF_h = list(map(lambda x: "EpiMap_" + x, TF))
RM_h = list(map(lambda x: "RoadMap_" + x, RM))


def TDGenome_count(TDTissue_list):
    TDTissue_list = list(np.unique(TDTissue_list))

    try:
        Tissue_count = len(TDTissue_list.remove("."))
    except:
        Tissue_count = len(TDTissue_list)
        
    return(Tissue_count)


def score_change(score_list):
    for i in range(len(score_list)):
        try:
            score_list[i] = round(float(score_list[i]), 6)
        except:
            score_list[i] = "NA"

    return(score_list)


def score_change_reg(score_list):
    score_list_out = []
    for i in range(len(score_list)):
        try:
            score_list_out.append(round(float(score_list[i]), 6))
        except:
            pass
    
    if len(score_list_out) == 0:
        score_list_out = ["NA"]

    return(score_list_out)


def TDGenome_Pro(TDLine_in):
    TDLine_out = []

    for i in range(12):
        temp_col = TDLine_in[i].split(",")
        temp_num = TDGenome_count(temp_col)
        TDLine_out.append(temp_num)

    for i in range(12,16):
        temp_num = max(score_change_reg(TDLine_in[i].split(",")))
        TDLine_out.append(temp_num)

    return(TDLine_out)


def RMark_change(mark_list):
    for i in range(len(mark_list)):
        try:
            mark_list[i] = mark_list[i].split("-")[1]
        except:
            mark_list[i] = "."
            
    return(mark_list)


def RMark_Scoure(MarkLine_all, MarkLine_on, MarkLine_score):
    Mark_new = []
    Score_new = []

    for MarkType in list(np.unique(MarkLine_on)):
        temp_index = [i for i,x in enumerate(MarkLine_on) if x == MarkType]
        Mark_new.append(MarkType)
        Score_new.append(max([MarkLine_score[x] for x in temp_index]))

    Mark_num = ["NA"] * len(MarkLine_all)
        
    for MarkType in MarkLine_all:
        if MarkType in Mark_new:
            index_all = MarkLine_all.index(MarkType)
            index_on = Mark_new.index(MarkType)
            Mark_num[index_all] = Score_new[index_on]

    return(Mark_num)


def EpiData_Pro(EpiLine_in, CA, HM, TF, RM):
    
    HumanCAge_score = score_change(EpiLine_in[0].split(","))
    HumanCAge_mark = EpiLine_in[1].split(",")
    HumanHMge_score = score_change(EpiLine_in[2].split(","))
    HumanHMge_mark = EpiLine_in[3].split(",")
    HumanTFge_score = score_change(EpiLine_in[4].split(","))
    HumanTFge_mark = EpiLine_in[5].split(",")

    RoadMap_score = score_change(EpiLine_in[6].split(","))
    RoadMap_cellMark = RMark_change(EpiLine_in[7].split(","))

    CA_num = RMark_Scoure(CA, HumanCAge_mark, HumanCAge_score)
    HM_num = RMark_Scoure(HM, HumanHMge_mark, HumanHMge_score)
    TF_num = RMark_Scoure(TF, HumanTFge_mark, HumanTFge_score)
    RM_num = RMark_Scoure(RM, RoadMap_cellMark, RoadMap_score)

    EpiLine_out = CA_num + HM_num + TF_num + RM_num
    
    return(EpiLine_out)


def regBase_Pro(regLine_in, colnms_in):
    
    for ColNum in range(len(colnms_in)):
        regLine_in[ColNum] = max(score_change_reg(regLine_in[ColNum].split(",")))

    return(regLine_in)

### DATA PROCESSING ###
for line in stdin:
    line = line.strip("\n")
    # detect header line
    if line.startswith('#') or line.startswith('CHROM'):
        columnNames = line.strip('#').split('\t')
        CADD_colnm = columnNames[0:117]
        reg_colnm = columnNames[117:158] + columnNames[159:223] + columnNames[247:264]
        TDLine_colnm = columnNames[231:235] + columnNames[236:242] + columnNames[245:247] + \
                        [columnNames[235]] + columnNames[242:245]

        newColumnNames = CADD_colnm + CA_h + HM_h + TF_h + RM_h + reg_colnm + TDLine_colnm
        stdout.write('\t'.join(newColumnNames) + '\n')

        continue

    line = line.split('\t')

    CADDLine = line[0:117]
    regLine = line[117:158] + line[159:223] + line[247:264]
    EpiLine = line[223:231]
    TDLine_in = line[231:235] + line[236:242] + line[245:247] + [line[235]] + line[242:245]

    Epi_outl = EpiData_Pro(EpiLine, CA, HM, TF, RM)
    TD_outl = TDGenome_Pro(TDLine_in)
    reg_outl = regBase_Pro(regLine, reg_colnm)

    out_line = CADDLine + Epi_outl + reg_outl + TD_outl
    out_line = [str(x) for x in out_line]
    stdout.write('\t'.join(out_line) + '\n')
    
