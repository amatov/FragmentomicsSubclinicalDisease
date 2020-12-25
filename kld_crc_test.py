import math
import numpy as np
from numpy import array
import csv
import pandas as pd
from est_rel_entro_HJW import *

#sampP = pd.read_csv('delfi1_crc364.csv')
#sampP = pd.read_csv('delfi1_h364_23310.csv')
#sampP = pd.read_csv('delfi1_crc.csv')
#sampP = pd.read_csv('delfi1_ctl198.csv')
#sampP = pd.read_csv('delfi2_ctl3_71.csv') # of 131
#sampP = pd.read_csv('delfi2_ctl2_71.csv')
sampP = pd.read_csv('delfi2_ctl1_74.csv')

#sampQ = pd.read_csv('delfi1_crc198.csv')
#sampQ = pd.read_csv('delfi1_lcc.csv')
#sampQ = pd.read_csv('delfi1_dcc.csv')
#sampQ = pd.read_csv('delfi1_bcc.csv')
#sampQ = pd.read_csv('delfi1_gcc.csv')
#sampQ = pd.read_csv('delfi1_pcc.csv')
#sampQ = pd.read_csv('delfi1_ovc.csv')
#sampQ = pd.read_csv('delfi1_c364_14430.csv')
#sampQ = pd.read_csv('delfi1_ctl.csv') # to test the reverse divergence
#sampQ = pd.read_csv('delfi2_rec_all50.csv')
#sampQ = pd.read_csv('delfi2_col1_all8.csv')
sampQ = pd.read_csv('delfi2_col2_all30.csv')


est = est_rel_entro_HJW(sampP, sampQ)
print(est)
estP = pd.DataFrame(est)
estP.to_csv('KLdivergenceD2_COL2_30_ctl1.csv')
#estP.to_csv('KLdivergenceD2_COL1_8_ctl1.csv')
#estP.to_csv('KLdivergenceD2_REC50_ctl1.csv')
#estP.to_csv('KLdivergenceD2_REC50_ctl2.csv')
#estP.to_csv('KLdivergenceD2_REC50_ctl3.csv')
#estP.to_csv('KLdivergenceCRC_Reverse.csv')
#estP.to_csv('KLdivergenceCRC364.csv')
#estP.to_csv('KLdivergenceCRC.csv')
#estP.to_csv('KLdivergenceOVC.csv')
#estP.to_csv('KLdivergenceCRC_OVC.csv')
#estP.to_csv('KLdivergencePCC.csv')
#estP.to_csv('KLdivergenceCRC_PCC.csv')
#estP.to_csv('KLdivergenceGCC.csv')
#estP.to_csv('KLdivergenceCRC_GCC.csv')
#estP.to_csv('KLdivergenceBCC.csv')
#estP.to_csv('KLdivergenceCRC_BCC.csv')
#estP.to_csv('KLdivergenceDCC.csv')
#estP.to_csv('KLdivergenceCRC_DCC.csv')
#estP.to_csv('KLdivergenceLCC.csv')
#estP.to_csv('KLdivergenceCRC_LCC.csv')
