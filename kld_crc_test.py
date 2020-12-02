import math
import numpy as np
from numpy import array
import csv
import pandas as pd
from est_rel_entro_HJW import *

sampP = pd.read_csv('delfi1_h364_23310.csv')
#sampP = pd.read_csv('delfi1_crc.csv')
#sampP = pd.read_csv('delfi1_ctl198.csv')

#sampQ = pd.read_csv('delfi1_crc198.csv')
#sampQ = pd.read_csv('delfi1_lcc.csv')
#sampQ = pd.read_csv('delfi1_dcc.csv')
#sampQ = pd.read_csv('delfi1_bcc.csv')
#sampQ = pd.read_csv('delfi1_gcc.csv')
#sampQ = pd.read_csv('delfi1_pcc.csv')
#sampQ = pd.read_csv('delfi1_ovc.csv')
sampQ = pd.read_csv('delfi1_c364_14430.csv')

est = est_rel_entro_HJW(sampP, sampQ)
print(est)
estP = pd.DataFrame(est)
estP.to_csv('KLdivergenceCRC_fr364_s14430.csv')
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
