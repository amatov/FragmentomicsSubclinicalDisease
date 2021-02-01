
#!/usr/bin/env python

import math
import numpy as np
from numpy import array
import csv
import pandas as pd
from est_rel_entro_HJW import *
#sampP = pd.read_csv('umiseq_pon45_frl356no0.csv')
#sampP = pd.read_csv('umiseq_pon45_frl356.csv')
#sampP = pd.read_csv('umiseq_pon30_frl130.csv')
#sampP = pd.read_csv('delfi2_ctl1_74.csv')
#sampP = pd.read_csv('delfi1_crc364.csv')
#sampP = pd.read_csv('delfi2_ctl1_74_binned10frl.csv')
#sampP = pd.read_csv('delfi1_ctl86HG38.csv')
#sampP = pd.read_csv('delfi1_ctl1Bins.csv')
#sampP = pd.read_csv('delfi2_ctl1Bins.csv')
#sampP = pd.read_csv('delfi1_ctl198.csv')
#sampP = pd.read_csv('delfi2_ctl3_71.csv') # of 131
#sampP = pd.read_csv('delfi2_ctl2_71.csv')
#sampP = pd.read_csv('delfi2_ctl1_74.csv')
#sampP = pd.read_csv('delfi2_ctl1_73.csv')# no outlier
#sampP = pd.read_csv('delfi2_ctl1_fr195.csv')
#sampP = pd.read_csv('delfi2_ctl1_fr195.csv')
#sampP = pd.read_csv('delfi2_ctl73_fr195.csv')
#sampP = pd.read_csv('delfi2_ctl1_fr365.csv')
#sampP = pd.read_csv('delfi1_ctl43_frl499.csv')
#sampP = pd.read_csv('delfi1_ctl43_frl499_189.csv')
#sampP = pd.read_csv('delfi1_healthy43_555bins_fr189.csv')
#sampP = pd.read_csv('delfi2_ctl1_fr137.csv')
#sampP = pd.read_csv('delfi2_ctl1_fr205.csv')
#sampP = pd.read_csv('delfi2_ctl1_fr364.csv')
sampP = pd.read_csv('delfi2_ctl1_74_1M.csv')

#sampQ = pd.read_csv('delfi1_crc27HG38.csv')
#sampQ = pd.read_csv('delfi1_lcc.csv')
#sampQ = pd.read_csv('delfi1_dcc.csv')
#sampQ = pd.read_csv('delfi1_bcc.csv')
#sampQ = pd.read_csv('delfi1_gcc.csv')
#sampQ = pd.read_csv('delfi1_pcc.csv')
#sampQ = pd.read_csv('delfi1_ovc.csv')
#sampQ = pd.read_csv('delfi1_c364_14430.csv')
#sampQ = pd.read_csv('delfi1_ctl.csv') # to test the reverse divergence
#sampQ = pd.read_csv('delfi2_rec_all50.csv')
#sampQ = pd.read_csv('delfi2_col_all79.csv')
#sampQ = pd.read_csv('delfi2_crc_all129.csv')
#sampQ = pd.read_csv('delfi2_col1_alll8.csv')
#sampQ = pd.read_csv('delfi2_col2_all30.csv')
#sampQ = pd.read_csv('delfi2_col3_all18.csv')
#sampQ = pd.read_csv('delfi2_col4_all23.csv')
#sampQ = pd.read_csv('delfi2_col_fr195.csv')
sampQ = pd.read_csv('delfi2_col_all79_1M.csv')
#sampQ = pd.read_csv('delfi2_col79_bin555_189.csv')
#sampQ = pd.read_csv('delfi2_col79_binned10frl.csv')
#sampQ = pd.read_csv('umiseq_pre30_frl130.csv')
#sampQ = pd.read_csv('umiseq_pre56_frl356no0.csv')
#sampQ = pd.read_csv('umiseq_CRpre68_binned10frl.csv')

#list = [17, 18, 19]
#print(list)
#for i in list: #range(13,79):
#file_name = 'delfi2_col_adeH20.csv'
#kld_name = 'KLdivergenceD2col_adeH20_ctl1.csv'
#file_name = 'delfi2_col_adeL28.csv'
#kld_name = 'KLdivergenceD2col_adeL28_ctl1.csv'
#file_name = 'delfi2_rec_adeH11.csv'
#kld_name = 'KLdivergenceD2rec_adeH11_ctl1.csv'
#file_name = 'delfi2_rec_adeL7.csv'
#kld_name = 'KLdivergenceD2rec_adeL7_ctl1.csv'
#file_name = 'delfi2_col1_fr205.csv'
#kld_name = 'KLdivergenceD2col1_fr205_ctl1.csv'
#file_name = 'delfi2_col2_fr137.csv'
#kld_name = 'KLdivergenceD2col2_fr137_ctl1.csv'
#file_name = 'delfi2_col3_fr364.csv'
#kld_name = 'KLdivergenceD2col3_fr364_ctl1.csv'
#file_name = 'delfi2_col4_fr364.csv'
#kld_name = 'KLdivergenceD2col4_fr364_ctl1.csv'
#file_name = 'delfi1_crcBins.csv'
#kld_name = 'KLdivergenceCRC_Bins.csv'
#file_name = 'delfi2_colBins.csv'
#kld_name = 'KLdivergenceD2col_Bins.csv'
#file_name = 'delfi2_ctl1_individual%d.csv'%(i)
#kld_name = 'KLdivergenceD2ctl1_individual%d.csv'%(i)
#i=8
#for i in range(12, 43, 1):
#file_name = 'delfi1_crc4_individual%d.csv'%(i)
#kld_name = 'KLdivergenceD1crc4_individual%d.csv'%(i)
#for i in range(1,79):
#sampQ = pd.read_csv(file_name)
#sampQ = pandas.read_csv('delfi2_col_individual{0}.csv'.format(i))

est = est_rel_entro_HJW(sampP, sampQ)
print(est)
estP = pd.DataFrame(est)
estP.to_csv('KLdivergenceD2_COL79_CTL1_1M.csv')
#estP.to_csv('KLdivergenceD2_COL79_CTL1_binned10frl.csv')
#estP.to_csv('KLdivergenceUMIseqPre56Pon45_frl356no0.csv')
#estP.to_csv('KLdivergenceUMIseqPre56Pon45_frl356.csv')
#estP.to_csv('KLdivergenceUMIseqPre30Pon30_frl130.csv')
#estP.to_csv(kld_name)
#estP.to_csv('KLdivergenceD2_COL79.csv')
#estP.to_csv('KLdivergenceD2val210_D1ctl43_fr189.csv')
#estP.to_csv('KLdivergenceD2col79_D1ctl43_fr189.csv')
#estP.to_csv('KLdivergenceD2col79_D1ctl43.csv')
#estP.to_csv('KLdivergenceD2col79_D1ctl86HG38.csv')
#estP.to_csv('KLdivergenceD2_CRC129_ctl1.csv')
#estP.to_csv('KLdivergenceD2_COL4_23_ctl2.csv')
#estP.to_csv('KLdivergenceD2_COL3_18_ctl2.csv')
#estP.to_csv('KLdivergenceD2_COL2_30_ctl2.csv')
#estP.to_csv('KLdivergenceD2_COLfr195_ctl1.csv')
#estP.to_csv('KLdivergenceD2_COLfr195_ctl73.csv')
#estP.to_csv('KLdivergenceCRC_HG38.csv')
#estP.to_csv('KLdivergenceD2_COLfr365_ctl1.csv')
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
