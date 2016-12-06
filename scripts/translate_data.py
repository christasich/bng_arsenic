# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 14:46:34 2015

@author: Jonathan Gilligan
"""

import scipy.io as sio
import csv

mat_data = sio.loadmat('data.mat')
bamwsp = mat_data['data']

with open('bamwsp_full.csv', 'wb') as csvfile:
    bamwsp_writer = csv.writer(csvfile)
    bamwsp_writer.writerow(['geocode','depth_ft','arsenic_ppb'])
    for i in range(bamwsp.shape[0]):
        bamwsp_writer.writerow(bamwsp[i,])
print(i)
