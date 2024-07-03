#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, sys, glob

import matplotlib.pyplot as plt
import numpy as np
import h5py
import tables as tb
#import json,yaml

from sample_limits import setup_axes, calc_limits

data_path = os.getenv("DATA_PATH","/Users/douglste/data2/molusc_outputs")
output_dir = os.path.join(data_path,"molusc_outputs/")

h5file = "JS117.h5"

filename = os.path.join(output_dir,h5file)
#res = h5py.File(filename,"r+")
res0 = tb.open_file(filename, drive="H5FD_CORE")
res = res0.root

#print(res.keys())
#print(res["all"].keys())

star_mass = 0.593


"""
try:
    res["all"]['mass (M_Jup)'] = res["all"]['mass ratio'][:]*star_mass*1047.35
except:
    print("mass (M_Jup) array already exists!")
    raise

try:
    res["all"]['logPeriod'] = np.log10(res["all"]['period(days)'])
except:
    print("logPeriod array already exists!")
    
try:
    res["all"]['mass (M_Sun)'] = res["all"]['mass ratio'][:]*star_mass
except:
    print("mass (M_Sun) array already exists!")
#     raise
"""

def calc_limits2(mass,logperiod,select,percentile=95):
    # Need to split into logarithmic period bins 
    # and take the desired percentile of the masses in each bin
    bins = np.arange(0., 10., .5)
    # Create a subchart for each bin. Find the desired percentile mass in that bin, and the median period

    bin_idx = np.digitize(logperiod[:],bins=bins,right=False)
    
    nbins = len(bins)-1
    period, perc = np.zeros(nbins), np.zeros(nbins)
    
    for i in range(nbins):
        in_bin = np.where(bin_idx==i)[0]
        nbin = len(in_bin)
#        print(i,nbin)
        if nbin > 0:
            perc[i] = np.percentile(mass[in_bin], percentile)
            period[i] = np.median(res["all"]['period(days)'][select][in_bin])
#        print("done",i)

    return period, perc


#for val in [0,1]:
#    print(len(np.where(res["all"]["AO Rejected 1"][:]==val)[0]))



select = np.where((res["all"]["AO Rejected 1"][:]==0))[0]
ax = setup_axes()
for percentile in [60,80,90,95]:
    per_ao, perc_ao = calc_limits2(res["all"]['mass (M_Sun)'][:][select],
                                   res["all"]['logPeriod'][:][select],
                                   select,percentile=percentile)
    ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
plt.savefig(os.path.join(output_dir,h5file.replace(".h5",".png")))

res0.close()

