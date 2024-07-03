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

h5file = "JS187.h5"

filename = os.path.join(output_dir,h5file)
#res = h5py.File(filename,"r+")
res0 = tb.open_file(filename, drive="H5FD_CORE")
res = res0.root

#print(res.keys())
#print(res["all"].keys())

star_mass = 0.593

def calc_limits2(mass_ratio, star_mass, period, colnames,percentile=95):
    
    # Need to split into logarithmic period bins 
    # and take the desired percentile of the masses in each bin
    bins = np.arange(0., 10., .5)
    # Create a subchart for each bin. Find the desired percentile mass in that bin, and the median period

    # Select which stars to include. We want to keep things that were not rejected.
    if isinstance(colnames,str):
        select = np.where(res["all"][colnames][:]==0)[0]
    else:
        select0 = np.ones(len(mass_ratio),bool)
        for colname in colnames:
            select0 = select0 & (res["all"][colname][:]==0)
        select = np.where(select0)[0]

    mass_msun = mass_ratio[select]*star_mass
    logperiod = np.log10(period[select])

    bin_idx = np.digitize(logperiod,bins=bins,right=False)
    
    nbins = len(bins)-1
    per_out, perc = np.zeros(nbins), np.zeros(nbins)
    
    for i in range(nbins):
        in_bin = np.where(bin_idx==i)[0]
        nbin = len(in_bin)
#        print(i,nbin)
        if nbin > 0:
            perc[i] = np.percentile(mass_msun[in_bin], percentile)
            per_out[i] = np.median(period[select][in_bin])
#        print("done",i)

    return per_out, perc


#for val in [0,1]:
#    print(len(np.where(res["all"]["AO Rejected 1"][:]==val)[0]))



#select = np.where((res["all"]["AO Rejected 1"][:]==0))[0]
ax = setup_axes()
for percentile in [90]:
    per_ao, perc_ao = calc_limits2(res["all"]['mass ratio'][:],
                                   star_mass,
                                   res["all"]['period(days)'][:],
                                   ["AO Rejected 1"],percentile=percentile)
    ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
    per_ao, perc_ao = calc_limits2(res["all"]['mass ratio'][:],
                                   star_mass,
                                   res["all"]['period(days)'][:],
                                   ["AO Rejected 1","Gaia Rejected"],percentile=percentile)
    ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
    per_ao, perc_ao = calc_limits2(res["all"]['mass ratio'][:],
                                   star_mass,
                                   res["all"]['period(days)'][:],
                                   ["AO Rejected 1","Gaia Rejected","RUWE Rejected"],percentile=percentile)
    # WHY are they all just showing the same as AO??
    ax.plot(per_ao,np.asarray(perc_ao)*1047.35,'o')
plt.savefig(os.path.join(output_dir,h5file.replace(".h5",".png")))

res0.close()

