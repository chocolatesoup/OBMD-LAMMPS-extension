#!/usr/bin/env python3
import numpy as np
from sys import argv
# uncomment here and below to plot chart
#import matplotlib
#matplotlib.use('TkAgg')
#import matplotlib.pyplot as plt

Lx = 33.59462486002239
LbufferEnd = 5.0
RbufferStart = 27.0

# ----------------------------------------------------------------------------------------------------- #

def split(list,size):
    return [list[i:i+size] for i in range(0, len(list), size)]

# ----------------------------------------------------------------------------------------------------- #

def compute_ndenprof(file,bins,start,stop):
    data = (bins*[0])
    
    timestep, nds_all, dist_all, nds_avg_all = [], [], [], []
    for lines in file:
        splittedline = lines.split()

        if lines.startswith('#'):
            continue   

        if len(splittedline) < 4: # step & No. of bins
            timestep.append(int(splittedline[0]))

        if len(splittedline) == 4:
            dist = float(splittedline[1])
            nds = float(splittedline[3])
            dist_all.append(dist)
            nds_all.append(nds)
            if (dist > LbufferEnd and dist < RbufferStart):
                nds_avg_all.append(nds)
    
    distances = split(dist_all,bins)
    nds = split(nds_all,bins)

    distancesEDT = distances[start:stop] # SKIP some steps 
    ndsEDT = nds[start:stop]

    dist_avg = np.average(distancesEDT,axis=0)
    nds_avg = np.average(ndsEDT,axis=0)
    nds_err = np.std(ndsEDT,axis=0)
    
    # mean NDS should be around 3.0 (+-0.05) in between buffer regions
    mean_nds = np.mean(nds_avg_all)
    if (abs(mean_nds - 3.0) > 0.05):
        print("NDS is WRONG!")
    else:
        print("NDS is OK!")
     


    # print('np.shape(dist_avg) =', np.shape(dist_avg))
    # print('np.shape(nds_avg) =', np.shape(nds_avg))
    # print('np.shape(nds_err) =', np.shape(nds_err))

    # uncomment to plot chart
    #fig = plt.figure(figsize=(4,3))
    #ax = fig.add_subplot(111)
    #ax.errorbar(dist_avg,nds_avg,yerr=nds_err,capsize=1.0,elinewidth=0.5,ecolor='grey',c='darkblue')
    #ax.plot(dist_avg,nds_avg,linewidth=0.7,color='k')
    #ax.axhline(y=3.0,ls='--',linewidth=0.7,c='k')
    #ax.axvline(x=5.039193729003359,ls='-.',linewidth=0.7,c='darkred')
    #ax.axvline(x=(Lx-5.039193729003359),ls='-.',linewidth=0.7,c='darkred')
    #ax.set_xlim(dist_avg[0],dist_avg[-1])
    #ax.set_xlabel('$L_{x}\ [U_{L}]$')
    #ax.set_ylabel('$\\rho\ [U_{L}^{-3}]$')
    #plt.show()
    #plt.close()


# ----------------------------------------------------------------------------------------------------- #

if __name__ == "__main__":
    ndenprof_bins = int(argv[1]) 
    start = int(argv[2])
    stop = int(argv[3])

    file = open('nden_profile.out', 'r')
    compute_ndenprof(file,ndenprof_bins,start,stop)
