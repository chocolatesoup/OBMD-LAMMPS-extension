import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rc
from sys import argv
from os import listdir

# ----------------------------------------------------------------------------------------------------- #

def split(list,size):
    return [list[i:i+size] for i in range(0, len(list), size)]

# ----------------------------------------------------------------------------------------------------- #

def compute_momentum(file):
    timestep, px, py, pz = [], [], [], []
    totalp = []
    for lines in file:
        splittedline = lines.split()

        if lines.startswith('#'):
            continue   

        if len(splittedline) < 4:
            continue

        if len(splittedline) == 4:
            timestep.append(float(int(splittedline[0]))) # !
            px.append(float(splittedline[1]))
            py.append(float(splittedline[2]))
            pz.append(float(splittedline[3]))
            totalp.append(np.sqrt(float(splittedline[1])**2 + float(splittedline[2])**2 + float(splittedline[3])**2))

    # totalp_avg = np.average(totalp)
    # print('totalp_avg =', np.average(totalp))     
    
    # plot graph
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    ax.axhline(y=0,ls='--',linewidth=0.8,c='darkred')
    ax.plot(timestep[10:],px[10:],linewidth=1.0,c='k',label='$p_{x}$')
    ax.plot(timestep[10:],py[10:],linewidth=1.0,c='darkgreen',label='$p_{y}$')
    ax.plot(timestep[10:],pz[10:],linewidth=1.0,c='darkorange',label='$p_{z}$')
    ax.set_xlim(timestep[0],timestep[-1])
    ax.set_xlabel('Timestep')
    ax.set_ylabel('$p_{i}\ [U_{M} U_{L} / U_{T}]$')

    ax.legend(loc='best',fontsize='small')
    plt.show()
    plt.close()

    return 0 

# ----------------------------------------------------------------------------------------------------- #

if __name__ == "__main__":

    file = open('momentum_all.out', 'r')
    compute_momentum(file)
