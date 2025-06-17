# create and read input data 
# use reduced units
# system: DPD fluid

import os
import sys
import numpy as np
from random import randint

# -------------------------------------------------------- #
def compute_N(rho, Lx, Ly, Lz):
    N = int(rho*Lx*Ly*Lz)
    return N
# -------------------------------------------------------- #

# simulation box   
xlo = 0.0
ylo = 0.0
zlo = 0.0                                                                                                     
xhi = 33.594
yhi = 11.198
zhi = 11.198

# model details
Nm = 8 
rho = 3
rc = 1.0
gamma_dpd = 4.5
seed_dpd = randint(0,pow(10,4))
temp = 1.0
aij = 209.6 
NDPD = compute_N(rho, xhi, yhi, zhi)

# OBMD settings
usher_flag = 1
charge_flag = 0

seed = randint(0,pow(10,4))

pxx = 188.0
pxy = 0.0
pxz = 0.0
dpxx = 0.0
freq = 0.0
alpha = 0.7
tau = 0.005
nbuf = 1327
etarget = 31.03
ds0 = 1.0
dtheta = 0.02 # here not important
uovlp = pow(10,4)
dsolvp = 1.5
eps = 1.0
nattempt = 40

buffer_size = 0.15 * xhi
shear_size = 0.0
gfac = 0.25
offset = 0.0
step_parallel = 0
step_perp = 1
maxattempt = 1

# simulation details
skin = 0.4
dt = 0.001464
steps = int(20 * pow(10,5)) 
out = 1000

# --------------------------------------------------------------------------------------------- #
# write LAMMPS input script

in_file = "in.simulation"

# -------------------------------------------------------- #
def write_in(in_file):

    file_to_write = open(in_file, 'w')

    content_4in = f"""
# ----------------- Variables Section -----------------
units           lj
boundary        f p p
atom_style      atomic
comm_modify     vel yes
newton          on

region          leftB block {xlo} {buffer_size} {ylo} {yhi} {zlo} {zhi}
region          rightB block {xhi-buffer_size} {xhi} {ylo} {yhi} {zlo} {zhi}
region          leftshear block {0.0} {0.0} {0.0} {0.0} {0.0} {0.0} 
region          rightshear block {0.0} {0.0} {0.0} {0.0} {0.0} {0.0} 
region          leftBin block {offset} {buffer_size-offset} {offset} {yhi-offset} {offset} {zhi-offset}
region          rightBin block {xhi-buffer_size+offset} {xhi-offset} {offset} {yhi-offset} {offset} {zhi-offset}
region          roi block {buffer_size} {xhi-buffer_size} {ylo} {yhi} {zlo} {zhi}

# ----------------- Interaction Section -----------------
pair_style      dpd {temp} {rc} {seed_dpd}

# ----------------- Atom Definition Section -----------------
read_data       dpd_8map_obmd.data

# ----------------- Settings Section -----------------
pair_coeff      * * {aij} {gamma_dpd} {rc}

neighbor        {skin} bin
neigh_modify    delay 0 every 1

timestep        {dt}

fix             1 all nve
fix             2 all obmd 1 1 {seed} {pxx} {pxy} {pxz} {dpxx} {freq} {alpha} {tau} {nbuf} &
                region1 leftB region2 rightB region3 leftshear & 
                region4 rightshear region5 leftBin region6 rightBin &
                buffersize {buffer_size} gfac {gfac} stepparallel {step_parallel} stepperp {step_perp} &
                maxattempt {maxattempt} usher {usher_flag} {etarget} {ds0} {dtheta} {uovlp} {dsolvp} {eps} {nattempt} charged {charge_flag}

# ----------------- Compute Section -----------------

#

# ----------------- Output Section -----------------
thermo          {out}
thermo_style    custom step temp 
run             {steps}
    """
    file_to_write.write(content_4in)
    file_to_write.close()

# -------------------------------------------------------- #

# write input script
if __name__ == "__main__":
     write_in(in_file)
     print(f"file {in_file} is written")
