This document serves as a guide for users and developers working with the `obmd` extension in LAMMPS. Therefore, **general knowledge about LAMMPS and molecular dynamics simulations is assumed**. Additionally, basic familiarity with compiling and running C++ code is expected. 

To perform open-boundary molecular dynamics (OBMD) simulations using the `obmd` extension, we provide a **Python script** named `input.py`, which generates `in.simulation` input script for LAMMPS. LAMMPS is then started with this input as a parameter, for example,

```bash
lmp_mpi -in in.simulation
```

The script `input.py` is located [here](https://github.com/chocolatesoup/OBMD-LAMMPS-extension/tree/main/examples/OBMD_DPD).

##### INITIALIZATION SETTINGS

A code snippet of the **function** in the `input.py` Python script, which writes the input script for the LAMMPS OBMD simulation of a **DPD fluid**, is given and explained below.

First, the generated output file name, which serves as the input for LAMMPS (`in_file`), is set. The default value is `in.simulation`. Next, the `write_in(in_file)` function is defined, taking the name of the input file to be created as an argument. Using an f-string, the block of code that forms the input script is stored in `content_4in`. Finally, the content is written to a file using `file_to_write.write(content_4in)`.

```Python
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
region          leftshear block {xlo} {shear_size} {ylo} {yhi} {zlo} {zhi}
region          rightshear block {xhi-shear_size} {xhi} {ylo} {yhi} {zlo} {zhi}
region          leftBin block {offset} {buffer_size-offset} {offset} {yhi-offset} {offset} {zhi-offset}
region          rightBin block {xhi-buffer_size+offset} {xhi-offset} {offset} {yhi-offset} {offset}{zhi-offset}
region          roi block {buffer_size} {xhi-buffer_size} {ylo} {yhi} {zlo} {zhi}
```

OBMD simulations should be performed using **reduced units** (`units lj`) to avoid errors in the calculation of the external force and the style of boundaries (`boundary`) for the global simulation box in each dimension should be set to `f p p`. The latter applies the same style to the lower and upper sides of the simulation box, i.e., in the case of OBMD simulations, the non-periodic style `f` is applied in the $x$-direction, while the box is periodic in $y$- and $z$-directions by choosing the periodic style `p`.  

Following is the explanation of the inputs.
- `atom_style value` 
  The style of atoms in the system. The `obmd` extension was tested using `atomic`, `full`, and `molecular` atom styles.

When performing OBMD simulations, a **momentum-conserving thermostat** should be applied. Usually, the **DPD thermostat** is used, where velocities should be communicated between neighboring processors and stored as properties of ghost atoms using `comm_modify vel yes`. Importantly, Newton's third law for pairwise and bonded interactions should be enabled via `newton on`.

One should also define geometric regions. Here, **6 regions** (namely `leftB`, `rightB`, `leftshear`, `rightshear`, `leftBin`, and `rightBin`) must be specified. All regions must be defined using the following format: `xlo xhi ylo yhi zlo zhi` (for instructions on `block` style, readers are referred to the [LAMMPS documentation](https://docs.lammps.org/region.html)).
- `leftB block args`
  Region defining the **left buffer**, where `Region-ID = leftB`. `buffer_size` is the length of the buffer region, `xlo`, `ylo`, and `zlo` are the lower bounds of the simulation box in the $x$-, $y$-, and $z$-direction, respectively, while `xhi`, `yhi`, and `zhi` are the upper bounds.
- `rightB block args`
  Region defining the **right buffer**, where `Region-ID = rightB`.
- `leftshear block args`
  Region defining **part of the left buffer where shear forces are applied**. Here, `Region-ID = leftshear`. `shear_size` is the length of the region where shear flow is imposed. If only normal pressure forces are applied, e.g., in the case of equilibrium simulations, `shear_size`, `xlo`, `xhi`, `ylo`, `yhi`, `zlo`, and `zhi` should be set to $0.0$.
- `rightshear block args`
  Region defining **part of the right buffer where shear forces are applied**. Here, `Region-ID = rightshear`. If only normal pressure forces are applied, e.g., in the case of equilibrium simulations, `shear_size`, `xlo`, `xhi`, `ylo`, `yhi`, `zlo`, and `zhi` should be set to $0.0$.
- `leftBin block args`
  Region defining **part of the left buffer where new particles are inserted**.  Here, `Region-ID = leftBin`. `offset` is the value used to narrow the insertion region.
- `rightBin block args`
  Region defining **part of  the right buffer where new particles are inserted**. Here, `Region-ID = rightBin`.

##### SIMULATION SETTINGS AND SYSTEM DEFINITION

```Python
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
```
  
- `pair_style dpd` is used to perform simulations employing DPD force field. According to the [LAMMPS documentation](https://docs.lammps.org/pair_dpd.html) it takes temperature (`temp`), cutoff of the interaction (`rc`), and seed (`seed_dpd`) as input parameters. Alongside, the following coefficients (`pair_coeff`) must be given:  
	  - interaction parameter between atoms $i$ and $j$ (`aij`) 
	  - friction parameter (`gamma_dpd`)
	  - cutoff of the interaction (`rc`)
- `read_data` is used to read **pre-equilibrated configuration** of the DPD fluid given in the `dpd_8map_obmd.data` file. The latter contains information about the positions of the atoms in **reduced units**. It is important that the position of the atoms do not exceed $xlo$ or $xhi$ of the simulation box; if so, they should be manually "deleted" before conducting the OBMD simulations. 

```Python
fix             2 all obmd {ntype} {nfreq} {seed} {pxx} {pxy} {pxz} {dpxx} {freq} {alpha} {tau} {nbuf} region1 leftB region2 rightB region3 leftshear region4 rightshear region5 leftBin region6 rightBin buffersize {buffer_size} gfac {gfac} stepparallel {step_parallel} stepperp {step_perp} maxattempt {maxattempt} usher {usher_flag} {etarget} {ds0} {dtheta0} {uovlp} {dsolvp} {eps} {nattempt} charged {charge_flag}
```

The above fix **invokes** the newly added `obmd` extension in LAMMPS. The associated parameters are given as follows:
- FIX_ID
  Unique user-assigned name for the fix. Here, `FIX_ID = 2` . 
- `GROUP_ID` 
  ID of the group of atoms to which `fix obmd` is applied. The group used here must already exist. Here, `GROUP_ID = all`.
- `FIX_STYLE`
  Name of a `fix` style to be used. In the case of OBMD simulations, one must choose `obmd`. 
- `ntype`
  Type of inserted particles. Commonly, the `ntype` is chosen to correspond to the type of particles that are part of the `GROUP_ID` group, while for molecules, it should be set to $0$.
- `nfreq`
  Specifies the number of timesteps after which the `fix obmd` is executed. Typically it set to $1$, meaning that the `obmd` extension is invoked at every timestep.
- `seed` 
  A random number used to generate the initial position of newly inserted particles.
- `pxx`
  Pressure applied onto the buffer regions. It represents the **equilibrium pressure of a fluid** and can be derived from the equation of state, or obtained from an NVT simulation by monitoring the pressure to which the system relaxes. Can be provided as a real constant number or as a variable using  `v_pxx`.
- `pxy`
  Shear stress with forces along the $y$-axis. It can be given as a real constant number or as a variable using `v_pxy`.
- `pxz`
  Shear stress with forces along the $z$-direction. It can be given as a real constant number or as a variable using `v_pxz`.
- `dpxx value`
  Pressure amplitude added to the variable `pxx` and applied only onto the left buffer. It can be utilized to simulate mechanical pressure wave with a frequency `freq` and can be provided as a real constant number or as a variable using `v_dpxx`.
- `freq`
  Frequency of the mechanical pressure wave (if excited, otherwise should be set to $0.0$). It can be given as a real constant number or as a variable using `v_freq`.
- `alpha`
  Parameter to (further) reduce the desired number of particles in the buffers. 
- `tau` 
  Characteristic relaxation time of the buffers.
- `nbuf`
  Desired number of particles in the buffer. It can be computed using $nbuf = N * l_{b} / lx$, where $N$ is the total number of particles in the simulation box, $l_{b}$ is the buffer length, and $l_{x}$ is the length of the simulation box. 
- `region1 value`
  Region ID of the left buffer, where normal forces are applied. Here, region ID is `value = leftB`.
- `region2 value`
  Region ID of the right buffer, where normal forces are applied. Here, region ID is `value = rightB`.
- `region3 value`
  Region ID of the left buffer, where tangential shear forces are applied. Here, region ID is `value = leftshear`.
- `region4 value` `
  Region ID of the right buffer, where tangential shear forces are applied. Here, `value = rightshear`.
- `region5 value`
  Region ID of the left buffer, where new particles will be inserted. Here, `value = leftBin`.
- `region6 value`
  Region ID of the right buffer, where new particles will be inserted. Here, `value = rightBin`.
- `buffersize value`
  Length of the buffer region. Here, `value = buffer_size`.
- `gfac value`
  Smoothing length. Its value must be between $0.0$ and $1.0$. Here, `value = gfac`.
- `stepparallel value`
  Weighting function used to distribute the total external force acting in **normal direction**. For now only **"smooth" weighting function is implemented** and `value = 1` must be used. Here, `value = step_parallel = 1`.
- `stepperp value`
  Weighting function used to distribute the total external force acting in **tangential direction**. For now only **Heaviside step weighting function is implemented** and `value = 0` must be used. Here, `value = step_perp = 0`.
- `maxattempt value`
  The maximum number of attempts to perform the insertion by calling the `try_inserting()` member function. Here, `value = max_attempt`. If `maxattempt value` is not given, the default value is set to $1$.
- `usher options
  Insertion algorithm, where the `options` are:
	- `usher_flag` 
	  `usher_flag` must be set to $1$ if one aims to employ the USHER insertion algorithm.
	- `etarget`
	  Target energy of the newly inserted particle.
	- `ds0`
	  Displacement step used by the USHER algorithm, which should be small but not too small compared to the interaction range.
	- `dtheta0` 
	  Angular steps used when inserting a molecule.
	- `uovlp`
	  Very large energy characterizing the overlap position, which is commonly set to $10000$.
	- `dsolvp`
	  Parameter set to the value of the first RDF maximum.
	- `eps`
	  Parameter used in case of overlaps. 
	- `nattempt`
	  Maximum number of iterations performed by the USHER algorithm, which is typically set to $40$.

##### SIMULATION EXECUTION

```Python
# ----------------- Compute Section -----------------

# 

# ----------------- Output Section -----------------
run             {steps}
	"""
	file_to_write.write(content_4in)
	file_to_write.close()

# -------------------------------------------------------- #

# write input script
if __name__ = "__main__":
	write_in(in_file)
	print(f"file {in_file} is written")
```

The `run` option is used to execute the simulation, where `steps` represent the number of simulation steps to be performed. Finally, the `in.simulation` is written using `file_to_write.write(content_4in)`.

#### COMBINED INPUT SCRIPT TO RUN OBMD SIMULATIONS USING LAMMPS

Below are example values to generate a file for the OBMD simulation.

```Python
# ----------------- Variables Section -----------------
units           lj
boundary        f p p
atom_style      atomic
comm_modify     vel yes
newton          on

region          leftB block 0.0 5.0391 0.0 11.198 0.0 11.198
region          rightB block 28.5549 33.594 0.0 11.198 0.0 11.198
region          leftshear block 0.0 0.0 0.0 0.0 0.0 0.0 
region          rightshear block 0.0 0.0 0.0 0.0 0.0 0.0 
region          leftBin block 0.0 5.0391 0.0 11.198 0.0 11.198
region          rightBin block 28.5549 33.594 0.0 11.198 0.0 11.198
region          roi block 5.0391 28.5549 0.0 11.198 0.0 11.198

# ----------------- Interaction Section -----------------
pair_style      dpd 1.0 1.0 2616

# ----------------- Atom Definition Section -----------------
read_data       dpd_8map_obmd.data

# ----------------- Settings Section -----------------
pair_coeff      * * 209.6 4.5 1.0

neighbor        0.4 bin
neigh_modify    delay 0 every 1

timestep        0.001464

fix             1 all nve
fix             2 all obmd 1 1 6111 188.0 0.0 0.0 0.0 0.0 0.7 0.005 1327 &
                region1 leftB region2 rightB region3 leftshear & 
                region4 rightshear region5 leftBin region6 rightBin &
                buffersize 5.0391 gfac 0.25 stepparallel 0 stepperp 1 & 
                maxattempt 1 usher 1 31.03 1.0 0.02 10000 1.5 1.0 40 charged 0

# ----------------- Compute Section -----------------

#

# ----------------- Output Section -----------------
run             2000000
```


If you are using `fix obmd`, please **cite the following articles**:
- R. Delgado-Buscalioni et. al., *Eur. Phys. J. Special Topics* **224**, 2331-2349 (2015).
- J. Sablić et. al., *Soft Matter* **12**, 2416-2439 (2016).
-  L. Delle Site and M. Praprotnik, *Phys. Rep.* **693**, 1-56 (2017).
- P. Papež and M. Praprotnik, *J. Chem. Theory Comput.* **18**, 1227-1240 (2022).
- R. Delgado-Buscalioni and P. V. Coveney, *J. Chem. Phys.* **119**, 978-987 (2003).