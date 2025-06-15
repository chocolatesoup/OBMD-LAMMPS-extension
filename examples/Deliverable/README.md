# User guide
To perform OBMD simulations using the `obmd` extension, we provide a **python script** `input.py` that generates input script for LAMMPS named `in.simulation`.
A code snippet of the **function** in the `input.py`  script that writes an input script for OBMD simulation of a **DPD fluid** is given and explained below.

### Initialization settings
The function `write_in(in_file)` is defined and takes as an argument the name of the input file to be created.
Here, `in.simulation`.
Using f-string, the block of code that forms the input script is stored in `content_4in`.
The latter is ultimately written to a file using `file_to_write.write(content_4in)`.
```python
in_file = "in.simulation"

# ------------------------------------------------------- #
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
region          rightBin block {xhi-buffer_size+offset} {xhi-offset} {offset} {yhi-offset} {offset} {yhi-offset}
region          roi block {buffer_size} {xhi-buffer_size} {ylo} {yhi} {zlo} {zhi}
```
> [!IMPORTANT]
> OBMD simulations should be performed using **reduced units** (`units lj`) to avoid errors in the calculation of the external force and the style of boundaries (i.e., `boundary`) for the global simulation box in each dimension should be set to `f p p`. The latter applies the same style to the lower and upper sides of the simulation box, i.e., in the case of OBMD simulations, the non-periodic style `f` is applied in the $x$-direction, while the box is periodic in $y$- and $z$-directions by choosing the periodic style `p`.

> [!NOTE]
> `atom_style`: The style of atoms in the system. 
>>The OBMD extensions was tested using `atomic`, `full`, and `molecular` atom styles.

> [!IMPORTANT]
> When performing OBMD simulations, a momentum-conserving thermostat should be applied.
> Usually, the **DPD thermostat** is used, where velocities should be communicated between neighboring processors and stored as properties of ghost atoms using `comm_modify vel yes`.
> Importantly, the Newton’s third law for pairwise and bonded interactions should be enabled via `newton on`.

One should also define geometric regions of space. 
Here, **6 regions**, namely `leftB`, `rightB`, `leafshear`, `rightshear`, `leftBin`, and `rightBin`, should be given. 
All regions should be defined using the following format: $xlo$ $xhi$ $ylo$ $yhi$ $zlo$ $zhi$ (units of box).

> [!NOTE]
> `leftB`: Region defining the left buffer.
>>`buffer_size`: length of the buffer region, `xlo`, `ylo`, `zlo` are the lower bounds of the simulation box in the $x$-, $y$-, and $z$-direction, respectively, while `xhi`, `yhi`, `zhi` are the upper bounds.
> 
> `rightB`: Region defining the right buffer.
> 
> `leftshear`: Region defining part of the left buffer where shear forces are applied.
>>`shear_size`: length of the region where shear flow is imposed. If only normal pressure forces will be applied, e.g., equilibrium simulations, `shear_size`, `xlo`, `ylo`, `zlo`, `xhi`, `yhi`, and `zhi` should be set to 0.0.
>
> `rightshear`: Region defining part of the right buffer where shear forces are applied.
>
> `leftBin`: Region defining part of the left buffer where new particles are inserted.
>> `offset`: offset the value of the $xlo$, $xhi$, $ylo$, $yhi$, $zlo$, $zhi$ insertion region.
>
> `rightBin`: Region defining part of the right buffer where new particles are inserted.

### Simulation settings and system definition
```python
# ----------------- Interaction Section -----------------
pair_style      dpd {temp} {rc} {rnd_dpd}

# ----------------- Atom Definition Section -----------------
read_data       dpd_8map_obmd.data

# ----------------- Settings Section -----------------
pair_coeff      * * {aij} {gamma_dpd} {rc}

neighbor        {skin} bin
neigh_modify    delay 0 every 1

timestep        {dt}

fix             1 all nve
```
In this case, `pair_style dpd` is used to perform simulations employing the DPD force filed.
According to the [LAMMPS documentation](https://docs.lammps.org/pair_dpd.html) it takes temperature (`temp`), cutoff of the interaction (`rc`), and seed (`rnd_dpd`) as input parameters.
Alongside, the following [coefficients](https://docs.lammps.org/pair_dpd.html) (`pair_coeff`) must be given: interaction parameter between atoms $i$ and $j$ (`aij`), friction parameter (`gamma_dpd`), and cutoff of the interaction (`rc`).
> [!NOTE]
> `read_data` is used to read pre-equilibrated configuration of the DPD fluid given in the `dpd_8map_obmd.data` file.

> [!CAUTION]
> `dpd_8map_obmd.data` contains information about the positions of the atoms in **reduced units**.
> It is important that the position of the atoms do not exceed $xlo$ or $xhi$ of the simulation box; if so, they should be manually "deleted" before conducting the OBMD simulations.

```python
fix             2 all obmd {insert_type} {nevery} {rnd_obmd} {pxx} {pxy} {pxz} {dpxx} {freq} {buffer_alpha} &
                {buffer_tau} {n_buffer} region1 leftB  region2 rightB region3 leftshear region4 rightshear &
                region5 leftBin region6 rightBin buffersize {buffer_size} gfac {gfac} stepparallel {step_parallel} &
                stepperp {step_perp} maxattempt {max_attempt} usher {usher_flag} {etgt} {ds0} {dtheta} {uovlp} &
                {dsolvp} {eps} {max_attempt_usher} charged {charge_flag}
```

> [!NOTE]
> The above `fix` **invokes** the newly added `obmd` extension in LAMMPS.
> The associated parameters are given as follows:
>
>`FIX_ID`: Unique user-assigned name for the fix.
>> Here, `FIX_ID` is set to $2$.
>
>`FIX_GROUP`: ID of the group of atoms to which `fix` is applied (the group used here must already exist). 
>>Commonly `all` is selected.
>
>`FIX_STYLE`: Name of a `fix`style to be used. In case of OBMD, one must define `obmd`.
>
>`insert_type`: Type of the inserted particle.
>>Commonly, the `insert_type` equals to the type of particles that are part of the group `FIX_GROUP`, while for molecules `insert_type` is set to $0$.
>
>`nevery`: Specifies the number of timesteps after which the `fix` is executed.
>>Tipically `nevery` is set to $1$, meaning that `obmd` extension is invoked at every timestep.
>
>`rnd_obmd`: A random number used to generate the initial position of newly insetred particles.
>
>`pxx`: Pressure applied onto the buffer.
>> It represents the equilibirum pressure of the fluid and can be derived from the equation of state (for the DPD fluid), or obtained from NVT simulation by monitoring the pressure to which the system relaxes. Can be provided as a real constant or as a variable (`v_pxx`).
>
>`pxy`: Shear stress with forces along the $y$-axis.
>
>`pxz`: Shear stress with forces along the $z$-axis.
>>Parameters `pxy` and `pxz` are used to compute shear forces that are applied only in the shear region. Both can be provided as a real constant or as a variable (e.g., `v_pxy` and `v_pxz`).
>
>`dpxx`: Pressure amplitude added to the `pxx` and applied only onto the left buffer. 
>>Can be utilized to simulate acoustic wave with a frequency `freq`. Both can be provided as a real constant or as a variable (e.g., `v_dpxx` and `v_freq`).
>
>`buffer_alpha`: Parameter to (further) reduce the desired number of particles in the buffers.
>> If `buffer_alpha` is not provided, its default value is $0.75$.
>
>`buffer_tau`: Characteristic relaxation time of the buffers.
>>Usually is of the order of $100$ MD timesteps. If `buffer_tau` is not provided, its default value is $0.005$.
>
>`n_buffer`: The desired number of particles in the buffer.
>>It can be computed using $`n\_buffer`$ = $N$ * $`buffer\_size`$ / $lx$, where $N$ is the total number of particles in the simulation box and $lx$ is its length. If `n_buffer` is not provided, its default value is computed using $`n\_buffer= group->count(igroup)*buffer\_size/(xhi-xlo)`$.
>
> `region1 leftB`: Region ID of the left buffer, where normal pressure forces are applied.
>> Here, region ID is `leftB`.
>
> `region2 rightB`: Region ID of the right buffer, where normal pressure forces are applied.
>> Here, region ID is `rightB`.
>
> `region3 leftshear`: Region ID of the left buffer, where tangential shear forces are applied.
>> Here, region ID is `leftshear`.
>
> `region4 rightshear`: Region ID of the right buffer, where tangential shear forces are applied.
>> Here, region ID is `rightshear`.
>
> `region5 leftBin`: Region ID of the left buffer, where new particles will be inserted.
>> Here, region ID is `leftBin`.
>
> `region6 rightBin`: Region ID of the right buffer, where new particles will be inserted.
>> Here, region ID is `rightBin`.
>
> `buffersize buffer_size`: Length of the buffer region.
>> Here, length is set to the value of `buffer_size`.
>> If `buffersize`is not provided, its default is computed using $0.20(xhi-xlo)$.
>
> `gfac gfac`: Smoothing length having value between $0.0$ and $1.0$.
>> Here, smoothing length is set to the value of `gfac`.
>> If `gfac` is not provided, its default value is $0.25$.
>
> `stepparallel step_parallel`: Weighting function used to distribute the total external force acting in normal direction. 
>> Here, `stepparallel` is set to the value of `step_parallel`. For now only **smooth weighting function is applied**, therefore, `step_parallel` must be set to `0` (also default value).
>
> `stepperp step_perp`: Weighting function used to distribute the total external force acting in tangetial direction.
>> Here, `stepperp` is set to the value of `step_perp`. For now only **heaviside step weighting function is applied**, therefore, `step_perp` must be set to `1` (also default value).
>
> `maxattempt max_attempt`: The maximum number of attempts to perform the insertion by calling `try_inserting()` function.
>> Here, the maximum number of attempts is set to the value of `max_attempt`.
>> If `maxattempt` is not provided, its default value is $40$.
>
> `usher`: Insertion algorithm.
> 
>> `usher_flag`: Should be set to $1$ if one aims to employ USHER insertion algorithm.
>> 
>> `etgt`: Target energy of the inserted particle.
>> The default `etgt` value is set to $100.0$.
>> 
>> `ds0`: Displacement step, which should be small but not too small compared to the interaction range.
>> The default `ds0` value is set to $0.1$.
>> 
>> `dtheta`: Angular steps, which should be small but not too small. This paramater is relevant when inserting molecules.
>> The default `ds0` value is set to $0.1$.
>> 
>> `uovlp`: Very large energy representing an overlap position.
>> The default `uovlp` value is set to $10000.0$.
>> 
>> `dsolvp`: Parameter that should be set to about the first RDF maximum.
>> The default `dsolvp` value is set to $1.0$.
>> 
>> `eps`: In case of overlaps, the overlap between atoms must be computed.
>> The default `eps` value is set to $1.0$.
>> 
>> `max_attempt_usher`: Maximum number of iterations performed by USHER.
>> The default `max_attempt_usher` value is set to $40$.
>> 
>> `charged charge_flag`: To distinguish simulations of non-charged (`charge_flag` is $0$) and charged (`charge_flag`is $1$) atoms.

### Simulation execution
```python
# ----------------- Compute Section -----------------

#

# ----------------- Output Section -----------------
run             {steps}
    """
    file_to_write.write(content_4in)
    file_to_write.close()

# ------------------------------------------------------- #

# write input script
if __name__ == "__main__":
    write_in(in_file)
    print(f"file {in_file} is written")
```
The `run` option is used to run the simulation, where `steps` represents the number of simulation steps to be executed.
Finally, the `in.simulation` is written using `file_to_write.write(content_4in)`.

## Input script to perform OBMD simulation of DPD fluid in LAMMPS, created using `input.py`
```bash
# ----------------- Variables Section -----------------
units           lj
dimension       3
boundary        f p p
atom_style      atomic
comm_modify     vel yes

region          leftB block 0.0 5.039193729003359 0.0 11.198208286674133 0.0 11.198208286674133
region          rightB block 28.555431131019034 33.59462486002239 0.0 11.198208286674133 0.0 11.198208286674133
region          leftshear block 0.0 0.0 0.0 0.0 0.0 0.0
region          rightshear block 0.0 0.0 0.0 0.0 0.0 0.0
region          roi block 5.039193729003359 28.555431131019034 0.0 11.198208286674133 0.0 11.198208286674133

# ----------------- Init Section -----------------
pair_style      dpd 1.0 1.0 51949

# ----------------- Atom Definition Section -----------------
read_data       dpd_8map_obmd.data

# ----------------- Settings Section -----------------
pair_coeff      * * 209.6 4.5 1.0

neighbor        0.4 bin 
neigh_modify    every 1 delay 0 

timestep        0.0014641288433382138
fix             1 all nve

fix             2 all obmd 1 1 54152 188.0 0.0 0.0 0.0 0.0 0.7 0.005 1326 region1 leftB region2 rightB &
                region3 leftshear region4 rightshear region5 leftB region6 rightB &
                buffersize 5.039193729003359 shearsize 0.0 gfac 0.25 stepparallel 0 maxattempt 1 &
                usher 1 31.0272 1.0 0.017453292519943295 10000 1.5 1 40

# ----------------- Compute Section -----------------

#

# ----------------- Output Section -----------------
run             2000000
```

# Programmer's manual


# Architectural document
