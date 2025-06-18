
In the Programmer's manual, instructions on how to compile LAMMPS (**version 7 Feb 2024**) with the `obmd` extension are provided, along with explanation of the structure, classes, and function of the added `lj/cut/rf` pair style. The code introduced below can be downloaded [here](https://github.com/chocolatesoup/OBMD-LAMMPS-extension/tree/main). 

## OBMD PACKAGE

The OBMD package is added to the LAMMPS code as a subdirectory with the package name in capital letters. **OBMD package represents the extension that allows grand-canonical simulations in and out-of-equilibrium using particles**. To this end, it implements routines that take care of the deletion and insertion of particles, as well as computation and distribution of external forces, used to impose user-defined external boundary conditions (e.g., constant normal load, shear flow, mechanical pressure wave).

###### Build LAMMPS using CMake with OBMD extension enabled

To perform OBMD simulations, LAMMPS code needs to be compiled with the OBMD package enabled. Instructions on how to compile LAMMPS can be found [here](https://docs.lammps.org/Build.html), while the OBMD package can be included using `-D PKG_OBMD=on` when compiling.

```bash
#!/bin/sh

mkdir -p build_mpi

cd build_mpi

cmake ../cmake \
-D PKG_OBMD=on 

make -j 12
```

**Additions of the OBMD package**
- `src/OBMD`
- `fix obmd`
### FIX STYLE obmd

In this section, we explain the `obmd` extension, which is implemented as a new **fix style** within the `FixObmdMerged` class and is used for conducting OBMD simulations both in and out-of-equilibrium. **This fix manages the deletion and insertion of particles while concurrently monitoring the outgoing and incoming linear momentum of both deleted and newly inserted particles**, as implemented in the `fix_obmd_merged.cpp` and its corresponding `fix_obmd_merged.h` header file. Accurately measuring the outgoing and ingoing linear momentum is essential for **imposing the external boundary conditions through external forces**, which is also performed by the `obmd` fix style.

The LAMMPS input syntax is as follows:

```bash
fix FIX_ID GROUP_ID obmd ntype nfreq seed pxx pxy pxz dpxx freq alpha tau nbuf keyword value
```

The `FIX_ID` serves as the unique identifier for the fix, while `GROUP_ID` specifies the ID of the group to which the fix will be applied. The name of the `obmd` fix style is followed by **11 essential parameters** that are required for its execution. 
Namely:
- `ntype` 
  Atom type to assign to the newly inserted particles (offset when molecule is to be inserted).
- `nfreq`
  Specifies the number of timesteps after which the fix is executed.
- `seed` 
  A random number (positive integer).
- `pxx` 
  Pressure (real constant number or equal style variable).
- `pxy` 
  Shear stress with forces acting along $y$-axis (real constant number or equal style variable).
- `pxz`
  Shear stress with forces acting along $z$-axis (real constant number or equal style variable).
- `dpxx`
  Pressure amplitude (real constant number or equal style variable).
- `freq` 
  Frequency (real constant number or equal style variable).
- `alpha` 
  Parameter to (further) reduce the desired number of particles in the buffers (real constant number or equal style variable).
- `tau`
  Relaxation time of the buffers (real constant number or equal style variable).
- `nbuf`
  Desired number of particles in the buffer (real constant number or equal style variable).

Additionally, there are several (optional) parameters (`keyword value`) that can be included to optimize the performance, execution, and various factors within the simulation. 
Keyword:
- `region1 value`
  Region ID of the left buffer, where **normal forces** are applied.
- `region2 value`
  Region ID of the right buffer, where **normal forces** are applied.
- `region3 value`
  Region ID of the left buffer, where **tangential shear forces** are applied.
- `region4 value` 
  Region ID of the right buffer, where **tangential shear forces** are applied.
- `region5 value`
  Region ID of the left buffer, where new particles will be inserted.
-  `region6 value`
  Region ID of the right buffer, where new particles will be inserted. 
- `buffersize value`
  Length of the buffer region.
- `gfac value`
  Smoothing length (value between $0.0$ and $1.0$).
- `stepparallel value`
  Weighting function used to distribute the total external force acting in **normal direction**. `value = 1` represents the "smooth" weighting function.
- `stepperp value`
  Weighting function used to distribute the total external force acting in **tangential direction**. `value = 0` represents the Heaviside step weighting function.
- `maxattempt value`
  The maximum number of attempts to perform the insertion by calling the `try_inserting()` member function.
- `usher value`
  To select `usher` as the insertion algorithm, the `value` should be set to $1$, and followed by the parameters below.
	- `etarget` 
	  Target energy of the newly inserted particle. 
	- `ds0`
	  Displacement step used by the USHER algorithm.
	- `dtheta0`
	  Angular step used when inserting a molecule.
	- `uovlp`
	  Very large energy characterizing the overlap position.
	- `dsolvp`
	  Parameter set to the value of the first RDF maximum.
	- `eps`
	  Parameter used in case of overlaps.
	- `nattempt`
	  Maximum number of iterations performed by the USHER algorithm.
- `charged value`
  To simulate **charged** particles, the `value` should be set to $1$, otherwise to $0$.
- `mol value` 
  To insert molecule, the `value` should contain information about molecular template and the number of particles that need to be inserted. For more details, readers are referred to the [LAMMPS documentation](https://www.lammps.org/#gsc.tab=0).

Several other options are adopted from the LAMMPS code:
- `molfrac` 
- `rigid` 
- `shake`
- `global` 
- `local` 
- `near`
  Represents the insertion algorithm. It is not possible to select both `near` and `usher` at the same time.
- `rate`
- `vx`
- `vy`
- `vz`
- `orient`
- `units` 
- `gaussian` 
- `target`

and we recommend that readers are refer to the [LAMMPS documentation](https://www.lammps.org/#gsc.tab=0) for instructions on how to use them.

#### HEADER FILE

The first segment of the header file includes a copyright and license statement and attribution statement with a description of the code base and a brief description of the modifications and adjustments made.

```C++
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.


   Attribution:
   - Original code by Sandia Corporation.
   - OBMD extension uses portions of code derived from evaporate and deposit fix style:
   https://docs.lammps.org/fix_evaporate.html
   https://docs.lammps.org/fix_deposit.html

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
```

New fix style is registered with LAMMPS. The name of the fix is `obmd` and `FixObmdMerged` is the name of the class that implements it.

```C++
#ifdef FIX_CLASS
// clang-format off
FixStyle(obmd, FixObmdMerged)
// clang-format on
#else

// Content

#endif
```

The `// Content` part of the header file contains the complete definition of the `FixObmdMerged` class, including essential `setmask()` method. This method specifies functions that are invoked during the timestep.

```C++
#ifndef LMP_FIX_OBMD_MERGED_H
#define LMP_FIX_OBMD_MERGED_H

#include "fix.h"

namespace LAMMPS_NS {
class FixObmdMerged : public Fix {
 public:
  FixObmdMerged(class LAMMPS *, int, char **);
  ~FixObmdMerged();

  int setmask();
  void init();
  void pre_exchange();
  void post_force(int);
  void setup(int);
  void *extract(const char *, int &);
```

Several variables, arrays, and parameters are declared as `private` members of the derived class `FixObmdMerged`. Additionally, the following functions are implemented to perform the insertion and deletion of particles, as well as to apply external forces required for imposing boundary conditions:
- `try_deleting()`
- `try_inserting()`
- `usher()`
- `energy()`
- `energy_atomistic_obmd()`
- `reg_force()`
- `reg_force_perp()`
- `g_par_global_charged()`
- `g_par_local_charged()`
- `g_perp_global_charged()`

Each of the functions is explained in detail in the [[#IMPLEMENTATION FILE]] section.

#### IMPLEMENTATION FILE

The implementation of the `FixObmdMerged` class can be found in the `fix_obmd_merged.cpp`, which also contains copyright and license statement and attribution statement with a description of the code base and a brief description of the modifications and adjustments made. The `#include "fix_obmd_merged.h"` statement for the class header is listed first, followed by other essential include statements that provide access to other LAMMPS variables and functions, along with the `#include "pair_lj_cut_rf.h"`. The latter is needed to compute the potential energy of the charged particles to be inserted. 

```C++
#include "fix_obmd_merged.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "pair_lj_cut_rf.h"
#include "random_park.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include <cmath>
#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { ATOM, MOLECULE };
enum { DIST_UNIFORM, DIST_GAUSSIAN };
enum { NONE, CONSTANT, EQUAL };
enum { EXCHATOM, EXCHMOL };   
enum { MOVEATOM, MOVEMOL };   

#define EPSILON 1.0e-6
```

Constructor defined the `obmd` style that inherits from the `Fix` class. The arguments `FIX_ID`, `GROUP_ID`, and `FIX_STYLE = obmd` are read by the parent (i.e., `Fix`) class. To store several variables, the allocation code is added to the constructor.

```C++
FixObmdMerged::FixObmdMerged(LAMMPS *lmp, int narg, char **arg) :
	Fix(lmp, narg, arg), idregion(nullptr), idregion2(nullptr),idregion3(nullptr), idregion4(nullptr), idregion5(nullptr), idregion6(nullptr), idrigid(nullptr), idshake(nullptr), onemols(nullptr), molfrac(nullptr), coords(nullptr), imageflags(nullptr), fixrigid(nullptr), fixshake(nullptr), random(nullptr), list(nullptr), mark(nullptr)
```

In total there are **14 mandatory parameters** (including `FIX_ID`, `GROUP_ID`, and `FIX_STYLE = obmd`), which are specified in and read from the input script.

```C++
if (narg < 14) error->all(FLERR, "Illegal fix obmd command (check # of arguments)");
```

These parameters include:
- `ntype` (in the above e)
- `nfreq`
- `seed`
- `pxx`
- `pxy`
- `pxz`
- `dpxx`
- `freq`
- `alpha`
- `tau`
- `nbuf`

and several of them can be added as variables when preceded using `"v_"`.

```C++
ntype = utils::inumeric(FLERR, arg[3], false, lmp);
nfreq = utils::inumeric(FLERR, arg[4], false, lmp);
seed = utils::inumeric(FLERR, arg[5], false, lmp);
if (seed <= 0) error->all(FLERR, "Illegal fix obmd command");

// OBMD settings
if (strstr(arg[6], "v_") == arg[6]) {
	int n = strlen(&arg[6][2]) + 1;
    pxxstr = new char[n];
    strcpy(xstr, &arg[6][2]);
	} else {
    pxx = utils::numeric(FLERR, arg[6], false, lmp);
    pxxstyle = CONSTANT;
}

if (strstr(arg[7], "v_") == arg[7]) {
    int n = strlen(&arg[7][2]) + 1;
    pxystr = new char[n];
    strcpy(ystr, &arg[7][2]);
    } else {
    pxy = utils::numeric(FLERR, arg[7], false, lmp);
    pxystyle = CONSTANT;
}

if (strstr(arg[8], "v_") == arg[8]) {
    int n = strlen(&arg[8][2]) + 1;
    pxzstr = new char[n];
    strcpy(zstr, &arg[8][2]);
    } else {
    pxz = utils::numeric(FLERR, arg[8], false, lmp);
    pxzstyle = CONSTANT;
}

if (strstr(arg[9], "v_") == arg[9]) {
    int n = strlen(&arg[9][2]) + 1;
    dpstr = new char[n];
    strcpy(dpstr, &arg[9][2]);
    } else {
    dpxx = utils::numeric(FLERR, arg[9], false, lmp);
    dpstyle = CONSTANT;
}

if (strstr(arg[10], "v_") == arg[10]) {
    int n = strlen(&arg[10][2]) + 1;
    freqstr = new char[n];
    strcpy(freqstr, &arg[10][2]);
    } else {
    freq = utils::numeric(FLERR, arg[10], false, lmp);
    freqstyle = CONSTANT;
}

if (strstr(arg[11], "v_") == arg[11]) {
    int n = strlen(&arg[11][2]) + 1;
    alphastr = new char[n];
    strcpy(alphastr, &arg[11][2]);
    } else {
    alpha = utils::numeric(FLERR, arg[11], false, lmp);
    alphastyle = CONSTANT;
}

if (strstr(arg[12], "v_") == arg[12]) {
    int n = strlen(&arg[12][2]) + 1;
    taustr = new char[n];
    strcpy(taustr, &arg[12][2]);
    } else {
    tau = utils::numeric(FLERR, arg[12], false, lmp);
    taustyle = CONSTANT;
}

if (strstr(arg[13], "v_") == arg[13]) {
    int n = strlen(&arg[13][2]) + 1;
    nbufstr = new char[n];
    strcpy(nbufstr, &arg[13][2]);
    } else {
    nbuf = utils::numeric(FLERR, arg[13], false, lmp);
    nbufstyle = CONSTANT;
}
```

If `"v_"` is defined, the arguments are read as variables and their variable name is stored using the corresponding `str` string (i.e., `pxxstr`, `pxystr`, `pxzstr`, `dpstr`, `freqstr`, `alphastr`, `taustr`, and `nbufstr`). These are utilized in the `init()` method. The `init()` method defines the `pxxstyle`, `pxystyle`, `pxzstyle`, `dpstyle`, `freqstyle`, `alphastyle`, `taustyle`, and `nbufstyle` variables.

```C++
/* ---------------------------------------------------------------------- */
void FixObmdMerged::init()
{
	if (pxxstr) {
    pxxvar = input->variable->find(pxxstr);
    if (pxxvar < 0) error->all(FLERR, "Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(pxxvar))
      pxxstyle = EQUAL;
    else if (input->variable->atomstyle(pxxvar))
      pxxstyle = ATOM;
    else
      error->all(FLERR, "Variable for fix obmd is invalid style");
	}

	if (pxystr) {
    pxyvar = input->variable->find(pxystr);
    if (pxyvar < 0) error->all(FLERR, "Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(pxyvar))
      pxystyle = EQUAL;
    else if (input->variable->atomstyle(pxyvar))
      pxystyle = ATOM;
    else
      error->all(FLERR, "Variable for fix obmd is invalid style");
    }

    if (pxzstr) {
    pxzvar = input->variable->find(pxzstr);
    if (pxzvar < 0) error->all(FLERR, "Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(pxzvar))
      pxzstyle = EQUAL;
    else if (input->variable->atomstyle(pxzvar))
      pxzstyle = ATOM;
    else
      error->all(FLERR, "Variable for fix obmd is invalid style");
    }

    if (dpstr) {
    dpvar = input->variable->find(dpstr);
    if (dpvar < 0) error->all(FLERR, "Variable {} for fix obmd does not exist", dpstr);
    if (input->variable->equalstyle(dpvar))
      dpstyle = EQUAL;
    else if (input->variable->atomstyle(dpvar))
      dpstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix obmd is invalid", dpstr);
    }

    if (freqstr) {
    freqvar = input->variable->find(freqstr);
    if (freqvar < 0) error->all(FLERR, "Variable {} for fix obmd does not exist", freqstr);
    if (input->variable->equalstyle(freqvar))
      freqstyle = EQUAL;
    else if (input->variable->atomstyle(freqvar))
      freqstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix obmd is invalid", freqstr);
    }

    if (alphastr) {
    alphavar = input->variable->find(alphastr);
    if (alphavar < 0) error->all(FLERR, "Variable {} for fix obmd does not exist", alphastr);
    if (input->variable->equalstyle(alphavar))
      alphastyle = EQUAL;
    else if (input->variable->atomstyle(alphavar))
      alphastyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix obmd is invalid", alphastr);
    }

    if (taustr) {
    tauvar = input->variable->find(taustr);
    if (tauvar < 0) error->all(FLERR, "Variable {} for fix obmd does not exist", taustr);
    if (input->variable->equalstyle(tauvar))
      taustyle = EQUAL;
    else if (input->variable->atomstyle(tauvar))
      taustyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix obmd is invalid", taustr);
    }

    if (nbufstr) {
    nbufvar = input->variable->find(nbufstr);
    if (nbufvar < 0) error->all(FLERR, "Variable {} for fix obmd does not exist", nbufstr);
    if (input->variable->equalstyle(nbufvar))
      nbufstyle = EQUAL;
    else if (input->variable->atomstyle(nbufvar))
      nbufstyle = ATOM;
    else
      error->all(FLERR, "Variable {} for fix obmd is invalid", nbufstr);
    }
```

They can be set to `EQUAL` or `ATOM`, depending whether the variable is of `equal` or `atom` style, and are used in the `pre_exchange()` method.

```C++
if (input->variable->equalstyle(pxxvar)) pxxstyle = EQUAL;
else if (input->variable->atomstyle(pxxvar)) pxxstyle = ATOM;
```

If `"v_"` is not given, the arguments are treated as numeric constant values.

```C++
pxx = utils::numeric(FLERR, arg[6], false, lmp);
pxxstyle = CONSTANT;
```

In the `init()` method, checks are performed to ensure that the regions required for conducting OBMD simulations are defined, and error is raised if they are not. 
 
```C++
// get regions
iregion = domain->get_region_by_id(idregion);    // LEFT BUFFER
if (!iregion) error->all(FLERR, "Region ID for fix obmd does not exist");

iregion2 = domain->get_region_by_id(idregion2);    // RIGHT BUFFER
if (!iregion2) error->all(FLERR, "Region ID for fix obmd does not exist");

iregion3 = domain->get_region_by_id(idregion3);    // SHEAR FLOW LEFT
if (!iregion3) error->all(FLERR, "Region ID for fix obmd does not exist");

iregion4 = domain->get_region_by_id(idregion4);    // SHEAR FLOW RIGHT
if (!iregion4) error->all(FLERR, "Region ID for fix obmd does not exist");

iregion5 = domain->get_region_by_id(idregion5);    // INSERTION LEFT
if (!iregion5) error->all(FLERR, "Region ID for fix obmd does not exist");

iregion6 = domain->get_region_by_id(idregion6);    // INSERTION RIGHT
if (!iregion6) error->all(FLERR, "Region ID for fix obmd does not exist");
```

During the OBMD simulation, particles will be deleted from and inserted into the system. Ensuring that the IDs of particles remain unique, the `maxtag_all` (representing current maximal atom ID) or/and `maxmol_all` (representing current maximal molecule ID), is/are also determined (and distributed across all processes).

```C++
// find current max atom and molecule IDs if necessary
if (idnext) find_maxid();
```

```C++
/* ---------------------------------------------------------------------- */
void FixObmdMerged::find_maxid()
{
	tagint *tag = atom->tag;
	tagint *molecule = atom->molecule;
	int nlocal = atom->nlocal;

	tagint max = 0;
	for (int i = 0; i < nlocal; i++) max = MAX(max, tag[i]);
	MPI_Allreduce(&max, &maxtag_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);

	if (mode == MOLECULE && molecule) {
    max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max, molecule[i]);
    MPI_Allreduce(&max, &maxmol_all, 1, MPI_LMP_TAGINT, MPI_MAX, world);
	}
}
```

The allocated memory is freed in the destructor that follows the constructor.

```C++
/* ---------------------------------------------------------------------- */
FixObmdMerged::~FixObmdMerged()
{
	delete random;
	delete[] molfrac;
	delete[] idrigid;
	delete[] idshake;
	delete[] idregion;
	delete[] idregion2;
	delete[] idregion3;
	delete[] idregion4;
	delete[] idregion5;
	delete[] idregion6;
	memory->destroy(coords);
	memory->destroy(imageflags);
	memory->destroy(list);
	memory->destroy(mark);
}
```

The execution of the `obmd` fix during the timestep is controlled by the `setmask()` method. The latter invokes methods of the `FixObmdMerged` class at the `pre_exchange()` and `post_force()` stages of the timestep.

```C++
/* ---------------------------------------------------------------------- */
int FixObmdMerged::setmask()
{
	int mask = 0;
	mask |= PRE_EXCHANGE;
	mask |= POST_FORCE;
	return mask;
}
```

The **deletion and insertion of particles should be performed before the neighbor list is rebuilt**, therefore, these routines are implemented in the `pre_exchange()` method.  On the other hand, **external forces, which impose the external boundary conditions, should be applied after evaluating the bonded and non-bonded forces**. Therefore, these forces are applied to the particles in the buffers in the `post_force()` method.

Invoking the `pre_exchange()` method, call to the `try_deleting()` member function is made first. This function takes `Region_id`, `vnewl`, and `vnewr` as arguments. Here, `Region-ID` is set to `iregion` when considering the left buffer, and to `iregion2` when considering the right buffer.

```C++
/* ----------------------------------------------------------------------
   perform particle/molecule insertions/deletions
------------------------------------------------------------------------- */
void FixObmdMerged::pre_exchange()
{
	// delete particles that cross the open boundaries after first half of velocity-Verlet algorithm
	try_deleting(iregion, vnewl, vnewr);
	try_deleting(iregion2, vnewl, vnewr);
```

In the `try_deleting()` member function, a loop over all `nlocal` particles in made. If the particle's position in the open direction (i.e., `x[][0]`, corresponding to the position in the $x$-direction) is less than the lower limit (`boxl`) or greater than the upper limit (`boxh`) of the simulation box in the $x$-direction, the particle's local index (`i`) is added to the list of particles that need to be deleted (`list[ncount++] = i`). Each process stores the number of particles to be deleted in variable `ncount`. The total number of particles to be deleted is stored in value `nall`, representing summed value of `ncount` from all processes. The sum of `ncount` from all processes with lower rank (i.e., the processes before the current one) is stored in `nbefore`.

```C++
int ncount = 0;
for (i = 0; i < nlocal; i++) {
	if (x[i][0] < boxl || x[i][0] > boxh) {
	    std::cout << "Deleting x[i][0] = " << x[i][0] << " i = " << i << " type[i] = " << type[i]
                << " global id = " << atom->tag[i] << std::endl;
        list[ncount++] = i;
	}
}

int nall, nbefore;
MPI_Allreduce(&ncount, &nall, 1, MPI_INT, MPI_SUM, world);
MPI_Scan(&ncount, &nbefore, 1, MPI_INT, MPI_SUM, world);
nbefore -= ncount;
```

After, the deletion is performed. If `mode == ATOM`, the index of a particle to be deleted is randomly chosen across processors and marked for deletion (`mark[] = 1`). Variables `ncount`, `ndel`, and `nall` are updated accordingly. Using the post-decrement operator, the values of `ncount` and `nall` are decreased by $1$, while using post-increment operator, the value of variable `del` is increased by $1$. The latter stores the total number of deletions performed. 

```C++
// atomic deletions
// choose atoms randomly across all procs and mark them for deletion
// shrink eligible list as my atoms get marked
// keep ndel,ncount,nall,nbefore current after each atom deletion

if (mode == ATOM) {
    while (nall) {
	    iwhichglobal = static_cast<int>(nall * random->uniform());
	    if (iwhichglobal < nbefore) nbefore--;
	    else if (iwhichglobal < nbefore + ncount) {
	        iwhichlocal = iwhichglobal - nbefore;
	        mark[list[iwhichlocal]] = 1;
	        list[iwhichlocal] = list[ncount - 1];
	        ncount--;
	    }
    ndel++;
    nall--;
    }
}
```

If a molecule is to be deleted, a similar routine is executed, followed by the deletion of bonds, angles, proper dihedral angles, and improper dihedral angles (if any).

The vital part of the implemented OBMD method is monitoring the outgoing linear momentum. The positions (`x[][0]`) of all `nlocal` particles marked for deletion are extracted and checked to determine if they correspond to the crossing of the outer boundary of the left or right buffer region. If the particle has exited the left buffer, the variable `vnewl` is incremented by its linear momentum, while `vnewr` is incremented if the particle has exited the right buffer.

```C++
for (i = nlocal - 1; i >= 0; i--) {
    if (mark[i]) {
	      if (x[i][0] < 0.5 * (boxh + boxl)) {
	        vnewl[0] += mass[type[i]] * vel[i][0];
	        vnewl[1] += mass[type[i]] * vel[i][1];
	        vnewl[2] += mass[type[i]] * vel[i][2];
	      } else {
	        vnewr[0] += mass[type[i]] * vel[i][0];
	        vnewr[1] += mass[type[i]] * vel[i][1];
	        vnewr[2] += mass[type[i]] * vel[i][2];
	      }
      avec->copy(atom->nlocal - 1, i, 1);
      atom->nlocal--;
    }
}
```

The deletion of particles in the `pre_exchange()` method is followed by the computation of the number of particles to be inserted. The number of particles in the left and right buffer is stored in variables `cnt_left` and `cnt_right`, respectively, which are computed utilizing the `group->count(GROUP_ID, Region-ID)` implemented in the LAMMPS code. Equation

$$
\Delta N_{b} = \frac{\delta t}{\tau_{b}} \left( \langle \alpha N_{b} \rangle - N_{b}\right)
$$

gives the number of particles to be inserted into the left and right buffer, which are stored in variables `ninsert_left` and `ninsert_right`, respectively.

```C++
// count number of particles in left/right buffer
cnt_left = group->count(igroup, iregion);
cnt_right = group->count(igroup, iregion2);

// calculate number of particles needed for insertion
ninsert_left = -static_cast<int>((static_cast<double>(cnt_left) / mol_len - alpha * nbuf) * update->dt / tau);
ninsert_right = -static_cast<int>((static_cast<double>(cnt_right) / mol_len - alpha * nbuf) * update->dt / tau);
```

By calling `try_inserting()` member function (for each buffer region separately), the insertion procedure starts. The `try_inserting()` method takes `Region-ID`, `ninsert_left/right`, `vnewl`, and `vnewr` as arguments.

```C++
// tries inserting ninsert_left/right particles into left/right buffer
try_inserting(iregion5, ninsert_left, vnewl, vnewr);
try_inserting(iregion6, ninsert_right, vnewl, vnewr);
```

The `try_inserting()` member function supports both `near` and `usher` insertion algorithms, however, here we focus on the latter. Readers interested in the former are referred to the [LAMMPS documentation](https://www.lammps.org/#gsc.tab=0).

After selecting a random position within the buffer region for the particle to be inserted, its (potential) energy is computed by calling the `usher()` member function and stored in the variable `entmp`. 

```C++
else if (usherflag) {
    me = -1;
    entmp = usher(iregion_var, coords, etarget, natom, imol, iter);
    if (entmp < etarget + EPSILON) {
        if (comm->me == 0)
            std::cout << "USHER accepts at E = " << entmp << " in attempt No. " << attempt << " with " << iter << " iterations" << std::endl;
    } else {
        if (comm->me == 0)
            std::cout << "USHER denies at E = " << entmp << " at attempt No. " << attempt << std::endl;
            flag = 1;
    }
}
```

In the `usher()` method, call to `energy_atomistic_obmd()` or `energy()` member function is made. The former is used when simulating **charged particles** and takes the following parameters:
- `iregion_var`
  Region ID where particle to be inserted is located.
- `onemols[]->q[]`
  Charge of the particle to be inserted.
- `type_temp`
  Type of the particle to be inserted.
- `coords[]`
  Randomly selected position of the particle to be inserted (vector of length three).
- `fusher`
  Force acting on the particle to be inserted.

When simulating **non-charged particles**, the `energy()` member function is called and it takes
- `index` 
  (Not relevant when simulating DPD particles.)
- `type_temp`
- `coords[]`
- `fusher` 

as arguments.

```C++
if (chargeflag) {
    entmp +=
    energy_atomistic_obmd(iregion_var, onemols[imol]->q[m], type_temp, coords[m], fusher);
} else {
    entmp += energy(1, type_temp, coords[m],fusher);    // i - does not matter, itype inserted, coords[m], fusher
}
```
In the `energy_atomistic_obmd()` member function, the downcast conversion of the pointer to the pair potential class `PairLJCutRF` is made, so that newly added `single_atomistic_obmd()` method implemented in the `PairLJCutRF` class can be accessed.

```C++
// potential acting between atomistic particles
auto pair = dynamic_cast<PairLJCutRF *>(force->pair_match("lj/cut/rf", 1));
```

Looping over all `nlocal` particles, the distance between the particle to be inserted and the $j$-th particle of `nlocal` particles is calculated. If the distance (`rsq`) is less than a cutoff (`cutsq`), the energy and the force acting on the particle to be inserted are computed.

```C++
for (int j = 0; j < nlocal; j++) {
    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];

    domain->minimum_image(delx, dely, delz);
    rsq = delx * delx + dely * dely + delz * delz;

    jtype = type[j];

    if (rsq < cutsq[itype][jtype]) {
	    total_energy += pair->single_atomistic_obmd(qi, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

	    fusher[0] += fpair * delx;
	    fusher[1] += fpair * dely;
	    fusher[2] += fpair * delz;
    }
}
```

The `energy_atomistic_obmd()` member function returns the total energy of the particle to be inserted, which is stored in variable `total_energy`. Similar is performed when a call to the `energy()` function is made. The computed values of the total energy and force for each process are stored in variables `entmp` and `fusher`, respectively. The summed values of energy and force over all processes are stored in variables `entmp_all` and `fusher_all`, respectively.

```C++
MPI_Allreduce(fusher, fusher_all, 3, MPI_DOUBLE, MPI_SUM, world);
MPI_Allreduce(&entmp, &entmp_all, 1, MPI_DOUBLE, MPI_SUM, world);
```

The computed total energy (`entmp_all`) is first compared to the predefined threshold, and if it is lower, the position is accepted and the `usher()` member function is exited.

```C++
if (entmp_all < etarget + EPSILON) break;
```

Otherwise, the position is iteratively corrected according to the following update rule:

$$
\mathbf{r}^{n+1} = \mathbf{r}^{n} + \frac{\mathbf{f}^{n}}{|\mathbf{f}^{n}|} \delta s^{n},
$$
where displacement $\delta s^{n}$ depends on the local potential energy (stored in variable `entmp_all`) and force magnitude (stored in variable `fabs`). If the computed energy is larger than `uovlp`, which is chosen to be a very large energy representing the overlap position, $\delta s^{n}$ is computed using

$$
\delta s^{n} = \Delta s_{ovlp} = r_{\sigma} - \left( \frac{4\varepsilon}{U^{n}} \right)^{\frac{1}{2}},
$$

where particle position will be translated for a distance $r_{\sigma}$ away from the center of the overlapped particle.

```C++
else if (entmp_all > uovlp) {
    fabs = sqrt(fusher_all[0] * fusher_all[0] + fusher_all[1] * fusher_all[1] +
                  fusher_all[2] * fusher_all[2]);
    if (fabs < EPSILON) continue;
	ds = dsovlp - pow(4 * eps / entmp_all, 1.0 / 12.0);

    // loop over atoms of molecule
    for (m = 0; m < natom; m++) {
	    coords[m][0] += fusher_all[0] / fabs * ds;
        coords[m][1] += fusher_all[1] / fabs * ds;
        coords[m][2] += fusher_all[2] / fabs * ds;
	}
    int check = check_mol_region(iregion_var, coords, natom);
    if (check == 1) break;    // exit - particle is now within the buffer region
}
```
When the computed energy is lower than `uovlp`, $\delta s^{n}$ is determined using

$$
\delta s^{n} = \min \left(\Delta s, \frac{U^{n} - U^{0}}{|\mathbf{f}^{n}|} \right).
$$
Here, $U^{0}$ represents the target energy, which is  stored in variable `etarget`, while $\Delta s$ stands for the displacement step, which is stored in variable `ds0`.

```C++
else {
    fabs = sqrt(fusher_all[0] * fusher_all[0] + fusher_all[1] * fusher_all[1] + fusher_all[2] * fusher_all[2]);
	if (fabs < EPSILON) continue;
    ds = std::min((entmp_all - etarget) / fabs, ds0);

    if (mode == MOLECULE) {
        torqabs = sqrt(torq_all[0] * torq_all[0] + torq_all[1] * torq_all[1] + torq_all[2] * torq_all[2]);
        dtheta = std::min((entmp_all - etarget) / torqabs, dtheta0);
        MathExtra::norm3(torq_all);
        MathExtra::axisangle_to_quat(torq_all, dtheta, quat);
        MathExtra::quat_to_mat(quat, rotmat);
    }

    for (m = 0; m < natom; m++) {
        coords[m][0] += fusher_all[0] / fabs * ds;
        coords[m][1] += fusher_all[1] / fabs * ds;
        coords[m][2] += fusher_all[2] / fabs * ds;

        if (mode == MOLECULE) {
          MathExtra::matvec(rotmat, coords[m], coords_tmp);
          MathExtra::copy3(coords_tmp, coords[m]);
        }
    }
    int check = check_mol_region(iregion_var, coords, natom);
    if (check == 1) break;    // exit
}
```

In case of molecule insertion, rotation of the molecule about its center of mass is also performed. The position obtained by the iterative procedure is verified to remain within the buffer region by calling the `check_mol_region()` member function. If the iterative correction of the position leads to particle leaving the buffer region, the `usher()` method is exited, and the insertion procedure is terminated. On the other hand, if the iterative procedure is successful, a new particle/molecule is inserted. In case of molecule insertion, bonds, angles, proper dihedrals, and improper dihedrals (if any) are also created. As the newly created particles are inserted with zero velocity, the incoming linear momentum is also zero. Therefore, the computed change in the linear momentum comes only from the particles that have exited the outer boundary of the left and right buffer region. Values for the left and right buffer are for each process stored in variables `vnewl` and `vnewr`, respectively, while the corresponding summed values from all processes are stored in variables `vnewl_all` and `vnewr_all`, respectively.

```C++
MPI_Allreduce(vnewl, vnewl_all, 3, MPI_DOUBLE, MPI_SUM, world);
MPI_Allreduce(vnewr, vnewr_all, 3, MPI_DOUBLE, MPI_SUM, world);
```

The insertion procedure is followed by the computation of the external forces:

$$
\mathbf{F}^{ext} = \displaystyle{\sum_{i \in B}} \mathbf{f}_{i}^{ext} =  A \mathsf{J}^{p} \cdot \mathbf{n} - \frac{\displaystyle\sum_{i'\in B} \Delta (m_{i'} \mathbf{v_{i'}})}{\delta t},
$$

which are utilized in the `post_force()` method. The first term on the right-hand side of the above equation read as:

$$
A \mathsf{J}^{p} \cdot \mathbf{n} = p_{xx} A \mathbf{n}.
$$
when a system under a constant normal load applied from both sides of the simulation domain is simulated. To introduce a sinusoidal disturbance, it should be rewritten into

$$
A \mathsf{J}^{p} \cdot \mathbf{n} = p_{xx} A \mathbf{n} + \Delta p \sin(2\pi \nu t) A \mathbf{n}
$$
 while it is expressed as:

$$
A \mathsf{J}^{p} \cdot \mathbf{n} = p_{xx} A \mathbf{n} + p_{xy} A \mathbf{t}
$$
when a shear flow is introduced. 

```C++
// area of the buffer - ROI interface
double area = ly * lz;

simulation_time += update->dt;
double factor = pxx + dpxx * sin(2.0 * MY_PI * freq * simulation_time);

// calculates momentum forces on the left buffer
MathExtra::zero3(momentumForce_left);
momentumForce_left[0] = vnewl_all[0] / (update->dt) + factor * area;
momentumForce_left[1] = vnewl_all[1] / (update->dt);
momentumForce_left[2] = vnewl_all[2] / (update->dt);
shearForce_left[0] = 0.0;
shearForce_left[1] = pxy * area;
shearForce_left[2] = pxz * area;

// calculates momentum forces on the right buffer
MathExtra::zero3(momentumForce_right);
momentumForce_right[0] = vnewr_all[0] / (update->dt) - pxx * area;    // constant normal load
momentumForce_right[1] = vnewr_all[1] / (update->dt);
momentumForce_right[2] = vnewr_all[2] / (update->dt);
shearForce_right[0] = 0.0;
shearForce_right[1] = -pxy * area;
shearForce_right[2] = -pxz * area;
```

Variables `momentumForce[]_left/right` and `shearForce_left/right[]` represent the array of length three and carry information about the external boundary condition of the constant normal load (or mechanical pressure wave) and shear stress to be applied in the $x$-, $y$-, and $z$-directions, respectively. 

The `obmd` fix is also executed at the `post_force()` stage of the timestep, as defined in the `setmask()` method. At this stage, external forces computed in the `pre_exchange()` method are applied to the particles in the buffers. Forces acting in normal direction are distributed first by calling the `reg_force()` member function, which takes `Region-ID`, `momentumForce_left/right`, and `step_parallel` as arguments. Here, `Region-ID` denotes `iregion` corresponding to the left buffer region, where the total external force stored in variable `momentumForce_left` will be applied. The same applies to the imposition of the total external force (`momentumForce_right`) to the right buffer. The `step_parallel` parameter defines the weighting function used to distribute the total external force among particles within the buffers. 

```C++
/* ---------------------------------------------------------------------- */
void FixObmdMerged::post_force(int vflag)
{
	if (update->ntimestep % nevery) { return; }

	// parallel forces
	reg_force(vflag, iregion, momentumForce_left, step_parallel);
	reg_force(vflag, iregion2, momentumForce_right, step_parallel);
```

In the `reg_force()` method the total mass of particles in the corresponding buffer is computed by calling the `g_par_global_charged()` member function. The computed value in stored in variable `gtmp`. 

```C++
/* ---------------------------------------------------------------------- */
void FixObmdMerged::reg_force(int vflag, Region *region, double *momentumForce, int step)
{
  if (region) region->prematch();
  v_init(vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  int nlocal = atom->nlocal;
  double mass_tmp;
  int *type = atom->type;
  double gtmp = g_par_global_charged(region, step);
  double gloctmp;
  double unwrap[3];
  double xval, yval, zval;
```

For now, only **smooth weighting function** is enabled (i.e., `step_parallel == 0`), where parameter `gfac` is used to define the region in which the position dependent sigmoidal function is multiplied with the mass of a particle.

```C++
/* ---------------------------------------------------------------------- */

double FixObmdMerged::g_par_global_charged(Region *region, int step)
{
	if (region) region->prematch();

	int nlocal = atom->nlocal;
	tagint *molecule = atom->molecule;
	int *rep_atom = atom->rep_atom;
	double **cms = atom->cms_mol;
	double *mass = atom->mass;
	int *type = atom->type;
	double **x = atom->x;
	int imol;
	int *tag = atom->tag;
	double **f = atom->f;

	double upper_x = domain->boxhi[0];    // upper
	double lower_x = domain->boxlo[0];    // lower

	// init for distribution function
	double g_par_all = 0.0;
	double g_par_all_tmp = 0.0;

	// for now only smooth transition -> ADD for step = 1
	double carg;
	if (!step) {
	    for (int i = 0; i < nlocal; i++) {
		    // mass of a particle
		    double mass_temp = mass[type[i]];

		    // (looping over nlocals) check if in prescribed region, continue if not
		    if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

		    // if here ... molecule in buffer
		    if (x[i][0] < lower_x + buffer_size) {    // LEFT BUFFER
			    if (x[i][0] < (lower_x + (1.0 - g_fac) * buffer_size)) {    // only mass
		        g_par_all += mass_temp;
				} else {    // sigmoidal mass
			    carg = 1.0 / g_fac * MY_PI * (x[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI;    
		        g_par_all += 0.5 * (1.0 + cos(carg)) * mass_temp;
		        }
		    }
		    if (x[i][0] > upper_x - buffer_size) {    // RIGHT BUFFER
		        if (x[i][0] > (upper_x - (1.0 - g_fac) * buffer_size)) {
				    g_par_all += mass_temp;
		        } else {
			        carg = 1.0 / g_fac * MY_PI * (x[i][0] - upper_x + buffer_size) / (buffer_size) -MY_PI;    
				    g_par_all += 0.5 * (1.0 + cos(carg)) * mass_temp;
		        }
		      }
	    // OTHER IS NOT OF MY INTEREST ... ROI
	    }
	} else {
    error->all(FLERR, "For now implemented only for !step");
	}
  // sum and distribute
  MPI_Allreduce(&g_par_all, &g_par_all_tmp, 1, MPI_DOUBLE, MPI_SUM, world);

  return g_par_all_tmp;
}
```

The total mass computed by each process is stored in variable `g_par_all` , while the value returned by the `g_par_global_charged()` method is stored in variable `g_par_all_tmp`, which represents the summed value of `g_par_all` from all processes. 

In the `reg_force()` member function, the mass (or weighted mass) of a particle in the buffer region (stored in variable `gloctmp`) is divided by the total mass of the buffer region (`gtmp`), which is computed by calling the `g_par_global_charged()` method. Multiplying this ratio by the `momentumForce_left/right` yields the force acting on the particle within the corresponding buffer region.

```C++
  for (int i = 0; i < nlocal; i++) {
    if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
    domain->unmap(x[i], image[i], unwrap);
    mass_tmp = mass[type[i]];
    if (mode == MOLECULE) mass_tmp = mtot;

    gloctmp = g_par_local_charged(mass[type[i]], region, x[i][0], step);

    xval = momentumForce[0] * gloctmp / gtmp;
    yval = momentumForce[1] * gloctmp / gtmp;
    zval = momentumForce[2] * gloctmp / gtmp;
```

The latter is applied to each particle within the buffer region.

```C++
f[i][0] += xval;
f[i][1] += yval;
f[i][2] += zval;
```

To impose shear stress, the `reg_force_perp()` member function is called separately for each buffer region in the `post_force()` method. It takes `vflag`, `Region-ID`, `shearForce_left/right` , and `step_perp` as arguments. Here, the `step_perp` parameter defines the weighting function that is used to distribute shear forces (stored in `shearForce_left/right` representing the array of length three) among particles in the buffer regions.

```C++
 // tangential forces
reg_force_perp(vflag, iregion3, shearForce_left, step_perp);
reg_force_perp(vflag, iregion4, shearForce_right, step_perp);
```

In the `reg_force_perp()` member function the total mass of particles in the corresponding buffer is computed by summing their masses by calling the `g_perp_global_charged()` member function. The computed value is stored in variable `gtmp`. 

```C++
/* ---------------------------------------------------------------------- */

void FixObmdMerged::reg_force_perp(int vflag, Region *region, double *shearForce, int step)
{
	if (region) region->prematch();
	v_init(vflag);

	double **x = atom->x;
	double **f = atom->f;
	double *mass = atom->mass;
	int *mask = atom->mask;
	imageint *image = atom->image;
	double v[6];
	int nlocal = atom->nlocal;
	int *type = atom->type;
	double gtmp = g_perp_global_charged(region, step);
	double mass_tmp, gloctmp;
	double unwrap[3];
	double xval, yval, zval;
```

For now, only **Heaviside step function** is enabled (i.e., `step_perp == 1`).

```C++
/* ---------------------------------------------------------------------- */

double FixObmdMerged::g_perp_global_charged(Region *iregion_var, int step)
{
	if (iregion_var) iregion_var->prematch();

	int nlocal = atom->nlocal;
	tagint *molecule = atom->molecule;
	double *mass = atom->mass;
	int *type = atom->type;
	double **x = atom->x;
	int imol;
	int *tag = atom->tag;
	double **f = atom->f;

	double upper_x = domain->boxhi[0];
	double lower_x = domain->boxlo[0];

	double g_perp_all = 0.0;
	double g_perp_all_tmp = 0.0;

	if (step) {
	    for (int i = 0; i < nlocal; i++) {

	      // mass of a particle
	      double mass_temp = mass[type[i]];

	      if (iregion_var && !iregion_var->match(x[i][0], x[i][1], x[i][2])) continue;
	      g_perp_all += mass_temp;
	    }
	} else {
    error->all(FLERR, "For now implemented only for step");
	}

	// sum and distribute
	MPI_Allreduce(&g_perp_all, &g_perp_all_tmp, 1, MPI_DOUBLE, MPI_SUM, world);
	return g_perp_all_tmp;
}
```

In the same manner as before, the total mass computed by each process is stored in variable `g_perp_all` , while the value returned by the `g_perp_global_charged()` method is stored in the variable `g_perp_all_tmp`, representing the summed value of `g_perp_all` from all processes. Shear forces are distributed by calling the `reg_force_perp()` member function. This routine follows the same approach as the one presented above for calling `reg_force()`, except this time, `shearForce_left/right` are utilized. 

```C++
/* ---------------------------------------------------------------------- */
void FixObmdMerged::reg_force_perp(int vflag, Region *region, double *shearForce, int step)
{
	if (region) region->prematch();
	v_init(vflag);

	double **x = atom->x;
	double **f = atom->f;
	double *mass = atom->mass;
	int *mask = atom->mask;
	imageint *image = atom->image;
	double v[6];
	int nlocal = atom->nlocal;
	int *type = atom->type;
	double gtmp = g_perp_global_charged(region, step);
	double mass_tmp, gloctmp;
	double unwrap[3];
	double xval, yval, zval;

	for (int i = 0; i < nlocal; i++) {
	    if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
	    domain->unmap(x[i], image[i], unwrap);
	
		gloctmp = mass[type[i]];

	    xval = shearForce[0] * gloctmp / gtmp;
	    yval = shearForce[1] * gloctmp / gtmp;
	    zval = shearForce[2] * gloctmp / gtmp;

	    f[i][0] += xval;
	    f[i][1] += yval;
	    f[i][2] += zval;

	    if (evflag) {
		    v[0] = xval * unwrap[0];
		    v[1] = yval * unwrap[1];
		    v[2] = zval * unwrap[2];
		    v[3] = xval * unwrap[1];
		    v[4] = xval * unwrap[2];
		    v[5] = yval * unwrap[2];
		    v_tally(i, v);
	    }
	}
}
```

If you are using `fix obmd`, please **cite the following articles**:
- R. Delgado-Buscalioni et. al., *Eur. Phys. J. Special Topics* **224**, 2331-2349 (2015).
- J. Sablić et. al., *Soft Matter* **12**, 2416-2439 (2016).
-  L. Delle Site and M. Praprotnik, *Phys. Rep.* **693**, 1-56 (2017).
- P. Papež and M. Praprotnik, *J. Chem. Theory Comput.* **18**, 1227-1240 (2022).
- R. Delgado-Buscalioni and P. V. Coveney, *J. Chem. Phys.* **119**, 978-987 (2003).

### PAIR STYLE lj/cut/rf

Because in the open-boundary molecular dynamics (OBMD),  electrostatic interactions are treated using reaction filed (rf) method, the latter was implemented as a `lj/cut/rf` pair style and added to the `src`.

Pair style `lj/cut/rf` is derived from the `Pair` class, and it is **used to compute the Lennard-Jones (LJ) potential and treat electrostatic interactions using the rf method**.  

$$
U_{rf} (r_{i\alpha j\beta}) = \frac{q_{i\alpha} q_{j\beta}}{4 \pi
\varepsilon_0 r_{i\alpha j\beta}} \left[ 1 + \frac{\varepsilon_{rf} - 1}
{2\varepsilon_{rf} +1} \left(\frac{r_{i\alpha j\beta}}{r_c} \right)^3
\right] - \frac{q_{i\alpha} q_{j\beta}}{4 \pi \varepsilon_0 r_c}
\frac{3\varepsilon_{rf}}{2\varepsilon_{rf} +1}
$$
Here, $q$, $\varepsilon_{0}$, and $\varepsilon_{rf}$ stand for the electric charge, vacuum permittivity, and dielectric constant, respectively.

The input script syntax reads as

```bash
pair_style lj/cut/rf rc_lj rc_rf
pair_coeff TYPE1 TYPE2 epsilon sigma rc_lj rc_rf epsilon_rf
```

The newly implemented pair style takes the following arguments:
- `rc_lj`
  Cutoff of the LJ interaction.
- `rc_rf`
  Cutoff of the electrostatic interaction.

For each pair of atom types (defined by `TYPE1` and `TYPE2`), the following coefficients must be determined:
- `epsilon`
  (Energy units.)
- `sigma`
  (Distance units.)
- `rc_lj` 
  (Distance units and optional parameter.)
- `rc_rf`
  (Distance units and optional parameter.)
- `epsilon_rf`
  Dielectric permittivity.

The class name of the `lj/cut/rf` pair style is `PairLJCutRF` and source files are `pair_lj_cut_rf.cpp` and `pair_lj_cut_rf.h`.

#### HEADER FILE

The first segment of the header file contains the copyright and license statement. 

```C++
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

```

Inclusion of the following lines in the second part of the header file registers pair style `lj/cut/rf` with LAMMPS. This block is included by the `Force` class in `force.cpp` that creates an instance of the `PairLJCutRF` class and returns a pointer to it, and connects the name of the `lj/cut/rf` pair style to the `PairLJCutRF` class.

```C++
#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/rf,PairLJCutRF);
// clang-format on
#else

// Content

#endif
```

The `// Content` part of the header file represents an actual definition of the `PairLJCutRF` class. The `PairLJCutRF` class includes required (i.e., `compute`, `settings`, and `coeff`) and some additional (i.e., `init_one`, `single`, `single_atomistic_obmd`, and `extract`) member functions. Several variables, arrays, and interaction parameters are declared as `protected` members of the derived `PairLJCutRF` class.

```C++
#ifndef LMP_PAIR_LJ_CUT_RF_H
#define LMP_PAIR_LJ_CUT_RF_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJCutRF : public Pair {
	public:
		PairLJCutRF(class LAMMPS *);
		~PairLJCutRF() override;
		void compute(int, int) override;
		void settings(int, char **) override;
		void coeff(int, char **) override;
		void init_style() override;
		double init_one(int, int) override;
		void write_restart(FILE *) override;
		void read_restart(FILE *) override;
		void write_restart_settings(FILE *) override;
		void read_restart_settings(FILE *) override;
		void write_data(FILE *) override;
		void write_data_all(FILE *) override;
		double single(int, int, int, int, double, double, double, double &) override;
		double single_atomistic_obmd(double, int, int, int, double, double, double, double &); // PP
		void born_matrix(int, int, int, int, double, double, double, double &, double &) override;
		void *extract(const char *, int &) override;
	protected:
		double cut_lj_global, cut_coul_global;
		double **cut_lj, **cut_ljsq;
		double **cut_coul, **cut_coulsq;
		double **epsilon, **sigma;
		double **lj1, **lj2, **lj3, **lj4, **offset;

		// for reaction field
		double epsilon_rf_one;
		double **epsilon_rf;
		virtual void allocate();
};

}   // namespace LAMMPS_NS

#endif
#endif
```

#### IMPLEMENTATION FILE

In the `pair_lj_cut_rf.cpp`, the implementation of the `PairLJCutRF` class can be found. The `pair_lj.cut_rf.cpp` file also contains copyright and license statement with several include statements, where `#include "pair_lj_cut_rf.h"` for the class header is listed first. 

```C++
#include "pair_lj_cut_rf.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include <cmath>
#include <cstring>

#include <iostream>

using namespace LAMMPS_NS;
using namespace MathConst;
```

The first section of the implementation file is followed by constructor of the `PairLJCutRF` class, which takes the LAMMPS class instance pointer `*lmp` as an argument. 

```C++
/* ---------------------------------------------------------------------- */
PairLJCutRF::PairLJCutRF(LAMMPS *lmp) : Pair(lmp)
{
	born_matrix_enable = 1;
	writedata = 1;
}
```

In destructor the memory allocated to the objects of `lj/cut/rf` pair style is freed. 

```C++
/* ---------------------------------------------------------------------- */

PairLJCutRF::~PairLJCutRF()
{
	if (copymode) return;

	if (allocated) {
	    memory->destroy(setflag);
	    memory->destroy(cutsq);

	    memory->destroy(cut_lj);
	    memory->destroy(cut_ljsq);
	    memory->destroy(cut_coul);
	    memory->destroy(cut_coulsq);
	    memory->destroy(epsilon);
	    memory->destroy(sigma);
	    memory->destroy(lj1);
	    memory->destroy(lj2);
	    memory->destroy(lj3);
	    memory->destroy(lj4);
	    memory->destroy(offset);

	    memory->destroy(epsilon_rf);
	}
}
```

The arguments of the `pair_style lj/cut/rf` command are processed in the `settings()` member function. Here, global parameters (to define cutoffs) are set.

```C++
/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */
void PairLJCutRF::settings(int narg, char **arg)
{
	if (narg < 1 || narg > 2) error->all(FLERR, "Illegal pair_style command");

	cut_lj_global = utils::numeric(FLERR, arg[0], false, lmp);
	if (narg == 1) 
		cut_coul_global = cut_lj_global;
	else
	    cut_coul_global = utils::numeric(FLERR, arg[1], false, lmp);

	// reset cutoffs that have been explicitly set

	if (allocated) {
	    int i, j;
	    for (i = 1; i <= atom->ntypes; i++)
		    for (j = i; j <= atom->ntypes; j++)
			    if (setflag[i][j]) {
			        cut_lj[i][j] = cut_lj_global;
			        cut_coul[i][j] = cut_coul_global;
		        }
	}
}
```

The `allocate()` member function handles memory allocation.

```C++
/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairLJCutRF::allocate()
{
	allocated = 1;
	int np1 = atom->ntypes + 1;

	memory->create(setflag, np1, np1, "pair:setflag");
		for (int i = 1; i < np1; i++)
		    for (int j = i; j < np1; j++) setflag[i][j] = 0;

	memory->create(cutsq, np1, np1, "pair:cutsq");
	 memory->create(cut_lj, np1, np1, "pair:cut_lj");
	memory->create(cut_ljsq, np1, np1, "pair:cut_ljsq");
	memory->create(cut_coul, np1, np1, "pair:cut_coul");
	memory->create(cut_coulsq, np1, np1, "pair:cut_coulsq");
	memory->create(epsilon, np1, np1, "pair:epsilon");
	memory->create(sigma, np1, np1, "pair:sigma");
	memory->create(lj1, np1, np1, "pair:lj1");
	memory->create(lj2, np1, np1, "pair:lj2");
	memory->create(lj3, np1, np1, "pair:lj3");
	memory->create(lj4, np1, np1, "pair:lj4");
	memory->create(offset, np1, np1, "pair:offset");
	memory->create(epsilon_rf, np1, np1, "pair:epsilon_rf");
}
```

While the arguments of the `pair_style lj/cut/rf` command are processes in the `settings()` method, the arguments provided to the `pair_coeff` command are handled in the `coeff()` member function. A total of **5 or 7 arguments** should be provided to the `pair_coeff` command, where global cutoffs are set as default, i.e., `cut_lj_one = cut_lj_global` and `cut_coul_one = cut_coul_global`, unless the sixth and seventh parameter are specified. Therefore, if the `rc_lj` and `rc_rf` arguments are not provided, the global LJ and rf cutoffs specified in the `pair_style` command are utilized. If only one cutoff is given, it is applied for both interactions. If both cutoffs are provided, they are used as the LJ and rf cutoffs for this type pair. Ultimately, the 2D arrays, namely `epsilon`, `sigma`, `cut_lj`, `cut_coul`, and `epsilon_rf`, are created for each atom pair.

```C++
/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */
void PairLJCutRF::coeff(int narg, char **arg)
{
	if (narg < 5 || narg > 7) error->all(FLERR, "Incorrect args for pair coefficients");
	if (!allocated) allocate();

	int ilo, ihi, jlo, jhi;
	utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
	utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

	double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
	double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);

	double cut_lj_one = cut_lj_global;
	double cut_coul_one = cut_coul_global;
	epsilon_rf_one = utils::numeric(FLERR, arg[4], false, lmp);

	if (narg >= 6) {
	    cut_coul_one = cut_lj_one = utils::numeric(FLERR, arg[4], false, lmp);
	    epsilon_rf_one = utils::numeric(FLERR, arg[5], false, lmp);
    }
    if (narg == 7) {
	    cut_lj_one = utils::numeric(FLERR, arg[4], false, lmp);
	    cut_coul_one = utils::numeric(FLERR, arg[5], false, lmp);
	    epsilon_rf_one = utils::numeric(FLERR, arg[6], false, lmp);
    }

	int count = 0;
	for (int i = ilo; i <= ihi; i++) {
	    for (int j = MAX(jlo, i); j <= jhi; j++) {
	      epsilon[i][j] = epsilon_one;
	      sigma[i][j] = sigma_one;
	      cut_lj[i][j] = cut_lj_one;
	      cut_coul[i][j] = cut_coul_one;

	      // reaction field param
	      epsilon_rf[i][j] = epsilon_rf_one; 

	      setflag[i][j] = 1;
	      count++;
	    }
	}

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}
```

The completeness of the potential parameters is verified in the `init_one()` method. This method invokes the `Pair::mix_energy()` and `Pair:mix_distance()` member functions if the `setflag` for particle pairs is $0$ or if the command `pair_modify mix {geometric,arithemic,sixthpower}` is used. Additionally, the quantities `cut_ljsq`, `cut_coulsq`, `lj1`, `lj2`, `lj3`, and `lj4` are defined and stored as 2D arrays. The offset to shift the LJ potential value to zero at the cutoff distance is computed when the `pair_modify shift yes` command is employed. Finally, the potential parameters arrays are symmetrized. 

```C++
/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */
double PairLJCutRF::init_one(int i, int j)
{
	if (setflag[i][j] == 0) {
	    epsilon[i][j] = mix_energy(epsilon[i][i], epsilon[j][j], sigma[i][i], sigma[j][j]);
	    sigma[i][j] = mix_distance(sigma[i][i], sigma[j][j]);
	    cut_lj[i][j] = mix_distance(cut_lj[i][i], cut_lj[j][j]);
	    cut_coul[i][j] = mix_distance(cut_coul[i][i], cut_coul[j][j]);
	}

	double cut = MAX(cut_lj[i][j], cut_coul[i][j]);
	cut_ljsq[i][j] = cut_lj[i][j] * cut_lj[i][j];
    cut_coulsq[i][j] = cut_coul[i][j] * cut_coul[i][j];

	lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
	lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
	lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
	lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

	if (offset_flag && (cut_lj[i][j] > 0.0)) {
	    double ratio = sigma[i][j] / cut_lj[i][j];
	    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
    } else
	    offset[i][j] = 0.0;

	cut_ljsq[j][i] = cut_ljsq[i][j];
	cut_coul[j][i] = cut_coul[i][j]; 
	cut_coulsq[j][i] = cut_coulsq[i][j];
	lj1[j][i] = lj1[i][j];
	lj2[j][i] = lj2[i][j];
	lj3[j][i] = lj3[i][j];
	lj4[j][i] = lj4[i][j];

	epsilon_rf[j][i] = epsilon_rf[i][j];

	offset[j][i] = offset[i][j];

	// compute I,J contribution to long-range tail correction
	// count total # of atoms of type I and J via Allreduce

	if (tail_flag) {
	    int *type = atom->type;
	    int nlocal = atom->nlocal;

	    double count[2], all[2];
	    count[0] = count[1] = 0.0;
	    for (int k = 0; k < nlocal; k++) {
	      if (type[k] == i) count[0] += 1.0;
	      if (type[k] == j) count[1] += 1.0;
	    }
	    MPI_Allreduce(count, all, 2, MPI_DOUBLE, MPI_SUM, world);

	    double sig2 = sigma[i][j] * sigma[i][j];
	    double sig6 = sig2 * sig2 * sig2;
	    double rc3 = cut_lj[i][j] * cut_lj[i][j] * cut_lj[i][j];
	    double rc6 = rc3 * rc3;
	    double rc9 = rc3 * rc6;
	    etail_ij =
        8.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 * (sig6 - 3.0 * rc6) / (9.0 * rc9);
	    ptail_ij = 16.0 * MY_PI * all[0] * all[1] * epsilon[i][j] * sig6 * (2.0 * sig6 - 3.0 * rc6) /
        (9.0 * rc9);
    }

	return cut;
}
```

The pairwise forces and potential between atoms $i$ and $j$ are computed in the `compute()` method. The charge (`q`), position vector components (in $x$-, $y$-, $z$-directions, i.e., `x[][0]`, `x[][1]`, `x[][2]`, respectively), and neighbors of the central atom $i$ are extracted by applying a for loop over all central atoms. Applying an additional for loop over neighbors $j$ of atom $i$, the square of the distance between atoms $i$ and $j$ is computed and stored in variable `rsq`. The computed distance is compared with the cutoff distance, and if it is smaller, force between atoms $i$ and $j$ is evaluated. Focusing only on the computation of the rf contribution, the following variables are defined:
- `rf_fctr_0`
  $$ 
  rf\_fctr\_0 = \frac{1}{r_{c}^{3}}
  $$
- `rf_fctr_1`
  $$
  rf\_fctr\_1 = \varepsilon_{rf} - 1
  $$
- `rf_fctr_2`
  $$
  rf\_fctr\_2 = 1 + 2\varepsilon_{rf}
  $$

and the `forcecoul` variable is calculated using

$$
F_{rf} (r_{i\alpha j\beta}) = \frac{q_{i\alpha} q_{j\beta}}{4 \pi
\varepsilon_0} \left[ \frac{1}{r{_{i\alpha j\beta}}^{3}} - \frac{1}
{r{_c}^{3}} \frac{2(\varepsilon_{rf}-1)}{2\varepsilon_{rf} + 1} \right] %
\mathbf{r}_{i\alpha j\beta}.
$$
The obtained force is scaled by `factor_coul` and added to the `fpairlj`, defining the `fpair` variable. The `fpair` variable is further multiplied with the displacement vector components (`delx`, `dely`, and `delz`), and forces acting in $x$-, $y$-, and $z$-directions are assigned to atom $i$. If command `newton_pair on` is specified, forces acting in the opposite directions are applied to the atom $j$. When energy contributions are requested (i.e., when `eflag = 1`), the energy interaction between atoms $i$ and $j$ is also computed. 

```C++
/* ---------------------------------------------------------------------- */

void PairLJCutRF::compute(int eflag, int vflag)
{
	std::cout<<"compute"<<"\n";
	int i, j, ii, jj, inum, jnum, itype, jtype;
	double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
	double rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
	int *ilist, *jlist, *numneigh, **firstneigh;

	double fpairlj, fpaircoul;
	evdwl = ecoul = 0.0;
	ev_init(eflag, vflag);

	double **x = atom->x;
	double **f = atom->f;
	double *q = atom->q;
	int *type = atom->type;
	int nlocal = atom->nlocal;
	double *special_coul = force->special_coul;
	double *special_lj = force->special_lj;
	int newton_pair = force->newton_pair;
	double qqrd2e = force->qqrd2e;

	inum = list->inum;
	ilist = list->ilist;
	numneigh = list->numneigh;
	firstneigh = list->firstneigh;

	// loop over neighbors of my atoms

	for (ii = 0; ii < inum; ii++) {
	    i = ilist[ii];
	    qtmp = q[i];
	    xtmp = x[i][0];
	    ytmp = x[i][1];
	    ztmp = x[i][2];
	    itype = type[i];
	    jlist = firstneigh[i];
	    jnum = numneigh[i];

	    for (jj = 0; jj < jnum; jj++) {
		    j = jlist[jj];
		    factor_lj = special_lj[sbmask(j)];
		    factor_coul = special_coul[sbmask(j)];
		    j &= NEIGHMASK;

		    delx = xtmp - x[j][0];
		    dely = ytmp - x[j][1];
		    delz = ztmp - x[j][2];
			rsq = delx * delx + dely * dely + delz * delz;
		    jtype = type[j];

		    if (rsq < cutsq[itype][jtype]) {
		        r2inv = 1.0 / rsq;

		        // terms for reaction field 
		        double rf_fctr_0 = 1.0 / pow(sqrt(rsq),3.0);
		        double rf_fctr_1 = (epsilon_rf[itype][jtype] - 1.0);
		        double rf_fctr_2 = 1.0 + 2.0 * epsilon_rf[itype][jtype];

		        if (rsq < cut_coulsq[itype][jtype]) {
		          forcecoul = (qqrd2e * qtmp * q[j]) * ((r2inv*sqrt(r2inv)) - (1.0 / pow(cut_coul[itype][jtype],3.0) * (2.0 * rf_fctr_1 / rf_fctr_2)));
		        }
		        else {
			        forcecoul = 0.0;
		        }
		
		        if (rsq < cut_ljsq[itype][jtype]) {
			        r6inv = r2inv * r2inv * r2inv;
			        forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
		        } else
			        forcelj = 0.0;
		
		        fpairlj = factor_lj * forcelj * r2inv;
		        fpaircoul = factor_coul * forcecoul; 
		        fpair = fpairlj + fpaircoul;

		        f[i][0] += delx * fpair;
		        f[i][1] += dely * fpair;
		        f[i][2] += delz * fpair;
		        if (newton_pair || j < nlocal) {
		          f[j][0] -= delx * fpair;
		          f[j][1] -= dely * fpair;
		          f[j][2] -= delz * fpair;
		        }

		        if (f[i][0] > 10000 || f[i][1] > 10000 || f[i][2] > 10000) {
		          std::cout<<"PairLJCutRF::compute"<<"\n";
		          std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
		          std::cout<<"f[i][1]: "<<f[i][1]<<"\n";
		          std::cout<<"f[i][2]: "<<f[i][2]<<"\n";
		        }

		        if (eflag) {
			        if (rsq < cut_coulsq[itype][jtype]) {
		            ecoul = (qqrd2e * qtmp * q[j]) * sqrt(r2inv) * (1.0 + (rf_fctr_1 / rf_fctr_2) * (pow(sqrt(rsq)/cut_coul[itype][jtype],3.0)))
		            - (qqrd2e * qtmp * q[j]) * (1.0 / cut_coul[itype][jtype]) * (3.0 * epsilon_rf[itype][jtype] / rf_fctr_2); 
		            ecoul *= factor_coul;
		        }
		        else {
			        ecoul = 0.0;
		        }
		        if (rsq < cut_ljsq[itype][jtype]) {
		            evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
		            evdwl *= factor_lj;
		        } else
		            evdwl = 0.0;
		        }

		        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, ecoul, fpair, delx, dely, delz);
		    }
	    }
	}

if (vflag_fdotr) virial_fdotr_compute();
}
```

Other classes can access potential energy by calling `single()` member function. Again, focusing only on the rf method, the energy is computed as

$$
U_{rf} (r_{i\alpha j\beta}) = \frac{q_{i\alpha} q_{j\beta}}{4 \pi
\varepsilon_0 r_{i\alpha j\beta}} \left[ 1 + \frac{\varepsilon_{rf} - 1}
{2\varepsilon_{rf} +1} \left(\frac{r_{i\alpha j\beta}}{r_c} \right)^3
\right] - \frac{q_{i\alpha} q_{j\beta}}{4 \pi \varepsilon_0 r_c}
\frac{3\varepsilon_{rf}}{2\varepsilon_{rf} +1}
$$
and scaled by factor `factor_coul`. The `single()` method returns the interaction energy stored in variable `eng`, which represents a sum of the LJ and rf energy terms.

```C++
/* ---------------------------------------------------------------------- */
double PairLJCutRF::single(int i, int j, int itype, int jtype, double rsq, double factor_coul,
                                double factor_lj, double &fforce)
{
	double r2inv, r6inv, forcecoul, forcelj, phicoul, philj;
	double rf_fctr_1 = (epsilon_rf[itype][jtype] - 1.0);
    double rf_fctr_2 = 1.0 + 2.0 * epsilon_rf[itype][jtype];
	r2inv = 1.0 / rsq;

	if (rsq < cut_coulsq[itype][jtype]) {
	    forcecoul = (force->qqrd2e * atom->q[i] * atom->q[j]) 
	    * (r2inv*sqrt(r2inv) 
	    - (1.0 / pow(cut_coul[itype][jtype],3.0) * (2.0 * rf_fctr_1 / rf_fctr_2))); 
	}
	else {
	    forcecoul = 0.0;
	}
	if (rsq < cut_ljsq[itype][jtype]) {
	    r6inv = r2inv * r2inv * r2inv;
	    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
    } else
	    forcelj = 0.0;

	fforce = factor_coul * forcecoul + factor_lj * forcelj * r2inv;

	double eng = 0.0;
	if (rsq < cut_coulsq[itype][jtype]) {
	    phicoul = force->qqrd2e * atom->q[i] * atom->q[j] * sqrt(r2inv) * (1.0 + (rf_fctr_1 / rf_fctr_2) * (pow(sqrt(rsq)/cut_coul[itype][jtype],3.0))) 
	    - force->qqrd2e * atom->q[i] * atom->q[j] * (1.0 / cut_coul[itype][jtype]) * (3.0 * epsilon_rf[itype][jtype] / rf_fctr_2); 
	    eng += factor_coul * phicoul;
	}
	if (rsq < cut_ljsq[itype][jtype]) {
	    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
	    eng += factor_lj * philj;
	}

  return eng;
}
```

To compute the energy of a particle to be inserted during OBMD simulations, an additional method called `single_atomistic_obmd()` has been implemented. The key difference between `single_atomistic_obmd()` and the `single()` method lies in the arguments that these member functions accept. While the `single()` member function takes (local) indices of atoms $i$ and $j$, the `single_atomistic_obmd()` function requires the charge of particle $i$ (`qi`) to be inserted instead of index $i$. It is important to note that in this case, the particle "$i$" does not yet exist. 

```C++
/* ---------------------------------------------------------------------- */
double PairLJCutRF::single_atomistic_obmd(double qi, int j, int itype, int jtype, double rsq, double factor_coul,
                                double factor_lj, double &fforce)
{
	double r2inv, r6inv, forcecoul, forcelj, phicoul, philj;

	double rf_fctr_1 = (epsilon_rf[itype][jtype] - 1.0);
	double rf_fctr_2 = 1.0 + 2.0 * epsilon_rf[itype][jtype];

	 r2inv = 1.0 / rsq;
	if (rsq < cut_coulsq[itype][jtype]) {
	    forcecoul = (force->qqrd2e * qi * atom->q[j]) 
    * (r2inv*sqrt(r2inv) 
	    - (1.0 / pow(cut_coul[itype][jtype],3.0) * (2.0 * rf_fctr_1 / rf_fctr_2))); 
	}
	else {
	    forcecoul = 0.0;
    }
    if (rsq < cut_ljsq[itype][jtype]) {
	    r6inv = r2inv * r2inv * r2inv;
	    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
    } else
	    forcelj = 0.0;

	fforce = factor_coul * forcecoul + factor_lj * forcelj * r2inv;

	double eng = 0.0;
	if (rsq < cut_coulsq[itype][jtype]) {
		phicoul = force->qqrd2e * qi * atom->q[j] * sqrt(r2inv) * (1.0 + (rf_fctr_1 / rf_fctr_2) * (pow(sqrt(rsq)/cut_coul[itype][jtype],3.0))) 
	    - force->qqrd2e * qi * atom->q[j] * (1.0 / cut_coul[itype][jtype]) * (3.0 * epsilon_rf[itype][jtype] / rf_fctr_2); 
	    eng += factor_coul * phicoul;
    }
	if (rsq < cut_ljsq[itype][jtype]) {
	    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
	    eng += factor_lj * philj;
    }

  return eng;
}
```
