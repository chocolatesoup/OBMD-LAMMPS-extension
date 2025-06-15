/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_atomistic_obmd.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix.h"
#include "lattice.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "random_park.h"
#include "region.h"
#include "update.h"
#include "pair.h"
#include "group.h"
#include "force.h"
#include "input.h"
#include "variable.h"
#include "neighbor.h"
#include "compute.h"

#include "pair_lj_cut_rf.h" // PPapez COMMENT: using this potential (option ? ? ?)

#include <cmath>
#include <cstring> // char[]
#include <iostream> // standard output

#include <fenv.h> // PPapez COMMENT: used for debug

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{ATOM,MOLECULE};
enum{DIST_UNIFORM,DIST_GAUSSIAN};
enum{NONE,CONSTANT,EQUAL};

#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

FixATMOBMD::FixATMOBMD(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(nullptr), idregion2(nullptr), idregion3(nullptr),
  idregion4(nullptr), idregion5(nullptr), idregion6(nullptr), /* idregion7(nullptr),
  idregion8(nullptr),*/ idshake(nullptr), onemols(nullptr), molfrac(nullptr), 
  coords(nullptr), imageflags(nullptr), fixshake(nullptr), random(nullptr), list(nullptr), mark(nullptr)
{

  if (narg < 11) error->all(FLERR,"Illegal fix atomistic/obmd command (check # of arguments)"); // after options

  restart_global = 1; // 1 if Fix saves global state
  time_depend = 1; // 1 if requires continuous timestepping
  size_vector = 1; // length of global vector

  dpstr = nullptr;
  freqstr = nullptr;

  // other args
  ntype = utils::inumeric(FLERR, arg[3], false, lmp);
  nfreq = utils::inumeric(FLERR, arg[4], false, lmp); // how often insertion is performed
  seed = utils::inumeric(FLERR, arg[5], false, lmp);

  if (seed <= 0) error->all(FLERR, "Illegal seed in fix atomistic/obmd command");

  // obmd settings 
  pxx = utils::numeric(FLERR, arg[6], false, lmp);
  pxy = utils::numeric(FLERR, arg[7], false, lmp);
  pxz = utils::numeric(FLERR, arg[8], false, lmp);
  // can be given as variables
  // pressure amplitude
  if (strstr(arg[9], "v_") == arg[9]) {
    int n = strlen(&arg[9][2]) + 1; // v_name value
    dpstr = new char[n];
    strcpy(dpstr, &arg[9][2]); // second
  }
  else {
    dpxx = utils::numeric(FLERR, arg[9], false, lmp);
    dpstyle = CONSTANT;
  }
  // frequency
  if (strstr(arg[10], "v_") == arg[10]) {
    int n = strlen(&arg[10][2]) + 1;
    freqstr = new char[n];
    strcpy(freqstr, &arg[10][2]);
  }
  else {
    freq = utils::numeric(FLERR, arg[10], false, lmp);
    freqstyle = CONSTANT;
  }

  // options
  varflag = CONSTANT;
  options(narg-11,&arg[11]); // 11 nargs

  // continue
  /* ------------------------------------------------------------- */
  // only molecules
  /* if (mode1 == ATOM || mode2 == ATOM || mode3 == ATOM) {
    error->all(FLERR, "Fix at/obmd works only for mode == MOLECULE");
  } */
  /* ------------------------------------------------------------- */

  // error checks on iregion & iregion2 (buffers) and their extent being inside simulation box
  if (!iregion || !iregion2) error->all(FLERR,"Must specify region1 and region2 in fix atomistic/obmd");
  if (iregion->bboxflag == 0 || iregion2->bboxflag == 0) error->all(FLERR,"Fix atomistic/obmd region does not support a bounding box");
  if (iregion->dynamic_check() || iregion2->dynamic_check()) error->all(FLERR,"Fix atomistic/obmd regions cannot be dynamic"); // in region.cpp 1 is returned if region is dynamic

  xlo = iregion->extent_xlo;
  xhi = iregion->extent_xhi;
  ylo = iregion->extent_ylo;
  yhi = iregion->extent_yhi;
  zlo = iregion->extent_zlo;
  zhi = iregion->extent_zhi;

  // omitting triclinic geometry 
  /* std::cout<<"FixATOBMD(): boxlo[0]: "<<domain->boxlo[0]<<" boxlo[1]: "<<domain->boxlo[1]<<" boxlo[2]: "<<domain->boxlo[2]<<"\n";
  std::cout<<"FixATOBMD(): boxhi[0]: "<<domain->boxhi[0]<<" boxhi[1]: "<<domain->boxhi[1]<<" boxhi[2]: "<<domain->boxhi[2]<<"\n";
  std::cout<<"FixATOBMD(): xlo: "<<xlo<<" xhi: "<<xhi<<" ylo: "<<ylo<<" yhi: "<<yhi<<" zlo: "<<zlo<<" zhi: "<<zhi<<"\n"; */
  if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
      ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
      zlo < domain->boxlo[2] || zhi > domain->boxhi[2]) {
        error->all(FLERR,"Insertion region extends outside simulation box");
  }

  // error check and further setup for mode = MOLECULE
  if (atom->tag_enable == 0) error->all(FLERR,"Cannot use fix atomistic/obmd unless atoms have IDs");

  // this is just a check and computation of geometric center of mass
  // assuming mode == MOLECULE 
  for (int i = 0; i < nmol; i++) {
    if (onemols[i]->xflag == 0)
      error->all(FLERR,"Fix atomistic/obmd molecule must have coordinates");
    if (onemols[i]->typeflag == 0)
          error->all(FLERR,"Fix atomistic/obmd molecule must have atom types");
    if (ntype+onemols[i]->ntypes <= 0 ||
        ntype+onemols[i]->ntypes > atom->ntypes)
        error->all(FLERR,"Invalid atom type in fix atomistic/obmd molecule command");

    if (atom->molecular == Atom::TEMPLATE && onemols != atom->avec->onemols)
      error->all(FLERR,"Fix atomistic/obmd molecule template ID must be same " "as atom_style template ID");
    
    onemols[i]->check_attributes(); // old check_attributes(0)
    // fix atomistic/obmd uses geoemetric center of molecule for insertion 
    // yeah, why not 
    onemols[i]->compute_center(); // in molecule.cpp // mol i is taken
  }

  // setup of coords and imageflags array
  // for insertion and iteration of coordinates
  natom_max = 0; 
  // molecule
  if (molflag) {
    for (int i = 0; i < nmol; i++) {
      natom_max = MAX(natom_max,onemols[i]->natoms);
    }
  memory->create(coords,natom_max,3,"atomistic/obmd:coords"); 
  memory->create(imageflags,natom_max,"atomistic/obmd:imageflags");
  }
  
  // std::cout<<"FixATOBMD(): natom_max: "<<natom_max<<"\n";
  
  // omitting scaling 

  // random number generator, same for all procs
  // warm up the generator 30x to avoid correlations in first-particle
  // positions if runs are repeated with consecutive seeds
  random = new RanPark(lmp, seed);
  for (int ii = 0; ii < 30; ii++) random->uniform();

  // set up reneighboring
  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nfirst = next_reneighbor;
  ninserted = 0;
  
  ndeleted = 0;
  nmax = 0;
  list = nullptr;
  mark = nullptr;

  // PPapez COMMENT: track deletion of particles
  delete_left_file.open("./../out/data/delete_left.out");
  delete_right_file.open("./../out/data/delete_right.out");
  insert_left_file.open("./../out/data/insert_left.out");
  insert_right_file.open("./../out/data/insert_right.out");
}

/* ---------------------------------------------------------------------- */

FixATMOBMD::~FixATMOBMD()
{
  delete random;
  delete [] molfrac;
  delete [] idshake;
  delete [] idregion;
  delete [] idregion2;
  delete [] idregion3;
  delete [] idregion4;
  delete [] idregion5;
  delete [] idregion6;
  /* delete [] idregion7;
  delete [] idregion8; */
  memory->destroy(coords);
  memory->destroy(imageflags);
  memory->destroy(list);
  memory->destroy(mark);

  // close files
  delete_left_file.close();
  delete_right_file.close();
  insert_left_file.close();
  insert_right_file.close();
}

/* ---------------------------------------------------------------------- */

int FixATMOBMD::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::init()
{ 
  
  // check variables (dpstr and freqstr)
  if (dpstr) {
    dpvar = input->variable->find(dpstr);
    if (dpvar < 0) error->all(FLERR, "Variable {} for fix atomistic/obmd does not exist", dpstr);
    if (input->variable->equalstyle(dpvar)) dpstyle = EQUAL;
    else if (input->variable->atomstyle(dpvar)) dpstyle = ATOM;
    else error->all(FLERR, "Variable {} for fix atomistic/obmd is invalid", dpstr);
  }
  if (freqstr) {
    freqvar = input->variable->find(freqstr);
    if (freqvar < 0) error->all(FLERR, "Variable {} for fix atomistic/obmd does not exist", freqstr);
    if (input->variable->equalstyle(freqvar)) freqstyle = EQUAL;
    else if (input->variable->atomstyle(freqvar)) freqstyle = ATOM;
    else error->all(FLERR, "Variable {} for fix atomistic/obmd is invalid", freqstr);
  }
  
  if (dpstyle == ATOM || freqstyle == ATOM) {
    std::cout<<"FixATOBMD::init(): varflag = ATOM"<<"\n";
    varflag = ATOM;
  }
  else if (dpstyle == EQUAL || freqstyle == EQUAL) {
    std::cout<<"FixATOBMD::init(): varflag = EQUAL"<<"\n";
    varflag = EQUAL;
  }
  else {
    std::cout<<"FixATOBMD::init(): varflag = CONSTANT"<<"\n";
    varflag = CONSTANT;
  }

  // get regions
  iregion = domain->get_region_by_id(idregion); // LEFT BUFFER
  if (!iregion)
    error->all(FLERR, "Region ID for fix obmd does not exist");

  iregion2 = domain->get_region_by_id(idregion2); // RIGHT BUFFER
  if (!iregion2)
    error->all(FLERR, "Region ID for fix obmd does not exist");
    
  iregion3 = domain->get_region_by_id(idregion3); // SHEAR FLOW LEFT
  if (!iregion3)
    error->all(FLERR, "Region ID for fix obmd does not exist");

  iregion4 = domain->get_region_by_id(idregion4); // SHEAR FLOW RIGHT
  if (!iregion4)
    error->all(FLERR, "Region ID for fix obmd does not exist");

  iregion5 = domain->get_region_by_id(idregion5); // INSERTION LEFT
  if (!iregion5)
    error->all(FLERR, "Region ID for fix obmd does not exist");
 
  iregion6 = domain->get_region_by_id(idregion6); // INSERTION RIGHT
  if (!iregion6)
    error->all(FLERR, "Region ID for fix obmd does not exist");

  /* iregion7 = domain->get_region_by_id(idregion7); // SOUND LEFT
  if (!iregion7)
    error->all(FLERR, "Region ID for fix atomistic/obmd does not exist");

  iregion8 = domain->get_region_by_id(idregion8); // SOUND RIGHT
  if (!iregion8)
    error->all(FLERR, "Region ID for fix atomistic/obmd does not exist"); */
      
  // check that no deletable atoms are in atom->firstgroup
  // deleting such an atom would not leave firstgroup atoms first
  if (atom->firstgroup >= 0) {
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int firstgroupbit = group->bitmask[atom->firstgroup];

    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if ((mask[i] & groupbit) && (mask[i] && firstgroupbit)) flag = 1;

    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);

    if (flagall)
      error->all(FLERR,"Cannot delete atoms in atom_modify first group");
  }
  
  // if molflag not set, warn if any deletable atom has a mol ID
  if (molflag == 0 && atom->molecule_flag) {
    tagint *molecule = atom->molecule;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (molecule[i]) flag = 1;
    int flagall;
    MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
    if (flagall && comm->me == 0)
      error->warning(FLERR,"Fix atomistic/obmd may delete atom with non-zero molecule ID");
  }

  if (molflag && atom->molecule_flag == 0)
      error->all(FLERR,"Fix atomistic/obmd molecule requires atom attribute molecule"); 

  // PPapez COMMENT: using SHAKE instead of RATTLE
  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one
  fixshake = nullptr;
  if (shakeflag) {
    fixshake = modify->get_fix_by_id(idshake);
    if (!fixshake) error->all(FLERR,"Fix deposit shake fix ID {} does not exist", idshake);
    int tmp;
    if (onemols != (Molecule **) fixshake->extract("onemol",tmp))
      error->all(FLERR,"Fix atomistic/obmd and fix shake are not using the same molecule template ID");
    /* int ifix = modify->find_fix(idshake);
    if (ifix < 0) error->all(FLERR,"Fix atomistic/obmd shake fix does not exist");
    fixshake = modify->fix[ifix];
    int tmp;
    if (onemols != (Molecule **) fixshake->extract("onemol",tmp))
      error->all(FLERR,"Fix atomistic/obmd and fix shake not using " "same molecule template ID"); */
  }
}

/* ----------------------------------------------------------------------
   setup
------------------------------------------------------------------------- */

void FixATMOBMD::setup(int vflag)
{ 
  if (utils::strmatch(update->integrate_style,"^verlet"))
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   perform particle/molecule insertions/deletions
------------------------------------------------------------------------- */

void FixATMOBMD::pre_exchange()
{ 
  // /* ==================================================================== */
  // int *type = atom->type;
  // int *tag = atom->tag;
  // double **v = atom->v;
  // double **x = atom->x;
  // double *mass = atom->mass;

  // double linear_momentum_pe[3];
  // MathExtra::zero3(linear_momentum_pe);
  // int nlocal302 = atom->nlocal;
  // for (int i = 0; i < nlocal302; i++) {
  //   if (tag[i] == 7386) {
  //     std::cout<<"it is particle #7386"<<"\n";
  //   }
  //   linear_momentum_pe[0] += mass[type[i]] * v[i][0];
  //   linear_momentum_pe[1] += mass[type[i]] * v[i][1];
  //   linear_momentum_pe[2] += mass[type[i]] * v[i][2];
  // }

  // std::cout<<"pre_exchange ::: START :: before try_deleting and try_inserting"<<"\n";
  // std::cout<<"nlocal302: "<<nlocal302<<"\n";
  // std::cout<<"linear_momentum_pe[0]: "<<linear_momentum_pe[0]<<"\n";
  // std::cout<<"linear_momentum_pe[1]: "<<linear_momentum_pe[1]<<"\n";
  // std::cout<<"linear_momentum_pe[2]: "<<linear_momentum_pe[2]<<"\n"; 
  // /* ==================================================================== */

  int cnt_left, cnt_right, stev_left, stev_right;
  double ninsert_left, ninsert_right;

  // calculate box dimensions
  double lx = domain->boxhi[0] - domain->boxlo[0]; // length of x
  double ly = domain->boxhi[1] - domain->boxlo[1];
  double lz = domain->boxhi[2] - domain->boxlo[2];
  
  if (varflag != CONSTANT) {
    if (dpstyle == EQUAL) dpxx = input->variable->compute_equal(dpvar);
    if (freqstyle == EQUAL) freq = input->variable->compute_equal(freqvar);
  }  
  
  // initialize momentum (out \& in)
  MathExtra::zero3(vnewl);
  MathExtra::zero3(vnewr);
  /* std::cout<<"vnewl[0]: "<<vnewl[0]<<" vnewr[0]: "<<vnewr[0]<<"\n";
  std::cout<<"vnewl[1]: "<<vnewl[1]<<" vnewr[1]: "<<vnewr[1]<<"\n";
  std::cout<<"vnewl[2]: "<<vnewl[2]<<" vnewr[2]: "<<vnewr[2]<<"\n"; */

  // OBMD
  // DELETING after particle crosses open end of simulation box
  try_deleting(iregion, vnewl, vnewr); // whole LEFT buffer
  try_deleting(iregion2, vnewl, vnewr); // whole RIGHT buffer  

  // count number of particles in LEFT and RIGHT buffer \& compute # of particles needed for insertion
  cnt_left = group->count(igroup, iregion);
  cnt_right = group->count(igroup, iregion2);
  ninsert_left = -static_cast<int>((static_cast<double>(cnt_left) / mol_len - alpha * nbuf) * update->dt / tau);
  ninsert_right = -static_cast<int>((static_cast<double>(cnt_right) / mol_len - alpha * nbuf) * update->dt / tau);

  // std::cout<<"FixATOBMD::pre_exchange(): cnt_left/mol_len: "<<cnt_left/mol_len<<" mol_len: "<<mol_len<<" alpha: "<<alpha<<" nbuf: "<<nbuf<<" tau: "<<tau<<" ninsert_left: "<<ninsert_left<<"\n";
  // std::cout<<"FixATOBMD::pre_exchange(): cnt_right/mol_len: "<<cnt_right/mol_len<<" mol_len: "<<mol_len<<" alpha: "<<alpha<<" nbuf: "<<nbuf<<" tau: "<<tau<<" ninsert_right: "<<ninsert_right<<"\n";
  
  // INSERTION of particles
  if (ninsert_left > 0) {
    // std::cout<<"FixATOBMD::pre_exchange()"<<"\n";
    // std::cout<<"ninsert_left: "<<ninsert_right<<"\n";
    try_inserting(iregion5, ninsert_left, vnewl, vnewr); 
  }
  if (ninsert_right > 0) {
    // std::cout<<"FixATOBMD::pre_exchange()"<<"\n";
    // std::cout<<"ninsert_right: "<<ninsert_right<<"\n";
    try_inserting(iregion6, ninsert_right, vnewl, vnewr); 
  } 
   
  // DELETING
  // additional check

  try_deleting(iregion, vnewl, vnewr);
  try_deleting(iregion2, vnewl, vnewr);

  // sum \& distribute
  MPI_Allreduce(vnewl, vnewl_all, 3, MPI_DOUBLE, MPI_SUM, world); // P_left
  MPI_Allreduce(vnewr, vnewr_all, 3, MPI_DOUBLE, MPI_SUM, world); // P_right
  
  // area of buffer - ROI interface
  double area = ly * lz;
  
  // calculates momentum forces on the left buffer
  // substracting momentum of deleted \& inserted particles (EQUATION)
  simulation_time += update->dt;
  pressure_wave = pxx + dpxx * sin(2.0 * MY_PI * freq * simulation_time);
  // std::cout<<"pxx: "<<pxx<<" dpxx: "<<dpxx<<" freq: "<<freq<<" simulation_time: "<<simulation_time<<"\n";

  // 2390.91
  // 0.000015
  MathExtra::zero3(momentumForce_left);
  MathExtra::zero3(momentumForce_right);
  momentumForce_left[0] = vnewl_all[0] / (update->dt) + pressure_wave * area; // = -vnewl_all[0] / (update->dt) + pressure_wave * area; // (pxx + dpxx) * area; // substract ? ? ? YES ! ! ! // unit_vector = (1.0, 0.0, 0.0) ... toward the center of the simulation box (L->R)
  momentumForce_left[1] = vnewl_all[1] / (update->dt); // = -vnewl_all[1] / (update->dt);
  momentumForce_left[2] = vnewl_all[2] / (update->dt); // = -vnewl_all[2] / (update->dt);
  shearForce_left[0] = 0.0;
  shearForce_left[1] = pxy * area;
  shearForce_left[2] = pxz * area;

  /* std::cout<<"vnewl_all[0]: "<<vnewl_all[0]<<" vnewr_all[0]: "<<vnewr_all[0]<<"\n";
  std::cout<<"vnewl_all[1]: "<<vnewl_all[1]<<" vnewr_all[1]: "<<vnewr_all[1]<<"\n";
  std::cout<<"vnewl_all[2]: "<<vnewl_all[2]<<" vnewr_all[2]: "<<vnewr_all[2]<<"\n";
  std::cout<<"update->dt: "<<update->dt<<"\n";
  std::cout<<"pressure_wave: "<<pressure_wave<<"\n";
  std::cout<<"area: "<<area<<"\n";
  std::cout<<"momentumForce_left[0]: "<<momentumForce_left[0]<<"\n";
  std::cout<<"momentumForce_left[1]: "<<momentumForce_left[1]<<"\n";
  std::cout<<"momentumForce_left[2]: "<<momentumForce_left[2]<<"\n"; 
  std::cout<<"shearForce_left[0]: "<<shearForce_left[0]<<"\n";
  std::cout<<"shearForce_left[1]: "<<shearForce_left[1]<<"\n";
  std::cout<<"shearForce_left[2]: "<<shearForce_left[2]<<"\n"; */


  // calculates momentum forces on the right buffer
  // substracting momentum of deleted \& inserted particles (EQUATION)
  momentumForce_right[0] = vnewr_all[0] / (update->dt) - pxx * area; // = -vnewr_all[0] / (update->dt) - pxx * area;  // unit_vector = (-1.0, 0.0, 0.0) ... toward the center of the simulation box (R->L)
  momentumForce_right[1] = vnewr_all[1] / (update->dt); // = -vnewr_all[1] / (update->dt);
  momentumForce_right[2] = vnewr_all[2] / (update->dt); // = -vnewr_all[2] / (update->dt);
  shearForce_right[0] = 0.0;
  shearForce_right[1] = -pxy * area;
  shearForce_right[2] = -pxz * area;  

  /* std::cout<<"momentumForce_right[0]: "<<momentumForce_right[0]<<"\n";
  std::cout<<"momentumForce_right[1]: "<<momentumForce_right[1]<<"\n";
  std::cout<<"momentumForce_right[2]: "<<momentumForce_right[2]<<"\n";
  std::cout<<"shearForce_right[0]: "<<shearForce_right[0]<<"\n";
  std::cout<<"shearForce_right[1]: "<<shearForce_right[1]<<"\n";
  std::cout<<"shearForce_right[2]: "<<shearForce_right[2]<<"\n"; */
  
  // /* ==================================================================== */
  // double linear_momentum_pee[3];
  // MathExtra::zero3(linear_momentum_pee);
  // std::cout<<"linear_momentum_pee[0]: "<<linear_momentum_pee[0]<<"\n";
  // std::cout<<"linear_momentum_pee[1]: "<<linear_momentum_pee[1]<<"\n";
  // std::cout<<"linear_momentum_pee[2]: "<<linear_momentum_pee[2]<<"\n"; 
  // int nlocal312 = atom->nlocal;
  // double **vel = atom->v;
  // int type1 = 0;
  // int type2 = 0;
  // for (int i = 0; i < nlocal312; i++) {

  //   if (type[i] == 1) {
  //     type1 += 1;
  //   }
  //   if (type[i] == 2) {
  //     type2 += 1;
  //   }

  //   if (i == 7246) {
  //     std::cout<<"it is i #7246"<<"\n";
  //     std::cout<<"tag[i]: "<<tag[i]<<"\n";
  //     std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
  //   }
  //   if (tag[i] == 7386){ 
  //     std::cout<<"it is particle #7386"<<"\n";
  //   }
  //   if (tag[i] == 7385){ 
  //     std::cout<<"it is particle #7385"<<"\n";
  //   }
  //   if (tag[i] == 7384){ 
  //     std::cout<<"it is particle #7384"<<"\n";
  //   }
  //   linear_momentum_pee[0] += mass[type[i]] * vel[i][0];
  //   linear_momentum_pee[1] += mass[type[i]] * vel[i][1];
  //   linear_momentum_pee[2] += mass[type[i]] * vel[i][2];
  // }

  // std::cout<<"pre_exchange ::: END ::: after try_deleting"<<"\n";
  // std::cout<<"nlocal312: "<<nlocal312<<"\n";
  // std::cout<<"linear_momentum_pee[0]: "<<linear_momentum_pee[0]<<"\n";
  // std::cout<<"linear_momentum_pee[1]: "<<linear_momentum_pee[1]<<"\n";
  // std::cout<<"linear_momentum_pee[2]: "<<linear_momentum_pee[2]<<"\n"; 
  // std::cout<<"type1: "<<type1<<" type2: "<<type2<<"\n";
  // /* ==================================================================== */

  next_reneighbor += nfreq; 

  // /* ==================================================================== */
  // double linear_momentum_peee[3];
  // MathExtra::zero3(linear_momentum_peee);
  // int nlocal322 = atom->nlocal;
  // for (int i = 0; i < nlocal322; i++) {
    
  //   if (tag[i] == 7386) {
  //     std::cout<<"it is particle #7386"<<"\n";
  //   }
  //   linear_momentum_peee[0] += mass[type[i]] * v[i][0];
  //   linear_momentum_peee[1] += mass[type[i]] * v[i][1];
  //   linear_momentum_peee[2] += mass[type[i]] * v[i][2];
  // }

  // std::cout<<"pre_exchange ::: END"<<"\n";
  // std::cout<<"nlocal322: "<<nlocal322<<"\n";
  // std::cout<<"linear_momentum_peee[0]: "<<linear_momentum_peee[0]<<"\n";
  // std::cout<<"linear_momentum_peee[1]: "<<linear_momentum_peee[1]<<"\n";
  // std::cout<<"linear_momentum_peee[2]: "<<linear_momentum_peee[2]<<"\n"; 
  // /* ==================================================================== */

} // closing pre_exchange()

/* ---------------------------------------------------------------------- */

void FixATMOBMD::try_deleting(Region* region, double *vnewl, double *vnewr) 
{ 
  int i, j, m, iwhichglobal, iwhichlocal;
  int ndel, ndeltopo[4];
  
  // grow list and mark arrays if necessary
  if (atom->nmax > nmax) {
    memory->destroy(list);
    memory->destroy(mark);
    nmax = atom->nmax;
    memory->create(list,nmax,"atomistic/obmd:list");
    memory->create(mark,nmax,"atomistic/obmd:mark");
  }

  // ncount = # of deletable atoms in region that I own
  // nall = # on all procs
  // nbefore = # on procs before me
  // list[ncount] = list of local indices of atoms I can delete ... local ! ! !
  region->prematch();

  double **x = atom->x;
  double **vel = atom->v;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  double *mass = atom->mass;
  int *type = atom->type;
  double boxl = domain->boxlo[0];
  double boxh = domain->boxhi[0];
  double fusher[3];

  // std::cout<<"nlocal in try_deleting: "<<nlocal<<"\n";

  int ncount = 0;
  for (i = 0; i < nlocal; i++) {
      if (x[i][0] < boxl || x[i][0] > boxh) {  
		    std::cout<<"Deleting x[i][0]: "<<x[i][0]<<" i: "<<i<<" type[i]: "<<type[i]<<" tag[i]: "<<tag[i]<<"\n";
        list[ncount++] = i;	// add index
	  }
  }
  int nall, nbefore;
  MPI_Allreduce(&ncount, &nall, 1, MPI_INT, MPI_SUM, world);
  MPI_Scan(&ncount, &nbefore, 1, MPI_INT, MPI_SUM, world);

  nbefore -= ncount;

  // ndel = total # of atom deletions, in or out of region
  // ndeltopo[1,2,3,4] = ditto for bonds, angles, dihedrals, impropers
  // mark[] = 1 if deleted
  ndel = 0;
  for (i = 0; i < nlocal; i++) mark[i] = 0; // init

  // PPapez COMMENT: only MOLECULE
  // molecule deletions
  // bcast mol ID and delete all atoms in that molecule on any proc
  // update deletion count by total # of atoms in molecule
  // shrink list of eligible candidates as any of my atoms get marked
  // keep ndel,ndeltopo,ncount,nall,nbefore current after each mol deletions
  int me, proc, iatom, ndelone, ndelall, index;
  tagint imolecule;
  tagint *molecule = atom->molecule;
  int *molindex = atom->molindex;
  int *molatom = atom->molatom;
  int molecular = atom->molecular;
  Molecule **onemols = atom->avec->onemols;
	
  ndeltopo[0] = ndeltopo[1] = ndeltopo[2] = ndeltopo[3] = 0;

  while (nall) {
    // pick an iatom,imolecule on proc me to delete
    iwhichglobal = static_cast<int> (nall*random->uniform());
    if (iwhichglobal >= nbefore && iwhichglobal < nbefore + ncount) {
      iwhichlocal = iwhichglobal - nbefore;
      iatom = list[iwhichlocal];
      imolecule = molecule[iatom];
      me = comm->me;
    } 
    else {
      me = -1;
    }

    // bcast mol ID to delete all atoms from
    // if mol ID > 0, delete any atom in molecule and decrement counters
    // if mol ID == 0, delete single iatom
    // logic with ndeltopo is to count # of deleted bonds,angles,etc
    // for atom->molecular = Atom::MOLECULAR, do this for each deleted atom in molecule
    // for atom->molecular = Atom::TEMPLATE, use Molecule counts for just 1st atom in mol

    MPI_Allreduce(&me, &proc, 1, MPI_INT, MPI_MAX, world);
    MPI_Bcast(&imolecule, 1, MPI_LMP_TAGINT, proc, world);
    ndelone = 0;
    for (i = 0; i < nlocal; i++) {
      if (imolecule && molecule[i] == imolecule) {
        mark[i] = 1;
        ndelone++;

        if (molecular == Atom::MOLECULAR) {
          if (atom->avec->bonds_allow) { // bonds
            if (force->newton_bond) 
              ndeltopo[0] += atom->num_bond[i];
            else {
              for (j = 0; j < atom->num_bond[i]; j++) {
                if (tag[i] < atom->bond_atom[i][j]) ndeltopo[0]++;
              }
            }
          }
          if (atom->avec->angles_allow) { // angles
            if (force->newton_bond) 
              ndeltopo[1] += atom->num_angle[i];
            else {
              for (j = 0; j < atom->num_angle[i]; j++) {
                m = atom->map(atom->angle_atom2[i][j]);
                if (m >= 0 && m < nlocal) ndeltopo[1]++;
              }
            }
          }
          if (atom->avec->dihedrals_allow) { // dihedrals
            if (force->newton_bond) 
              ndeltopo[2] += atom->num_dihedral[i];
            else {
              for (j = 0; j < atom->num_dihedral[i]; j++) {
                m = atom->map(atom->dihedral_atom2[i][j]);
                if (m >= 0 && m < nlocal) ndeltopo[2]++;
              }
            }
          }
          if (atom->avec->impropers_allow) { // improper dihedrals
            if (force->newton_bond) 
              ndeltopo[3] += atom->num_improper[i];
            else {
              for (j = 0; j < atom->num_improper[i]; j++) {
                m = atom->map(atom->improper_atom2[i][j]);
                if (m >= 0 && m < nlocal) ndeltopo[3]++;
              }
            }
          }

        } else if (molecular == Atom::TEMPLATE) {
            if (molatom[i] == 0) {
              index = molindex[i];
              ndeltopo[0] += onemols[index]->nbonds;
              ndeltopo[1] += onemols[index]->nangles;
              ndeltopo[2] += onemols[index]->ndihedrals;
              ndeltopo[3] += onemols[index]->nimpropers;
            }
        }

      } else if (me == proc && i == iatom) {
        mark[i] = 1;
        ndelone++;
      }
    }

    // remove any atoms marked for deletion from my eligible list
    i = 0;
    while (i < ncount) {
      if (mark[list[i]]) {
        list[i] = list[ncount-1];
        ncount--;
      } else i++;
    }

    // update ndel,ncount,nall,nbefore
    // ndelall is total atoms deleted on this iteration
    // ncount is already correct, so resum to get nall and nbefore
    MPI_Allreduce(&ndelone, &ndelall, 1, MPI_INT, MPI_SUM, world);
    ndel += ndelall;
    MPI_Allreduce(&ncount, &nall, 1, MPI_INT, MPI_SUM, world);
    MPI_Scan(&ncount, &nbefore, 1, MPI_INT, MPI_SUM, world);
    nbefore -= ncount;
  }
  
  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms
  AtomVec *avec = atom->avec;
  
  for (i = nlocal - 1; i >= 0; i--) {
    if (mark[i]) {
      if(x[i][0] < 0.5*(boxh+boxl)) {
        std::cout<<"DELETE: type[i]: "<<type[i]<<" mass[type[i]]: "<<mass[type[i]]<<" vel[i][0]: "<<vel[i][0]<<" vel[i][1]: "<<vel[i][1]<<" vel[i][2]: "<<vel[i][2]<<"\n";
        vnewl[0] += mass[type[i]] * vel[i][0]; // -= mass[type[i]] * vel[i][0]; // outgoing momentum LEFT // particles
        vnewl[1] += mass[type[i]] * vel[i][1]; // -= mass[type[i]] * vel[i][1];
        vnewl[2] += mass[type[i]] * vel[i][2]; // -= mass[type[i]] * vel[i][2];
        // delete_left_file<<type[i]<<"\t"<<x[i][0]<<"\t"<<vel[i][0]<<"\t"<<"\n";
        // delete_left_file.flush();
      } 
      else {
        std::cout<<"DELETE: type[i]: "<<type[i]<<" mass[type[i]]: "<<mass[type[i]]<<" vel[i][0]: "<<vel[i][0]<<" vel[i][1]: "<<vel[i][1]<<" vel[i][2]: "<<vel[i][2]<<"\n";
        vnewr[0] += mass[type[i]] * vel[i][0]; // -= mass[type[i]] * vel[i][0]; // outgoing momentum RIGHT // particles
        vnewr[1] += mass[type[i]] * vel[i][1]; // -= mass[type[i]] * vel[i][1];
        vnewr[2] += mass[type[i]] * vel[i][2]; // -= mass[type[i]] * vel[i][2];
        // delete_right_file<<type[i]<<"\t"<<x[i][0]<<"\t"<<vel[i][0]<<"\t"<<"\n";
        // delete_right_file.flush();
      }
      avec->copy(atom->nlocal - 1, i, 1);
      atom->nlocal--;
    }
  }

  // reset global natoms and bonds, angles, etc
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts
  atom->natoms -= ndel;

  if (molflag) {
    int all[4];
    MPI_Allreduce(ndeltopo, all, 4, MPI_INT ,MPI_SUM, world);
    atom->nbonds -= all[0];
    atom->nangles -= all[1];
    atom->ndihedrals -= all[2];
    atom->nimpropers -= all[3];
  }

  if (ndel && (atom->map_style != Atom::MAP_NONE)) {
    atom->nghost = 0;
    atom->map_init();
    atom->map_set();
  }

  // statistics
  ndeleted += ndel;	
  // std::cout<<"ndeleted: "<<ndeleted<<"\n";
  // std::cout<<"vnewl[0]: "<<vnewl[0]<<" vnewl[1]: "<<vnewl[1]<<" vnewl[2]: "<<vnewl[2]<<"\n";
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::try_inserting(Region* iregion_var, int stev, double *vnewl, double *vnewr) 
{  
  // std::cout<<"try_inserting"<<"\n";
  int i, m, n, nlocalprev, imol, natom, flag, flagall, me, sum_all, ninsert;
  double coord[3], lamda[3], delx, dely, delz, rsq;
  double r[3], vnew[3], rotmat[3][3], quat[4];
  double *newcoord;
  double entmp, mtmp;
  int iter;
  double *sublo, *subhi; 
  
  double *mass = atom->mass;
  int dimension = domain->dimension;
  double boxl = domain->boxlo[0];
  double boxh = domain->boxhi[0];
    
  xlo = iregion_var->extent_xlo;
  xhi = iregion_var->extent_xhi;
  ylo = iregion_var->extent_ylo;
  yhi = iregion_var->extent_yhi;
  zlo = iregion_var->extent_zlo;
  zhi = iregion_var->extent_zhi;
  
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }

  // clear ghost count (and atom map) and any ghost bonus data
  // internal to AtomVec
  // same logic as beginning of Comm::exchange()
  // do it now b/c inserting atoms will overwrite ghost atoms

  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();
  
  if (!idnext) find_maxid(); // idnext is set to zero
      
  if(stev > 0) {
    for(ninsert = 0; ninsert != stev; ninsert++) { 
      int success = 0;
      int attempt = 0;      
      while (attempt < maxattempt) {

        // PPapez COMMENT: for debug
        // feenableexcept(FE_ALL_EXCEPT);

        attempt++;

        // choose random position for new particle within region
        // PPapez COMMENT: assuming DIST_UNIFORM

        do {
          coord[0] = xlo + random->uniform() * (xhi-xlo);
          coord[1] = ylo + random->uniform() * (yhi-ylo);
          coord[2] = zlo + random->uniform() * (zhi-zlo);
        } 
        while (iregion_var->match(coord[0],coord[1],coord[2]) == 0);
        
        // PPapez COMMENT: only MOLECULE insertions
        // coords = coords of all atoms
        // for molecule, perform random rotation around center pt
        // apply PBC so final coords are inside box
        // also modify image flags due to PBC

        double rng = random->uniform();
        imol = 0;
        while (rng > molfrac[imol]) imol++;
        natom = onemols[imol]->natoms;
        if (dimension == 3) {
          r[0] = random->uniform() - 0.5;
          r[1] = random->uniform() - 0.5;
          r[2] = random->uniform() - 0.5; 
        } 
        else {
          r[0] = r[1] = 0.0;
          r[2] = 1.0;
        }
        double theta = random->uniform() * MY_2PI;
        MathExtra::norm3(r); // unit vector 
        MathExtra::axisangle_to_quat(r, theta, quat); // quaternion
        MathExtra::quat_to_mat(quat, rotmat); // rotmat
        for (i = 0; i < natom; i++) {
          MathExtra::matvec(rotmat, onemols[imol]->dx[i], coords[i]); // getting coords after rotation of dx[i] using rotmat
          coords[i][0] += coord[0]; // rotated coords // vrstica - stolpec
          coords[i][1] += coord[1]; // rotated coords
          coords[i][2] += coord[2]; // rotated coords
          // PPapez COMMENT: due to this routine select / limit insertion region (reasonable dist from its edge)
          imageflags[i] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
          domain->remap(coords[i], imageflags[i]); // REMAP
        }
      
        // PPapez COMMENT: omitting near ! ! !
        // check distance between any existing atom and any inserted atom
        // if less than near, try again
        // use minimum_image() to account for PBC
      
        double **x = atom->x;
        int nlocal = atom->nlocal;
        flag = 0;

        // PPapez COMMENT: only USHER ! ! ! 
	      me = -1;
	      entmp = usher(iregion_var, coords, etgt, natom, imol, iter); 
      
      	if(entmp < etgt+EPSILON) {
	       	if(comm->me==0) std::cout << "USHER accepts at E = " << entmp << " in attempt No. " << attempt << " with " << iter << " iterations" << std::endl;
        } 
        else { 
	       	if(comm->me==0) std::cout << "USHER denies at E = " << entmp << " at attempt No. " << attempt << std::endl; 
	       	flag=1;
	      }
	      
        
        MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
        if (flagall) continue;
              
        // proceed with insertion
        nlocalprev = atom->nlocal;
        
        // choose random velocity for new particle
        // used for every atom in molecule
        // PPapez COMMENT: for now I'm inserting with zero velocity

        vnew[0] = 0.0; // vxlo + random->uniform() * (vxhi-vxlo);
        vnew[1] = 0.0; // vylo + random->uniform() * (vyhi-vylo);
        vnew[2] = 0.0; // vzlo + random->uniform() * (vzhi-vzlo);  
        
        // this loop first just checks if an entire molecule can be inserted
	      sum_all = check_mol_proc(coords, newcoord, lamda, sublo, subhi, dimension, natom, me);
        if(sum_all != natom) {
          if(comm->me==0) std::cout << "Can not insert the particle/molecule at timestep = " << update->ntimestep << " me = " << me << std::endl;
	        continue;
        }
          
        // happens rarely, e.g. if usher moves particle/molecule outside of the (open) boundaries
        int check;
        check = check_mol_region(iregion_var, coords, natom);
        if(check == 1) continue;       
        
        // check if new atoms are in my sub-box or above it if I am highest proc
        // if so, add atom to my list via create_atom()
        // initialize additional info about the atoms
        // set group mask to "all" plus fix group

        // NOT OPTIMAL SOLUTION ... TESTING FOR NOW
        // store all ilocal ids of newly inserted particles of a molecule
        // std::vector<int> locals;
        for (m = 0; m < natom; m++) {   
          if (domain->triclinic) {
            domain->x2lamda(coords[m],lamda);
            newcoord = lamda;
          } else newcoord = coords[m];
          flag = check_proc(coords, newcoord, sublo, subhi, dimension);             				
          if (flag) {
	          mtmp = mtot;
            atom->avec->create_atom(ntype+onemols[imol]->type[m],coords[m]); 
	                
            n = atom->nlocal - 1;
            atom->tag[n] = maxtag_all + m + 1;

            if (atom->molecule_flag) {
              if (onemols[imol]->moleculeflag) {
                atom->molecule[n] = maxmol_all + onemols[imol]->molecule[m];
              } 
              else {
                atom->molecule[n] = maxmol_all + 1;
              }
            }
            if (atom->molecular == Atom::TEMPLATE) {
              atom->molindex[n] = 0;
              atom->molatom[n] = m;
            }
            
            atom->mask[n] = 1 | groupbit;
            atom->image[n] = imageflags[m];
            atom->v[n][0] = vnew[0];
            atom->v[n][1] = vnew[1];
            atom->v[n][2] = vnew[2];

            onemols[imol]->quat_external = quat;
            atom->add_molecule_atom(onemols[imol], m, n, maxtag_all); // creating atom
            // std::cout<<"inserted particle: atom->tag[n]: "<<atom->tag[n]<<" atom->x[n][0]: "<<atom->x[n][0]<<"\n";
            // std::cout<<"atom->type[n]: "<<atom->type[n]<<"\n";
            // std::cout<<"atom->x[n][0]: "<<atom->x[n][0]<<" atom->x[n][1]: "<<atom->x[n][1]<<" atom->x[n][2]: "<<atom->x[n][2]<<"\n";
            // std::cout<<"atom->q[n]: "<<atom->q[n]<<"\n";
  
            modify->create_attribute(n);

          } // closing flag 
        } // closing while()
        
        if (shakeflag) {
          // std::cout<<"shakeflag in fix_atomistic_obmd"<<"\n";
          // std::cout<<"imol: "<<imol<<"\n";
          fixshake->set_molecule(nlocalprev,maxtag_all,imol,coord,vnew,quat);
        }

        success = 1;
        break;
      }

      if (!success && comm->me == 0) {
        // error->warning(FLERR,"Particle/molecule insertion was unsuccessful");
        ;
      }
      
      // reset global natoms,nbonds,etc
      // increment maxtag_all and maxmol_all if necessary
      // if global map exists, reset it now instead of waiting for comm
      // since other pre-exchange fixes may use it
      // invoke map_init() b/c atom count has grown
        
      if(success) {
        /// std::cout<<"FixATMOBMD::try_inserting()"<<"\n";
        // std::cout<<"if(success)"<<"\n";
        // std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
        // std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
        // std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 
        if(coords[0][0] < 0.5*(boxh + boxl)) {
          // std::cout<<"mtmp: "<<mtmp<<"\n";
	        vnewl[0] -= 0.0; // mtmp * vnew[0]; // momentum in
          vnewl[1] -= 0.0; // mtmp * vnew[1];
          vnewl[2] -= 0.0; // mtmp * vnew[2];
          insert_left_file<<"inserted"<<"\n";
          insert_left_file.flush();
        } else {
          vnewr[0] -= 0.0; // mtmp * vnew[0];
          vnewr[1] -= 0.0; // mtmp * vnew[1];
          vnewr[2] -= 0.0; // mtmp * vnew[2];		
          insert_right_file<<"inserted"<<"\n";
          insert_right_file.flush();
        }
        atom->natoms += natom;
        
        if (atom->natoms < 0)
          error->all(FLERR,"Too many total atoms");
        if (mode == MOLECULE) {
          atom->nbonds += onemols[imol]->nbonds;
          atom->nangles += onemols[imol]->nangles;
          atom->ndihedrals += onemols[imol]->ndihedrals;
          atom->nimpropers += onemols[imol]->nimpropers;
        }
        maxtag_all += natom;
        if (maxtag_all >= MAXTAGINT)
          error->all(FLERR,"New atom IDs exceed maximum allowed ID");
        if (mode == MOLECULE && atom->molecule_flag) {
          if (onemols[imol]->moleculeflag) {
            maxmol_all += onemols[imol]->nmolecules;
          } else {
            maxmol_all++;
          }
        } 
      }
      
      // rebuild atom map
      if (atom->map_style != Atom::MAP_NONE) {
        if (success) atom->map_init(); // update ! ! !
        atom->map_set();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::post_force(int vflag)
{
  if (update->ntimestep % nevery) { return ;}

  /* double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int *tag = atom->tag;
  
  for (int i = 0; i < nlocal; i++) {
    if (atom->tag[i] == 226) {
        std::cout<<"obmd::post_force"<<"\n";
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
    }
  } */

  // std::cout<<"FixATMOBMD::post_force(): stepflag: "<<stepflag<<"\n";
  // using parallel and perp
  // parallel for normal forces (and sound)
  // perp for shear
  // std::cout<<"FixATMOBMD::post_force(): step_parallel: "<<step_parallel<<"\n"; // 0 for smooth
  // std::cout<<"FixATMOBMD::post_force(): step_perp: "<<step_perp<<"\n"; // 1 for heaviside
  // std::cout<<"momentumForce_left[0]: "<<momentumForce_left[0]<<" momentumForce_left[1]: "<<momentumForce_left[1]<<" momentumForce_left[2]: "<<momentumForce_left[2]<<"\n";
  // std::cout<<"momentumForce_right[0]: "<<momentumForce_right[0]<<" momentumForce_right[1]: "<<momentumForce_right[1]<<" momentumForce_right[2]: "<<momentumForce_right[2]<<"\n";
  
  // parallel
  // std::cout<<"before reg_force for LEFT BUFFER"<<"\n";
  reg_force(vflag, iregion, momentumForce_left, step_parallel); // stepflag);
  // std::cout<<"before reg_force for RIGHT BUFFER"<<"\n";
  reg_force(vflag, iregion2, momentumForce_right, step_parallel); // stepflag);

  // tangential
  // std::cout<<"before reg_force_perp for LEFT BUFFER"<<"\n";
  // reg_force_perp(vflag, iregion3, shearForce_left, step_perp); // stepflag);
  // std::cout<<"before reg_force_perp for RIGHT BUFFER"<<"\n";
  // reg_force_perp(vflag, iregion4, shearForce_right, step_perp); // stepflag);  
}

/* --------------------------------------------------------------------------------------------- */

double FixATMOBMD::g_par_global(Region* region, int step)
{
  /* double **x = atom->x; */
  // if (region) region->prematch();
  
  /* int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double xh = domain->boxhi[0];
  double xl = domain->boxlo[0];
  double g_par_glob = 0.0;
  double g_par_glob_tmp = 0.0; 
  double y;
  double *mass = atom->mass;
  int  *type = atom->type;
  
  if(!step) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
        if(x[i][0] < xl+0.5*(xh+xl)) { //tp 23.10.2023
		  if(x[i][0] < xl+(1.0-g_fac)*buffer_size) g_par_glob += mass[type[i]]; // 1.0; //tp 23.10.2023
		  else {
		    y = g_fac_inv*MY_PI*(x[i][0]-buffer_size-xl)/(0.0-buffer_size)-MY_PI;//tp 23.10.2023
			g_par_glob += mass[type[i]] * 0.5*(1.0+cos(y));
		  }	
		} else {
		  if(x[i][0] > xh-(1.0-g_fac)*buffer_size) g_par_glob += mass[type[i]]; // 1.0;
		  else { 
		    y = g_fac_inv*MY_PI*(x[i][0]-(xh-buffer_size))/(buffer_size)-MY_PI;
			g_par_glob += mass[type[i]] * 0.5*(1.0+cos(y));
		  }  
		}
      }
    }

  MPI_Allreduce(&g_par_glob,&g_par_glob_tmp,1,MPI_DOUBLE,MPI_SUM,world);
  } else g_par_glob_tmp = group->count(igroup,region);
  
  return g_par_glob_tmp; */

  // using molecule's center-of-mass
  // std::cout<<"FixATOBMD::g_par_global()"<<"\n";
  if (region) region->prematch();
  
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  int *rep_atom = atom->rep_atom;
  double **cms = atom->cms_mol; 
  double *mass = atom->mass;
  int  *type = atom->type;
  double **x = atom->x;
  int imol; 
  int *tag = atom->tag;
  double **f = atom->f;

  double upper_x = domain->boxhi[0]; // upper
  double lower_x = domain->boxlo[0]; // lower
  // std::cout<<"upper_x: "<<upper_x<<" lower_x: "<<lower_x<<"\n";
  // init for distribution function
  double g_par_all = 0.0;
  double g_par_all_tmp = 0.0; 

  /* int counting_global = 0;
  // testing step function
  for (int i = 0; i < nlocal; i++) {
    if (rep_atom[i] == 0) { // PPapez COMMENT: only center-of-mass
      continue;
    }

    if (region && !region->match(cms[i][0], cms[i][1], cms[i][2])) continue; // PPapez COMMENT: is it here?

    imol = molecule[i];
    if (molmap_tmp) imol = molmap_tmp[imol-idlo];
    else imol--;
    double molmass = masstotal_tmp[imol];

    if (cms[i][0] < lower_x + buffer_size) {
      g_par_all += molmass;
      counting_global += 1;
    }
    else {
      g_par_all += molmass;
      counting_global += 1;
    }

  }
  std::cout<<"counting_global: "<<counting_global<<"\n"; */

  // for now only smooth transition -> ADD for step = 1
  double carg;
  if (!step) {
    for (int i = 0; i < nlocal; i++) {
      
      // mass of a particle
      double mass_temp = mass[type[i]];

      // (looping over nlocals) check if in prescribed region, continue if not 
      if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;

      // if here ... molecule in buffer
      if (x[i][0] < lower_x + buffer_size) { // LEFT BUFFER
        // std::cout<<"LEFT BUFFER"<<"\n";
        // std::cout<<"cms[i][0]: "<<cms[i][0]<<" lower_x + buffer_size: "<<lower_x + buffer_size<<"\n"; 
        if (x[i][0] < (lower_x + (1.0 - g_fac) * buffer_size)) { // only mass
          // std::cout<<"LEFT if"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<"(lower_x + (1.0 - g_fac) * buffer_size): "<<(lower_x + (1.0 - g_fac) * buffer_size)<<"\n"; 
          
          g_par_all += mass_temp;
        }
        else { // sigmoidal mass          
          carg = 1.0 / g_fac * MY_PI * (x[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI; // OK (also for -lower_x)
          g_par_all += 0.5 * (1.0 + cos(carg)) * mass_temp;

          // std::cout<<"LEFT else"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<" 1.0 / g_fac * MY_PI * (cms[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI: "<<1.0 / g_fac * MY_PI * (cms[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI<<"\n";
          // std::cout<<"0.5 * (1.0 + cos(carg)) * molmass: "<<0.5 * (1.0 + cos(carg)) * molmass<<"\n"; 
        }
      }
      // test LEFT buffer ... OK ! ! ! 
      if (x[i][0] > upper_x - buffer_size) { // RIGHT BUFFER
        // std::cout<<"RIGHT BUFFER"<<"\n";
        // std::cout<<"cms[i][0]: "<<cms[i][0]<<" upper_x - buffer_size: "<<upper_x - buffer_size<<"\n";
        if (x[i][0] > (upper_x - (1.0 - g_fac) * buffer_size)) {
          // std::cout<<"RIGHT if"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<"(upper_x - (1.0 - g_fac) * buffer_size): "<<(upper_x - (1.0 - g_fac) * buffer_size)<<"\n"; 

          g_par_all += mass_temp;
        }
        else {
          carg = 1.0 / g_fac * MY_PI * (x[i][0] - upper_x + buffer_size) / (buffer_size) - MY_PI; // OK (also for -lower_x)
          g_par_all += 0.5 * (1.0 + cos(carg)) * mass_temp;

          // std::cout<<"RIGHT else"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<" 1.0 / g_fac * MY_PI * (cms[i][0] - upper_x + buffer_size) / (buffer_size) - MY_PI: "<<1.0 / g_fac * MY_PI * (cms[i][0] - upper_x + buffer_size) / (buffer_size) - MY_PI<<"\n";
          // std::cout<<"0.5 * (1.0 + cos(carg)) * molmass: "<<0.5 * (1.0 + cos(carg)) * molmass<<"\n"; s
        }
      }
      // OTHER IS NOT OF MY INTEREST ... ROI
    }
  } 
  else {
    error->all(FLERR, "For now implemented only for !step");
  } 

  // sum and distribute
  MPI_Allreduce(&g_par_all, &g_par_all_tmp, 1, MPI_DOUBLE, MPI_SUM, world);

  return g_par_all_tmp; 

}

/* ---------------------------------------------------------------------- */

double FixATMOBMD::g_par_local(double masss, Region* region, double xv, int step)
{
  /* double xh = domain->boxhi[0];
  double xl = domain->boxlo[0];
  double y;

  std::cout<<"mass: "<<mass<<"\n";
  if(!step) {
    if(xv < xl+(1.0-g_fac)*buffer_size) return mass; // 1.0; //tp 23.10.2023
    else if(xv < xl+buffer_size) { //tp 23.10.2023
	  y = g_fac_inv*MY_PI*(xv-buffer_size-xl)/(0.0-buffer_size)-MY_PI;//tp 23.10.2023
	  return mass * 0.5 * (1.0+cos(y));
	} else if(xv > xh-(1.0-g_fac)*buffer_size) return mass; // 1.0;
	else if(xv > xh-buffer_size) {
	  y = g_fac_inv*MY_PI*(xv-(xh-buffer_size))/(buffer_size)-MY_PI;
	  return mass * 0.5 * (1.0+cos(y));  
    } else return 0.0;
  } else return 1.0; */
  // using molecule's center-of-mass
  double upper_x = domain->boxhi[0]; // upper
  double lower_x = domain->boxlo[0]; // lower

  /* // testing step function
  if (cms_x < lower_x + buffer_size) {
    return molmass;
  }
  else {
    return molmass;
  } */
  double carg = 0.0; 
  if (!step) {
    if (xv < lower_x + buffer_size) { // LEFT BUFFER
      if (xv < (lower_x + (1.0 - g_fac) * buffer_size)) {
        return masss;
      }
      else { // sigmodial mass
        carg = 1.0 / g_fac * MY_PI * (xv - buffer_size - lower_x) / (-buffer_size) - MY_PI;
        return 0.5 * (1.0 + cos(carg)) * masss;
      }
    }
    if (xv > upper_x - buffer_size) { // RIGHT BUFFER
      if (xv > (upper_x - (1.0 - g_fac) * buffer_size)) {
        return masss;
      }
      else { // sigmodial mass
        carg = 1.0 / g_fac * MY_PI * (xv - upper_x + buffer_size) / (buffer_size) - MY_PI;
        return 0.5 * (1.0 + cos(carg)) * masss;
      }
    }
    // OTHER IS NOT OF MY INTEREST ... ROI
  }
  else {
    error->all(FLERR, "For now implemented only for !step");
  } 
  
}

/* ---------------------------------------------------------------------- */

double FixATMOBMD::g_perp_global(Region* iregion_var, int step)
{
  // return 1.0*group->count(igroup,iregion_var);

  /* std::cout<<"FixATMOBMD::g_perp_global()"<<"\n";
  std::cout<<"FixATMOBMD::g_perp_global(): step: "<<step<<"\n"; */
  if (iregion_var) iregion_var->prematch();

  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  double *mass = atom->mass;
  int  *type = atom->type;
  double **x = atom->x;
  int imol; 
  int *tag = atom->tag;
  double **f = atom->f;

  double upper_x = domain->boxhi[0]; 
  double lower_x = domain->boxlo[0];

  double g_perp_all = 0.0;
  double g_perp_all_tmp = 0.0;

  // int count = 0;

  if (step) {
    // std::cout<<"step"<<"\n";

    for (int i = 0; i < nlocal; i++) {
      
      // mass of a particle
      double mass_temp = mass[type[i]];

      if (iregion_var && !iregion_var->match(x[i][0], x[i][1], x[i][2])) continue;
      g_perp_all += mass_temp;
      // count += 1;
    }
  }
  else {
    error->all(FLERR, "For now implemented only for step");
  }

  // sum and distribute
  MPI_Allreduce(&g_perp_all, &g_perp_all_tmp, 1, MPI_DOUBLE, MPI_SUM, world);

  /* std::cout<<"count: "<<count<<"\n";
  std::cout<<"g_perp_all_tmp: "<<g_perp_all_tmp<<"\n"; */
  return g_perp_all_tmp;

}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::reg_force(int vflag, Region* region, double *momentumForce, int step)
{ 
  // std::cout<<"FixATMOBMD::reg_force()"<<"\n";
  if (region) region->prematch(); 
  v_init(vflag);
  
  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  int *mask = atom->mask;
  imageint *image = atom->image;
  tagint *tag = atom->tag;
  double v[6];
  int nlocal = atom->nlocal;
  force_flag = 0;
  double mass_tmp;
  int *type = atom->type;
  double gtmp = g_par_global(region, step);
  double gloctmp;
  double unwrap[3];
  double xval, yval, zval;

  // create test_force to check obmd
  // double test_force[3];
  // test_force[0] = test_force[1] = test_force[2] = 0.0;
  // double sestevek = 0.0;

  // double test_force_all[3];
  // test_force_all[0] = test_force_all[2] = test_force_all[3] = 0.0;

  // std::cout<<"nlocal in reg_force: "<<nlocal<<"\n";
  // std::cout<<"momentumForce[0]: "<<momentumForce[0]<<" momentumForce[1]: "<<momentumForce[1]<<" momentumForce[2]: "<<momentumForce[2]<<"\n"; 
  
  // int allnu = 0;

  for (int i = 0; i < nlocal; i++){
      // if (mask[i] & groupbit) { // removing because force should be distributed among all particles (i.e., ions and water) within buffer (ions are NOT deleted)
      
      /* if (x[i][0] < 158.0 && x[i][0] > 0.0) {
        allnu += 1;
      } */

      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
      domain->unmap(x[i],image[i],unwrap);
      // mass_tmp = mass[type[i]];
      // gloctmp = g_par_local(region,x[i][0],step)/gtmp;

      if(mode == MOLECULE) {
        // std::cout<<"molecule"<<"\n";
        mass_tmp = mtot;
        gloctmp = g_par_local(mass[type[i]], region, x[i][0], step); // gloctmp = g_par_local(region,x[i][0],step)/gtmp;  // tp 24.7.2023  : TODO: cms of molecule not x[i][0] ... does not matter if g_par = step function 

        /* if (type[i] == 49) {
          std::cout<<"HERE: mass[type[i]]: "<<mass[type[i]]<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        }
        if (atom->tag[i] == 71299) {
          std::cout<<"HERE: mass[type[i]]: "<<mass[type[i]]<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        }
        if (atom->tag[i] == 71300) {
          std::cout<<"HERE: mass[type[i]]: "<<mass[type[i]]<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        }
        if (atom->tag[i] == 71301) {
          std::cout<<"HERE: mass[type[i]]: "<<mass[type[i]]<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        } */
      }
      
      xval = momentumForce[0] * gloctmp / gtmp; // mass[type[i]] * momentumForce[0] * gloctmp / mass_tmp;
      yval = momentumForce[1] * gloctmp / gtmp; // mass[type[i]] * momentumForce[1] * gloctmp / mass_tmp;
      zval = momentumForce[2] * gloctmp / gtmp; // mass[type[i]] * momentumForce[2] * gloctmp / mass_tmp;

      // std::cout<<"gloctmp: "<<gloctmp<<" gtmp: "<<gtmp<<" mass[type[i]]: "<<mass[type[i]]<<" mass_tmp: "<<mass_tmp<<"\n";
            
      /* if (atom->tag[i] == 71291) {
        std::cout<<"obmd: first: 1"<<"\n";
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      } 

      if (atom->type[i] == 49) {
        std::cout<<"obmd: first: 2"<<"\n";
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"xval: "<<xval<<" yval: "<<yval<<" zval: "<<zval<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      } */

      f[i][0] += xval; 
      f[i][1] += yval; 
      f[i][2] += zval; 

      // std::cout<<"xval: "<<xval<<" yval: "<<yval<<" zval: "<<zval<<"\n";
      
      /* if (atom->tag[i] == 71291) {
        std::cout<<"obmd: first: 2"<<"\n";
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"xval: "<<xval<<" yval: "<<yval<<" zval: "<<zval<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      } 

      if (atom->type[i] == 49) {
        std::cout<<"obmd: first: 2"<<"\n";
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"xval: "<<xval<<" yval: "<<yval<<" zval: "<<zval<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      } */


      /* if (atom->tag[i] == 226) {
        std::cout<<"obmd: second"<<"\n";
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
      } */

      // test_force[0] += momentumForce[0] * gloctmp / gtmp; // mass[type[i]] * momentumForce[0] * gloctmp / mass_tmp;
      // test_force[1] += momentumForce[1] * gloctmp / gtmp; // mass[type[i]] * momentumForce[1] * gloctmp / mass_tmp;
      // test_force[2] += momentumForce[2] * gloctmp / gtmp; // mass[type[i]] * momentumForce[2] * gloctmp / mass_tmp; 

      // sestevek += gloctmp;
      
      if (evflag) {
        v[0] = xval * unwrap[0];
        v[1] = yval * unwrap[1];
        v[2] = zval * unwrap[2];
        v[3] = xval * unwrap[1];
        v[4] = xval * unwrap[2];
        v[5] = yval * unwrap[2];
        v_tally(i,v);
      }
    }
 
  // output test_force to see if sum is equal to the prescribed external force 
  // testing only equilibrium case (deleteing)
  /* std::cout<<"FixATMOBMD::reg_force(): test_force[0]: "<<test_force[0]<<" momentumForce[0]: "<<momentumForce[0]<<"\n";
  std::cout<<"FixATMOBMD::reg_force(): test_force[1]: "<<test_force[1]<<" momentumForce[1]: "<<momentumForce[1]<<"\n";
  std::cout<<"FixATMOBMD::reg_force(): test_force[2]: "<<test_force[2]<<" momentumForce[2]: "<<momentumForce[2]<<"\n";
  std::cout<<"sestevek: "<<sestevek<<" gtmp: "<<gtmp<<"\n";
  std::cout<<"allnu: "<<allnu<<"\n";

  MPI_Allreduce(test_force, test_force_all, 3, MPI_DOUBLE, MPI_SUM, world); 
  std::cout<<"FixATMOBMD::reg_force(): test_force_all[0]: "<<test_force_all[0]<<" momentumForce[0]: "<<momentumForce[0]<<"\n";
  std::cout<<"FixATMOBMD::reg_force(): test_force_all[1]: "<<test_force_all[1]<<" momentumForce[1]: "<<momentumForce[1]<<"\n";
  std::cout<<"FixATMOBMD::reg_force(): test_force_all[2]: "<<test_force_all[2]<<" momentumForce[2]: "<<momentumForce[2]<<"\n"; */
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::reg_force_perp(int vflag, Region* region, double *shearForce, int step)
{ 
  std::cout<<"FixATOBMD::reg_force_perp()"<<"\n";
  if (region) region->prematch();
  v_init(vflag);
  
  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  int nlocal = atom->nlocal;
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;
  int *type = atom->type;
  double gtmp = g_perp_global(region, step);
  double mass_tmp,gloctmp;
  double unwrap[3];
  double xval, yval, zval;

  // create test_force to check obmd
  double test_perp_force[3];
  test_perp_force[0] = test_perp_force[1] = test_perp_force[2] = 0.0;
  double sestevek_perp = 0.0;
  
  double test_perp_force_all[3];
  test_perp_force_all[0] = test_perp_force_all[1] = test_perp_force_all[2] = 0.0;

  /* std::cout<<"nlocal in reg_perp_force: "<<nlocal<<"\n";
  std::cout<<"shearForce_left[0]: "<<shearForce_left[0]<<" shearForce_left[1]: "<<shearForce_left[1]<<" shearForce_left[2]: "<<shearForce_left[2]<<"\n"; */

  for (int i = 0; i < nlocal; i++){
    // std::cout<<"nlocal: "<<nlocal<<"\n";
    // if (mask[i] & groupbit) {
    if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue; // is it in region ?
    domain->unmap(x[i],image[i],unwrap);
      
    /* if (mode == MOLECULE) {
      // std::cout<<"mode==MOLECULE"<<"\n";
      mass_tmp = mtot; // mass[type[i]];
      gloctmp = mass[type[i]]; //g_perp_local(mass[type[i]], region, x[i][0], step); 
    } */

    gloctmp = mass[type[i]]; // heaviside (only mass[type[i]])

    xval = shearForce[0] * gloctmp / gtmp; // mass[type[i]] * shearForce[0] * gloctmp / mass_tmp;
    yval = shearForce[1] * gloctmp / gtmp; // mass[type[i]] * shearForce[1] * gloctmp / mass_tmp;
    zval = shearForce[2] * gloctmp / gtmp; // mass[type[i]] * shearForce[2] * gloctmp / mass_tmp;

    f[i][0] += xval; 
    f[i][1] += yval; 
    f[i][2] += zval; 

    /* std::cout<<"gloctmp: "<<gloctmp<<" gtmp: "<<gtmp<<" mass[type[i]]: "<<mass[type[i]]<<" mass_tmp: "<<mass_tmp<<"\n";
    std::cout<<"xval: "<<xval<<" yval: "<<yval<<" zval: "<<zval<<"\n"; */

    // sestevek_perp += gloctmp;

    /* test_perp_force[0] += shearForce[0] * gloctmp / gtmp;
    test_perp_force[1] += shearForce[1] * gloctmp / gtmp;
    test_perp_force[2] += shearForce[2] * gloctmp / gtmp; */
      
    if (evflag) {
      v[0] = xval * unwrap[0];
      v[1] = yval * unwrap[1];
      v[2] = zval * unwrap[2];
      v[3] = xval * unwrap[1];
      v[4] = xval * unwrap[2];
      v[5] = yval * unwrap[2];
      v_tally(i,v);
    }
  }

  /* std::cout<<"FixATMOBMD::reg_force_perp(): test_perp_force[0]: "<<test_perp_force[0]<<" shearForce[0]: "<<shearForce[0]<<"\n";
  std::cout<<"FixATMOBMD::reg_force_perp(): test_perp_force[1]: "<<test_perp_force[1]<<" shearForce[1]: "<<shearForce[1]<<"\n";
  std::cout<<"FixATMOBMD::reg_force_perp(): test_perp_force[2]: "<<test_perp_force[2]<<" shearForce[2]: "<<shearForce[2]<<"\n";
  std::cout<<"sestevek_perp: "<<sestevek_perp<<" gtmp: "<<gtmp<<"\n";

  MPI_Allreduce(test_perp_force, test_perp_force_all, 3, MPI_DOUBLE, MPI_SUM, world); 
  std::cout<<"FixATMOBMD::reg_perp_force(): test_perp_force_all[0]: "<<test_perp_force_all[0]<<" shearForce[0]: "<<shearForce[0]<<"\n";
  std::cout<<"FixATMOBMD::reg_perp_force(): test_perp_force_all[1]: "<<test_perp_force_all[1]<<" shearForce[1]: "<<shearForce[1]<<"\n";
  std::cout<<"FixATMOBMD::reg_perp_force(): test_perp_force_all[2]: "<<test_perp_force_all[2]<<" shearForce[2]: "<<shearForce[2]<<"\n"; */
}

/* ----------------------------------------------------------------------- */

double FixATMOBMD::usher(Region *iregion_var, double **coords, double etarget, int natom, int imol, int &iter)
{ 
  int i, m, type_temp;
  double entmp, torqabs, dtheta, ds, fabs;
  double quat[4], rotmat[3][3], fusher_tmp[3], fusher_all[3], xcom[3];
  double torq[3], torq_all[3], torq_tmp[3], coords_tmp[3];

  double entmp_all;
  MathExtra::zero3(fusher_tmp);
  
  double *mass = atom->mass; // PPapez COMMENT: to distribute fusher 

  i = 0;
  while(i < nattempt) {
    /* std::cout<<"FixATOBMD::usher()"<<"\n";
    std::cout<<"nattempt: "<<i<<"\n";
    std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
    std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
    std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; */ 
    
    entmp = 0.0;
    MathExtra::zero3(fusher);
    MathExtra::zero3(torq);
    
    mol_center_of_mass(natom, imol, coords, xcom); // PPapez COMMENT: find mol's cms \& store it in xcom
    // std::cout<<"xcom[0]: "<<xcom[0]<<" xcom[1]: "<<xcom[1]<<" xcom[2]: "<<xcom[2]<<"\n";

    for (m = 0; m < natom; m++) { // PPapez COMMENT: looping over all atoms of my molecule 
      type_temp = ntype+onemols[imol]->type[m];

      entmp += energy_atomistic_obmd(iregion_var, onemols[imol]->q[m], type_temp, coords[m], fusher); 
      
      // torque
      MathExtra::copy3(fusher, fusher_tmp);
      calc_torque(natom, coords, xcom, fusher_tmp, torq_tmp);
      MathExtra::add3(torq, torq_tmp, torq);
    }

    MPI_Allreduce(fusher, fusher_all, 3, MPI_DOUBLE, MPI_SUM, world); 
    MPI_Allreduce(&entmp, &entmp_all, 1, MPI_DOUBLE, MPI_SUM, world); 
    MPI_Allreduce(&torq, &torq_all, 3, MPI_DOUBLE, MPI_SUM, world);
    // std::cout<<"FixATOBMD::usher(): entmp: "<<entmp<<"\n";

    // PPapez COMMENT: USHER routine adopted for mode == MOLECULE (carefully with center-of-mass property \& pbc)
    if(entmp_all < etarget+EPSILON) 
    {
      break;
    }
    else if(entmp_all > uovlp) {
      fabs = sqrt(fusher_all[0] * fusher_all[0] + fusher_all[1] * fusher_all[1] + fusher_all[2] * fusher_all[2]);
      if(fabs < EPSILON) continue;
      ds = dsovlp - pow(4 * eps / entmp_all, 1.0/12.0); 

      for(m = 0; m < natom; m++) { // PPapez COMMENT: loop over atoms of my molecule

        /* std::cout<<"HERE 1"<<"\n";
        std::cout<<"m: "<<m<<"\n";
        std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n"; */
        coords[m][0] += fusher_all[0] / fabs * ds; 
        coords[m][1] += fusher_all[1] / fabs * ds; 
        coords[m][2] += fusher_all[2] / fabs * ds; 
        /* std::cout<<"HERE 1"<<"\n";
        std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
        std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
        std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; */
        
        // PPapez COMMENT: rather checking region and correction position
        // imageflags[m] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        // domain->remap(coords[m],imageflags[m]);
      }
      
      // PPapez COMMENT: check region ... returns 0 / 1 if yes or not
      int check = check_mol_region(iregion_var, coords, natom); 
      if (check == 1) {
        // std::cout<<"in first check_mol_region"<<"\n";
        /* for(m = 0; m < natom; m++) {     
          std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
        } */
        break; // exit 
        // PPapez COMMENT: try with new position
        // try_new_position(iregion_var, coords, imol);
      }
    } 
    else {
      fabs = sqrt(fusher_all[0] * fusher_all[0] + fusher_all[1] * fusher_all[1] + fusher_all[2] * fusher_all[2]);
      if(fabs < EPSILON) continue;
      ds = std::min((entmp_all - etarget) / fabs, ds0);

      torqabs = sqrt(torq_all[0] * torq_all[0] + torq_all[1] * torq_all[1] + torq_all[2] * torq_all[2]);
      dtheta == std::min((entmp_all - etarget) / torqabs, dtheta0);
      MathExtra::norm3(torq_all);
      MathExtra::axisangle_to_quat(torq_all, dtheta, quat);
      MathExtra::quat_to_mat(quat, rotmat);

      for(m = 0; m < natom; m++) {
        /* std::cout<<"HERE 2"<<"\n";
        std::cout<<"m: "<<m<<"\n";
        std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n"; */
        coords[m][0] += fusher_all[0] / fabs * ds; 
        coords[m][1] += fusher_all[1] / fabs * ds; 
        coords[m][2] += fusher_all[2] / fabs * ds; 
        /* std::cout<<"HERE 2"<<"\n";
        std::cout<<"fusher_all[0]: "<<fusher_all[0]<<" fusher_all[1]: "<<fusher_all[1]<<" fusher_all[2]: "<<fusher_all[2]<<"\n";
        std::cout<<"fabs: "<<fabs<<"\n";
        std::cout<<"ds: "<<ds<<"\n";
        std::cout<<"fusher_all[0]/fabs: "<<fusher_all[0]/fabs<<"\n";
        std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";  */

        // imageflags[m] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        // domain->remap(coords[m],imageflags[m]);

        MathExtra::matvec(rotmat, coords[m], coords_tmp);
        MathExtra::copy3(coords_tmp, coords[m]);
        // imageflags[m] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        // domain->remap(coords[m], imageflags[m]);
        // std::cout<<"HERE 2"<<"\n";
        // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
      }

      // PPapez COMMENT: check region ... returns 0 / 1 if yes or not
      int check = check_mol_region(iregion_var, coords, natom);
      if (check == 1) {
        // std::cout<<"in second check_mol_region"<<"\n";
        /* for(m = 0; m < natom; m++) {     
          std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
        } */
        break; // exit
        // PPapez COMMENT: try with new position
        // try_new_position(iregion_var, coords, imol);
      }  
    }
    i+=1;
  }
  iter = i;
  entmp = entmp_all;
  
  return entmp;
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::try_new_position(Region *iregion_var, double **coords, int imol)
{ 
  /*
  PPapez COMMENT
  try to find new position if the previous one is outside "iregion_var"
  adopted for mode == MOLECULE
  checking positions
  */ 

  int i, m, n, natom;
  double coord[3];
  double r[3], rotmat[3][3], quat[4];
  int insert_value = 0; // PPapez COMMENT: control

  std::cout<<"FixAdResSObmd::try_new_position"<<"\n";
  std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
  std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
  std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 

  keep_trying:
  do {
    coord[0] = xlo + random->uniform() * (xhi-xlo);
    coord[1] = ylo + random->uniform() * (yhi-ylo);
    coord[2] = zlo + random->uniform() * (zhi-zlo);
    std::cout<<"coord[0]: "<<coord[0]<<" coord[1]: "<<coord[1]<<" coord[2]: "<<coord[2]<<"\n";
  } 
  while (iregion_var->match(coord[0],coord[1],coord[2]) == 0);

  double rng = random->uniform();
  natom = onemols[imol]->natoms;
  if (orientflag) {
    r[0] = rx;
    r[1] = ry;
    r[2] = rz;
  } 
  else {
    r[0] = random->uniform() - 0.5;
    r[1] = random->uniform() - 0.5;
    r[2] = random->uniform() - 0.5;
  }

  double theta = random->uniform() * MY_2PI;
  MathExtra::norm3(r);
  MathExtra::axisangle_to_quat(r,theta,quat);
  MathExtra::quat_to_mat(quat,rotmat);
  for (i = 0; i < natom; i++) {
    MathExtra::matvec(rotmat,onemols[imol]->dx[i],coords[i]);
    std::cout<<"coords[i][0]: "<<coords[i][0]<<"\n";
    coords[i][0] += coord[0];
    coords[i][1] += coord[1];
    coords[i][2] += coord[2];
  }

  std::cout<<"FixATMOBMD::try_new_position"<<"\n";
  std::cout<<"before check_mol_region"<<"\n";
  std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
  std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
  std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 
      
  int check = check_mol_region(iregion_var,coords,3);
  if (check == 1) {
    insert_value = 1;
    std::cout<<"in check_mol_region"<<"\n";
    for(m = 0; m < 3; m++) {     
      std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
    }
  }
      
  if (insert_value == 1) { // PPapez COMMENT: should "goto" be limited ? ? ? usually takes two or three iterations ... I will leave it for now
    std::cout<<"if (insert_value == 1)"<<"\n";

    insert_value = 0; // set to 0 \& repeat 
    goto keep_trying;      
  }

  std::cout<<"in try_new_position after correction"<<"\n";
  std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
  std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
  std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 
}

/* ---------------------------------------------------------------------- */

int FixATMOBMD::check_proc(double **coords, double *newcoord, double *sublo, double*subhi, int dimension)
{
  int flag = 0;
  if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
      newcoord[1] >= sublo[1] && newcoord[1] < subhi[1] &&
      newcoord[2] >= sublo[2] && newcoord[2] < subhi[2]) flag = 1;
  else if (dimension == 3 && newcoord[2] >= domain->boxhi[2]) {     /// for 3D
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[2] == comm->procgrid[2]-1 &&
          newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
          newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
    } else {
      if (comm->mysplit[2][1] == 1.0 &&
          newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
          newcoord[1] >= sublo[1] && newcoord[1] < subhi[1]) flag = 1;
    }
  } else if (dimension == 2 && newcoord[1] >= domain->boxhi[1]) {    /// for 2D
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[1] == comm->procgrid[1]-1 &&
          newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
    } else {
      if (comm->mysplit[1][1] == 1.0 &&
          newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
    }
  }
  return flag;
}

/* ---------------------------------------------------------------------- */

int FixATMOBMD::check_mol_proc(double **coords, double *newcoord, double *lamda, double *sublo, double*subhi, int dimension, int natom, int &me)
{
  int m;
  int sum_all_tmp = 0;
  int sum_all_tmpall,flag;
  
  for (m = 0; m < natom; m++) {
    if (domain->triclinic) {
      domain->x2lamda(coords[m],lamda);
      newcoord = lamda;
    } else newcoord = coords[m];
    flag = check_proc(coords,newcoord,sublo,subhi,dimension);
    sum_all_tmp += flag;
  }
  MPI_Allreduce(&sum_all_tmp,&sum_all_tmpall,1,MPI_INT,MPI_SUM,world); 
  if(sum_all_tmp == natom) me = comm->me;  ///ce je cela molekula na enem, potem posebej oznaci me = comm->me .. zaenkrat brez posebne funkcije
  
  /*if(sum_all_tmpall == natom) {
	for (m = 0; m < natom; m++) {
		newcoord = coords[m];
		flag = check_proc(coords,newcoord,sublo,subhi,dimension);
		if(flag == 1) std::cout << "Atom m = " << m << " on proc " << comm->me 
		<< ". coords[m][0] = " << coords[m][0] << ". coords[m][1] = " << coords[m][1] << ". coords[m][2] = " << coords[m][2] << std::endl;
	}  
  }*/

  return sum_all_tmpall;
  //return sum_all_tmp; //samo testiram tp 9.12.2023 //ce dam to potem je mpiallreduce error ?? zakaj?
}

/* ---------------------------------------------------------------------- */

int FixATMOBMD::check_mol_region(Region* region, double **coords, int natom)
{
  int m,flag,flagall;
  flag = 0;
  for (m = 0; m < natom; m++) {
    if (region && !region->match(coords[m][0],coords[m][1],coords[m][2])) flag = 1;
  }
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
  if(comm->me==0 && flagall==1) std::cout << "USHER moved the particle/molecule too much." << std::endl;
  
  return flagall;
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::mol_center_of_mass(int natom, int imol, double **coords, double *xcom)
{
  double *mass = atom->mass;
  int m;
  double massmol = 0.0;
  MathExtra::zero3(xcom);
  
  for(m = 0; m < natom; m++) {
    xcom[0] += mass[onemols[imol]->type[m]] * coords[m][0];
    xcom[1] += mass[onemols[imol]->type[m]] * coords[m][1];
    xcom[2] += mass[onemols[imol]->type[m]] * coords[m][2];

    massmol += mass[onemols[imol]->type[m]];

  }
  // std::cout<<"massmol: "<<massmol<<"\n";

  xcom[0] = xcom[0] / massmol;
  xcom[1] = xcom[1] / massmol;
  xcom[2] = xcom[2] / massmol;
  // std::cout<<"xcm[0]: "<<xcm[0]<<" xcm[1]: "<<xcm[1]<<" xcm[2]: "<<xcm[2]<<"\n";
}

/* ---------------------------------------------------------------------- */

void FixATMOBMD::calc_torque(int natom, double **coords, double *xcom, double *force_tmp, double *torq)
{
  int m;
  double xrel[3];
  
  for(m = 0; m < natom; m++) {
    xrel[0] = coords[m][0] - xcom[0];
    xrel[1] = coords[m][1] - xcom[1];
    xrel[2] = coords[m][2] - xcom[2];
  }   
  
  MathExtra::cross3(xrel,force_tmp,torq); // HERE ! ! ! 
}

/* ---------------------------------------------------------------------- */

double FixATMOBMD::energy_atomistic_obmd(Region *iregion_var, double qi, int itype, double *coord, double *fusher)
{  
  double delx, dely, delz, rsq;
  int jtype;
 
  double **x = atom->x;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  auto pair = dynamic_cast<PairLJCutRF*>(force->pair_match("lj/cut/rf", 1)); // potential acting between atomistic particles

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;
  double total_energy = 0.0;

  // std::cout<<"FixATOBMD::energy_atomistic_obmd(): itype: "<<itype<<"\n";

  for (int j = 0; j < nlocal; j++) { //nall is probably not needed, could just be nlocal ...  PPapez COMMENT: USING THIS
    // PPapez COMMENT: match region
    // if (iregion_var->match(x[j][0], x[j][1], x[j][2])) {
      delx = coord[0] - x[j][0];
      dely = coord[1] - x[j][1];
      delz = coord[2] - x[j][2];

      domain->minimum_image(delx, dely, delz);
      rsq = delx * delx + dely * dely + delz * delz;

      jtype = type[j];

      if (rsq < cutsq[itype][jtype]){
        total_energy += pair->single_atomistic_obmd(qi, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

        fusher[0] += fpair * delx;
        fusher[1] += fpair * dely;
        fusher[2] += fpair * delz;

         /* if (fpair > 10000) {
          std::cout<<"FixATOBMD::energy_atomistic_obmd"<<"\n";
          std::cout<<"fpair: "<<fpair<<"\n";
        } */
      }
    // }
  }

  // std::cout<<"total_energy: "<<total_energy<<"\n";
  return total_energy;
} 

/* ---------------------------------------------------------------------- */

/* double FixATobmd::energy(int i, int itype, double *coord, double *fusher)
{  
  // std::cout<<"in energy function"<<"\n";
  // std::cout<<"i: "<<i<<" itype: "<<itype<<" coord[0]: "<<coord[0]<<" coord[1]: "<<coord[1]<<" coord[2]: "<<coord[2]<<"\n";

  double delx,dely,delz,rsq;
  int jtype;
  double total_energy = 0.0;
  double **x = atom->x;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  pair = force->pair;
  cutsq = force->pair->cutsq;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;

  //PairLJCut *pstyle = dynamic_cast<PairLJCut *>(force->pair_match("lj/cut",1)); // not needed
  /// does pair_hybrid_overlay sum all energies and forces with single() ??
  /// answer: yes, if they contain single() functions, otherwise it throws an error
  /// should this error be kept or replaced by a warning?
  /// answer: probably best to keep it as an error or even better to add an "ignore flag" in the fix obmd command if the user knows what he is doing
  
  //std::cout << " inside usher on proc " << comm->me << ": atom->nghost = " << atom->nghost << " atom->nlocal = " << atom->nlocal << "." << std::endl;
  // my understanding so far:
  // let say you want lj/cut + dpd thermostat
  // in this case you can not use pair_style hybrid/overlay lj/cut lj_cutoff dpd/ext/tstat T rc seed
  // but you can use hybrid/overlay with two force fields that have a single() function
  // solution: use the logic in pair_hybrid which checks if pair has a single function (there it throws an error), here we could just add a warning once
  double temp_energy;
  for (int j = 0; j < nall; j++) { //nall is probably not needed, could just be nlocal
    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];
    domain->minimum_image(delx,dely,delz); //added tp 20.10.2023
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    if (rsq < cutsq[itype][jtype]) {
      // std::cout<<"atom->tag[j]: "<<atom->tag[j]<<"\n"; */

      /* if (atom->tag[j] == 3763) {
        std::cout<<"atom->tag[j]: "<<atom->tag[j]<<"\n";
        std::cout<<"coord[0]: "<<coord[0]<<" x[j][0]: "<<x[j][0]<<"\n";
        std::cout<<"coord[1]: "<<coord[1]<<" x[j][1]: "<<x[j][1]<<"\n";
        std::cout<<"coord[2]: "<<coord[2]<<" x[j][2]: "<<x[j][2]<<"\n"; 
      }
      if (atom->tag[j] == 3764) {
        std::cout<<"atom->tag[j]: "<<atom->tag[j]<<"\n";
        std::cout<<"coord[0]: "<<coord[0]<<" x[j][0]: "<<x[j][0]<<"\n";
        std::cout<<"coord[1]: "<<coord[1]<<" x[j][1]: "<<x[j][1]<<"\n";
        std::cout<<"coord[2]: "<<coord[2]<<" x[j][2]: "<<x[j][2]<<"\n"; 
      }
      if (atom->tag[j] == 3765) {
        std::cout<<"atom->tag[j]: "<<atom->tag[j]<<"\n";
        std::cout<<"coord[0]: "<<coord[0]<<" x[j][0]: "<<x[j][0]<<"\n";
        std::cout<<"coord[1]: "<<coord[1]<<" x[j][1]: "<<x[j][1]<<"\n";
        std::cout<<"coord[2]: "<<coord[2]<<" x[j][2]: "<<x[j][2]<<"\n"; 
      } */

      // std::cout<<"before single"<<"\n";
      // temp_energy = pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair); // fpair returns energy ? ? ? 

      // avoiding table/adress/obmd/cg
      /* auto pair = dynamic_cast<PairLJCutRFAdResSAT*>(force->pair_match("lj/cut/rf/adress/at",1));
      total_energy += pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

      // total_energy += pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);
      // std::cout<<"pair->single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fpair): "<<temp_energy<<" rsq: "<<rsq<<"\n";
      // std::cout<<"atom->tag[j]: "<<atom->tag[j]<<" jtype: "<<jtype<<" itype: "<<itype<<"\n";
      // std::cout<<"after single"<<"\n";
      // std::cout<<"i: "<<i<<" itype: "<<itype<<" fpair: "<<fpair<<"\n";
      // std::cout<<"j: "<<j<<" jtype: "<<jtype<<" fpair: "<<fpair<<"\n";

      fusher[0] += fpair * delx;
      fusher[1] += fpair * dely;
      fusher[2] += fpair * delz;
	  }
  }

  // std::cout<<"total_energy: "<<total_energy<<"\n";
  return total_energy;
} */

/* ---------------------------------------------------------------------- */

/* double FixATMOBMD::compute_scalar()
{
  return 1.0;
} */

/* ---------------------------------------------------------------------- */

/* double FixATMOBMD::compute_vector(int n)
{
  return 1.0;
} */

/* ---------------------------------------------------------------------- */

void FixATMOBMD::find_maxid()
{
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  tagint max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,tag[i]);
  MPI_Allreduce(&max,&maxtag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

  if (mode == MOLECULE && molecule) {
    max = 0;
    for (int i = 0; i < nlocal; i++) max = MAX(max,molecule[i]);
    MPI_Allreduce(&max,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
  }
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

/* void FixATMOBMD::min_post_force(int vflag)
{
  post_force(vflag);
} */

/* ---------------------------------------------------------------------- */

void FixATMOBMD::options(int narg, char **arg)
{
  // default settings
  // regions
  iregion = nullptr; // region1
  idregion = nullptr;
  iregion2 = nullptr; // region2
  idregion2 = nullptr;
  iregion3 = nullptr; // region3
  idregion3 = nullptr;
  iregion4 = nullptr; // region4
  idregion4 = nullptr;  
  iregion5 = nullptr; // PPapez COMMENT: region5 for insertion in LEFT BUFFER
  idregion5 = nullptr;
  iregion6 = nullptr; // PPapez COMMENT: region5 for insertion in RIGHT BUFFER
  idregion6 = nullptr; 
  /* iregion7 = nullptr;
  idregion7 = nullptr;
  iregion8 = nullptr;
  idregion8 = nullptr; */

  buffer_size = 0.30 * (domain->boxhi[0] - domain->boxlo[0]); 
  // constraints
  shakeflag = 0;
  idshake = nullptr;
  stepflag = 0;
  step_parallel = 0; // normal
  step_perp = 1; // tangential
  g_fac = 0.25;
  g_fac_inv = 1.0 / g_fac;
  // usher params
  alpha = 1.0;
  tau = 1.0;
  // nbuf = group->count(igroup,iregion); // how many of them are in left buffer
  nbuf = nbuf = group->count(igroup)*buffer_size/(domain->boxhi[0]-domain->boxlo[0]);
  etgt = 3.6;
  ds0 = 1.0;
  dtheta0 = 0.35;
  uovlp = 10000.0;
  dsovlp = 3.0;
  eps = 0.15;
  nattempt = 40;
  // other params
  maxattempt = 1;
  // distflag = DIST_UNIFORM; // for init positions ... assuming DIST_UNIFORM
  sigma = 1.0; // distribution
  xmid = ymid = zmid = 0.0; // distribution
  molfrac = nullptr;
  // mode1 = mode2 = mode3 = MOLECULE; // atom if not molecule // only MOLECULE
  idnext = 0;
  mol_len = -1;
  imol = -1; // to check whether mol template is given
  iarg = 0;
  rx = ry = rz = 0.0; // by default each molecule is inserted at random orienttaion ... perhaps not needed (only for molecules)
  
  // std::cout<<"fix atomistic/obmd::options: narg: "<<narg<<"\n";
  // loop
  while (iarg < narg) {
    // std::cout<<"iarg: "<<iarg<<"\n";
    if (strcmp(arg[iarg], "region1") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion = domain->get_region_by_id(arg[iarg+1]); // id = utils::strdup(arg[0]) in region.cpp // in our case leftB
      if (!iregion) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n]; // array & new char
      strcpy(idregion, arg[iarg+1]);
      iarg += 2; // length of two
    }
    else if (strcmp(arg[iarg], "region2") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion2 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion2) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion2 = new char[n];
      strcpy(idregion2, arg[iarg+1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "region3") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion3 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion3) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion3 = new char[n];
      strcpy(idregion3, arg[iarg+1]);
      iarg += 2; 
    }
    else if (strcmp(arg[iarg], "region4") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion4 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion4) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion4 = new char[n];
      strcpy(idregion4, arg[iarg+1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "region5") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion5 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion5) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion5 = new char[n];
      strcpy(idregion5, arg[iarg+1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "region6") == 0) {
      // std::cout<<"region6"<<"\n";
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion6 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion6) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion6 = new char[n];
      strcpy(idregion6, arg[iarg+1]);
      iarg += 2;
    }
    /* else if (strcmp(arg[iarg], "region7") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion7 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion7) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion7 = new char[n];
      strcpy(idregion7, arg[iarg+1]);
      iarg += 2;
    }
    else if (strcmp(arg[iarg+2], "region8") == 0) {
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      iregion8 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion8) error->all(FLERR, "Region ID {} for fix atomistic/obmd does not exist", arg[iarg+1]);
      int n = strlen(arg[iarg+1]) + 1;
      idregion8 = new char[n];
      strcpy(idregion8, arg[iarg+1]);
      iarg += 2;
    } */
    else if (strcmp(arg[iarg], "buffersize") == 0) {
      // std::cout<<"buffersize"<<"\n";
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      buffer_size = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      if (buffer_size <= 0.0) { // check that is not less than zero or equal to zero ... division y zero leads to NaN (distribution of external forces)
        error->all(FLERR, "Illegal fix atomistic/obmd command. Parameter buffersize should be larger than 0.");
      }
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "gfac") == 0) {
      // std::cout<<"gfac: iarg: "<<iarg<<"\n";
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      g_fac = utils::numeric(FLERR, arg[iarg+1], false, lmp);
      g_fac_inv = 1.0 / g_fac;

      // some checks
      if (g_fac < 0.0 || g_fac > 1.0) {
        error->all(FLERR,"Illegal fix atomistic/obmd command. Parameter gfac should be between 0 and 1");
      }
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "stepparallel") == 0) { // step function for parallel forces
      // std::cout<<"stepparallel: iarg: "<<iarg<<"\n";
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      step_parallel = utils::inumeric(FLERR, arg[iarg+1], false, lmp); // 1 / 0 
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "stepperp") == 0) { // step function for perpendicular forces
      // std::cout<<"stepperp: iarg: "<<iarg<<"\n";
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      step_perp = utils::inumeric(FLERR, arg[iarg+1], false, lmp); // 1 / 0
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "maxattempt") == 0) { // max # of trying to insert molecule (not USHER attempt)
      // std::cout<<"maxattempt: iarg: "<<iarg<<"\n";
      if (iarg+2 > narg) error->all(FLERR, "Illegal fix at/obmd command");
      maxattempt = utils::inumeric(FLERR, arg[iarg+1], false, lmp); 
      iarg += 2;
    }
    else if (strcmp(arg[iarg], "usher") == 0) { // usher settings for first group
      // std::cout<<"usher: iarg: "<<iarg<<"\n";
      if (iarg+12 > narg) error->all(FLERR, "Illegal fix at/obmd command");
      usherflag = utils::inumeric(FLERR, arg[iarg+1], false, lmp);

      // options
      alpha = utils::numeric(FLERR, arg[iarg+2], false, lmp);
      tau = utils::numeric(FLERR, arg[iarg+3], false, lmp);
      nbuf = utils::inumeric(FLERR, arg[iarg+4], false, lmp);
      etgt = utils::numeric(FLERR, arg[iarg+5], false, lmp);
      ds0 = utils::numeric(FLERR, arg[iarg+6], false, lmp);
      dtheta0 = utils::numeric(FLERR, arg[iarg+7], false, lmp);
      uovlp = utils::numeric(FLERR, arg[iarg+8], false, lmp);
      dsovlp = utils::numeric(FLERR, arg[iarg+9], false, lmp);
      eps = utils::numeric(FLERR, arg[iarg+10], false, lmp);
      nattempt = utils::inumeric(FLERR, arg[iarg+11], false, lmp);
      // std::cout<<"nattempt: "<<nattempt<<"\n";
      // some checks
      if (alpha < 0.0 || tau < 0.0 || nbuf < 0) {
        error->all(FLERR, "Illegal fix atomistic/obmd settings");
      }
      if (nattempt < 1) {
        error->all(FLERR, "Illegal fix atomistic/obmd settings");
      }

      iarg += 12; // 12 params for usher ... mol info is separated
      // std::cout<<"iarg: "<<iarg<<"\n";
    }
    else if (strcmp(arg[iarg], "mol") == 0) {
      // std::cout<<"mol: iarg: "<<iarg<<"\n";
      if (iarg+3 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      imol = atom->find_molecule(arg[iarg+1]); // imol for instance water (find_molecule(const char *id))
      if (imol == -1) { // molecule in set with template ID does not exist
        error->all(FLERR, "Molecule template ID for fix atomistic/obmd does not exist");
      }
      mol_len = utils::inumeric(FLERR, arg[iarg+2], false, lmp); // computing number
      if (mol_len == -1) { // provide mol1
        error->all(FLERR, "Molecule template ID for fix atomistic/obmd does not exist");
      }
      molflag = 1;
      mode = MOLECULE; // PPapez COMMENT: only MOLECULE 
      onemols = &atom->molecules[imol]; // def ?
      onemols[imol]->compute_mass();
      mtot = onemols[imol]->masstotal; // masstotal += atom->mass[type[i]] in molecule.cpp
      nmol = onemols[0]->nset; // # of molecules in list in atom_vec.h
      delete [] molfrac;
      molfrac = new double[nmol];
      molfrac[0] = 1.0 / nmol;
      for (int i = 1; i < nmol-1; i++) molfrac[i] = molfrac[i-1] + 1.0/nmol;
      molfrac[nmol-1] = 1.0;

      iarg += 3;
    }
    else if (strcmp(arg[iarg], "orient") == 0) {
      if (iarg+4 > narg) error->all(FLERR, "Illegal fix atomistic/obmd command");
      orientflag = 1;
      rx = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ry = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      rz = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (domain->dimension == 2 && (rx != 0.0 || ry != 0.0))
        error->all(FLERR,"Illegal fix atomistic/obmd orient settings");
      if (rx == 0.0 && ry == 0.0 && rz == 0.0)
        error->all(FLERR,"Illegal fix atomistic/obmd orient settings");
      iarg += 4;
    } 
    else if (strcmp(arg[iarg],"shake") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix atomistic/obmd command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idshake;
      idshake = new char[n];
      strcpy(idshake,arg[iarg+1]);
      shakeflag = 1;
      iarg += 2;
   }   
  } // closing while

  // USHER in combination with molecular template
  // only molecules
  // NO ATOM ... use fix obmd
  // NO ADRESS ... use fix adress/obmd
} 

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

/* void FixATMOBMD::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = random->state();
  list[n++] = ninserted;
  list[n++] = nfirst;
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
} */

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

/* void FixATMOBMD::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  ninserted = static_cast<int> (list[n++]);
  nfirst = static_cast<int> (list[n++]);
  next_reneighbor = (bigint) ubuf(list[n++]).i;

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR,"Must not reset timestep when restarting this fix");

  random->reset(seed);
} */

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

/* void *FixATMOBMD::extract(const char *str, int &itype)
{
  if (strcmp(str,"radius") == 0) {
    if (mode == ATOM) {
      if (itype == ntype) oneradius = 0.5;
      else oneradius = 0.0;

    } else {

      // loop over onemols molecules
      // skip a molecule with no atoms as large as itype

      oneradius = 0.0;
      for (int i = 0; i < nmol; i++) {
        if (itype > ntype+onemols[i]->ntypes) continue;
        double *radius = onemols[i]->radius;
        int *type = onemols[i]->type;
        int natoms = onemols[i]->natoms;

        // check radii of atoms in Molecule with matching types
        // default to 0.5, if radii not defined in Molecule
        //   same as atom->avec->create_atom(), invoked in pre_exchange()

        for (int i = 0; i < natoms; i++)
          if (type[i]+ntype == itype) {
            if (radius) oneradius = MAX(oneradius,radius[i]);
            else oneradius = MAX(oneradius,0.5);
          }
      }
    }
    itype = 0;
    return &oneradius;
  }

  return nullptr;
} */