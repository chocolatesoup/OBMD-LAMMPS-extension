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
#include "fix_adress_obmd.h"

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

#include "pair_lj_cut_rf_adress_at.h"
#include "pair_table_adress_obmd_cg.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <algorithm> // PPapez COMMENT: to call find

#include <fenv.h> // PP DEBUG

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{ATOM,MOLECULE};
enum{DIST_UNIFORM,DIST_GAUSSIAN};
enum{NONE,CONSTANT,EQUAL};
enum{EXCHATOM,EXCHMOL}; // exchmode
enum{MOVEATOM,MOVEMOL}; // movemode

#define EPSILON 1.0e-6
#define BIG MAXTAGINT

/* ---------------------------------------------------------------------- */

FixAdResSObmd::FixAdResSObmd(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), idregion(nullptr), idregion2(nullptr), idregion3(nullptr), idregion4(nullptr), 
  idregion5(nullptr), idregion6(nullptr), idrigid(nullptr), idshake(nullptr), onemols(nullptr), molfrac(nullptr), 
  coords(nullptr), imageflags(nullptr), fixrigid(nullptr), fixshake(nullptr), random(nullptr), list(nullptr), mark(nullptr)
{
  if (narg < 13) error->all(FLERR,"Illegal fix obmd command");

  restart_global = 1;
  time_depend = 1;

  xstr = nullptr;
  ystr = nullptr;
  zstr = nullptr;
  qstr = nullptr;
  t0_left_str = nullptr;
  t0_right_str = nullptr;
  id_press = nullptr;
  
  size_vector = 1; 
  vector_flag = 1;
  scalar_flag = 1;
  ntype = utils::inumeric(FLERR,arg[3],false,lmp);
  nfreq = utils::inumeric(FLERR,arg[4],false,lmp); //insert every m timesteps
  seed = utils::inumeric(FLERR,arg[5],false,lmp);
 
  if (strstr(arg[6],"v_") == arg[6]) {
    int n = strlen(&arg[6][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[6][2]);
  } else {
    pxx = utils::numeric(FLERR,arg[6],false,lmp);
    xstyle = CONSTANT;
  }
  
  if (strstr(arg[7],"v_") == arg[7]) {
    int n = strlen(&arg[7][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[7][2]);
  } else {
    pxy = utils::numeric(FLERR,arg[7],false,lmp);
    ystyle = CONSTANT;
  }
  
  if (strstr(arg[8],"v_") == arg[8]) {
    int n = strlen(&arg[8][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[8][2]);
  } else {
    pxz = utils::numeric(FLERR,arg[8],false,lmp);
    zstyle = CONSTANT;
  } 
  
  if (strstr(arg[9],"v_") == arg[9]) {
    int n = strlen(&arg[9][2]) + 1;
    qstr = new char[n];
    strcpy(qstr,&arg[9][2]);
  } else {
    dpxx = utils::numeric(FLERR,arg[9],false,lmp);
    qstyle = CONSTANT;
  }   
  
  if (strstr(arg[10],"v_") == arg[10]) {
    int n = strlen(&arg[10][2]) + 1;
    t0_left_str = new char[n];
    strcpy(t0_left_str,&arg[10][2]);
  } else {
    t0_left = utils::numeric(FLERR,arg[10],false,lmp);
    t0_left_style = CONSTANT;
  }  

  if (strstr(arg[11],"v_") == arg[11]) {
    int n = strlen(&arg[11][2]) + 1;
    t0_right_str = new char[n];
    strcpy(t0_right_str,&arg[11][2]);
  } else {
    t0_right = utils::numeric(FLERR,arg[11],false,lmp);
    t0_right_style = CONSTANT;
  }
  
  lambda = utils::numeric(FLERR,arg[12],false,lmp);  
    
  if (seed <= 0) error->all(FLERR,"Illegal fix obmd command");

  //id_press = utils::strdup(arg[10]);

  // read options from end of input line
  varflag = CONSTANT;
  options(narg-13,&arg[13]);

  // error check on type

  if (mode == ATOM && (ntype <= 0 || ntype > atom->ntypes))
    error->all(FLERR,"Invalid atom type in fix obmd command");

  // error checks on region and its extent being inside simulation box

  if (!iregion || !iregion2) error->all(FLERR,"Must specify a region in fix obmd");
  if (iregion->bboxflag == 0 || iregion2->bboxflag == 0)
    error->all(FLERR,"Fix obmd region does not support a bounding box");
  if (iregion->dynamic_check() || iregion2->dynamic_check())
    error->all(FLERR,"Fix obmd region cannot be dynamic");

  xlo = iregion->extent_xlo;
  xhi = iregion->extent_xhi;
  ylo = iregion->extent_ylo;
  yhi = iregion->extent_yhi;
  zlo = iregion->extent_zlo;
  zhi = iregion->extent_zhi;

  if (domain->triclinic == 0) {
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] ||
        ylo < domain->boxlo[1] || yhi > domain->boxhi[1] ||
        zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR,"Insertion region extends outside simulation box");
  } else {
    if (xlo < domain->boxlo_bound[0] || xhi > domain->boxhi_bound[0] ||
        ylo < domain->boxlo_bound[1] || yhi > domain->boxhi_bound[1] ||
        zlo < domain->boxlo_bound[2] || zhi > domain->boxhi_bound[2])
      error->all(FLERR,"Insertion region extends outside simulation box");
  }

  // error check and further setup for mode = MOLECULE

  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use fix_obmd unless atoms have IDs");

  if (mode == MOLECULE) {
    for (int i = 0; i < nmol; i++) {
      if (onemols[i]->xflag == 0)
        error->all(FLERR,"Fix obmd molecule must have coordinates");
      if (onemols[i]->typeflag == 0)
        error->all(FLERR,"Fix obmd molecule must have atom types");
      if (ntype+onemols[i]->ntypes <= 0 ||
          ntype+onemols[i]->ntypes > atom->ntypes)
        error->all(FLERR,"Invalid atom type in fix obmd mol command");

      if (atom->molecular == Atom::TEMPLATE && onemols != atom->avec->onemols)
        error->all(FLERR,"Fix obmd molecule template ID must be same " "as atom_style template ID");
      onemols[i]->check_attributes(); // old check_attributes(0)

      // fix obmd uses geoemetric center of molecule for insertion

      onemols[i]->compute_center();
    }
  }

  if (rigidflag && mode == ATOM)
    error->all(FLERR,"Cannot use fix obmd rigid and not molecule");
  if (shakeflag && mode == ATOM)
    error->all(FLERR,"Cannot use fix obmd shake and not molecule");
  if (rigidflag && shakeflag)
    error->all(FLERR,"Cannot use fix obmd rigid and shake");

  // setup of coords and imageflags array

  if (mode == ATOM) natom_max = 1;
  else {
    natom_max = 0;
    for (int i = 0; i < nmol; i++)
      natom_max = MAX(natom_max,onemols[i]->natoms);
  }
  memory->create(coords,natom_max,3,"adress/obmd:coords");
  memory->create(imageflags,natom_max,"adress/obmd:imageflags");
  
  // setup scaling

  double xscale,yscale,zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  // apply scaling to all input parameters with dist/vel units

  if (domain->dimension == 2) {
    lo *= yscale;
    hi *= yscale;
    rate *= yscale;
  } else {
    lo *= zscale;
    hi *= zscale;
    rate *= zscale;
  }
  deltasq *= xscale*xscale;
  nearsq *= xscale*xscale;
  vxlo *= xscale;
  vxhi *= xscale;
  vylo *= yscale;
  vyhi *= yscale;
  vzlo *= zscale;
  vzhi *= zscale;
  xmid *= xscale;
  ymid *= yscale;
  zmid *= zscale;
  sigma *= xscale; // same as in region sphere
  tx *= xscale;
  ty *= yscale;
  tz *= zscale;

  // find current max atom and molecule IDs if necessary

  if (idnext) find_maxid();

  // random number generator, same for all procs
  // warm up the generator 30x to avoid correlations in first-particle
  // positions if runs are repeated with consecutive seeds

  random = new RanPark(lmp,seed);
  for (int ii=0; ii < 30; ii++) random->uniform();

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

FixAdResSObmd::~FixAdResSObmd()
{
  delete random;
  delete [] molfrac;
  delete [] idrigid;
  delete [] idshake;
  delete [] idregion;
  delete [] idregion2;
  delete [] idregion3;
  delete [] idregion4;
  delete [] idregion5;
  delete [] idregion6;
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

int FixAdResSObmd::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::init()
{
  
  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0)
      error->all(FLERR,"Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix obmd is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0)
      error->all(FLERR,"Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix obmd is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0)
      error->all(FLERR,"Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix obmd is invalid style");
  }
  if (qstr) {
    qvar = input->variable->find(qstr);
    if (qvar < 0)
      error->all(FLERR,"Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(qvar)) qstyle = EQUAL;
    else if (input->variable->atomstyle(qvar)) qstyle = ATOM;
    else error->all(FLERR,"Variable for fix obmd is invalid style");
  }
  if (t0_left_str) {
    t0_left_var = input->variable->find(t0_left_str);
    if (t0_left_var < 0)
      error->all(FLERR,"Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(t0_left_var)) t0_left_style = EQUAL;
    else if (input->variable->atomstyle(t0_left_var)) t0_left_style = ATOM;
    else error->all(FLERR,"Variable for fix obmd is invalid style");
  }   
  if (t0_right_str) {
    t0_right_var = input->variable->find(t0_right_str);
    if (t0_right_var < 0)
      error->all(FLERR,"Variable name for fix obmd does not exist");
    if (input->variable->equalstyle(t0_right_var)) t0_right_style = EQUAL;
    else if (input->variable->atomstyle(t0_right_var)) t0_right_style = ATOM;
    else error->all(FLERR,"Variable for fix obmd is invalid style");
  }   

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM || qstyle == ATOM || t0_left_style == ATOM || t0_right_style == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL || t0_left_style == EQUAL || t0_right_style == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT; 
 
  iregion = domain->get_region_by_id(idregion);
  if (!iregion)
    error->all(FLERR,"Region ID for fix obmd does not exist");

  iregion2 = domain->get_region_by_id(idregion2);
  if (!iregion2)
    error->all(FLERR,"Region ID for fix obmd does not exist");
    
  iregion3 = domain->get_region_by_id(idregion3);
  if (!iregion3)
    error->all(FLERR,"Region ID for fix obmd does not exist");

  iregion4 = domain->get_region_by_id(idregion4);
  if (!iregion4)
    error->all(FLERR,"Region ID for fix obmd does not exist");

  iregion5 = domain->get_region_by_id(idregion5);
  if (!iregion5)
    error->all(FLERR,"Region ID for fix obmd does not exist");
 
  iregion6 = domain->get_region_by_id(idregion6);
  if (!iregion6)
    error->all(FLERR,"Region ID for fix obmd does not exist");
      
  // if rigidflag defined, check for rigid/small fix
  // its molecule template must be same as this one

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
      error->all(FLERR,"Cannot evaporate atoms in atom_modify first group");
  }
  
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
      error->warning(FLERR,"Fix obmd may delete atom with non-zero molecule ID");
  }

  if (molflag && atom->molecule_flag == 0)
      error->all(FLERR,"Fix obmd molecule requires atom attribute molecule"); 

  fixrigid = nullptr;
  if (rigidflag) {
    int ifix = modify->find_fix(idrigid);
    if (ifix < 0) error->all(FLERR,"Fix obmd rigid fix does not exist");
    fixrigid = modify->fix[ifix];
    int tmp;
    if (onemols != (Molecule **) fixrigid->extract("onemol",tmp))
      error->all(FLERR, "Fix obmd and fix rigid/small not using " "same molecule template ID");
  }

  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one

  fixshake = nullptr;
  if (shakeflag) {
    int ifix = modify->find_fix(idshake);
    if (ifix < 0) error->all(FLERR,"Fix obmd shake fix does not exist");
    fixshake = modify->fix[ifix];
    int tmp;
    if (onemols != (Molecule **) fixshake->extract("onemol",tmp))
      error->all(FLERR,"Fix obmd and fix shake not using " "same molecule template ID");
  }

  // for finite size spherical particles:
  // warn if near < 2 * maxrad of existing and inserted particles
  // since may lead to overlaps
  // if inserted molecule does not define diameters,
  // use AtomVecSphere::create_atom() default radius = 0.5

  if (atom->radius_flag) {
    double *radius = atom->radius;
    int nlocal = atom->nlocal;

    double maxrad = 0.0;
    for (int i = 0; i < nlocal; i++)
      maxrad = MAX(maxrad,radius[i]);

    double maxradall;
    MPI_Allreduce(&maxrad,&maxradall,1,MPI_DOUBLE,MPI_MAX,world);

    double maxradinsert = 0.0;
    if (mode == MOLECULE) {
      for (int i = 0; i < nmol; i++) {
        if (onemols[i]->radiusflag)
          maxradinsert = MAX(maxradinsert,onemols[i]->maxradius);
        else maxradinsert = MAX(maxradinsert,0.5);
      }
    } else maxradinsert = 0.5;

    double separation = MAX(2.0*maxradinsert,maxradall+maxradinsert);
    if (sqrt(nearsq) < separation && comm->me == 0)
      error->warning(FLERR,fmt::format("Fix obmd near setting < possible " "overlap separation {}",separation));  
  }
}

void FixAdResSObmd::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style,"^verlet"))
  post_force(vflag);
}

void FixAdResSObmd::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   perform particle/molecule insertions/deletions
------------------------------------------------------------------------- */
void FixAdResSObmd::pre_exchange()
{ 
  // std::cout<<"pre_exchange: "<<"\n";
  int cnt_left, cnt_right, stev_left, stev_right, i;
  double ninsert_left, ninsert_right;
  double masstotal_left,masstotal_right;
  double vcm_internal_left_sq, vcm_internal_right_sq;
  
  double temp_left;
  double temp_right;
  
  // calculate the box dimensions
  double lx = domain->boxhi[0] - domain->boxlo[0];
  double ly = domain->boxhi[1] - domain->boxlo[1];
  double lz = domain->boxhi[2] - domain->boxlo[2];
  
  // fetch the curret value of momentum/energy flux if it is not CONSTANT
  if (varflag != CONSTANT) {
    if (xstyle == EQUAL) pxx = input->variable->compute_equal(xvar);
    if (ystyle == EQUAL) pxy = input->variable->compute_equal(yvar);
    if (zstyle == EQUAL) pxz = input->variable->compute_equal(zvar);
    if (qstyle == EQUAL) dpxx = input->variable->compute_equal(qvar);
    if (t0_left_style == EQUAL) t0_left = input->variable->compute_equal(t0_left_var);
    if (t0_right_style == EQUAL) t0_right = input->variable->compute_equal(t0_right_var);
  }  
  
  // zeroes the momentum, energies
  MathExtra::zero3(vnewl);
  MathExtra::zero3(vnewr);
  // enl = 0.0;
  // enr = 0.0;

  // deletes particles that went across open boundaries after first half of velocity-Verlet algorithm
  try_deleting(iregion, vnewl, vnewr); // , &enl, &enr);
  try_deleting(iregion2, vnewl, vnewr); // , &enl, &enr);

  std::cout<<"after try_deleting ::: 1"<<"\n";

  // counts number of particles in left/right buffer
  cnt_left = group->count(igroup, iregion);
  cnt_right = group->count(igroup, iregion2);
  
  // calculates number of needed insertions
  ninsert_left = -static_cast<int>((static_cast<double>(cnt_left) / mol_len - alpha * nbuf) / tau);
  ninsert_right = -static_cast<int>((static_cast<double>(cnt_right) / mol_len - alpha * nbuf) / tau);

  std::cout<<"cnt_left/mol_len: "<<cnt_left/mol_len<<" mol_len: "<<mol_len<<" alpha: "<<alpha<<" nbuf: "<<nbuf<<" tau: "<<tau<<" ninsert_left: "<<ninsert_left<<"\n";
  std::cout<<"cnt_right/mol_len: "<<cnt_right/mol_len<<" mol_len: "<<mol_len<<" alpha: "<<alpha<<" nbuf: "<<nbuf<<" tau: "<<tau<<" ninsert_right: "<<ninsert_right<<"\n";
  
  // tries inserting ninsert_left/right particles into left/right buffer
  std::cout<<"FixAdResSObmd::pre_exchange(): before try_inserting(iregion5)"<<"\n";
  try_inserting(iregion5, ninsert_left, vnewl, vnewr); // , &enl, &enr);
  // testing only left buffer
  try_inserting(iregion6, ninsert_right, vnewl, vnewr); // , &enl, &enr);
  
  std::cout<<"after try_inserting ::: 1"<<"\n";
  
  // deletes again, just in case part of newly inserted molecule is not inside box
  // is it necessary? TODO: check
  try_deleting(iregion, vnewl, vnewr); // , &enl, &enr);
  try_deleting(iregion2, vnewl, vnewr); // , &enl, &enr);

  // std::cout<<"after try_deleting ::: 2"<<"\n";
  
  // calculates total mass in left/right buffer
  /* masstotal_left = group->mass(igroup, iregion);
  masstotal_right = group->mass(igroup, iregion2); */
  
  // calculates average velocity in left/right buffer
  /* group->vcm(igroup, masstotal_left, vcml, iregion);
  group->vcm(igroup, masstotal_right, vcmr, iregion2); */
  
  // calculates sum of squared velocities in left/right buffer
  /* vcm_internal_sq(igroup, masstotal_left, &vcm_internal_left_sq, &temp_left, vcml, iregion);
  vcm_internal_sq(igroup, masstotal_left, &vcm_internal_right_sq, &temp_right, vcmr, iregion2); */
    
  // sums the outgoing momentum
  MPI_Allreduce(vnewl, vnewl_all, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(vnewr, vnewr_all, 3, MPI_DOUBLE, MPI_SUM, world);
  
  //std::cout << " temp_left = " << temp_left << " temp_left = " << temp_left << std::endl;
  
  // sums the inserted energy
  // MPI_Allreduce(&enl, &enl_all, 1, MPI_DOUBLE, MPI_SUM, world);
  // MPI_Allreduce(&enr, &enr_all, 1, MPI_DOUBLE, MPI_SUM, world);

  // area of buffer - ROI interface
  double area = ly * lz;
  
  // calculate energy currents
  // double energyCurrent_left = -lambda * (temp_left - t0_left);
  // double energyCurrent_right = -lambda * (temp_right - t0_right);
  
  // calculates momentum forces on the left buffer
  // substracting momentum of deleted particles
  std::cout<<"FixAdResSObmd::pre_exchange()"<<"\n";
  std::cout<<"vnewl_all[0]: "<<vnewl_all[0]<<" vnewr_all[0]: "<<vnewr_all[0]<<"\n";
  std::cout<<"vnewl_all[1]: "<<vnewl_all[1]<<" vnewr_all[1]: "<<vnewr_all[1]<<"\n";
  std::cout<<"vnewl_all[2]: "<<vnewl_all[2]<<" vnewr_all[2]: "<<vnewr_all[2]<<"\n";
  std::cout<<"pxx: "<<pxx<<"\n";
  std::cout<<"area: "<<area<<"\n";

  momentumForce_left[0] = vnewl_all[0] / (update->dt) + pxx * area; // -vnewl_all[0] / (update->dt) + pxx * area; // substract ? ? ? // unit_vector = (1.0, 0.0, 0.0) ... toward the center of the simulation box (L->R)
  momentumForce_left[1] = vnewl_all[1] / (update->dt); // -vnewl_all[1] / (update->dt);
  momentumForce_left[2] = vnewl_all[2] / (update->dt); // -vnewl_all[2] / (update->dt);
  shearForce_left[0] = 0.0;
  shearForce_left[1] = pxy * area;
  shearForce_left[2] = pxz * area;

  std::cout<<"FixAdResSObmd::pre_exchange()"<<"\n";
  std::cout<<"momentumForce_left[0]: "<<momentumForce_left[0]<<"\n";
  std::cout<<"momentumForce_left[1]: "<<momentumForce_left[1]<<"\n";
  std::cout<<"momentumForce_left[2]: "<<momentumForce_left[2]<<"\n";

  // calculates momentum forces on the right buffer
  momentumForce_right[0] = vnewr_all[0] / (update->dt) - pxx * area; // -vnewr_all[0] / (update->dt) - pxx * area;  // unit_vector = (-1.0, 0.0, 0.0) ... toward the center of the simulation box (R->L)
  momentumForce_right[1] = vnewr_all[1] / (update->dt); // -vnewr_all[1] / (update->dt);
  momentumForce_right[2] = vnewr_all[2] / (update->dt); // -vnewr_all[2] / (update->dt);
  shearForce_right[0] = 0.0;
  shearForce_right[1] = -pxy * area;
  shearForce_right[2] = -pxz * area;

  std::cout<<"FixAdResSObmd::pre_exchange()"<<"\n"; 
  std::cout<<"momentumForce_right[0]: "<<momentumForce_right[0]<<"\n";
  std::cout<<"momentumForce_right[1]: "<<momentumForce_right[1]<<"\n";
  std::cout<<"momentumForce_right[2]: "<<momentumForce_right[2]<<"\n";
  
  // calculate energy force factors
  // 0 for now: need to do test in the "energy" branch
  // *energyForce_left = 0.0; //area / vcm_internal_left_sq * (energyCurrent_left - enl_all / area / (update->dt) - MathExtra::dot3(momentumForce_left, vcml));
  // *energyForce_right = 0.0; //area / vcm_internal_right_sq * (energyCurrent_right - enr_all / area / (update->dt) - MathExtra::dot3(momentumForce_right, vcmr));

  //int icompute = modify->find_compute(id_press);
  //pressure_tensor = modify->compute[icompute];
  
  //Compute *compute = modify->compute[icompute];
  //int n = 
  //compute->compute_vector();
  //double *compute_vector = compute->vector; //>computepressure_tensor->vector_atom;
  
  //for(i=1;i<2;i++) {
  //  std::cout << "compute_vector["<<i<<"] = " << compute_vector[i] << std::endl;
  //}
  //double a;
  //for(i=0;i<n;i++)
  //a = pressure_tensor->compute_vector(i);
  //std::cout << "pressure_tensor->compute_vector("<<i<<") = " << a << std::endl;
  
  // std::cout<<"before next_reneighbor"<<"\n";
  next_reneighbor += nfreq; //v fix deposit vemo da se bo v pre_exchange nekaj zgodilo... tukaj se lahko izognemo
  // std::cout<<"after next_reneighbor"<<"\n";
}
/* ---------------------------------------------------------------------- */
void FixAdResSObmd::try_deleting(Region* region, double *vnewl, double *vnewr) // ,  double *enl, double *enr)
{
	
  // std::cout<<"try_deleting"<<"\n";
  int i,j,m,iwhichglobal,iwhichlocal;
  int ndel,ndeltopo[4];
  
  // grow list and mark arrays if necessary

  if (atom->nmax > nmax) {
    memory->destroy(list);
    memory->destroy(mark);
    nmax = atom->nmax;
    memory->create(list,nmax,"obmd:list");
    memory->create(mark,nmax,"obmd:mark");
  }

  // ncount = # of deletable atoms in region that I own
  // nall = # on all procs
  // nbefore = # on procs before me
  // list[ncount] = list of local indices of atoms I can delete

  region->prematch();

  double **x = atom->x;
  double **vel = atom->v;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  double boxl = domain->boxlo[0];
  double boxh = domain->boxhi[0];
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *type = atom->type;
  tagint *molecule = atom->molecule;
  int *rep_atom = atom->rep_atom;
  double **cms = atom->cms_mol;
  double **vcms = atom->vcms_mol;
  double fusher[3];
  double massone;
  int imol;
  
  int ncount = 0;
  for (i = 0; i < nlocal; i++) {
    //if (mask[i] & groupbit) // commented out 17.11.2023 -> a molecule that is not a part of the group can also be deleted
      if (x[i][0] < boxl || x[i][0] > boxh) {  
		    std::cout<<"Deleting x[i][0]: "<<x[i][0]<<" i: "<<i<<" type[i]: "<<type[i]<<" tag[i]: "<<tag[i]<<"\n";
        list[ncount++] = i;	  
	  }
  }
  int nall,nbefore;
  MPI_Allreduce(&ncount,&nall,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&ncount,&nbefore,1,MPI_INT,MPI_SUM,world);
  // std::cout<<"nbefore: "<<nbefore<<" ncount: "<<ncount<<" nall: "<<nall<<"\n";
  nbefore -= ncount;
  // std::cout<<"nbefore: "<<nbefore<<" ncount: "<<ncount<<" nall: "<<nall<<"\n";

  // ndel = total # of atom deletions, in or out of region
  // ndeltopo[1,2,3,4] = ditto for bonds, angles, dihedrals, impropers
  // mark[] = 1 if deleted

  ndel = 0;
  for (i = 0; i < nlocal; i++) mark[i] = 0; 

  /* -------------------------------------------------------------------------------------------------- */
  // PPapez COMMENT: testing one idea ... probably it needs to be improved
  // create a list of molecule's ids that are marked for deletion
  // append only once 
  // if id is not on list
  // use it to compute velocity of the center-of-mass 
  // needed for outgoing momentum
  /* std::vector<int> molecules_to_delete;
  std::vector<int> atoms_to_delete;
  // find number of molecules
  nmolecules= find_mols(idlo, idhi);
  maxmol_all = idhi;
  std::cout<<"FixAdResSObmd::try_deleting(): nmolecules: "<<nmolecules<<" maxmol_all: "<<maxmol_all<<"\n";

  // init array with masses
  double *massproc_tmp, *masstotal_tmp;
  memory->create(massproc_tmp,maxmol_all,"adress/obmd:massproc_tmp");
  memory->create(masstotal_tmp,maxmol_all,"adress/obmd:masstotal_tmp");
  // fill with zeros
  for (int i = 0; i < maxmol_all; i++) {
        massproc_tmp[i] = 0.0;
  }
  // compute masses
  for (int i = 0; i < nlocal; i++) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else {
        massone = mass[type[i]];
      }

      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;
      massproc_tmp[imol] += massone;

    }
  }
  // sum and distribute
  MPI_Allreduce(massproc_tmp, masstotal_tmp, maxmol_all, MPI_DOUBLE, MPI_SUM, world); */
  /* -------------------------------------------------------------------------------------------------- */

  // atomic deletions
  // choose atoms randomly across all procs and mark them for deletion
  // shrink eligible list as my atoms get marked
  // keep ndel,ncount,nall,nbefore current after each atom deletion

  if (molflag == 0) {
    // std::cout<<"molflag == 0"<<"\n";
    while (nall) {
      iwhichglobal = static_cast<int> (nall*random->uniform());
      if (iwhichglobal < nbefore) nbefore--;
      else if (iwhichglobal < nbefore + ncount) {
        iwhichlocal = iwhichglobal - nbefore;
        mark[list[iwhichlocal]] = 1;
        list[iwhichlocal] = list[ncount-1];
        ncount--;
      }
      ndel++;
      nall--;
    }

  // molecule deletions
  // bcast mol ID and delete all atoms in that molecule on any proc
  // update deletion count by total # of atoms in molecule
  // shrink list of eligible candidates as any of my atoms get marked
  // keep ndel,ndeltopo,ncount,nall,nbefore current after each mol deletion

  } else {
    // std::cout<<"molflag != 0"<<"\n"; // should be here if deleting molecules
    int me,proc,iatom,ndelone,ndelall,index;
    tagint imolecule;
    tagint *molecule = atom->molecule;
    int *molindex = atom->molindex;
    int *molatom = atom->molatom;
    int molecular = atom->molecular;
    Molecule **onemols = atom->avec->onemols;
	
    ndeltopo[0] = ndeltopo[1] = ndeltopo[2] = ndeltopo[3] = 0;

    while (nall) {
      // std::cout<<"nall: "<<nall<<"\n";
      // pick an iatom,imolecule on proc me to delete

      iwhichglobal = static_cast<int> (nall*random->uniform());
      if (iwhichglobal >= nbefore && iwhichglobal < nbefore + ncount) {
        iwhichlocal = iwhichglobal - nbefore;
        iatom = list[iwhichlocal];
        imolecule = molecule[iatom];
        me = comm->me;
        // std::cout<<"iatom: "<<iatom<<" imolecule: "<<imolecule<<"\n";
        // std::cout<<"tag[iatom]: "<<tag[iatom]<<"\n";
      } else me = -1;

      // bcast mol ID to delete all atoms from
      // if mol ID > 0, delete any atom in molecule and decrement counters
      // if mol ID == 0, delete single iatom
      // logic with ndeltopo is to count # of deleted bonds,angles,etc
      // for atom->molecular = Atom::MOLECULAR, do this for each deleted atom in molecule
      // for atom->molecular = Atom::TEMPLATE, use Molecule counts for just 1st atom in mol

      MPI_Allreduce(&me,&proc,1,MPI_INT,MPI_MAX,world);
      MPI_Bcast(&imolecule,1,MPI_LMP_TAGINT,proc,world);
      ndelone = 0;
      for (i = 0; i < nlocal; i++) {
        if (imolecule && molecule[i] == imolecule) {
          std::cout<<"imolecule: "<<imolecule<<"\n";
          mark[i] = 1;
          ndelone++;
          // std::cout<<"imolecule: "<<imolecule<<" tag[i]: "<<tag[i]<<" type[i]: "<<type[i]<<"\n";
          if (molecular == Atom::MOLECULAR) {
            // std::cout<<"molecular == Atom::MOLECULAR"<<"\n";
            if (atom->avec->bonds_allow) {
              if (force->newton_bond) ndeltopo[0] += atom->num_bond[i];
              else {
                for (j = 0; j < atom->num_bond[i]; j++) {
                  if (tag[i] < atom->bond_atom[i][j]) ndeltopo[0]++;
                }
              }
            }
            if (atom->avec->angles_allow) {
              if (force->newton_bond) ndeltopo[1] += atom->num_angle[i];
              else {
                for (j = 0; j < atom->num_angle[i]; j++) {
                  m = atom->map(atom->angle_atom2[i][j]);
                  if (m >= 0 && m < nlocal) ndeltopo[1]++;
                }
              }
            }
            if (atom->avec->dihedrals_allow) {
              if (force->newton_bond) ndeltopo[2] += atom->num_dihedral[i];
              else {
                for (j = 0; j < atom->num_dihedral[i]; j++) {
                  m = atom->map(atom->dihedral_atom2[i][j]);
                  if (m >= 0 && m < nlocal) ndeltopo[2]++;
                }
              }
            }
            if (atom->avec->impropers_allow) {
              if (force->newton_bond) ndeltopo[3] += atom->num_improper[i];
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
      // std::cout<<"ncount: "<<ncount<<"\n";
      while (i < ncount) {
        if (mark[list[i]]) {
          list[i] = list[ncount-1];
          ncount--;
        } else i++;
      }

      // update ndel,ncount,nall,nbefore
      // ndelall is total atoms deleted on this iteration
      // ncount is already correct, so resum to get nall and nbefore

      MPI_Allreduce(&ndelone,&ndelall,1,MPI_INT,MPI_SUM,world);
      ndel += ndelall;
      MPI_Allreduce(&ncount,&nall,1,MPI_INT,MPI_SUM,world);
      MPI_Scan(&ncount,&nbefore,1,MPI_INT,MPI_SUM,world);
      nbefore -= ncount;
      // std::cout<<"ndel: "<<ndel<<"\n";
    }
  }
  
  // PPapez COMMENT: compute outgoing momentum 

  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms

  AtomVec *avec = atom->avec;
  int count_deleted_atoms = 0;
  for (i = nlocal-1; i >= 0; i--) {
    if (mark[i]) {
        std::cout<<"atom->tag[i]: "<<atom->tag[i]<<"\n";
        std::cout<<"vel[i][0]: "<<vel[i][0]<<" vel[i][1]: "<<vel[i][1]<<" vel[i][2]: "<<vel[i][2]<<"\n";
        std::cout<<"vcms[i][0]: "<<vcms[i][0]<<" vcms[i][1]: "<<vcms[i][1]<<" vcms[i][2]: "<<vcms[i][2]<<"\n";
        std::cout<<"mtot: "<<mtot<<"\n";

        // PPapez COMMENT: use rep_atom and center-of-mass
        if (rep_atom[i] == 0) {
          std::cout<<"rep_atom[i]: "<<rep_atom[i]<<"\n";
          std::cout<<"type[i]: "<<type[i]<<"\n";
          /* std::cout<<"rep_atom != 1"<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<" type[i]: "<<type[i]<<"\n";
          std::cout<<"vel[i][0]: "<<vel[i][0]<<" vel[i][1]: "<<vel[i][1]<<" vel[i][2]: "<<vel[i][2]<<"\n";
          std::cout<<"vcms[i][0]: "<<vcms[i][0]<<" vcms[i][1]: "<<vcms[i][1]<<" vcms[i][2]: "<<vcms[i][2]<<"\n"; */
          count_deleted_atoms += 1;
          avec->copy(atom->nlocal-1,i,1);
          atom->nlocal--;
          continue;
        }
        // std::cout<<"rep_atom[i]: "<<rep_atom[i]<<"\n";
        if (cms[i][0] < 0.5 * (boxh + boxl)) {
          std::cout<<"rep_atom == 1"<<"\n";
          std::cout<<"rep_atom[i]: "<<rep_atom[i]<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<" type[i]: "<<type[i]<<"\n";
          std::cout<<"vel[i][0]: "<<vel[i][0]<<" vel[i][1]: "<<vel[i][1]<<" vel[i][2]: "<<vel[i][2]<<"\n";
          std::cout<<"vcms[i][0]: "<<vcms[i][0]<<" vcms[i][1]: "<<vcms[i][1]<<" vcms[i][2]: "<<vcms[i][2]<<"\n"; 
          delete_left_file<<"deleted"<<"\n";
          delete_left_file.flush();
          count_deleted_atoms += 1;
          vnewl[0] += mtot * vcms[i][0];
          vnewl[1] += mtot * vcms[i][1];
          vnewl[2] += mtot * vcms[i][2];
          avec->copy(atom->nlocal-1,i,1);
          atom->nlocal--;
        }
        else {
          std::cout<<"rep_atom == 1"<<"\n";
          std::cout<<"rep_atom[i]: "<<rep_atom[i]<<"\n";
          std::cout<<"atom->tag[i]: "<<atom->tag[i]<<" type[i]: "<<type[i]<<"\n";
          std::cout<<"vel[i][0]: "<<vel[i][0]<<" vel[i][1]: "<<vel[i][1]<<" vel[i][2]: "<<vel[i][2]<<"\n";
          std::cout<<"vcms[i][0]: "<<vcms[i][0]<<" vcms[i][1]: "<<vcms[i][1]<<" vcms[i][2]: "<<vcms[i][2]<<"\n"; 
          delete_right_file<<"deleted"<<"\n";
          delete_right_file.flush();
          count_deleted_atoms += 1;
          vnewr[0] += mtot * vcms[i][0];
          vnewr[1] += mtot * vcms[i][1];
          vnewr[2] += mtot * vcms[i][2];
          avec->copy(atom->nlocal-1,i,1);
          atom->nlocal--;
        } 
      
		/* if(x[i][0] < 0.5*(boxh+boxl)) {
       count_deleted_atoms += 1;
       std::cout<<"rep_atom[i]: "<<rep_atom[i]<<"\n";
		  // *enl -= energy(i, type[i], x[i], fusher);
	      vnewl[0] += mass[type[i]] * vel[i][0]; //  / mtot; // -= mass[type[i]]*vel[i][0] / mtot;
	      vnewl[1] += mass[type[i]] * vel[i][1]; // / mtot; // -= mass[type[i]]*vel[i][1] / mtot;
	      vnewl[2] += mass[type[i]] * vel[i][2]; // / mtot; // -= mass[type[i]]*vel[i][2] / mtot;
        } else {
        count_deleted_atoms += 1;
	     //  *enr -= energy(i, type[i], x[i], fusher);
	      vnewr[0] += mass[type[i]] * vel[i][0]; //  / mtot; // -= mass[type[i]]*vel[i][0] / mtot;
	      vnewr[1] += mass[type[i]] * vel[i][1]; // / mtot; // -= mass[type[i]]*vel[i][1] / mtot;
	      vnewr[2] += mass[type[i]] * vel[i][2]; // / mtot; // -= mass[type[i]]*vel[i][2] / mtot;
        } */
      // avec->copy(atom->nlocal-1,i,1);
      // atom->nlocal--;
    }
  }
  std::cout<<"count_deleted_atoms: "<<count_deleted_atoms<<"\n";
  // reset global natoms and bonds, angles, etc
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  atom->natoms -= ndel;
  // std::cout<<"atom->natoms: "<<atom->natoms<<"\n";
  if (molflag) {
    int all[4];
    MPI_Allreduce(ndeltopo,all,4,MPI_INT,MPI_SUM,world);
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
  //next_reneighbor = update->ntimestep + nevery; // update every timestep ? ? ?
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::try_inserting(Region* iregion_var, int stev, double *vnewl, double *vnewr) // , double *enl, double *enr)
{ 
  // std::cout<<"try_inserting"<<"\n";
  // std::cout<<"stev: "<<stev<<"\n";

  int i,m,n,nlocalprev,imol,natom,flag,flagall,me,sum_all,ninsert;
  double coord[3],lamda[3],delx,dely,delz,rsq;
  double r[3],vnew[3],rotmat[3][3],quat[4];
  double *newcoord;
  double offset = 0.0;
  if (rateflag) offset = (update->ntimestep - nfirst) * update->dt * rate;
  double *sublo,*subhi;
  int iter;
  
  double boxl = domain->boxlo[0];
  double boxh = domain->boxhi[0];

  double *mass = atom->mass;
  double mtmp;
    
  xlo = iregion_var->extent_xlo;
  xhi = iregion_var->extent_xhi;
  ylo = iregion_var->extent_ylo;
  yhi = iregion_var->extent_yhi;
  zlo = iregion_var->extent_zlo;
  zhi = iregion_var->extent_zhi;
  
  double entmp;
  
  if (domain->triclinic == 0) {
    sublo = domain->sublo;
    subhi = domain->subhi;
  } else {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  }
  
  if (atom->map_style != Atom::MAP_NONE) atom->map_clear();
  atom->nghost = 0;
  atom->avec->clear_bonus();
  
  if (!idnext) find_maxid();  
    
  int dimension = domain->dimension;
  
  if(stev>0) {
    for(ninsert=0; ninsert!=stev; ninsert++) { 
      int success = 0;
      int attempt = 0;      
      while (attempt < maxattempt) {
        // PP DEBUG
        // feenableexcept(FE_ALL_EXCEPT);
        attempt++;
        // choose random position for new particle within region
        if (distflag == DIST_UNIFORM) {
          do {
            coord[0] = xlo + random->uniform() * (xhi-xlo);
            coord[1] = ylo + random->uniform() * (yhi-ylo);
            coord[2] = zlo + random->uniform() * (zhi-zlo);
          } while (iregion_var->match(coord[0],coord[1],coord[2]) == 0);
        } else if (distflag == DIST_GAUSSIAN) {
          do {
            coord[0] = xmid + random->gaussian() * sigma;
            coord[1] = ymid + random->gaussian() * sigma;
            coord[2] = zmid + random->gaussian() * sigma;
          } while (iregion_var->match(coord[0],coord[1],coord[2]) == 0);
        } else error->all(FLERR,"Unknown particle distribution in fix obmd");
        
        // adjust vertical coord by offset
        if (dimension == 2) coord[1] += offset;
        else coord[2] += offset;
        
        // if global, reset vertical coord to be lo-hi above highest atom
        // if local, reset vertical coord to be lo-hi above highest "nearby" atom
        // local computation computes lateral distance between 2 particles w/ PBC
        // when done, have final coord of atom or center pt of molecule
        
        if (globalflag || localflag) {
          int dim;
          double max,maxall,delx,dely,delz,rsq;
        
          if (dimension == 2) {
            dim = 1;
            max = domain->boxlo[1];
          } else {
            dim = 2;
            max = domain->boxlo[2];
          }
        
          double **x = atom->x;
          int nlocal = atom->nlocal;
          for (i = 0; i < nlocal; i++) {
            if (localflag) {
              delx = coord[0] - x[i][0];
              dely = coord[1] - x[i][1];
              delz = 0.0;
              domain->minimum_image(delx,dely,delz);
              if (dimension == 2) rsq = delx*delx;
              else rsq = delx*delx + dely*dely;
              if (rsq > deltasq) continue;
            }
            if (x[i][dim] > max) max = x[i][dim];
          }
        
          MPI_Allreduce(&max,&maxall,1,MPI_DOUBLE,MPI_MAX,world);
          if (dimension == 2)
            coord[1] = maxall + lo + random->uniform()*(hi-lo);
          else
            coord[2] = maxall + lo + random->uniform()*(hi-lo);
        }
        
        // coords = coords of all atoms
        // for molecule, perform random rotation around center pt
        // apply PBC so final coords are inside box
        // also modify image flags due to PBC

        if (mode == ATOM) {
          natom = 1;
          coords[0][0] = coord[0];
          coords[0][1] = coord[1];
          coords[0][2] = coord[2];
          imageflags[0] = ((imageint) IMGMAX << IMG2BITS) |
            ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        } else { // MOLECULE ... working 
          double rng = random->uniform();
          imol = 0;
          while (rng > molfrac[imol]) imol++;
          natom = onemols[imol]->natoms;
          if (dimension == 3) {
            if (orientflag) {
              r[0] = rx;
              r[1] = ry;
              r[2] = rz;
            } else {
              r[0] = random->uniform() - 0.5;
              r[1] = random->uniform() - 0.5;
              r[2] = random->uniform() - 0.5;
            }
          } else {
            r[0] = r[1] = 0.0;
            r[2] = 1.0;
          }
          double theta = random->uniform() * MY_2PI;
          MathExtra::norm3(r); // unit vector 
          MathExtra::axisangle_to_quat(r,theta,quat); // quaternion
          MathExtra::quat_to_mat(quat,rotmat); // rotmat
          for (i = 0; i < natom; i++) {
            MathExtra::matvec(rotmat,onemols[imol]->dx[i],coords[i]); // getting coords after rotation of dx[i] using rotmat
            coords[i][0] += coord[0]; // rotated coords // vrstica - stolpec
            coords[i][1] += coord[1]; // rotated coords
            coords[i][2] += coord[2]; // rotated coords
            // PPapez COMMENT: due to this routine select / limit insertion region (reasonable dist from its edge)
            imageflags[i] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
            domain->remap(coords[i],imageflags[i]);
          }
        }

        // PP ... OMITTING NEAR
        // check distance between any existing atom and any inserted atom
        // if less than near, try again
        // use minimum_image() to account for PBC
      
        double **x = atom->x;
        int nlocal = atom->nlocal;
        flag = 0;
        if(nearflag) {
          for (m = 0; m < natom; m++) {
            ; // PP no near option
            /* for (i = 0; i < nlocal; i++) {
              delx = coords[m][0] - x[i][0];
              dely = coords[m][1] - x[i][1];
              delz = coords[m][2] - x[i][2];
              domain->minimum_image(delx,dely,delz);
              rsq = delx*delx + dely*dely + delz*delz;
              if (rsq < nearsq) {
                // if(comm->me==0) std::cout << "NEAR denies in attempt No. " << attempt << "." << std::endl;
                flag = 1;
              }
            } */
          }
        } else if(usherflag) { // USING THIS ! ! ! 
          /* ============================================================================================================================================ */
          // USHER working
          std::cout<<"try_inserting -> usher(coords, etarget, natom, imol, iter)"<<"\n";
          std::cout<<"P1: coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
          std::cout<<"P2: coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
          std::cout<<"P3: coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 
          
	        me = -1;
	        //sum_all = check_mol_proc(coords,newcoord,lamda,sublo,subhi,dimension,natom,me);
	        entmp = usher(iregion_var, coords, etarget, natom, imol, iter); 

          std::cout<<"after usher"<<"\n";
          std::cout<<"P1: coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
          std::cout<<"P2: coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
          std::cout<<"P3: coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 
        
	        /// preveri še nghost ... za usher nall, nlocal, nghost
	        if(entmp < etarget+EPSILON) {
	       	  if(comm->me==0) std::cout << "USHER accepts at E = " << entmp << " in attempt No. " << attempt << " with " << iter << " iterations" << std::endl;
            ;
          } else { 
	       	  if(comm->me==0) std::cout << "USHER denies at E = " << entmp << " at attempt No. " << attempt << std::endl; 
	       	  flag=1;
	        }
	      }
        
        MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);
        if (flagall) continue;
              
        if(nearflag) {
          entmp = near_energy(coords, natom, imol);
          //if(comm->me==0) std::cout << "NEAR inserts at E = " << entmp << " in attempt No. " << attempt << "." <<  std::endl;
	      }
	      
        // proceed with insertion
        nlocalprev = atom->nlocal;
        
        // choose random velocity for new particle
        // used for every atom in molecule
        // PPapez COMMENT: for now I'm inserting with zero velocity
        vnew[0] = 0.0; // vxlo + random->uniform() * (vxhi-vxlo);
        vnew[1] = 0.0; // vylo + random->uniform() * (vyhi-vylo);
        vnew[2] = 0.0; // vzlo + random->uniform() * (vzhi-vzlo);
	      
        // if target specified, change velocity vector accordingly
        if (targetflag) {
          double vel = sqrt(vnew[0]*vnew[0] + vnew[1]*vnew[1] + vnew[2]*vnew[2]);
          delx = tx - coord[0];
          dely = ty - coord[1];
          delz = tz - coord[2];
          double rsq = delx*delx + dely*dely + delz*delz;
          if (rsq > 0.0) {
            double rinv = sqrt(1.0/rsq);
            vnew[0] = delx*rinv*vel;
            vnew[1] = dely*rinv*vel;
            vnew[2] = delz*rinv*vel;
          }
        }
        
        // this loop first just checks if an entire molecule can be inserted
	      sum_all = check_mol_proc(coords, newcoord, lamda, sublo, subhi, dimension, natom, me);
           
        if(sum_all != natom) {
          if(comm->me==0) std::cout << "Can not insert the particle/molecule at timestep = " << update->ntimestep << " me = " << me << std::endl;
	        continue;
        }
          
        // happens rarely, e.g. if usher moves the atom/molecule outside of the (open) boundaries
        int check;
        check = check_mol_region(iregion_var,coords,natom);
        if(check == 1) continue;       
        
        // check if new atoms are in my sub-box or above it if I am highest proc
        // if so, add atom to my list via create_atom()
        // initialize additional info about the atoms
        // set group mask to "all" plus fix group

        // NOT OPTIMAL SOLUTION ... TESTING FOR NOW
        // store all ilocal ids of newly inserted particles of a molecule
        std::vector<int> locals;
        for (m = 0; m < natom; m++) {   
          if (domain->triclinic) {
            domain->x2lamda(coords[m],lamda);
            newcoord = lamda;
          } else newcoord = coords[m];
          flag = check_proc(coords, newcoord, sublo, subhi, dimension);             				
          if (flag) {
            if (mode == ATOM) {
	          mtmp = mass[ntype];
	          atom->avec->create_atom(ntype,coords[m]);
	        } else { 
	          mtmp = mtot;
            atom->avec->create_atom(ntype+onemols[imol]->type[m],coords[m]); // ,onemols[imol]->rep_atom[m]); 
	        }
            
          n = atom->nlocal-1;
          atom->tag[n] = maxtag_all+m+1;
          if (mode == MOLECULE) {
            if (atom->molecule_flag) {
              if (onemols[imol]->moleculeflag) {
                atom->molecule[n] = maxmol_all + onemols[imol]->molecule[m];
              } else {
                atom->molecule[n] = maxmol_all+1;
              }
            }
            if (atom->molecular == Atom::TEMPLATE) {
              atom->molindex[n] = 0;
              atom->molatom[n] = m;
            }
          }

          atom->mask[n] = 1 | groupbit;
          atom->image[n] = imageflags[m];
          atom->v[n][0] = vnew[0];
          atom->v[n][1] = vnew[1];
          atom->v[n][2] = vnew[2];

          // PP 
          // inserting in CG region (lambda is zero) ... fix_adress to pre_exchange()
          // heterogeneus buffer ... still inserting in CG region (unwritten rule ... insert in the outer half of CG region)
          atom->lambdaF[n] = 0.0; 

          if (mode == MOLECULE) {
            onemols[imol]->quat_external = quat;
            atom->add_molecule_atom(onemols[imol],m,n,maxtag_all); // creating atom
            std::cout<<"inserted particle: atom->tag[n]: "<<atom->tag[n]<<" atom->x[n][0]: "<<atom->x[n][0]<<" lambda[n]: "<<atom->lambdaF[n]<<"\n";
            std::cout<<"atom->type[n]: "<<atom->type[n]<<" atom->rep_atom[n]: "<<atom->rep_atom[n]<<"\n";
            std::cout<<"atom->x[n][0]: "<<atom->x[n][0]<<" atom->x[n][1]: "<<atom->x[n][1]<<" atom->x[n][2]: "<<atom->x[n][2]<<"\n";
            std::cout<<"atom->q[n]: "<<atom->q[n]<<"\n";
          }
            
          modify->create_attribute(n);

          // fill nlocals
          locals.push_back(n);

          } // closing if(flag)
        }

        // loop over nlocals and compute center-of-mass
        /// double com_tmp[3];
        com_tmp[0] = com_tmp[1] = com_tmp[2] = 0.0;
        for (int lngth = 0; lngth < locals.size(); lngth++) {         
          com_tmp[0] += mass[atom->type[locals[lngth]]] * atom->x[locals[lngth]][0] / mtot;
          com_tmp[1] += mass[atom->type[locals[lngth]]] * atom->x[locals[lngth]][1] / mtot;
          com_tmp[2] += mass[atom->type[locals[lngth]]] * atom->x[locals[lngth]][2] / mtot;
        } 

        // std::cout<<"com_tmp[0]: "<<com_tmp[0]<<" com_tmp[1]: "<<com_tmp[1]<<" com_tmp[2]: "<<com_tmp[2]<<"\n";
        MPI_Allreduce(&com_tmp,&com_tmp_all,3,MPI_DOUBLE,MPI_SUM,world);

        // loop again to set atom's center of mass property 
        int locals_size = locals.size();
        for (int lngth = 0; lngth < locals_size; lngth++) { // locals.size(); lngth++) {
          atom->cms_mol[locals[lngth]][0] = com_tmp_all[0]; // com_tmp[0];
          atom->cms_mol[locals[lngth]][1] = com_tmp_all[1]; // com_tmp[1];
          atom->cms_mol[locals[lngth]][2] = com_tmp_all[2]; // com_tmp[2];

          /* std::cout<<"atom->cms_mol[locals[lngth]][0]: "<<atom->cms_mol[locals[lngth]][0]<<"\n";
          std::cout<<"atom->cms_mol[locals[lngth]][1]: "<<atom->cms_mol[locals[lngth]][1]<<"\n";
          std::cout<<"atom->cms_mol[locals[lngth]][2]: "<<atom->cms_mol[locals[lngth]][2]<<"\n"; */
        }     

        // FixRigidSmall::set_molecule stores rigid body attributes
        // coord is new position of geometric center of mol, not COM
        // FixShake::set_molecule stores shake info for molecule
        if (mode == MOLECULE) {    
          // std::cout<<"mode == MOLECULE ::: 2"<<"\n";
          if (rigidflag)
            fixrigid->set_molecule(nlocalprev,maxtag_all,imol,coord,vnew,quat);
          else if (shakeflag)
            fixshake->set_molecule(nlocalprev,maxtag_all,imol,coord,vnew,quat);
	      }
	      
        success = 1;
        break;
      }

      if (!success && comm->me == 0) error->warning(FLERR,"Particle/molecule insertion was unsuccessful");
      
      // reset global natoms,nbonds,etc
      // increment maxtag_all and maxmol_all if necessary
      // if global map exists, reset it now instead of waiting for comm
      // since other pre-exchange fixes may use it
      // invoke map_init() b/c atom count has grown
        
      if(success) {
        std::cout<<"if(success)"<<"\n";
        std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
        std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
        std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 
        if(coords[0][0] < 0.5*(boxh+boxl)) {
          insert_left_file<<"inserted"<<"\n";
          insert_left_file.flush();
	        vnewl[0] += 0.0; // mtmp*vnew[0]; // momentum in
          vnewl[1] += 0.0; // mtmp*vnew[1];
          vnewl[2] += 0.0; // mtmp*vnew[2];
        } else {
          insert_right_file<<"inserted"<<"\n";
          insert_right_file.flush();
          vnewr[0] += 0.0; //  mtmp*vnew[0];
          vnewr[1] += 0.0; // mtmp*vnew[1];
          vnewr[2] += 0.0; // mtmp*vnew[2];		
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

void FixAdResSObmd::post_force(int vflag)
{
  if (update->ntimestep % nevery) { return ;}

  std::cout<<"FixAdResSObmd::post_force(): stepflag: "<<stepflag<<"\n";
  std::cout<<"momentumForce_left[0]: "<<momentumForce_left[0]<<" momentumForce_left[1]: "<<momentumForce_left[1]<<" momentumForce_left[2]: "<<momentumForce_left[2]<<"\n";
  std::cout<<"momentumForce_right[0]: "<<momentumForce_right[0]<<" momentumForce_right[1]: "<<momentumForce_right[1]<<" momentumForce_right[2]: "<<momentumForce_right[2]<<"\n";

  std::cout<<"before reg_force for LEFT BUFFER"<<"\n";
  reg_force(vflag, iregion, momentumForce_left, stepflag);
  std::cout<<"before reg_force for RIGHT BUFFER"<<"\n";
  reg_force(vflag, iregion2, momentumForce_right, stepflag);
  // reg_force_perp(vflag, iregion3, shearForce_left, stepflag);
  // reg_force_perp(vflag, iregion4, shearForce_right, stepflag);  
  // add_energy_forces(vflag, iregion, energyForce_left, vcml);
  // add_energy_forces(vflag, iregion2, energyForce_right, vcmr);  
}

/* ---------------------------------------------------------------------- */

/* void FixObmd::add_energy_forces(int vflag, Region* region, double *fac, double *vcm)
{
  if (region) region->prematch(); 
  v_init(vflag);
  
  double **x = atom->x;
  double **vel = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  imageint *image = atom->image;
  double v[6];
  int nlocal = atom->nlocal;
  double unwrap[3];
  
  double xval, yval, zval;
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (region && !region->match(x[i][0],x[i][1],x[i][2])) continue;
      domain->unmap(x[i],image[i],unwrap);
      
      xval = *fac * (vel[i][0] - vcm[0]);
      yval = *fac * (vel[i][1] - vcm[1]);
      zval = *fac * (vel[i][2] - vcm[2]);
      
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
        v_tally(i,v);
      }
    }
} */

/* --------------------------------------------------------------------------------------------- */

double FixAdResSObmd::g_par_global(Region* region, int step, double* masstotal_tmp)
{  
  // using molecule's center-of-mass
  std::cout<<"FixAdResSObmd::g_par_global()"<<"\n";
  if (region) region->prematch();
  
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  int *rep_atom = atom->rep_atom;
  double **cms = atom->cms_mol; 
  int imol; 

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
  if (!step) {
    for (int i = 0; i < nlocal; i++) {

      // check if representative atom
      if (rep_atom[i] == 0) {
        continue;
      }

      // find molecule
      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;

      // mass of a molecule
      double molmass = masstotal_tmp[imol];

      // (looping over nlocals) check if in prescribed region, continue of not 
      if (region && !region->match(cms[i][0], cms[i][1], cms[i][2])) continue;

      // if here ... molecule in buffer
      if (cms[i][0] < lower_x + buffer_size) { // LEFT BUFFER
        // std::cout<<"LEFT BUFFER"<<"\n";
        // std::cout<<"cms[i][0]: "<<cms[i][0]<<" lower_x + buffer_size: "<<lower_x + buffer_size<<"\n"; 
        if (cms[i][0] < (lower_x + (1.0 - g_fac) * buffer_size)) { // only mass
          // std::cout<<"LEFT if"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<"(lower_x + (1.0 - g_fac) * buffer_size): "<<(lower_x + (1.0 - g_fac) * buffer_size)<<"\n"; 
          
          g_par_all += molmass;
        }
        else { // sigmoidal mass          
          carg = 1.0 / g_fac * MY_PI * (cms[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI; // OK (also for -lower_x)
          g_par_all += 0.5 * (1.0 + cos(carg)) * molmass;

          // std::cout<<"LEFT else"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<" 1.0 / g_fac * MY_PI * (cms[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI: "<<1.0 / g_fac * MY_PI * (cms[i][0] - buffer_size - lower_x) / (-buffer_size) - MY_PI<<"\n";
          // std::cout<<"0.5 * (1.0 + cos(carg)) * molmass: "<<0.5 * (1.0 + cos(carg)) * molmass<<"\n"; 
        }
      }
      // test LEFT buffer ... OK ! ! ! 
      if (cms[i][0] > upper_x - buffer_size) { // RIGHT BUFFER
        // std::cout<<"RIGHT BUFFER"<<"\n";
        // std::cout<<"cms[i][0]: "<<cms[i][0]<<" upper_x - buffer_size: "<<upper_x - buffer_size<<"\n";
        if (cms[i][0] > (upper_x - (1.0 - g_fac) * buffer_size)) {
          // std::cout<<"RIGHT if"<<"\n";
          // std::cout<<"molmass: "<<molmass<<" cms[i][0]: "<<cms[i][0]<<"(upper_x - (1.0 - g_fac) * buffer_size): "<<(upper_x - (1.0 - g_fac) * buffer_size)<<"\n"; 

          g_par_all += molmass;
        }
        else {
          carg = 1.0 / g_fac * MY_PI * (cms[i][0] - upper_x + buffer_size) / (buffer_size) - MY_PI; // OK (also for -lower_x)
          g_par_all += 0.5 * (1.0 + cos(carg)) * molmass;

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

double FixAdResSObmd::g_par_local(double cms_x, double molmass, int step)
{
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

  if (!step) {
    if (cms_x < lower_x + buffer_size) { // LEFT BUFFER
      if (cms_x < (lower_x + (1.0 - g_fac) * buffer_size)) {
        return molmass;
      }
      else { // sigmodial mass
        carg = 1.0 / g_fac * MY_PI * (cms_x - buffer_size - lower_x) / (-buffer_size) - MY_PI;
        return 0.5 * (1.0 + cos(carg)) * molmass;
      }
    }
    if (cms_x > upper_x - buffer_size) { // RIGHT BUFFER
      if (cms_x > (upper_x - (1.0 - g_fac) * buffer_size)) {
        return molmass;
      }
      else { // sigmodial mass
        carg = 1.0 / g_fac * MY_PI * (cms_x - upper_x + buffer_size) / (buffer_size) - MY_PI;
        return 0.5 * (1.0 + cos(carg)) * molmass;
      }
    }
    // OTHER IS NOT OF MY INTEREST ... ROI
  }
  else {
    error->all(FLERR, "For now implemented only for !step");
  } 
}

/* ---------------------------------------------------------------------- */

double FixAdResSObmd::g_perp_global(Region* iregion_var, int step)
{
  return 1.0*group->count(igroup,iregion_var);
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::reg_force(int vflag, Region* region, double *momentumForce, int step)
{ 
  std::cout<<"FixAdResSObmd::reg_force()"<<"\n";
  if (region) region->prematch(); 
  v_init(vflag);
  
  double **f = atom->f;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int *type = atom->type;
  // external force acts on molecules in buffer
  tagint *molecule = atom->molecule;
  int *rep_atom = atom->rep_atom;
  double **cms = atom->cms_mol;
  double massone;
  tagint imol;
  // create test_force to check obmd
  double test_force[3];
  test_force[0] = test_force[1] = test_force[2] = 0.0;

  // find number of molecules
  nmolecules= find_mols(idlo, idhi);
  maxmol_all = idhi;
  std::cout<<"FixAdResSObmd::reg_force(): nmolecules: "<<nmolecules<<" maxmol_all: "<<maxmol_all<<"\n";

  // init array with masses
  double *massproc_tmp, *masstotal_tmp;
  memory->create(massproc_tmp,maxmol_all,"adress/obmd:massproc_tmp");
  memory->create(masstotal_tmp,maxmol_all,"adress/obmd:masstotal_tmp");
  // fill with zeros
  for (int i = 0; i < maxmol_all; i++) {
        massproc_tmp[i] = 0.0;
  }
  // compute masses
  for (int i = 0; i < nlocal; i++) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else {
        massone = mass[type[i]];
      }

      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;
      massproc_tmp[imol] += massone;

    }
  }
  // sum and distribute
  MPI_Allreduce(massproc_tmp, masstotal_tmp, maxmol_all, MPI_DOUBLE, MPI_SUM, world);

  // compute "global" distribution function 
  // contribution of all particles within buffer
  // added masstotal_tmp info ... different molecules have different masses (for now only "one type" of a molecule) ? ? ? 
  std::cout<<"FixAdResSObmd::reg_force(): before g_par_global"<<"\n";
  double g_all = g_par_global(region, step, masstotal_tmp); // OK 

  // loop over all nlocal atoms
  // init external forces
  double **mol_f_ext, **mol_f_ext_all;
  memory->create(mol_f_ext, maxmol_all, 3, "adress/obmd:mol_f_ext");
  memory->create(mol_f_ext_all, maxmol_all, 3, "adress/obmd:mol_f_ext_all");
  // fill with zeros
  for (int i = 0; i < maxmol_all; i++) {
    mol_f_ext[i][0] = mol_f_ext[i][1] = mol_f_ext[i][2] = 0.0;
    mol_f_ext_all[i][0] = mol_f_ext_all[i][1] = mol_f_ext_all[i][2] = 0.0;
  }

  // loop over nlocals
  int counting_distribute = 0;
  double sestevek = 0.0;
  for (int i = 0; i < nlocal; i++) {
    // check if representative particle
    if (rep_atom[i] == 0) {
      continue;
    }

    // if (mask[i] & groupbit) { // not necessary in group ... consider future applications
      // std::cout<<"1: cms[i][0]: "<<cms[i][0]<<" cms[i][1]: "<<cms[i][1]<<" cms[i][2]: "<<cms[i][2]<<"\n";
      // check if in region
      if (region && !region->match(cms[i][0],cms[i][1],cms[i][2])) continue; // REGION ! ! ! 
      // std::cout<<"2: cms[i][0]: "<<cms[i][0]<<" cms[i][1]: "<<cms[i][1]<<" cms[i][2]: "<<cms[i][2]<<"\n";

      // find molecule
      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;

      // mass of a molecule
      double molmass = masstotal_tmp[imol];

      // compute "local" distribution function
      // removed first arg region 
      double g_one = g_par_local(cms[i][0], molmass, step);
      // std::cout<<"g_one: "<<g_one<<" g_all: "<<g_all<<"\n";

      // compute external force ! ! !
      mol_f_ext[imol][0] += momentumForce[0] * g_one / g_all;
      mol_f_ext[imol][1] += momentumForce[1] * g_one / g_all;
      mol_f_ext[imol][2] += momentumForce[2] * g_one / g_all;

      test_force[0] += momentumForce[0] * g_one / g_all;
      test_force[1] += momentumForce[1] * g_one / g_all;
      test_force[2] += momentumForce[2] * g_one / g_all;

      counting_distribute += 1;
      sestevek += g_one;
    // }
  }

  // sum and distribute
  MPI_Allreduce(&mol_f_ext[0][0], &mol_f_ext_all[0][0], 3 * maxmol_all, MPI_DOUBLE, MPI_SUM, world);
  
  std::cout<<"counting_distribute: "<<counting_distribute<<"\n";
  // output test_force to see if sum is equal to the prescribed external force 
  // testing only equilibrium case (deleteing)
  std::cout<<"FixAdResSObmd::reg_force(): test_force[0]: "<<test_force[0]<<" momentumForce[0]: "<<momentumForce[0]<<"\n";
  std::cout<<"FixAdResSObmd::reg_force(): test_force[1]: "<<test_force[1]<<" momentumForce[1]: "<<momentumForce[1]<<"\n";
  std::cout<<"FixAdResSObmd::reg_force(): test_force[2]: "<<test_force[2]<<" momentumForce[2]: "<<momentumForce[2]<<"\n";

  // distribute external force among particles
  // create test_force_at to check obmd
  double test_force_at[3];
  test_force_at[0] = test_force_at[1] = test_force_at[2] = 0.0;
  for (int i = 0; i < nlocal; i++) {
    imol = molecule[i];
    if (molmap_tmp) imol = molmap_tmp[imol-idlo]; // important ! ! !
    else imol--;

    double molmass = masstotal_tmp[imol];
    double mass_frac = mass[type[i]] / molmass;

    f[i][0] += mass_frac * mol_f_ext_all[imol][0];
    f[i][1] += mass_frac * mol_f_ext_all[imol][1];
    f[i][2] += mass_frac * mol_f_ext_all[imol][2];

    test_force_at[0] += mass_frac * mol_f_ext_all[imol][0];
    test_force_at[1] += mass_frac * mol_f_ext_all[imol][1];
    test_force_at[2] += mass_frac * mol_f_ext_all[imol][2];
  }

  std::cout<<"FixAdResSObmd::reg_force(): test_force_at[0]: "<<test_force_at[0]<<" momentumForce[0]: "<<momentumForce[0]<<"\n";
  std::cout<<"FixAdResSObmd::reg_force(): test_force_at[1]: "<<test_force_at[1]<<" momentumForce[1]: "<<momentumForce[1]<<"\n";
  std::cout<<"FixAdResSObmd::reg_force(): test_force_at[2]: "<<test_force_at[2]<<" momentumForce[2]: "<<momentumForce[2]<<"\n";
  std::cout<<"sestevek: "<<sestevek<<" g_all: "<<g_all<<"\n";

  memory->destroy(massproc_tmp);
  memory->destroy(masstotal_tmp);
  memory->destroy(mol_f_ext);
  memory->destroy(mol_f_ext_all);
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::reg_force_perp(int vflag, Region* region, double *shearForce, int step)
{ 
  std::cout<<"FixAdResSObmd::reg_force_perp"<<"\n";
  if (region) region->prematch();
  v_init(vflag); // ? ? ? peratom virial 
  
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
  double gtmp = g_perp_global(region,step);
  double mass_tmp,gloctmp;
  double unwrap[3];
  double xval, yval, zval;

  // PPapez ... apply on the mol's center-of-mass
  int *rep_atom = atom->rep_atom;
  double **cms = atom->cms_mol;
  double xtmp, ytmp, ztmp;
  
  for (int i = 0; i < nlocal; i++) {
      if (rep_atom[i] == 0) {
        continue;
      }

      if (mask[i] & groupbit) {
        xtmp = cms[i][0];
        ytmp = cms[i][1];
        ztmp = cms[i][2];
      
        if (region && !region->match(xtmp,ytmp,ztmp)) continue;
        // domain->unmap(x[i],image[i],unwrap);
        
        mass_tmp = mass[type[i]];
        gloctmp = 1.0/gtmp;
        if(mode == MOLECULE) {
          mass_tmp = mtot;    //onemols[imol]->masstotal;
          gloctmp = 1.0/gtmp; //tp 24.7.2023  : TODO: cms of molecule not x[i][0] ... does not matter if g_perp = step function
        }
        
        xval = mass[type[i]] * shearForce[0] * gloctmp / mass_tmp;
        yval = mass[type[i]] * shearForce[1] * gloctmp / mass_tmp;
        zval = mass[type[i]] * shearForce[2] * gloctmp / mass_tmp;
        
        foriginal[0] -= xval*unwrap[0] + yval*unwrap[1] + zval*unwrap[2];
        foriginal[1] += f[i][0];
        foriginal[2] += f[i][1];
        foriginal[3] += f[i][2];
        
        f[i][0] += xval; // tp 24.7.2023 
        f[i][1] += yval; // tp 24.7.2023
        f[i][2] += zval; // tp 24.7.2023
        
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
}
}

/* ----------------------------------------------------------------------- */

double FixAdResSObmd::usher(Region *iregion_var, double **coords, double etarget, int natom, int imol, int &iter)
{  
  /* PPapez COMMENT
  USHER for molecule (no rotation)
  omitting imageflags[m] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX \& domain->remap(coords[m],imageflags[m])
  */ 
  double *mass = atom->mass; // to distribute fusher 

  double entmp,torqabs,dtheta,ds,fabs;
  double fusher_tmp[3],fusher_all[3],xcom[3];
  int i,m,tip;
  double entmp_all;
  MathExtra::zero3(fusher_tmp);
  
  i = 0;
  while(i < nattempt) {
    std::cout<<"FixAdResSObmd::usher()"<<"\n";
    std::cout<<"nattempt: "<<i<<"\n";
    std::cout<<"coords[0][0]: "<<coords[0][0]<<" coords[0][1]: "<<coords[0][1]<<" coords[0][2]: "<<coords[0][2]<<"\n";
    std::cout<<"coords[1][0]: "<<coords[1][0]<<" coords[1][1]: "<<coords[1][1]<<" coords[1][2]: "<<coords[1][2]<<"\n";
    std::cout<<"coords[2][0]: "<<coords[2][0]<<" coords[2][1]: "<<coords[2][1]<<" coords[2][2]: "<<coords[2][2]<<"\n"; 

    /* int check = check_mol_region(iregion_var,coords,natom);
    if (check == 1) {
      std::cout<<"in usher brake"<<"\n";
      for(m = 0; m < natom; m++) {     
        std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
      }
    break; // exit
    } */
    
    // PPapez COMMENT: compute center-of-mass
    mol_center_of_mass(natom,imol,coords,xcom); // find mol's cms -> stored in xcom
    std::cout<<"xcom[0]: "<<xcom[0]<<" xcom[1]: "<<xcom[1]<<" xcom[2]: "<<xcom[2]<<"\n";
    entmp = 0.0;
    MathExtra::zero3(fusher);

    for (m = 0; m < natom; m++) {
      // PPapez COMMENT: find representative atom \& its "itype" to compute energy (of a coarse-grained particle in buffer)
      int type_temp = onemols[imol]->type[m];
      int rep_atom_temp = onemols[imol]->rep_atom[m];

      if (rep_atom_temp == 0) { // not representative atom
        continue;
      }

      // PPapez COMMENT: keeping "iregion" \& "iregion2"
      entmp += energy_cg_adress_obmd(iregion,iregion2,type_temp,xcom,fusher); 
    }

    // sum \& distribute
    MPI_Allreduce(fusher,fusher_all,3,MPI_DOUBLE,MPI_SUM,world); // force
    MPI_Allreduce(&entmp,&entmp_all,1,MPI_DOUBLE,MPI_SUM,world); // energy

    // PPapez COMMENT: USHER routine adopted for mode == MOLECULE (carefully with center-of-mass property \& pbc)
    if(entmp_all < etarget+EPSILON) 
    {
      break; // exit because energy is to high
    }
    else if(entmp_all > uovlp) { // corrections \& fun
        fabs = sqrt(fusher_all[0]*fusher_all[0]+fusher_all[1]*fusher_all[1]+fusher_all[2]*fusher_all[2]);
        if(fabs < EPSILON) continue;
        ds = dsovlp -pow(4*eps/entmp_all,1.0/12.0); 
        
        for(m = 0; m < natom; m++) {     
          // std::cout<<"else if(entmp_all > uovlp)"<<"\n";
          // std::cout<<"m: "<<m<<"\n";
          // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
          coords[m][0] += fusher_all[0] / fabs * ds; 
          coords[m][1] += fusher_all[1] / fabs * ds; 
          coords[m][2] += fusher_all[2] / fabs * ds; 
          // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
          
          // PPapez COMMENT: rather checking region and correction position
          // imageflags[m] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
          // domain->remap(coords[m],imageflags[m]);
          // std::cout<<"HERE 1 -> after domain->remap"<<"\n";
          // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
        }
        
        // PPapez COMMENT: check region ... returns 0 / 1 if yes or not
        int check = check_mol_region(iregion_var,coords,natom); 
        if (check == 1) {
          /* std::cout<<"in first check_mol_region"<<"\n";
          for(m = 0; m < natom; m++) {     
            std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
          } */
          break; // exit 
          // PPapez COMMENT: try with new position
          // try_new_position(iregion_var,coords,imol);
        }
    } 
    else {
      fabs = sqrt(fusher_all[0]*fusher_all[0]+fusher_all[1]*fusher_all[1]+fusher_all[2]*fusher_all[2]);
      if(fabs < EPSILON) continue;
      ds = std::min((entmp_all-etarget)/fabs,ds0);

      // std::cout<<"ds: "<<ds<<"\n";
      // std::cout<<"entmp_all: "<<entmp_all<<" etarget: "<<etarget<<" fabs: "<<fabs<<" ds0: "<<ds0<<"\n";

      for(m = 0; m < natom; m++) {
        // std::cout<<"else"<<"\n";
        // std::cout<<"m: "<<m<<"\n";
        // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
        coords[m][0] += fusher_all[0] / fabs * ds; 
        coords[m][1] += fusher_all[1] / fabs * ds; 
        coords[m][2] += fusher_all[2] / fabs * ds; 
        // std::cout<<"fusher_all[0]: "<<fusher_all[0]<<" fusher_all[1]: "<<fusher_all[1]<<" fusher_all[2]: "<<fusher_all[2]<<"\n";
        // std::cout<<"fabs: "<<fabs<<"\n";
        // std::cout<<"ds: "<<ds<<"\n";
        // std::cout<<"fusher_all[0]/fabs: "<<fusher_all[0]/fabs<<"\n";
        // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";

        // PPapez COMMENT: rather checking region and correction position
        // imageflags[m] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        // domain->remap(coords[m],imageflags[m]);
        // std::cout<<"HERE 2 -> after domain->remap"<<"\n";
        // std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
      } 

      // PPapez COMMENT: check region ... returns 0 / 1 if yes or not
      int check = check_mol_region(iregion_var,coords,natom);
      if (check == 1) {
        /* std::cout<<"in second check_mol_region"<<"\n";
        for(m = 0; m < natom; m++) {     
          std::cout<<"coords[m][0]: "<<coords[m][0]<<" coords[m][1]: "<<coords[m][1]<<" coords[m][2]: "<<coords[m][2]<<"\n";
        } */
        break; // exit
        // PPapez COMMENT: try with new position
        // try_new_position(iregion_var,coords,imol);
      }  
    }
    i+=1;
  }
  iter = i;
  entmp = entmp_all;
  
  return entmp;
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::try_new_position(Region *iregion_var, double **coords, int imol)
{ 
  /*
  PPapez COMMENT
  try to find new position if the previous one is outside "iregion_var"
  adopted for mode == MOLECULE
  checking positions
  */ 

  int i,m,n,natom;
  double coord[3];
  double r[3],rotmat[3][3],quat[4];
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

  std::cout<<"FixAdResSObmd::try_new_position"<<"\n";
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

double FixAdResSObmd::near_energy(double **coords, int natom, int imol)
{
  double entmp,entmp_all;
  double fusher[3];
  int m,tip;
  entmp = 0.0;
  tip = ntype;
  for (m = 0; m < natom; m++) {
    if(mode == MOLECULE) tip = ntype+onemols[imol]->type[m];
    entmp += energy(1, tip, coords[m], fusher); // i - does not matter, itype inserted, coords, fusher
  }
  MPI_Allreduce(&entmp,&entmp_all,1,MPI_DOUBLE,MPI_SUM,world); 
  return entmp_all;
}

/* ---------------------------------------------------------------------- */

int FixAdResSObmd::check_proc(double **coords, double *newcoord, double *sublo, double*subhi, int dimension)
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

int FixAdResSObmd::check_mol_proc(double **coords, double *newcoord, double *lamda, double *sublo, double*subhi, int dimension, int natom, int &me)
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

int FixAdResSObmd::check_mol_region(Region* region, double **coords, int natom)
{
  int m,flag,flagall;
  flag = 0;

  for (m = 0; m < natom; m++) {
    if (region && !region->match(coords[m][0],coords[m][1],coords[m][2])) flag = 1;
  }

  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_MAX,world);

  if(comm->me==0 && flagall==1) {
    std::cout<<"USHER moved the particle/molecule too much"<<"\n";
  }
  
  return flagall;
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::mol_center_of_mass(int natom, int imol, double **coords, double *xcm)
{
  double *mass = atom->mass;
  int m;
  double massmol = 0.0;
  MathExtra::zero3(xcm);
  
  for(m = 0; m < natom; m++) {
    xcm[0] += mass[onemols[imol]->type[m]] * coords[m][0];
    xcm[1] += mass[onemols[imol]->type[m]] * coords[m][1];
    xcm[2] += mass[onemols[imol]->type[m]] * coords[m][2];

    massmol += mass[onemols[imol]->type[m]];

  }
  // std::cout<<"massmol: "<<massmol<<"\n";

  xcm[0] = xcm[0] / massmol;
  xcm[1] = xcm[1] / massmol;
  xcm[2] = xcm[2] / massmol;
  // std::cout<<"xcm[0]: "<<xcm[0]<<" xcm[1]: "<<xcm[1]<<" xcm[2]: "<<xcm[2]<<"\n";
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::center_of_mass(int natom, double **coords, double *xcm)
{
  int m;
  MathExtra::zero3(xcm);
  
  for(m = 0; m < natom; m++) {
    xcm[0] += coords[m][0];
    xcm[1] += coords[m][1]; // ERR
    xcm[2] += coords[m][1];
  }
  
  xcm[0] = xcm[0] / natom;
  xcm[1] = xcm[1] / natom;
  xcm[2] = xcm[2] / natom;
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::calc_torque(int natom, double **coords, double *xcm, double *force_tmp, double *torq)
{
  int m;
  double xrel[3];
  
  for(m = 0; m < natom; m++) {
    xrel[0] = coords[m][0] - xcm[0];
    xrel[1] = coords[m][1] - xcm[1];
    xrel[2] = coords[m][2] - xcm[2];
  }   
  
  MathExtra::cross3(xrel,force_tmp,torq);  
}

/* ---------------------------------------------------------------------- */

double FixAdResSObmd::energy_cg_adress_obmd(Region *leftB, Region *rightB, int itype, double *xcom, double *fusher)
{  
  std::cout<<"energy_cg_adress_obmd"<<"\n";
  std::cout<<"itype: "<<itype<<"\n";
  // std::cout<<"xcom[0]: "<<xcom[0]<<" xcom[1]: "<<xcom[1]<<" xcom[2]: "<<xcom[2]<<"\n";

  int *rep_atom = atom->rep_atom; // get cms using rep_atom of mol
  double **cms = atom->cms_mol; // cms of existing mols
  int nlocal = atom->nlocal; // over nlocal 
  int *type = atom->type;
  // int *mask = atom->mask;
  // double *special_lj = force->special_lj;
  cutsq = force->pair->cutsq;
  auto pair = dynamic_cast<PairTableAdResSOBMDCG*>(force->pair_match("table/adress/obmd/cg",1)); // tabulated potential given by user

  int j;
  int jtype;
  double xtmp, ytmp, ztmp;
  double delx, dely, delz; 
  double rsq, dist; 
  double box_center_x;
  box_center_x = (domain->boxlo[0] + domain->boxhi[0])/2;

  double total_energy = 0.0;
  double fpair = 0.0;
  double factor_lj = 1.0;
  int count_lb, count_rb, count_roi;
  count_lb = count_rb = count_roi = 0;

  for (j = 0; j < nlocal; j++) {
    if (rep_atom[j] == 0) {
      continue;
    }

    jtype = type[j];
    // std::cout<<"jtype: "<<jtype<<"\n";
   
    xtmp = cms[j][0];
    ytmp = cms[j][1];
    ztmp = cms[j][2];

    if (xcom[0] < box_center_x && leftB->match(xtmp,ytmp,ztmp)) { // in which region i of nlocal falls ? ? ? match to region
      // std::cout<<"LEFT BUFFER"<<"\n";
      // std::cout<<"xcom[0] < box_center_x && leftB->match(xtmp,ytmp,ztmp)"<<"\n";
      // std::cout<<"xcm[0]: "<<xcom[0]<<" xcom[1]: "<<xcom[1]<<" xcom[2]: "<<xcom[2]<<"\n";
      // std::cout<<"xtmp: "<<xtmp<<" ytmp: "<<ytmp<<" ztmp: "<<ztmp<<"\n";
      
      // compute energy for left buffer
      count_lb += 1;
      // distance between particles
      delx = xcom[0] - xtmp;
      dely = xcom[1] - ytmp;
      delz = xcom[2] - ztmp;
      domain->minimum_image(delx,dely,delz);
      rsq = delx * delx + dely * dely + delz * delz;

      // std::cout<<"rsq: "<<rsq<<" cutsq[itype][jtype]: "<<cutsq[itype][jtype]<<" itype: "<<itype<<" jtype: "<<jtype<<"\n";

      if (rsq < cutsq[itype][jtype]) {
        
        total_energy += pair->single(0,j,itype,jtype,rsq,0,factor_lj,fpair); 

        if (pair->single(0,j,itype,jtype,rsq,0,factor_lj,fpair) == pow(10,6)) {
          std::cout<<pair->single(0,j,itype,jtype,rsq,0,factor_lj,fpair)<<"\n";
          std::cout<<"fpair: "<<fpair<<"\n";
        }
        
          // force on mol
        fusher[0] += fpair * delx; // mol
        fusher[1] += fpair * dely; // mol
        fusher[2] += fpair * delz; // mol
        // single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fforce) ... setting i in factor_coul to 0 (not needed) ! ! ! 
      }

      
      // std::cout<<"total_energy: "<<total_energy<<"\n";
      // std::cout<<"count_lb: "<<count_lb<<"\n";

    } else if (xcom[0] > box_center_x && rightB->match(xtmp,ytmp,ztmp)) { // in which region i of nlocal falls ? ? ? match to region
      // std::cout<<"RIGHT BUFFER"<<"\n";
      // std::cout<<"xcom[0] < box_center_x && rightB->match(xtmp,ytmp,ztmp)"<<"\n";
      // std::cout<<"xcm[0]: "<<xcom[0]<<" xcom[1]: "<<xcom[1]<<" xcom[2]: "<<xcom[2]<<"\n";
      // std::cout<<"xtmp: "<<xtmp<<" ytmp: "<<ytmp<<" ztmp: "<<ztmp<<"\n";
      
      // compute energy for left buffer
      count_rb += 1;
      // distance between particles
      delx = xcom[0] - xtmp;
      dely = xcom[1] - ytmp;
      delz = xcom[2] - ztmp;
      // domain->minimum_image(delx,dely,delz);
      rsq = delx * delx + dely * dely + delz * delz;

      // std::cout<<"rsq: "<<rsq<<" cutsq[itype][jtype]: "<<cutsq[itype][jtype]<<" itype: "<<itype<<" jtype: "<<jtype<<"\n";

      if (rsq < cutsq[itype][jtype]) {
        total_energy += pair->single(0,j,itype,jtype,rsq,0,factor_lj,fpair); 
        // force on mol
        fusher[0] += fpair * delx; // mol
        fusher[1] += fpair * dely; // mol
        fusher[2] += fpair * delz; // mol
        // single(i,j,itype,jtype,rsq,factor_coul,factor_lj,fforce) ... setting i in factor_coul to 0 (not needed) ! ! ! 
      }

      // std::cout<<"total_energy: "<<total_energy<<"\n";
      // std::cout<<"count_rb: "<<count_lb<<"\n";

    } else {
      // std::cout<<"ROI region and not the buffer of interest"<<"\n";
      // std::cout<<"xcm[0]: "<<xcom[0]<<" xcom[1]: "<<xcom[1]<<" xcom[2]: "<<xcom[2]<<"\n";
      // std::cout<<"xtmp: "<<xtmp<<" ytmp: "<<ytmp<<" ztmp: "<<ztmp<<"\n";
      count_roi += 1;

      // std::cout<<"total_energy: "<<total_energy<<"\n";
      // std::cout<<"count_roi: "<<count_roi<<"\n";
      continue;

    } 
  }

  std::cout<<"count_lb: "<<count_lb<<" count_rb: "<<count_rb<<" count_roi: "<<count_roi<<"\n";
  std::cout<<"total_energy: "<<total_energy<<"\n";
  return total_energy;

}

/* ---------------------------------------------------------------------- */


/* double FixAdResSObmd::energy_adress_obmd(double qi, int itype, double *coord, double *fusher)
{  
  // std::cout<<"in energy_adress_obmd"<<"\n";
  // std::cout<<"qi: "<<qi<<" itype: "<<itype<<" coord[0]: "<<coord[0]<<" coord[1]: "<<coord[1]<<" coord[2]: "<<coord[2]<<"\n";

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

  for (int j = 0; j < nall; j++) { //nall is probably not needed, could just be nlocal
    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];
    domain->minimum_image(delx,dely,delz); //added tp 20.10.2023
    rsq = delx*delx + dely*dely + delz*delz;
    jtype = type[j];

    if (rsq < cutsq[itype][jtype]) {
            // avoiding table/adress/obmd/cg
      auto pair = dynamic_cast<PairLJCutRFAdResSAT*>(force->pair_match("lj/cut/rf/adress/at",1));
      total_energy += pair->single_adress_obmd(qi,j,itype,jtype,rsq,factor_coul,factor_lj,fpair);

      fusher[0] += fpair * delx;
      fusher[1] += fpair * dely;
      fusher[2] += fpair * delz;
	  }
  }

  // std::cout<<"total_energy: "<<total_energy<<"\n";
  return total_energy;
} */
/* ---------------------------------------------------------------------- */

double FixAdResSObmd::energy(int i, int itype, double *coord, double *fusher)
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
      // std::cout<<"atom->tag[j]: "<<atom->tag[j]<<"\n";

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
      auto pair = dynamic_cast<PairLJCutRFAdResSAT*>(force->pair_match("lj/cut/rf/adress/at",1));
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
}

/* ----------------------------------------------------------------------
   compute the center-of-mass velocity of group of atoms in region
   masstotal = total mass
   return center-of-mass velocity in cm[]
------------------------------------------------------------------------- */

void FixAdResSObmd::vcm_internal_sq(int igroup, double masstotal, double *cm, double *temp, double *vcm, Region *region)
{
  int groupbit = group->bitmask[igroup];
  region->prematch();

  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double p, q;
  p = q = 0.0;
  
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit && region->match(x[i][0],x[i][1],x[i][2])) {
      p += (v[i][0] - vcm[0]) * (v[i][0] - vcm[0]);
      p += (v[i][1] - vcm[1]) * (v[i][1] - vcm[1]);
      p += (v[i][2] - vcm[2]) * (v[i][2] - vcm[2]);
      q += mass[type[i]] * (v[i][0] - vcm[0]) * (v[i][0] - vcm[0]);
      q += mass[type[i]] * (v[i][1] - vcm[1]) * (v[i][1] - vcm[1]);
      q += mass[type[i]] * (v[i][2] - vcm[2]) * (v[i][2] - vcm[2]);   
    }

  MPI_Allreduce(&p,cm,1,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&q,temp,1,MPI_DOUBLE,MPI_SUM,world);
  
  *temp = *temp / 3 / masstotal;
  
}

/* ---------------------------------------------------------------------- */

double FixAdResSObmd::compute_scalar()
{
  return 1.0;
}

/* ---------------------------------------------------------------------- */

double FixAdResSObmd::compute_vector(int n)
{
  return 1.0;
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::find_maxid()
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

/* ---------------------------------------------------------------------- */

// find molecules
int FixAdResSObmd::find_mols(tagint &idlo, tagint &idhi)
{
  int i; 
  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal; 
   
  tagint lo = BIG;
  tagint hi = -BIG;

  for (i = 0; i < nlocal; i++) {
    if (mask[i]) {
      lo = MIN(lo,molecule[i]);
      hi = MAX(hi,molecule[i]); 
    }
  }

  MPI_Allreduce(&lo,&idlo,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_LMP_TAGINT,MPI_MAX,world);
  
  if (idlo == BIG) return 0;
  
  int nlen = (int) idhi-idlo+1;
  
  memory->create(molmap_tmp,nlen,"tdforce:molmap");
  for (i = 0; i < nlen; i++) molmap_tmp[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i])
      molmap_tmp[molecule[i]-idlo] = 1;

  int *molmapall;
  memory->create(molmapall,nlen,"tdforce:molmapall");
  MPI_Allreduce(molmap_tmp,molmapall,nlen,MPI_INT,MPI_MAX,world);

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap_tmp[i] = nmolecules++;
    else molmap_tmp[i] = -1;
    
  memory->destroy(molmapall);  
  memory->destroy(molmap_tmp);  
  
  return nmolecules;
} 
// PP using again

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixAdResSObmd::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixAdResSObmd::options(int narg, char **arg)
{
  // defaults
  nearflag = 0;
  usherflag = 0;
  iregion = nullptr;
  idregion = nullptr;
  iregion2 = nullptr;
  idregion2 = nullptr;
  iregion3 = nullptr;
  idregion3 = nullptr;
  iregion4 = nullptr;
  idregion4 = nullptr;  
  iregion5 = nullptr;
  idregion5 = nullptr;
  iregion6 = nullptr;
  idregion6 = nullptr;  
  mode = ATOM;
  molfrac = nullptr;
  rigidflag = 0;
  idrigid = nullptr;
  shakeflag = 0;
  idshake = nullptr;
  idnext = 0;
  globalflag = localflag = 0;
  lo = hi = deltasq = 0.0;
  nearsq = 0.0;
  maxattempt = 10;
  rateflag = 0;
  vxlo = vxhi = vylo = vyhi = vzlo = vzhi = 0.0;
  distflag = DIST_UNIFORM;
  sigma = 1.0;
  xmid = ymid = zmid = 0.0;
  scaleflag = 1;
  targetflag = 0;
  orientflag = 0;
  rx = 0.0;
  ry = 0.0;
  rz = 0.0;
  print_insert_flag = 0;
  print_delete_flag = 0;
  
  g_fac = 0.25;
  g_fac_inv = 1.0/g_fac;
  mol_len=1.0;
  alpha = 0.75;
  tau = 100.0;
  buffer_size = 0.2*(domain->boxhi[0]-domain->boxlo[0]);
  shear_size = 0.125*(domain->boxhi[0]-domain->boxlo[0]);
  nbuf = group->count(igroup)*buffer_size/(domain->boxhi[0]-domain->boxlo[0]);
  etarget = 100.0;
  uovlp = 5000;
  dsovlp = 1.0;
  eps = 1.0;
  ds0 = 0.1;
  dtheta0 = 0.1;
  nattempt = 20;
  
  stepflag = 1;
  molflag = 0;
  int imol = -1;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"region1") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      iregion = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion)
        error->all(FLERR,"Region ID for fix obmd does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion = new char[n];
      strcpy(idregion,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"region2") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      iregion2 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion2)
        error->all(FLERR,"Region ID for fix obmd does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion2 = new char[n];
      strcpy(idregion2,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"region3") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      iregion3 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion3)
        error->all(FLERR,"Region ID for fix obmd does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion3 = new char[n];
      strcpy(idregion3,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"region4") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      iregion4 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion4)
        error->all(FLERR,"Region ID for fix obmd does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion4 = new char[n];
      strcpy(idregion4,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"region5") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      iregion5 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion5)
        error->all(FLERR,"Region ID for fix obmd does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion5 = new char[n];
      strcpy(idregion5,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"region6") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      iregion6 = domain->get_region_by_id(arg[iarg+1]);
      if (!iregion6)
        error->all(FLERR,"Region ID for fix obmd does not exist");
      int n = strlen(arg[iarg+1]) + 1;
      idregion6 = new char[n];
      strcpy(idregion6,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"usher") == 0) {
      if (iarg+9 > narg) error->all(FLERR,"Illegal fix obmd command");
      usherflag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      if (usherflag ==1 && nearflag==1)
        error->all(FLERR,"You can not have both usher and near");
      etarget = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      ds0 = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      dtheta0 = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      uovlp = utils::numeric(FLERR,arg[iarg+5],false,lmp);
      dsovlp = utils::numeric(FLERR,arg[iarg+6],false,lmp);
      eps = utils::numeric(FLERR,arg[iarg+7],false,lmp);
      nattempt = utils::inumeric(FLERR,arg[iarg+8],false,lmp);
      iarg += 9;
    } else if (strcmp(arg[iarg],"printinsert") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      print_insert_flag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"printdelete") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      print_delete_flag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"alpha") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      alpha = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (alpha < 0.0)
        error->all(FLERR,"alpha should be bigger than 0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"gfac") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      g_fac = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      g_fac_inv = 1.0/g_fac;
      if (g_fac <= 0.0 || g_fac >= 1.0)
        error->all(FLERR,"gfac should be between 0 and 1");
      iarg += 2;
    } else if (strcmp(arg[iarg],"tau") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      tau = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tau < 0.0)
        error->all(FLERR,"tau should be bigger than 0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nbuf") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      nbuf = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (nbuf < 0.0)
        error->all(FLERR,"nbuf should be bigger than 0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"step") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix evaporate command");
      stepflag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"buffersize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      buffer_size = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (buffer_size < 0.0)
        error->all(FLERR,"buffer size should be bigger than 0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"shearsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      shear_size = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (shear_size < 0.0)
        error->all(FLERR,"shear size should be bigger than 0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"mol") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix obmd command");//
      imol = atom->find_molecule(arg[iarg+1]);
      mol_len = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      molflag = 1;
      if (imol == -1)
        error->all(FLERR,"Molecule template ID for fix obmd does not exist");
      mode = MOLECULE;
      onemols = &atom->molecules[imol];
      onemols[imol]->compute_mass();
      mtot = onemols[imol]->masstotal;
      // std::cout<<"in mol ::: mtot: "<<mtot<<"\n";
      nmol = onemols[0]->nset;
      delete [] molfrac;
      molfrac = new double[nmol];
      molfrac[0] = 1.0/nmol;
      for (int i = 1; i < nmol-1; i++) molfrac[i] = molfrac[i-1] + 1.0/nmol;
      molfrac[nmol-1] = 1.0;
      iarg += 3;//
    } else if (strcmp(arg[iarg],"molfrac") == 0) {
      if (mode != MOLECULE) error->all(FLERR,"Illegal fix obmd command");
      if (iarg+nmol+1 > narg) error->all(FLERR,"Illegal fix obmd command");
      molfrac[0] = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      for (int i = 1; i < nmol; i++)
        molfrac[i] = molfrac[i-1] + utils::numeric(FLERR,arg[iarg+i+1],false,lmp);
      if (molfrac[nmol-1] < 1.0-EPSILON || molfrac[nmol-1] > 1.0+EPSILON)
        error->all(FLERR,"Illegal fix obmd command");
      molfrac[nmol-1] = 1.0;
      iarg += nmol+1;
    } else if (strcmp(arg[iarg],"rigid") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idrigid;
      idrigid = new char[n];
      strcpy(idrigid,arg[iarg+1]);
      rigidflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"shake") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      int n = strlen(arg[iarg+1]) + 1;
      delete [] idshake;
      idshake = new char[n];
      strcpy(idshake,arg[iarg+1]);
      shakeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"id") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      if (strcmp(arg[iarg+1],"max") == 0) idnext = 0;
      else if (strcmp(arg[iarg+1],"next") == 0) idnext = 1;
      else error->all(FLERR,"Illegal fix obmd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"global") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix obmd command");
      globalflag = 1;
      localflag = 0;
      lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg],"local") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix obmd command");
      localflag = 1;
      globalflag = 0;
      lo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      hi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      deltasq = utils::numeric(FLERR,arg[iarg+3],false,lmp) *
        utils::numeric(FLERR,arg[iarg+3],false,lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg],"near") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix obmd command");
      nearflag = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      nearsq = utils::numeric(FLERR,arg[iarg+2],false,lmp) *
        utils::numeric(FLERR,arg[iarg+2],false,lmp);
      if (usherflag ==1 && nearflag==1)
        error->all(FLERR,"You can not have both usher and near");      
      iarg += 3;
    } else if (strcmp(arg[iarg],"attempt") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      maxattempt = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      rateflag = 1;
      rate = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg],"vx") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix obmd command");
      vxlo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      vxhi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vy") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix obmd command");
      vylo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      vyhi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg],"vz") == 0) {
      if (iarg+3 > narg) error->all(FLERR,"Illegal fix obmd command");
      vzlo = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      vzhi = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg],"orient") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix obmd command");
      orientflag = 1;
      rx = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ry = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      rz = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      if (domain->dimension == 2 && (rx != 0.0 || ry != 0.0))
        error->all(FLERR,"Illegal fix obmd orient settings");
      if (rx == 0.0 && ry == 0.0 && rz == 0.0)
        error->all(FLERR,"Illegal fix obmd orient settings");
      iarg += 4;
    } else if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix obmd command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all(FLERR,"Illegal fix obmd command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"gaussian") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal fix obmd command");
      xmid = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ymid = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      zmid = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      sigma = utils::numeric(FLERR,arg[iarg+4],false,lmp);
      distflag = DIST_GAUSSIAN;
      iarg += 5;
    } else if (strcmp(arg[iarg],"target") == 0) {
      if (iarg+4 > narg) error->all(FLERR,"Illegal fix obmd command");
      tx = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      ty = utils::numeric(FLERR,arg[iarg+2],false,lmp);
      tz = utils::numeric(FLERR,arg[iarg+3],false,lmp);
      targetflag = 1;
      iarg += 4;
    } else error->all(FLERR,"Illegal fix obmd command");
  }
  if (usherflag == 0 && nearflag == 0)
        error->all(FLERR,"You must pick either usher or near");
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixAdResSObmd::write_restart(FILE *fp)
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
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixAdResSObmd::restart(char *buf)
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
}

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

void *FixAdResSObmd::extract(const char *str, int &itype)
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
}

