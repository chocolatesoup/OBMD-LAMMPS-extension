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
#include "random_park.h"
#include "region.h"
#include "update.h"
#include "variable.h"

#include "pair_lj_cut_rf.h"

#include <cmath>
#include <cstring>
#include <iostream>

#include <fenv.h>    // used for debug

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum { ATOM, MOLECULE };
enum { DIST_UNIFORM, DIST_GAUSSIAN };
enum { NONE, CONSTANT, EQUAL };
enum { EXCHATOM, EXCHMOL };    // exchmode
enum { MOVEATOM, MOVEMOL };    // movemode

#define EPSILON 1.0e-6

/* ---------------------------------------------------------------------- */

FixObmdMerged::FixObmdMerged(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), idregion(nullptr), idregion2(nullptr), idregion3(nullptr),
    idregion4(nullptr), idregion5(nullptr), idregion6(nullptr), idrigid(nullptr), idshake(nullptr),
    onemols(nullptr), molfrac(nullptr), coords(nullptr), imageflags(nullptr), fixrigid(nullptr),
    fixshake(nullptr), random(nullptr), list(nullptr), mark(nullptr)
{
  if (narg < 14) error->all(FLERR, "Illegal fix obmd command (check # of arguments)");

  restart_global = 1;
  time_depend = 1;
  size_vector = 1;

  pxxstr = nullptr;
  pxystr = nullptr;
  pxzstr = nullptr;
  dpstr = nullptr;       // pressure amplitude
  freqstr = nullptr;     // frequency
  alphastr = nullptr;    // alpha for feedback algorithm
  taustr = nullptr;      // tau for feedback algorithm
  nbufstr = nullptr;     // desired number of particles in the buffer

  // other args
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

  // read options from end of input line
  varflag = CONSTANT;
  options(narg - 14, &arg[14]);

  // error check on type
  if (mode == ATOM && (ntype <= 0 || ntype > atom->ntypes))
    error->all(FLERR, "Invalid atom type in fix obmd command");

  // error checks on region and its extent being inside simulation box
  if (!iregion || !iregion2) error->all(FLERR, "Must specify a region in fix obmd");
  if (iregion->bboxflag == 0 || iregion2->bboxflag == 0)
    error->all(FLERR, "Fix obmd region does not support a bounding box");
  if (iregion->dynamic_check() || iregion2->dynamic_check())
    error->all(FLERR, "Fix obmd region cannot be dynamic");

  xlo = iregion->extent_xlo;
  xhi = iregion->extent_xhi;
  ylo = iregion->extent_ylo;
  yhi = iregion->extent_yhi;
  zlo = iregion->extent_zlo;
  zhi = iregion->extent_zhi;

  if (domain->triclinic == 0) {
    if (xlo < domain->boxlo[0] || xhi > domain->boxhi[0] || ylo < domain->boxlo[1] ||
        yhi > domain->boxhi[1] || zlo < domain->boxlo[2] || zhi > domain->boxhi[2])
      error->all(FLERR, "Insertion region extends outside simulation box");
  } else {
    if (xlo < domain->boxlo_bound[0] || xhi > domain->boxhi_bound[0] ||
        ylo < domain->boxlo_bound[1] || yhi > domain->boxhi_bound[1] ||
        zlo < domain->boxlo_bound[2] || zhi > domain->boxhi_bound[2])
      error->all(FLERR, "Insertion region extends outside simulation box");
  }

  // error check and further setup for mode = MOLECULE
  if (atom->tag_enable == 0) error->all(FLERR, "Cannot use fix obmd unless atoms have IDs");

  if (mode == MOLECULE) {
    for (int i = 0; i < nmol; i++) {
      if (onemols[i]->xflag == 0) error->all(FLERR, "Fix obmd molecule must have coordinates");
      if (onemols[i]->typeflag == 0) error->all(FLERR, "Fix obmd molecule must have atom types");
      if (ntype + onemols[i]->ntypes <= 0 || ntype + onemols[i]->ntypes > atom->ntypes)
        error->all(FLERR, "Invalid atom type in fix obmd mol command");

      if (atom->molecular == Atom::TEMPLATE && onemols != atom->avec->onemols)
        error->all(FLERR,
                   "Fix obmd molecule template ID must be same "
                   "as atom_style template ID");
      onemols[i]->check_attributes();

      // fix obmd uses geoemetric center of molecule for insertion
      onemols[i]->compute_center();
    }
  }

  if (rigidflag && mode == ATOM) error->all(FLERR, "Cannot use fix obmd rigid and not molecule");
  if (shakeflag && mode == ATOM) error->all(FLERR, "Cannot use fix obmd shake and not molecule");
  if (rigidflag && shakeflag) error->all(FLERR, "Cannot use fix obmd rigid and shake");

  // setup of coords and imageflags array
  if (mode == ATOM)
    natom_max = 1;
  else {
    natom_max = 0;
    for (int i = 0; i < nmol; i++) natom_max = MAX(natom_max, onemols[i]->natoms);
  }
  memory->create(coords, natom_max, 3, "obmd:coords");
  memory->create(imageflags, natom_max, "obmd:imageflags");

  // setup scaling
  double xscale, yscale, zscale;
  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  } else
    xscale = yscale = zscale = 1.0;

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
  deltasq *= xscale * xscale;
  nearsq *= xscale * xscale;
  vxlo *= xscale;
  vxhi *= xscale;
  vylo *= yscale;
  vyhi *= yscale;
  vzlo *= zscale;
  vzhi *= zscale;
  xmid *= xscale;
  ymid *= yscale;
  zmid *= zscale;
  sigma *= xscale;    // same as in region sphere
  tx *= xscale;
  ty *= yscale;
  tz *= zscale;

  // find current max atom and molecule IDs if necessary
  if (idnext) find_maxid();

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

  delete_left_file.open("./../out/data/delete_left.out");
  delete_right_file.open("./../out/data/delete_right.out");
  insert_left_file.open("./../out/data/insert_left.out");
  insert_right_file.open("./../out/data/insert_right.out");
}

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

  delete_left_file.close();
  delete_right_file.close();
  insert_left_file.close();
  insert_right_file.close();
}

/* ---------------------------------------------------------------------- */

int FixObmdMerged::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= POST_FORCE;
  return mask;
}

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

  if (pxxstyle == ATOM || pxystyle == ATOM || pxzstyle == ATOM || dpstyle == ATOM ||
      freqstyle == ATOM || alphastyle == ATOM || taustyle == ATOM || nbufstyle == ATOM)
    varflag = ATOM;
  else if (pxxstyle == EQUAL || pxystyle == EQUAL || pxzstyle == EQUAL || dpstyle == EQUAL ||
           freqstyle == EQUAL || alphastyle == EQUAL || taustyle == EQUAL || nbufstyle == EQUAL)
    varflag = EQUAL;
  else
    varflag = CONSTANT;

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
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);

    if (flagall) error->all(FLERR, "Cannot delete atoms in atom_modify first group");
  }

  // if MOLECULE not set, warn if any deletable atom has a mol ID
  if (mode != MOLECULE && atom->molecule_flag) {
    tagint *molecule = atom->molecule;
    int *mask = atom->mask;
    int nlocal = atom->nlocal;
    int flag = 0;
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit)
        if (molecule[i]) flag = 1;
    int flagall;
    MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, world);
    if (flagall && comm->me == 0)
      error->warning(FLERR, "Fix obmd may delete atom with non-zero molecule ID");
  }

  if (mode == MOLECULE && atom->molecule_flag == 0)
    error->all(FLERR, "Fix obmd molecule requires atom attribute molecule");

  // if rigidflag defined, check for rigid/small fix
  // its molecule template must be same as this one
  fixrigid = nullptr;
  if (rigidflag) {
    int ifix = modify->find_fix(idrigid);
    if (ifix < 0) error->all(FLERR, "Fix obmd rigid fix does not exist");
    fixrigid = modify->fix[ifix];
    int tmp;
    if (onemols != (Molecule **) fixrigid->extract("onemol", tmp))
      error->all(FLERR,
                 "Fix obmd and fix rigid/small not using "
                 "same molecule template ID");
  }

  // if shakeflag defined, check for SHAKE fix
  // its molecule template must be same as this one
  fixshake = nullptr;
  if (shakeflag) {
    fixshake = modify->get_fix_by_id(idshake);
    if (!fixshake) error->all(FLERR, "Fix deposit shake fix ID {} does not exist", idshake);

    int tmp;
    if (onemols != (Molecule **) fixshake->extract("onemol", tmp))
      error->all(FLERR,
                 "Fix obmd and fix shake not using "
                 "same molecule template ID");
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
    for (int i = 0; i < nlocal; i++) maxrad = MAX(maxrad, radius[i]);

    double maxradall;
    MPI_Allreduce(&maxrad, &maxradall, 1, MPI_DOUBLE, MPI_MAX, world);

    double maxradinsert = 0.0;
    if (mode == MOLECULE) {
      for (int i = 0; i < nmol; i++) {
        if (onemols[i]->radiusflag)
          maxradinsert = MAX(maxradinsert, onemols[i]->maxradius);
        else
          maxradinsert = MAX(maxradinsert, 0.5);
      }
    } else
      maxradinsert = 0.5;

    double separation = MAX(2.0 * maxradinsert, maxradall + maxradinsert);
    if (sqrt(nearsq) < separation && comm->me == 0)
      error->warning(FLERR,
                     fmt::format("Fix obmd near setting < possible "
                                 "overlap separation {}",
                                 separation));
  }
}

/* ----------------------------------------------------------------------
   setup
------------------------------------------------------------------------- */

void FixObmdMerged::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet")) post_force(vflag);
}

/* ----------------------------------------------------------------------
   perform particle/molecule insertions/deletions
------------------------------------------------------------------------- */
void FixObmdMerged::pre_exchange()
{
  std::cout << "pre_exchange" << "\n";
  int cnt_left, cnt_right, stev_left, stev_right, i;
  double ninsert_left, ninsert_right;
  double masstotal_left, masstotal_right;

  // calculate box dimensions
  double lx = domain->boxhi[0] - domain->boxlo[0];
  double ly = domain->boxhi[1] - domain->boxlo[1];
  double lz = domain->boxhi[2] - domain->boxlo[2];

  // fetch the curret value of momentum/energy flux if it is not CONSTANT
  if (varflag != CONSTANT) {
    if (pxxstyle == EQUAL) pxx = input->variable->compute_equal(pxxvar);
    if (pxystyle == EQUAL) pxy = input->variable->compute_equal(pxyvar);
    if (pxzstyle == EQUAL) pxz = input->variable->compute_equal(pxzvar);
    if (dpstyle == EQUAL) dpxx = input->variable->compute_equal(dpvar);
    if (freqstyle == EQUAL) freq = input->variable->compute_equal(freqvar);
    if (alphastyle == EQUAL) alpha = input->variable->compute_equal(alphavar);
    if (taustyle == EQUAL) tau = input->variable->compute_equal(tauvar);
    if (nbufstyle == EQUAL) nbuf = input->variable->compute_equal(nbufvar);
  }

  MathExtra::zero3(vnewl);
  MathExtra::zero3(vnewr);

  // delete particles that cross the open boundaries after first half of velocity-Verlet algorithm
  try_deleting(iregion, vnewl, vnewr);
  try_deleting(iregion2, vnewl, vnewr);

  // count number of particles in left/right buffer
  cnt_left = group->count(igroup, iregion);
  cnt_right = group->count(igroup, iregion2);

  // calculate number of particles needed for insertion
  ninsert_left = -static_cast<int>((static_cast<double>(cnt_left) / mol_len - alpha * nbuf) *
                                   update->dt / tau);
  ninsert_right = -static_cast<int>((static_cast<double>(cnt_right) / mol_len - alpha * nbuf) *
                                    update->dt / tau);

  // tries inserting ninsert_left/right particles into left/right buffer
  try_inserting(iregion5, ninsert_left, vnewl, vnewr);
  try_inserting(iregion6, ninsert_right, vnewl, vnewr);

  // deletes again, just in case part of newly inserted molecule is not inside box
  try_deleting(iregion, vnewl, vnewr);
  try_deleting(iregion2, vnewl, vnewr);

  // calculates total mass in left/right buffer
  masstotal_left = group->mass(igroup, iregion);
  masstotal_right = group->mass(igroup, iregion2);

  // sum \& distribute
  MPI_Allreduce(vnewl, vnewl_all, 3, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(vnewr, vnewr_all, 3, MPI_DOUBLE, MPI_SUM, world);

  // area of the buffer - ROI interface
  double area = ly * lz;

  simulation_time += update->dt;
  double factor = pxx + dpxx * sin(2.0 * MY_PI * freq * simulation_time);

  // calculates momentum forces on the left buffer
  // substracting momentum of deleted particles
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

  next_reneighbor += nfreq;
}

/* ----------------------------------------------------------------------
   maxtag_all = current max atom ID for all atoms
   maxmol_all = current max molecule ID for all atoms
------------------------------------------------------------------------- */
void FixObmdMerged::try_deleting(Region *region, double *vnewl, double *vnewr)
{

  int i, j, m, iwhichglobal, iwhichlocal;
  int ndel, ndeltopo[4];

  // grow list and mark arrays if necessary
  if (atom->nmax > nmax) {
    memory->destroy(list);
    memory->destroy(mark);
    nmax = atom->nmax;
    memory->create(list, nmax, "obmd:list");
    memory->create(mark, nmax, "obmd:mark");
  }

  region->prematch();

  double **x = atom->x;
  double **vel = atom->v;
  int *mask = atom->mask;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;
  double boxl = domain->boxlo[0];
  double boxh = domain->boxhi[0];
  double *mass = atom->mass;
  int *type = atom->type;
  double fusher[3];

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

  ndel = 0;
  for (i = 0; i < nlocal; i++) mark[i] = 0;

  // atomic deletions
  // choose atoms randomly across all procs and mark them for deletion
  // shrink eligible list as my atoms get marked
  // keep ndel,ncount,nall,nbefore current after each atom deletion
  if (mode == ATOM) {
    while (nall) {
      iwhichglobal = static_cast<int>(nall * random->uniform());
      if (iwhichglobal < nbefore)
        nbefore--;
      else if (iwhichglobal < nbefore + ncount) {
        iwhichlocal = iwhichglobal - nbefore;
        mark[list[iwhichlocal]] = 1;
        list[iwhichlocal] = list[ncount - 1];
        ncount--;
      }
      ndel++;
      nall--;
    }

    // molecule deletions
    // choose one atom in one molecule randomly across all procs
    // bcast mol ID and delete all atoms in that molecule on any proc
    // update deletion count by total # of atoms in molecule
    // shrink list of eligible candidates as any of my atoms get marked
    // keep ndel,ndeltopo,ncount,nall,nbefore current after each mol deletion
  } else {
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
      iwhichglobal = static_cast<int>(nall * random->uniform());
      if (iwhichglobal >= nbefore && iwhichglobal < nbefore + ncount) {
        iwhichlocal = iwhichglobal - nbefore;
        iatom = list[iwhichlocal];
        imolecule = molecule[iatom];
        me = comm->me;
      } else
        me = -1;

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
            if (atom->avec->bonds_allow) {
              if (force->newton_bond)
                ndeltopo[0] += atom->num_bond[i];
              else {
                for (j = 0; j < atom->num_bond[i]; j++) {
                  if (tag[i] < atom->bond_atom[i][j]) ndeltopo[0]++;
                }
              }
            }
            if (atom->avec->angles_allow) {
              if (force->newton_bond)
                ndeltopo[1] += atom->num_angle[i];
              else {
                for (j = 0; j < atom->num_angle[i]; j++) {
                  m = atom->map(atom->angle_atom2[i][j]);
                  if (m >= 0 && m < nlocal) ndeltopo[1]++;
                }
              }
            }
            if (atom->avec->dihedrals_allow) {
              if (force->newton_bond)
                ndeltopo[2] += atom->num_dihedral[i];
              else {
                for (j = 0; j < atom->num_dihedral[i]; j++) {
                  m = atom->map(atom->dihedral_atom2[i][j]);
                  if (m >= 0 && m < nlocal) ndeltopo[2]++;
                }
              }
            }
            if (atom->avec->impropers_allow) {
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
          list[i] = list[ncount - 1];
          ncount--;
        } else
          i++;
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
  }

  // delete my marked atoms
  // loop in reverse order to avoid copying marked atoms
  AtomVec *avec = atom->avec;

  for (i = nlocal - 1; i >= 0; i--) {
    if (mark[i]) {
      if (x[i][0] < 0.5 * (boxh + boxl)) {
        vnewl[0] += mass[type[i]] * vel[i][0];
        vnewl[1] += mass[type[i]] * vel[i][1];
        vnewl[2] += mass[type[i]] * vel[i][2];
        delete_left_file << "deleted" << "\n";
        delete_left_file.flush();
      } else {
        vnewr[0] += mass[type[i]] * vel[i][0];
        vnewr[1] += mass[type[i]] * vel[i][1];
        vnewr[2] += mass[type[i]] * vel[i][2];
        delete_right_file << "deleted" << "\n";
        delete_right_file.flush();
      }
      avec->copy(atom->nlocal - 1, i, 1);
      atom->nlocal--;
    }
  }

  // reset global natoms and bonds, angles, etc
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts
  atom->natoms -= ndel;

  if (mode == MOLECULE) {
    int all[4];
    MPI_Allreduce(ndeltopo, all, 4, MPI_INT, MPI_SUM, world);
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
}

/* ---------------------------------------------------------------------- */

void FixObmdMerged::try_inserting(Region *iregion_var, int stev, double *vnewl, double *vnewr)
{
  int i, m, n, nlocalprev, imol, natom, flag, flagall, me, sum_all, ninsert;
  double coord[3], lamda[3], delx, dely, delz, rsq;
  double r[3], vnew[3], rotmat[3][3], quat[4];
  double *newcoord;
  double entmp, mtmp;
  double *sublo, *subhi;
  int iter;

  double offset = 0.0;
  if (rateflag) offset = (update->ntimestep - nfirst) * update->dt * rate;

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

  if (!idnext) find_maxid();    // idnext is set to zero

  if (stev > 0) {
    for (ninsert = 0; ninsert != stev; ninsert++) {
      int success = 0;
      int attempt = 0;

      while (attempt < maxattempt) {
        attempt++;

        // choose random position for new particle within region
        if (distflag == DIST_UNIFORM) {
          do {
            coord[0] = xlo + random->uniform() * (xhi - xlo);
            coord[1] = ylo + random->uniform() * (yhi - ylo);
            coord[2] = zlo + random->uniform() * (zhi - zlo);
          } while (iregion_var->match(coord[0], coord[1], coord[2]) == 0);
        } else if (distflag == DIST_GAUSSIAN) {
          do {
            coord[0] = xmid + random->gaussian() * sigma;
            coord[1] = ymid + random->gaussian() * sigma;
            coord[2] = zmid + random->gaussian() * sigma;
          } while (iregion_var->match(coord[0], coord[1], coord[2]) == 0);
        } else
          error->all(FLERR, "Unknown particle distribution in fix obmd");

        // adjust vertical coord by offset
        if (dimension == 2)
          coord[1] += offset;
        else
          coord[2] += offset;

        // if global, reset vertical coord to be lo-hi above highest atom
        // if local, reset vertical coord to be lo-hi above highest "nearby" atom
        // local computation computes lateral distance between 2 particles w/ PBC
        // when done, have final coord of atom or center pt of molecule
        if (globalflag || localflag) {
          int dim;
          double max, maxall, delx, dely, delz, rsq;

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
              domain->minimum_image(delx, dely, delz);
              if (dimension == 2)
                rsq = delx * delx;
              else
                rsq = delx * delx + dely * dely;
              if (rsq > deltasq) continue;
            }
            if (x[i][dim] > max) max = x[i][dim];
          }

          MPI_Allreduce(&max, &maxall, 1, MPI_DOUBLE, MPI_MAX, world);
          if (dimension == 2)
            coord[1] = maxall + lo + random->uniform() * (hi - lo);
          else
            coord[2] = maxall + lo + random->uniform() * (hi - lo);
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
          imageflags[0] = ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
        } else {

          double rng = random->uniform();
          imol = 0;
          while (rng > molfrac[imol])
            imol++;    // it chooses imol randomly out of all defined mymols
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
          MathExtra::norm3(r);                             // unit vector
          MathExtra::axisangle_to_quat(r, theta, quat);    // quaternion
          MathExtra::quat_to_mat(quat, rotmat);            // rotmat
          for (i = 0; i < natom; i++) {
            MathExtra::matvec(rotmat, onemols[imol]->dx[i], coords[i]);
            coords[i][0] += coord[0];
            coords[i][1] += coord[1];
            coords[i][2] += coord[2];
            imageflags[i] =
                ((imageint) IMGMAX << IMG2BITS) | ((imageint) IMGMAX << IMGBITS) | IMGMAX;
            domain->remap(coords[i], imageflags[i]);
          }
        }

        // check distance between any existing atom and any inserted atom
        // if less than near, try again
        // use minimum_image() to account for PBC
        double **x = atom->x;
        int nlocal = atom->nlocal;
        flag = 0;
        if (nearflag) {
          for (m = 0; m < natom; m++) {
            for (i = 0; i < nlocal; i++) {
              delx = coords[m][0] - x[i][0];
              dely = coords[m][1] - x[i][1];
              delz = coords[m][2] - x[i][2];
              domain->minimum_image(delx, dely, delz);
              rsq = delx * delx + dely * dely + delz * delz;
              if (rsq < nearsq) {
                if (comm->me == 0)
                  std::cout << "NEAR denies in attempt No. " << attempt << "." << std::endl;
                flag = 1;
              }
            }
          }
        } else if (usherflag) {
          me = -1;
          entmp = usher(iregion_var, coords, etarget, natom, imol, iter);
          if (entmp < etarget + EPSILON) {
            if (comm->me == 0)
              std::cout << "USHER accepts at E = " << entmp << " in attempt No. " << attempt
                        << " with " << iter << " iterations" << std::endl;
          } else {
            if (comm->me == 0)
              std::cout << "USHER denies at E = " << entmp << " at attempt No. " << attempt
                        << std::endl;
            flag = 1;
          }
        }

        MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
        if (flagall) continue;

        if (nearflag) { entmp = near_energy(iregion_var, coords, natom, imol); }

        // proceed with insertion
        nlocalprev = atom->nlocal;

        // choose random velocity for new particle
        // for now instering with zero velocity
        // used for every atom in molecule
        vnew[0] = 0.0;
        vnew[1] = 0.0;
        vnew[2] = 0.0;

        // if target specified, change velocity vector accordingly
        if (targetflag) {
          double vel = sqrt(vnew[0] * vnew[0] + vnew[1] * vnew[1] + vnew[2] * vnew[2]);
          delx = tx - coord[0];
          dely = ty - coord[1];
          delz = tz - coord[2];
          double rsq = delx * delx + dely * dely + delz * delz;
          if (rsq > 0.0) {
            double rinv = sqrt(1.0 / rsq);
            vnew[0] = delx * rinv * vel;
            vnew[1] = dely * rinv * vel;
            vnew[2] = delz * rinv * vel;
          }
        }

        // this loop first just checks if an entire molecule can be inserted
        sum_all = check_mol_proc(coords, newcoord, lamda, sublo, subhi, dimension, natom, me);

        if (sum_all != natom) {
          if (comm->me == 0)
            std::cout << "Can not insert particle/molecule at timestep = " << update->ntimestep
                      << " me = " << me << std::endl;
          continue;
        }

        // happens rarely, e.g. if usher moves the atom/molecule outside of the (open) boundaries
        int check;
        check = check_mol_region(iregion_var, coords, natom);
        if (check == 1) continue;

        // check if new atoms are in my sub-box or above it if I am highest proc
        // if so, add atom to my list via create_atom()
        // initialize additional info about the atoms
        // set group mask to "all" plus fix group
        for (m = 0; m < natom; m++) {
          if (domain->triclinic) {
            domain->x2lamda(coords[m], lamda);
            newcoord = lamda;
          } else
            newcoord = coords[m];
          flag = check_proc(coords, newcoord, sublo, subhi, dimension);
          if (flag) {
            if (mode == ATOM) {
              mtmp = mass[ntype];
              atom->avec->create_atom(ntype, coords[m]);
            } else {
              mtmp = mtot;
              atom->avec->create_atom(ntype + onemols[imol]->type[m], coords[m]);
            }
            n = atom->nlocal - 1;
            atom->tag[n] = maxtag_all + m + 1;

            if (mode == MOLECULE) {
              if (atom->molecule_flag) {
                if (onemols[imol]->moleculeflag) {
                  atom->molecule[n] = maxmol_all + onemols[imol]->molecule[m];
                } else {
                  atom->molecule[n] = maxmol_all + 1;
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

            if (mode == MOLECULE) {
              onemols[imol]->quat_external = quat;
              atom->add_molecule_atom(onemols[imol], m, n, maxtag_all);
            }
            modify->create_attribute(n);
          }
        }

        // FixRigidSmall::set_molecule stores rigid body attributes
        // coord is new position of geometric center of mol, not COM
        // FixShake::set_molecule stores shake info for molecule
        if (mode == MOLECULE) {
          if (rigidflag)
            fixrigid->set_molecule(nlocalprev, maxtag_all, imol, coord, vnew, quat);
          else if (shakeflag)
            fixshake->set_molecule(nlocalprev, maxtag_all, imol, coord, vnew, quat);
        }

        success = 1;
        break;
      }

      if (!success && comm->me == 0)
        error->warning(FLERR, "Particle/molecule insertion was unsuccessful");

      // reset global natoms,nbonds,etc
      // increment maxtag_all and maxmol_all if necessary
      // if global map exists, reset it now instead of waiting for comm
      // since other pre-exchange fixes may use it
      // invoke map_init() b/c atom count has grown
      if (success) {
        if (coords[0][0] < 0.5 * (boxh + boxl)) {
          vnewl[0] +=
              0.0;    // add zero as you insert with zero velocity (if some velocity will be chosen, change this accordingly)
          vnewl[1] +=
              0.0;    // add zero as you insert with zero velocity (if some velocity will be chosen, change this accordingly)
          vnewl[2] +=
              0.0;    // add zero as you insert with zero velocity (if some velocity will be chosen, change this accordingly)
          insert_left_file << "inserted" << "\n";
          insert_left_file.flush();
        } else {
          vnewr[0] +=
              0.0;    // add zero as you insert with zero velocity (if some velocity will be chosen, change this accordingly)
          vnewr[1] +=
              0.0;    // add zero as you insert with zero velocity (if some velocity will be chosen, change this accordingly)
          vnewr[2] +=
              0.0;    // add zero as you insert with zero velocity	(if some velocity will be chosen, change this accordingly)
          insert_right_file << "inserted" << "\n";
          insert_right_file.flush();
        }
        atom->natoms += natom;

        if (atom->natoms < 0) error->all(FLERR, "Too many total atoms");
        if (mode == MOLECULE) {
          atom->nbonds += onemols[imol]->nbonds;
          atom->nangles += onemols[imol]->nangles;
          atom->ndihedrals += onemols[imol]->ndihedrals;
          atom->nimpropers += onemols[imol]->nimpropers;
        }
        maxtag_all += natom;
        if (maxtag_all >= MAXTAGINT) error->all(FLERR, "New atom IDs exceed maximum allowed ID");
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
        if (success) atom->map_init();
        atom->map_set();
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixObmdMerged::post_force(int vflag)
{
  if (update->ntimestep % nevery) { return; }

  // parallel forces
  reg_force(vflag, iregion, momentumForce_left, step_parallel);
  reg_force(vflag, iregion2, momentumForce_right, step_parallel);

  // tangential forces
  reg_force_perp(vflag, iregion3, shearForce_left, step_perp);
  reg_force_perp(vflag, iregion4, shearForce_right, step_perp);
}

/* ---------------------------------------------------------------------- */

double FixObmdMerged::g_par_global_charged(Region *region, int step)
{

  // using molecule's center-of-mass
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
      if (x[i][0] < lower_x + buffer_size) {                        // LEFT BUFFER
        if (x[i][0] < (lower_x + (1.0 - g_fac) * buffer_size)) {    // only mass
          g_par_all += mass_temp;
        } else {    // sigmoidal mass
          carg = 1.0 / g_fac * MY_PI * (x[i][0] - buffer_size - lower_x) / (-buffer_size) -
              MY_PI;    // OK (also for -lower_x)
          g_par_all += 0.5 * (1.0 + cos(carg)) * mass_temp;
        }
      }
      if (x[i][0] > upper_x - buffer_size) {    // RIGHT BUFFER
        if (x[i][0] > (upper_x - (1.0 - g_fac) * buffer_size)) {
          g_par_all += mass_temp;
        } else {
          carg = 1.0 / g_fac * MY_PI * (x[i][0] - upper_x + buffer_size) /
              (buffer_size) -MY_PI;    // OK (also for -lower_x)
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

/* ---------------------------------------------------------------------- */

double FixObmdMerged::g_par_local_charged(double mass, Region *region, double xv, int step)
{
  // using molecule's center-of-mass
  double upper_x = domain->boxhi[0];    // upper
  double lower_x = domain->boxlo[0];    // lower

  double carg = 0.0;
  if (!step) {
    if (xv < lower_x + buffer_size) {    // LEFT BUFFER
      if (xv < (lower_x + (1.0 - g_fac) * buffer_size)) {
        return mass;
      } else {    // sigmodial mass
        carg = 1.0 / g_fac * MY_PI * (xv - buffer_size - lower_x) / (-buffer_size) - MY_PI;
        return 0.5 * (1.0 + cos(carg)) * mass;
      }
    }
    if (xv > upper_x - buffer_size) {    // RIGHT BUFFER
      if (xv > (upper_x - (1.0 - g_fac) * buffer_size)) {
        return mass;
      } else {    // sigmodial mass
        carg = 1.0 / g_fac * MY_PI * (xv - upper_x + buffer_size) / (buffer_size) -MY_PI;
        return 0.5 * (1.0 + cos(carg)) * mass;
      }
    }
    // OTHER IS NOT OF MY INTEREST ... ROI
  } else {
    error->all(FLERR, "For now implemented only for !step");
  }
}

/* ---------------------------------------------------------------------- */

double FixObmdMerged::g_perp_global_charged(Region *iregion_var, int step)
{
  std::cout << "g_perp_global()" << "\n";
  std::cout << "g_perp_global(): step: " << step << "\n";
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

  // create test_force to check obmd
  double test_force[3];
  test_force[0] = test_force[1] = test_force[2] = 0.0;
  double sestevek = 0.0;

  double test_force_all[3];
  test_force_all[0] = test_force_all[2] = test_force_all[3] = 0.0;

  int allnu = 0;

  for (int i = 0; i < nlocal; i++) {
    if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
    domain->unmap(x[i], image[i], unwrap);
    mass_tmp = mass[type[i]];
    if (mode == MOLECULE) mass_tmp = mtot;

    gloctmp = g_par_local_charged(mass[type[i]], region, x[i][0], step);

    xval = momentumForce[0] * gloctmp / gtmp;
    yval = momentumForce[1] * gloctmp / gtmp;
    zval = momentumForce[2] * gloctmp / gtmp;

    f[i][0] += xval;
    f[i][1] += yval;
    f[i][2] += zval;

    test_force[0] += momentumForce[0] * gloctmp / gtmp;
    test_force[1] += momentumForce[1] * gloctmp / gtmp;
    test_force[2] += momentumForce[2] * gloctmp / gtmp;

    sestevek += gloctmp;

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

  MPI_Allreduce(test_force, test_force_all, 3, MPI_DOUBLE, MPI_SUM, world);
}

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
  foriginal[0] = foriginal[1] = foriginal[2] = foriginal[3] = 0.0;
  force_flag = 0;
  int *type = atom->type;
  double gtmp = g_perp_global_charged(region, step);
  double mass_tmp, gloctmp;
  double unwrap[3];
  double xval, yval, zval;

  // create test_force to check obmd
  double test_perp_force[3];
  test_perp_force[0] = test_perp_force[1] = test_perp_force[2] = 0.0;
  double sestevek_perp = 0.0;

  double test_perp_force_all[3];
  test_perp_force_all[0] = test_perp_force_all[1] = test_perp_force_all[2] = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (region && !region->match(x[i][0], x[i][1], x[i][2])) continue;
    domain->unmap(x[i], image[i], unwrap);
    mass_tmp = mass[type[i]];
    if (mode == MOLECULE) mass_tmp = mtot;

    gloctmp = mass[type[i]];

    xval = shearForce[0] * gloctmp / gtmp;
    yval = shearForce[1] * gloctmp / gtmp;
    zval = shearForce[2] * gloctmp / gtmp;

    f[i][0] += xval;
    f[i][1] += yval;
    f[i][2] += zval;

    test_perp_force[0] += shearForce[0] * gloctmp / gtmp;
    test_perp_force[1] += shearForce[1] * gloctmp / gtmp;
    test_perp_force[2] += shearForce[2] * gloctmp / gtmp;

    sestevek_perp += gloctmp;

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

  MPI_Allreduce(test_perp_force, test_perp_force_all, 3, MPI_DOUBLE, MPI_SUM, world);
}

/* ---------------------------------------------------------------------- */

double FixObmdMerged::usher(Region *iregion_var, double **coords, double etarget, int natom,
                            int imol, int &iter)
{

  int i, m, type_temp;
  double entmp, torqabs, dtheta, ds, fabs;
  double quat[4], rotmat[3][3], fusher_tmp[3], fusher_all[3], xcom[3];
  double torq[3], torq_all[3], torq_tmp[3], coords_tmp[3];

  double entmp_all;
  MathExtra::zero3(fusher_tmp);

  i = 0;
  while (i < nattempt) {

    entmp = 0.0;
    type_temp = ntype;
    MathExtra::zero3(fusher);
    MathExtra::zero3(torq);

    if (chargeflag) {
      mol_center_of_mass(natom, imol, coords, xcom);
    } else {
      center_of_mass(natom, coords, xcom);
    }

    for (m = 0; m < natom; m++) {
      if (mode == MOLECULE) type_temp = ntype + onemols[imol]->type[m];

      if (chargeflag) {
        entmp +=
            energy_atomistic_obmd(iregion_var, onemols[imol]->q[m], type_temp, coords[m], fusher);
      } else {
        entmp += energy(1, type_temp, coords[m],
                        fusher);    // i - does not matter, itype inserted, coords[m], fusher
      }

      // torque
      MathExtra::copy3(fusher, fusher_tmp);
      calc_torque(natom, coords, xcom, fusher_tmp, torq_tmp);
      MathExtra::add3(torq, torq_tmp, torq);
    }

    MPI_Allreduce(fusher, fusher_all, 3, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&entmp, &entmp_all, 1, MPI_DOUBLE, MPI_SUM, world);
    MPI_Allreduce(&torq, &torq_all, 3, MPI_DOUBLE, MPI_SUM, world);

    if (entmp_all < etarget + EPSILON)
      break;
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
      if (check == 1) break;    // exit

    } else {
      fabs = sqrt(fusher_all[0] * fusher_all[0] + fusher_all[1] * fusher_all[1] +
                  fusher_all[2] * fusher_all[2]);
      if (fabs < EPSILON) continue;
      ds = std::min((entmp_all - etarget) / fabs, ds0);

      if (mode == MOLECULE) {
        torqabs =
            sqrt(torq_all[0] * torq_all[0] + torq_all[1] * torq_all[1] + torq_all[2] * torq_all[2]);
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
    i += 1;
  }
  iter = i;
  entmp = entmp_all;

  return entmp;
}

/* ---------------------------------------------------------------------- */

double FixObmdMerged::near_energy(Region *iregion_var, double **coords, int natom, int imol)
{
  double entmp, entmp_all;
  double fusher[3];
  int m, type_temp;
  entmp = 0.0;
  type_temp = ntype;

  for (m = 0; m < natom; m++) {
    if (mode == MOLECULE) type_temp = ntype + onemols[imol]->type[m];

    if (chargeflag) {
      entmp +=
          energy_atomistic_obmd(iregion_var, onemols[imol]->q[m], type_temp, coords[m], fusher);
    } else {
      entmp += energy(1, type_temp, coords[m], fusher);
    }
  }

  MPI_Allreduce(&entmp, &entmp_all, 1, MPI_DOUBLE, MPI_SUM, world);
  return entmp_all;
}

/* ---------------------------------------------------------------------- */

int FixObmdMerged::check_proc(double **coords, double *newcoord, double *sublo, double *subhi,
                              int dimension)
{
  int flag = 0;
  if (newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] && newcoord[1] >= sublo[1] &&
      newcoord[1] < subhi[1] && newcoord[2] >= sublo[2] && newcoord[2] < subhi[2])
    flag = 1;
  else if (dimension == 3 && newcoord[2] >= domain->boxhi[2]) {    /// for 3D
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[2] == comm->procgrid[2] - 1 && newcoord[0] >= sublo[0] &&
          newcoord[0] < subhi[0] && newcoord[1] >= sublo[1] && newcoord[1] < subhi[1])
        flag = 1;
    } else {
      if (comm->mysplit[2][1] == 1.0 && newcoord[0] >= sublo[0] && newcoord[0] < subhi[0] &&
          newcoord[1] >= sublo[1] && newcoord[1] < subhi[1])
        flag = 1;
    }
  } else if (dimension == 2 && newcoord[1] >= domain->boxhi[1]) {    /// for 2D
    if (comm->layout != Comm::LAYOUT_TILED) {
      if (comm->myloc[1] == comm->procgrid[1] - 1 && newcoord[0] >= sublo[0] &&
          newcoord[0] < subhi[0])
        flag = 1;
    } else {
      if (comm->mysplit[1][1] == 1.0 && newcoord[0] >= sublo[0] && newcoord[0] < subhi[0]) flag = 1;
    }
  }
  return flag;
}

/* ---------------------------------------------------------------------- */

int FixObmdMerged::check_mol_proc(double **coords, double *newcoord, double *lamda, double *sublo,
                                  double *subhi, int dimension, int natom, int &me)
{
  int m;
  int sum_all_tmp = 0;
  int sum_all_tmpall, flag;

  for (m = 0; m < natom; m++) {
    if (domain->triclinic) {
      domain->x2lamda(coords[m], lamda);
      newcoord = lamda;
    } else
      newcoord = coords[m];
    flag = check_proc(coords, newcoord, sublo, subhi, dimension);
    sum_all_tmp += flag;
  }
  MPI_Allreduce(&sum_all_tmp, &sum_all_tmpall, 1, MPI_INT, MPI_SUM, world);
  if (sum_all_tmp == natom) me = comm->me;

  return sum_all_tmpall;
}

/* ---------------------------------------------------------------------- */

int FixObmdMerged::check_mol_region(Region *region, double **coords, int natom)
{
  int m, flag, flagall;
  flag = 0;
  for (m = 0; m < natom; m++) {
    if (region && !region->match(coords[m][0], coords[m][1], coords[m][2])) flag = 1;
  }
  MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_MAX, world);
  if (comm->me == 0 && flagall == 1)
    std::cout << "USHER has moved the particle/molecule too much" << std::endl;

  return flagall;
}

/* ---------------------------------------------------------------------- */

void FixObmdMerged::center_of_mass(int natom, double **coords, double *xcom)
{
  int m;
  MathExtra::zero3(xcom);

  for (m = 0; m < natom; m++) {
    xcom[0] += coords[m][0];
    xcom[1] += coords[m][2];
    xcom[2] += coords[m][1];
  }

  xcom[0] = xcom[0] / natom;
  xcom[1] = xcom[1] / natom;
  xcom[2] = xcom[2] / natom;
}

/* ---------------------------------------------------------------------- */

void FixObmdMerged::mol_center_of_mass(int natom, int imol, double **coords, double *xcom)
{
  double *mass = atom->mass;
  int m;
  double massmol = 0.0;
  MathExtra::zero3(xcom);

  for (m = 0; m < natom; m++) {
    xcom[0] += mass[onemols[imol]->type[m]] * coords[m][0];
    xcom[1] += mass[onemols[imol]->type[m]] * coords[m][1];
    xcom[2] += mass[onemols[imol]->type[m]] * coords[m][2];

    massmol += mass[onemols[imol]->type[m]];
  }

  xcom[0] = xcom[0] / massmol;
  xcom[1] = xcom[1] / massmol;
  xcom[2] = xcom[2] / massmol;
}

/* ---------------------------------------------------------------------- */

void FixObmdMerged::calc_torque(int natom, double **coords, double *xcom, double *force_tmp,
                                double *torq)
{
  int m;
  double xrel[3];

  for (m = 0; m < natom; m++) {
    xrel[0] = coords[m][0] - xcom[0];
    xrel[1] = coords[m][1] - xcom[1];
    xrel[2] = coords[m][2] - xcom[2];
  }
  dtheta0 = 0.1;
  dtheta0 = 0.1;
  MathExtra::cross3(xrel, force_tmp, torq);
}

/* ---------------------------------------------------------------------- */

double FixObmdMerged::energy(int i, int itype, double *coord, double *fusher)
{
  double delx, dely, delz, rsq;
  int jtype;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  pair = force->pair;
  cutsq = force->pair->cutsq;

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;
  double total_energy = 0.0;

  for (int j = 0; j < nlocal; j++) {
    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];

    domain->minimum_image(delx, dely, delz);
    rsq = delx * delx + dely * dely + delz * delz;

    jtype = type[j];

    if (rsq < cutsq[itype][jtype]) {
      total_energy += pair->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

      fusher[0] += fpair * delx;
      fusher[1] += fpair * dely;
      fusher[2] += fpair * delz;
    }
  }

  return total_energy;
}

/* ---------------------------------------------------------------------- */

double FixObmdMerged::energy_atomistic_obmd(Region *iregion_var, double qi, int itype,
                                            double *coord, double *fusher)
{
  double delx, dely, delz, rsq;
  int jtype;

  double **x = atom->x;
  int *type = atom->type;
  int nall = atom->nlocal + atom->nghost;
  cutsq = force->pair->cutsq;
  tagint *tag = atom->tag;
  int nlocal = atom->nlocal;

  // potential acting between atomistic particles
  auto pair = dynamic_cast<PairLJCutRF *>(force->pair_match("lj/cut/rf", 1));

  double fpair = 0.0;
  double factor_coul = 1.0;
  double factor_lj = 1.0;
  double total_energy = 0.0;

  for (int j = 0; j < nlocal; j++) {
    delx = coord[0] - x[j][0];
    dely = coord[1] - x[j][1];
    delz = coord[2] - x[j][2];

    domain->minimum_image(delx, dely, delz);
    rsq = delx * delx + dely * dely + delz * delz;

    jtype = type[j];

    if (rsq < cutsq[itype][jtype]) {
      total_energy +=
          pair->single_atomistic_obmd(qi, j, itype, jtype, rsq, factor_coul, factor_lj, fpair);

      fusher[0] += fpair * delx;
      fusher[1] += fpair * dely;
      fusher[2] += fpair * delz;
    }
  }

  return total_energy;
}

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

/* ----------------------------------------------------------------------
   parse optional parameters at end of input line
------------------------------------------------------------------------- */

void FixObmdMerged::options(int narg, char **arg)
{
  // defaults
  nearflag = 0;      // near
  usherflag = 0;     // usher
  chargeflag = 0;    // particle

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

  g_fac = 0.25;
  g_fac_inv = 1.0 / g_fac;
  mol_len = 1.0;    // particle
  buffer_size = 0.30 * (domain->boxhi[0] - domain->boxlo[0]);

  etarget = 3.6;
  ds0 = 0.1;
  dtheta0 = 0.35;
  uovlp = 10000.0;
  dsovlp = 3.0;
  eps = 0.15;
  nattempt = 40;
  maxattempt = 1;

  globalflag = localflag = 0;
  lo = hi = deltasq = 0.0;
  nearsq = 0.0;
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

  step_parallel = 0;    // normal > smooth distribution function
  step_perp = 1;        // tangential > heaviside step distribution function

  int imol = -1;
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "region1") == 0) {
      iregion = domain->get_region_by_id(arg[iarg + 1]);
      if (!iregion) error->all(FLERR, "Region ID {} for fix obmd does not exist", arg[iarg + 1]);
      int n = strlen(arg[iarg + 1]) + 1;
      idregion = new char[n];
      strcpy(idregion, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "region2") == 0) {
      iregion2 = domain->get_region_by_id(arg[iarg + 1]);
      if (!iregion2) error->all(FLERR, "Region ID {} for fix obmd does not exist", arg[iarg + 1]);
      int n = strlen(arg[iarg + 1]) + 1;
      idregion2 = new char[n];
      strcpy(idregion2, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "region3") == 0) {
      iregion3 = domain->get_region_by_id(arg[iarg + 1]);
      if (!iregion3) error->all(FLERR, "Region ID {} for fix obmd does not exist", arg[iarg + 1]);
      int n = strlen(arg[iarg + 1]) + 1;
      idregion3 = new char[n];
      strcpy(idregion3, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "region4") == 0) {
      iregion4 = domain->get_region_by_id(arg[iarg + 1]);
      if (!iregion4) error->all(FLERR, "Region ID {} for fix obmd does not exist", arg[iarg + 1]);
      int n = strlen(arg[iarg + 1]) + 1;
      idregion4 = new char[n];
      strcpy(idregion4, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "region5") == 0) {
      iregion5 = domain->get_region_by_id(arg[iarg + 1]);
      if (!iregion5) error->all(FLERR, "Region ID {} for fix obmd does not exist", arg[iarg + 1]);
      int n = strlen(arg[iarg + 1]) + 1;
      idregion5 = new char[n];
      strcpy(idregion5, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "region6") == 0) {
      iregion6 = domain->get_region_by_id(arg[iarg + 1]);
      if (!iregion6) error->all(FLERR, "Region ID {} for fix obmd does not exist", arg[iarg + 1]);
      int n = strlen(arg[iarg + 1]) + 1;
      idregion6 = new char[n];
      strcpy(idregion6, arg[iarg + 1]);
      iarg += 2;
    } else if (strcmp(arg[iarg], "buffersize") == 0) {
      buffer_size = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (buffer_size <= 0.0)
        error->all(FLERR,
                   "Illegal fix obmd command. Parameter buffersize should be larger than 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "gfac") == 0) {
      g_fac = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      g_fac_inv = 1.0 / g_fac;
      if (g_fac < 0.0 || g_fac > 1.0)
        error->all(FLERR, "Illegal fix obmd command. Parameter gfac should be between 0.0 and 1.0");
      iarg += 2;
    } /* else if (strcmp(arg[iarg],"alpha") == 0) {
      alpha = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (alpha < 0.0) error->all(FLERR,"Parameter alpha should be bigger than 0.0");
      iarg += 2;
    }  else if (strcmp(arg[iarg],"tau") == 0) {
      tau = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (tau < 0.0) error->all(FLERR,"Parameter tau should be bigger than 0.0");
      iarg += 2;
    } else if (strcmp(arg[iarg],"nbuf") == 0) {
      nbuf = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      if (nbuf < 0.0) error->all(FLERR,"Parameter nbuf should be bigger than 0");
      iarg += 2;
    } */
    else if (strcmp(arg[iarg], "stepparallel") == 0) {    // step function for parallel forces
      step_parallel = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);    // 1 / 0
      if (step_parallel != 0)
        error->all(FLERR,
                   "For now, only the smooth transition distribution function is implemented.");
      iarg += 2;
    } else if (strcmp(arg[iarg], "stepperp") == 0) {    // step function for perpendicular forces
      step_perp = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);    // 1 / 0
      if (step_perp != 1)
        error->all(FLERR, "For now, only the heaviside step distribution function is implemented.");
      iarg += 2;
    } else if (strcmp(arg[iarg], "maxattempt") == 0) {    // not USHER attempts
      maxattempt = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "usher") == 0) {
      if (iarg + 9 > narg) error->all(FLERR, "Illegal fix obmd command");
      usherflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (usherflag == 1 && nearflag == 1)
        error->all(FLERR, "You can not have both usher and near");
      // options
      etarget = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      ds0 = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      dtheta0 = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      uovlp = utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      dsovlp = utils::numeric(FLERR, arg[iarg + 6], false, lmp);
      eps = utils::numeric(FLERR, arg[iarg + 7], false, lmp);
      nattempt = utils::inumeric(FLERR, arg[iarg + 8], false, lmp);
      iarg += 9;
    } else if (strcmp(arg[iarg], "mol") == 0) {
      if (iarg + 3 > narg) error->all(FLERR, "Illegal fix obmd command");
      imol = atom->find_molecule(arg[iarg + 1]);
      if (imol == -1) error->all(FLERR, "Molecule template ID for fix obmd does not exist");
      mol_len = utils::inumeric(FLERR, arg[iarg + 2], false, lmp);
      mode = MOLECULE;
      onemols = &atom->molecules[imol];
      onemols[imol]->compute_mass();
      mtot = onemols[imol]->masstotal;
      nmol = onemols[0]->nset;
      delete[] molfrac;
      molfrac = new double[nmol];
      molfrac[0] = 1.0 / nmol;
      for (int i = 1; i < nmol - 1; i++) molfrac[i] = molfrac[i - 1] + 1.0 / nmol;
      molfrac[nmol - 1] = 1.0;
      iarg += 3;
    } else if (strcmp(arg[iarg], "molfrac") == 0) {
      if (mode != MOLECULE) error->all(FLERR, "You can not use molfrac without MOLECULE mode");
      if (iarg + nmol + 1 > narg) error->all(FLERR, "Illegal fix obmd command");
      molfrac[0] = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      for (int i = 1; i < nmol; i++)
        molfrac[i] = molfrac[i - 1] + utils::numeric(FLERR, arg[iarg + i + 1], false, lmp);
      if (molfrac[nmol - 1] < 1.0 - EPSILON || molfrac[nmol - 1] > 1.0 + EPSILON)
        error->all(FLERR, "Illegal fix obmd command");
      molfrac[nmol - 1] = 1.0;
      iarg += nmol + 1;
    } else if (strcmp(arg[iarg], "rigid") == 0) {
      int n = strlen(arg[iarg + 1]) + 1;
      delete[] idrigid;
      idrigid = new char[n];
      strcpy(idrigid, arg[iarg + 1]);
      rigidflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "shake") == 0) {
      int n = strlen(arg[iarg + 1]) + 1;
      delete[] idshake;
      idshake = new char[n];
      strcpy(idshake, arg[iarg + 1]);
      shakeflag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg], "id") == 0) {
      if (strcmp(arg[iarg + 1], "max") == 0)
        idnext = 0;
      else if (strcmp(arg[iarg + 1], "next") == 0)
        idnext = 1;
      else
        error->all(FLERR, "Illegal fix obmd command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "global") == 0) {
      globalflag = 1;
      localflag = 0;
      lo = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      hi = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "local") == 0) {
      localflag = 1;
      globalflag = 0;
      lo = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      hi = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      deltasq = utils::numeric(FLERR, arg[iarg + 3], false, lmp) *
          utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      iarg += 4;
    } else if (strcmp(arg[iarg], "near") == 0) {
      nearflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      nearsq = utils::numeric(FLERR, arg[iarg + 2], false, lmp) *
          utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (usherflag == 1 && nearflag == 1)
        error->all(FLERR, "You can not have both usher and near");
      iarg += 3;
    } else if (strcmp(arg[iarg], "charged") == 0) {
      chargeflag = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (chargeflag == 1 && mode != MOLECULE)
        error->all(FLERR, "You can not use charged without MOLECULE mode");
      iarg += 2;
    } else if (strcmp(arg[iarg], "rate") == 0) {
      rateflag = 1;
      rate = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "vx") == 0) {
      vxlo = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vxhi = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "vy") == 0) {
      vylo = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vyhi = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "vz") == 0) {
      vzlo = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      vzhi = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      iarg += 3;
    } else if (strcmp(arg[iarg], "orient") == 0) {
      orientflag = 1;
      rx = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      ry = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      rz = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (domain->dimension == 2 && (rx != 0.0 || ry != 0.0))
        error->all(FLERR, "Illegal fix obmd orient settings");
      if (rx == 0.0 && ry == 0.0 && rz == 0.0)
        error->all(FLERR, "Illegal fix obmd orient settings");
      iarg += 4;
    } else if (strcmp(arg[iarg], "units") == 0) {
      if (strcmp(arg[iarg + 1], "box") == 0)
        scaleflag = 0;
      else if (strcmp(arg[iarg + 1], "lattice") == 0)
        scaleflag = 1;
      else
        error->all(FLERR, "Illegal fix obmd command");
      iarg += 2;
    } else if (strcmp(arg[iarg], "gaussian") == 0) {
      xmid = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      ymid = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      zmid = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      sigma = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      distflag = DIST_GAUSSIAN;
      iarg += 5;
    } else if (strcmp(arg[iarg], "target") == 0) {
      tx = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      ty = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      tz = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      targetflag = 1;
    } else
      error->all(FLERR, "Illegal fix obmd command");
  }

  if (usherflag == 0 && nearflag == 0) error->all(FLERR, "You must pick either usher or near");
}

/* ----------------------------------------------------------------------
   extract particle radius for atom type = itype
------------------------------------------------------------------------- */

void *FixObmdMerged::extract(const char *str, int &itype)
{
  if (strcmp(str, "radius") == 0) {
    if (mode == ATOM) {
      if (itype == ntype)
        oneradius = 0.5;
      else
        oneradius = 0.0;

    } else {

      // loop over onemols molecules
      // skip a molecule with no atoms as large as itype

      oneradius = 0.0;
      for (int i = 0; i < nmol; i++) {
        if (itype > ntype + onemols[i]->ntypes) continue;
        double *radius = onemols[i]->radius;
        int *type = onemols[i]->type;
        int natoms = onemols[i]->natoms;

        // check radii of atoms in Molecule with matching types
        // default to 0.5, if radii not defined in Molecule
        //   same as atom->avec->create_atom(), invoked in pre_exchange()

        for (int i = 0; i < natoms; i++)
          if (type[i] + ntype == itype) {
            if (radius)
              oneradius = MAX(oneradius, radius[i]);
            else
              oneradius = MAX(oneradius, 0.5);
          }
      }
    }
    itype = 0;
    return &oneradius;
  }

  return nullptr;
}
