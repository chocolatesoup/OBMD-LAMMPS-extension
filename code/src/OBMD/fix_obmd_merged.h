/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
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

#ifdef FIX_CLASS

FixStyle(obmd, FixObmdMerged)

#else

#ifndef LMP_FIX_OBMD_MERGED_H
#define LMP_FIX_OBMD_MERGED_H

#include "fix.h"
#include <fstream>

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

 private:
  int ninsert, ntype, nfreq, seed;
  class Region *iregion, *iregion2, *iregion3, *iregion4, *iregion5, *iregion6;
  int globalflag, localflag, maxattempt, rateflag, scaleflag, targetflag;
  int mode, rigidflag, shakeflag, idnext, distflag, orientflag, usherflag, nearflag, chargeflag;
  double lo, hi, deltasq, nearsq, rate, sigma;
  double vxlo, vxhi, vylo, vyhi, vzlo, vzhi;
  double xlo, xhi, ylo, yhi, zlo, zhi, xmid, ymid, zmid;
  double rx, ry, rz, tx, ty, tz;

  int imol;
  int iarg;    // args

  char *idregion;
  char *idregion2;
  char *idregion3;
  char *idregion4;
  char *idregion5;
  char *idregion6;
  char *idrigid, *idshake;

  int alphastyle, taustyle, nbufstyle;
  int alphavar, tauvar, nbufvar;
  char *alphastr, *taustr, *nbufstr;

  int pxxstyle, pxystyle, pxzstyle, dpstyle, freqstyle;
  int pxxvar, pxyvar, pxzvar, dpvar, freqvar;
  char *pxxstr, *pxystr, *pxzstr, *dpstr, *freqstr;

  double tau, alpha, buffer_size, shear_size, nbuf;
  double pxx, pxy, pxz, dpxx, t0_left, t0_right, lambda, freq;
  double mol_len;

  void check_ghosts();
  bigint lastcheck;

  double vcml[3];
  double vcmr[3];

  // for inserting and deleting of particles/molecules
  double vnewl[3];
  double vnewr[3];
  double vnewl_all[3];
  double vnewr_all[3];

  int me;
  class Molecule **onemols;
  int nmol, natom_max;
  double *molfrac;
  double **coords;
  imageint *imageflags;
  class Fix *fixrigid, *fixshake;
  double oneradius;

  int xstyle, xvar;
  char *xstr;
  int ystyle, yvar;
  char *ystr;
  int zstyle, zvar;
  char *zstr;
  int t0_left_style, t0_left_var;
  char *t0_left_str;
  int t0_right_style, t0_right_var;
  char *t0_right_str;

  int nfirst, ninserted;
  tagint maxtag_all, maxmol_all;
  class RanPark *random;

  void find_maxid();
  void options(int, char **);
  double momentumForce_left[3];
  double momentumForce_right[3];
  int *list, *mark;
  int ndeleted, ncount2, nmax, nmax2;

  void try_inserting(Region *, int, double *vnewl, double *vnewr);
  void try_deleting(Region *, double *vnewl, double *vnewr);
  void reg_force(int, Region *, double *, int);
  void reg_force_perp(int, Region *, double *, int);
  double g_par_global(Region *, int);
  double g_par_global_charged(Region *, int);
  double g_perp_global(Region *, int);
  double g_perp_global_charged(Region *, int);
  double g_par_local(Region *, double, int);
  double g_par_local_charged(double, Region *, double, int);
  double energy(int, int, double *, double *);
  double usher(Region *, double **, double, int, int, int &);
  double energy_atomistic_obmd(Region *, double, int, double *,
                               double *);    // charge qi, type, coords, fusher
  double near_energy(Region *, double **, int, int);
  void center_of_mass(int, double **, double *);
  void mol_center_of_mass(int, int, double **, double *);
  void calc_torque(int, double **, double *, double *, double *);
  int check_proc(double **, double *, double *, double *, int);
  int check_mol_proc(double **, double *, double *, double *, double *, int, int, int &);
  int check_mol_region(Region *, double **, int);
  void vcm_internal_sq(int, double, double *, double *, double *, Region *);
  double xvalue, yvalue, zvalue;
  int varflag;

  double foriginal[4], foriginal_all[4];
  int force_flag;
  int print_delete_flag;
  int print_insert_flag;
  int stepflag;
  int step_parallel, step_perp;
  double g_fac, g_fac_inv;
  int maxatom;
  double **sforce;

  class Compute *pressure_tensor;

  double mtot;

  double fusher[3];
  double shearForce_left[3];
  double shearForce_right[3];
  double simulation_time;
  double pressure_wave;

  double **cutsq;
  class Pair *pair;
  double masstotal;
  int nattempt;    // usher
  double uovlp, dsovlp, eps, ds0, dtheta0, etarget;
};

}    // namespace LAMMPS_NS

#endif
#endif
