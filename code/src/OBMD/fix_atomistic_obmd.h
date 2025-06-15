/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(atomistic/obmd,FixATMOBMD) // all-atom simulations 

#else

#ifndef LMP_FIX_ATOMISTIC_OBMD_H
#define LMP_FIX_ATOMISTIC_OBMD_H

#include "fix.h"
#include <fstream>

namespace LAMMPS_NS {

class FixATMOBMD : public Fix {

 public:
  FixATMOBMD(class LAMMPS *, int, char **);
  ~FixATMOBMD();

  int setmask();
  void init();
  void pre_exchange();
  void post_force(int);
  void setup(int);
  // PPapez COMMENT: not writing or reading restart \& no minimization
  // void min_setup(int);
  // void min_post_force(int);
  // void write_restart(FILE *);
  // void restart(char *);
  //  PPapez COMMENT: only molecules (using molecular template)
  // void *extract(const char *, int &);
  // virtual double compute_scalar();
  // virtual double compute_vector(int);

 private:
  
  char *idregion; // PPapez COMMENT: initialize idregion
  char *idregion2;
  char *idregion3;
  char *idregion4;
  char *idregion5;
  char *idregion6;
  /* char *idregion7;
  char *idregion8; */
  // char *idrigid;
  char *idshake;
  class Molecule **onemols;
  double *molfrac;
  double **coords;
  imageint *imageflags;
  class Fix *fixshake;
  class RanPark *random;
  int *list;
  int *mark;

  int ntype, nfreq, seed;
  int dpstyle, freqstyle;
  int dpvar, freqvar;
  char *dpstr, *freqstr;
  double pxx, pxy, pxz, dpxx, freq;

  void options(int, char **);
  class Region *iregion, *iregion2, *iregion3, *iregion4, *iregion5, *iregion6; // *iregion7, *iregion8;
  int shakeflag;
  double buffer_size;
  int step_parallel, step_perp;
  double g_fac, g_fac_inv;
  // one group
  double alpha;
  double tau;
  int nbuf;
  double etgt;
  double ds0;
  double dtheta0;
  double uovlp;
  double dsovlp;
  double eps;
  int nattempt; // usher
  int maxattempt;
  int distflag;
  double sigma, xmid, ymid, zmid;
  int idnext; 
  int imol;
  int mol_len;
  int usherflag; // PPapez COMMENT: omitting near 
  int molflag; // PPapez COMMENT: using mol template <-> usher
  int mode; // PPapez COMMENT: packed in mode == MOLECULE
  int iarg; // args
  double mtot;
  int nmol;
  int orientflag;
  double rx, ry, rz; 
  int natom_max;

  // for inserting and deleting of particles/molecules
  double vnewl[3];
  double vnewr[3];
  double vnewl_all[3];
  double vnewr_all[3];

  double momentumForce_left[3];
  double momentumForce_right[3];
  double shearForce_left[3];
  double shearForce_right[3];
  /* double pressureForce_left[3]; // only pressure contribution
  double pressureForce_right[3]; // only pressure contribution */
  double simulation_time;
  double pressure_wave;

  int nmax;
  int ndeleted;
  
  void find_maxid();
  void try_deleting(Region *, double *, double *);
  void try_inserting(Region *, int, double *, double *); 
  void mol_center_of_mass(int, int, double **, double *);
  double usher(Region *, double **, double, int, int,int&);
  double energy_atomistic_obmd(Region *, double, int, double *, double *); // charge qi, type, coords, fusher
  void calc_torque(int, double**, double*, double*, double*);
  int check_mol_region(Region*, double **, int); 
  void try_new_position(Region *, double **, int);
  int check_proc(double**, double*, double *, double*, int);
  int check_mol_proc(double**, double*, double*, double*, double*, int, int,int&);
    
  double fusher[3];
  double xlo, xhi, ylo, yhi, zlo, zhi;

  std::ofstream delete_left_file;
  std::ofstream delete_right_file;
  std::ofstream insert_left_file;
  std::ofstream insert_right_file;

  void reg_force(int, Region*,double*,int);
  void reg_force_perp(int, Region*,double*,int);
  double g_par_global(Region*, int);
  double g_perp_global(Region*, int);  
  double g_par_local(double, Region*, double, int); 
 
  /* ============================================================ */ 
 
  void check_ghosts();
  bigint lastcheck;

  double oneradius;

  int nfirst,ninserted;
  tagint maxtag_all,maxmol_all;
  int varflag;

  double foriginal[4], foriginal_all[4];
  int force_flag;
  int stepflag;
  int maxatom;
  double **sforce;
  
  char *id_press;
  class Compute *pressure_tensor;
  double **cutsq;
  class Pair *pair;
  double masstotal;
  double com_tmp[3], com_tmp_all[3]; // PP testing
  
};
}

#endif
#endif