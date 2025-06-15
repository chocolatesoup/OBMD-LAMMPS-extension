// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_nve.h"

#include "atom.h"
#include "force.h"
#include "respa.h"
#include "update.h"

#include <iostream>
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixNVE::FixNVE(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (!utils::strmatch(style,"^nve/sphere") && narg < 3)
    utils::missing_cmd_args(FLERR, "fix nve", error);

  dynamic_group_allow = 1;
  time_integrate = 1;

  px_file.open("./../out/data/px_nve.out");
  py_file.open("./../out/data/py_nve.out");
  pz_file.open("./../out/data/pz_nve.out");
}

/* ---------------------------------------------------------------------- */

FixNVE::~FixNVE() {
  px_file.close();
  py_file.close();
  pz_file.close();
}

/* ---------------------------------------------------------------------- */

int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= FINAL_INTEGRATE;
  mask |= INITIAL_INTEGRATE_RESPA;
  mask |= FINAL_INTEGRATE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixNVE::init()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;

  if (utils::strmatch(update->integrate_style,"^respa"))
    step_respa = (dynamic_cast<Respa *>(update->integrate))->step;

}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixNVE::initial_integrate(int /*vflag*/)
{
  std::cout<<"FixNVE::initial_integrate"<<"\n";
  double dtfm;

  // update v and x of atoms in group

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;

  double linear_momentum_iii[3];
  double linear_momentum_all_iii[3];
  MathExtra::zero3(linear_momentum_iii);
  MathExtra::zero3(linear_momentum_all_iii);
  for (int i = 0; i < nlocal; i++) {
    if (atom->tag[i] == 226) {
        std::cout<<"atom->tag[j]: "<<atom->tag[i]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
      }
    linear_momentum_iii[0] += mass[type[i]] * v[i][0];
    linear_momentum_iii[1] += mass[type[i]] * v[i][1];
    linear_momentum_iii[2] += mass[type[i]] * v[i][2];
  }

  std::cout<<"FixNVE::initial_integrate ::BEFORE update of velocities by a half-step and coordinates by a full step"<<"\n";
  std::cout<<"nlocal: "<<nlocal<<"\n";
  std::cout<<"linear_momentum_iii[0]: "<<linear_momentum_iii[0]<<"\n";
  std::cout<<"linear_momentum_iii[1]: "<<linear_momentum_iii[1]<<"\n";
  std::cout<<"linear_momentum_iii[2]: "<<linear_momentum_iii[2]<<"\n"; 
  MPI_Allreduce(linear_momentum_iii, linear_momentum_all_iii, 3, MPI_DOUBLE, MPI_SUM, world); 
  std::cout<<"linear_momentum_all_iii[0]: "<<linear_momentum_all_iii[0]<<"\n";
  std::cout<<"linear_momentum_all_iii[1]: "<<linear_momentum_all_iii[1]<<"\n";
  std::cout<<"linear_momentum_all_iii[2]: "<<linear_momentum_all_iii[2]<<"\n"; 

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }

  } else {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
        x[i][0] += dtv * v[i][0];
        x[i][1] += dtv * v[i][1];
        x[i][2] += dtv * v[i][2];
      }
  }

  double linear_momentum_ii[3];
  double linear_momentum_all_ii[3];
  MathExtra::zero3(linear_momentum_ii);
  MathExtra::zero3(linear_momentum_all_ii);
  for (int i = 0; i < nlocal; i++) {
    if (atom->tag[i] == 226) {
        std::cout<<"atom->tag[j]: "<<atom->tag[i]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
      }
    linear_momentum_ii[0] += mass[type[i]] * v[i][0];
    linear_momentum_ii[1] += mass[type[i]] * v[i][1];
    linear_momentum_ii[2] += mass[type[i]] * v[i][2];
  }

  std::cout<<"FixNVE::initial_integrate :: AFTER update of velocities by a half-step and coordinates by a full step"<<"\n";
  std::cout<<"nlocal: "<<nlocal<<"\n";
  std::cout<<"linear_momentum_ii[0]: "<<linear_momentum_ii[0]<<"\n";
  std::cout<<"linear_momentum_ii[1]: "<<linear_momentum_ii[1]<<"\n";
  std::cout<<"linear_momentum_ii[2]: "<<linear_momentum_ii[2]<<"\n"; 
  MPI_Allreduce(linear_momentum_ii, linear_momentum_all_ii, 3, MPI_DOUBLE, MPI_SUM, world); 
  std::cout<<"linear_momentum_all_ii[0]: "<<linear_momentum_all_ii[0]<<"\n";
  std::cout<<"linear_momentum_all_ii[1]: "<<linear_momentum_all_ii[1]<<"\n";
  std::cout<<"linear_momentum_all_ii[2]: "<<linear_momentum_all_ii[2]<<"\n"; 

  // PPapez COMMENT: write to file
  px_file<<linear_momentum_all_ii[0]<<"\n";
  px_file.flush();
  py_file<<linear_momentum_all_ii[1]<<"\n";
  py_file.flush();
  pz_file<<linear_momentum_all_ii[2]<<"\n";
  pz_file.flush();

  std::cout<<"FixNVE::initial_integrate END"<<"\n";
}

/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate()
{ 
  std::cout<<"FixNVE::final_integrate"<<"\n";
  std::cout<<"dtf: "<<dtf<<"\n";

  double **v = atom->v;
  double **f = atom->f;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  if (igroup == atom->firstgroup) nlocal = atom->nfirst;  

  double linear_momentum_fii[3];
  MathExtra::zero3(linear_momentum_fii);
  for (int i = 0; i < nlocal; i++) {
    if (tag[i] == 226) {
      std::cout<<"v[i][0]: "<<v[i][0]<<"\n";
      std::cout<<"v[i][1]: "<<v[i][1]<<"\n";
      std::cout<<"v[i][2]: "<<v[i][2]<<"\n";
      std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
    }
    linear_momentum_fii[0] += mass[type[i]] * v[i][0];
    linear_momentum_fii[1] += mass[type[i]] * v[i][1];
    linear_momentum_fii[2] += mass[type[i]] * v[i][2];
  }

  std::cout<<"FixNVE::final_integrate :: BEFORE half-step update of the velocities"<<"\n";
  std::cout<<"nlocal: "<<nlocal<<"\n";
  std::cout<<"linear_momentum_fii[0]: "<<linear_momentum_fii[0]<<"\n";
  std::cout<<"linear_momentum_fii[1]: "<<linear_momentum_fii[1]<<"\n";
  std::cout<<"linear_momentum_fii[2]: "<<linear_momentum_fii[2]<<"\n"; 

  double dtfm;

  // update v of atoms in group

  if (rmass) {
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        dtfm = dtf / rmass[i];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }

  } else {
    // std::cout<<"in fi ::: else"<<"\n";
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {
        // std::cout<<"in fi ::: else :: if"<<"\n";
        // std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
        dtfm = dtf / mass[type[i]];
        v[i][0] += dtfm * f[i][0];
        v[i][1] += dtfm * f[i][1];
        v[i][2] += dtfm * f[i][2];
      }
  }

  double linear_momentum_fi[3];
  MathExtra::zero3(linear_momentum_fi);

  double sumforces[3];
  MathExtra::zero3(sumforces);

  for (int i = 0; i < nlocal; i++) {
    if (tag[i] == 226) {
      std::cout<<"v[i][0]: "<<v[i][0]<<"\n";
      std::cout<<"v[i][1]: "<<v[i][1]<<"\n";
      std::cout<<"v[i][2]: "<<v[i][2]<<"\n";
      std::cout<<"f[i][0]: "<<f[i][0]<<"\n";
    }

    sumforces[0] += f[i][0];
    sumforces[1] += f[i][1];
    sumforces[2] += f[i][2];

    linear_momentum_fi[0] += mass[type[i]] * v[i][0];
    linear_momentum_fi[1] += mass[type[i]] * v[i][1];
    linear_momentum_fi[2] += mass[type[i]] * v[i][2];
  }

  std::cout<<"FixNVE::final_integrate :: AFTER half-step update of the velocities"<<"\n";
  std::cout<<"nlocal: "<<nlocal<<"\n";
  std::cout<<"linear_momentum_fi[0]: "<<linear_momentum_fi[0]<<"\n";
  std::cout<<"linear_momentum_fi[1]: "<<linear_momentum_fi[1]<<"\n";
  std::cout<<"linear_momentum_fi[2]: "<<linear_momentum_fi[2]<<"\n"; 
  std::cout<<"sumforces[0]: "<<sumforces[0]<<"\n";
  std::cout<<"sumforces[1]: "<<sumforces[1]<<"\n";
  std::cout<<"sumforces[2]: "<<sumforces[2]<<"\n";

  std::cout<<"FixNVE::final_integrate END"<<"\n";
}

/* ---------------------------------------------------------------------- */

void FixNVE::initial_integrate_respa(int vflag, int ilevel, int /*iloop*/)
{
  dtv = step_respa[ilevel];
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;

  // innermost level - NVE update of v and x
  // all other levels - NVE update of v

  if (ilevel == 0) initial_integrate(vflag);
  else final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVE::final_integrate_respa(int ilevel, int /*iloop*/)
{
  dtf = 0.5 * step_respa[ilevel] * force->ftm2v;
  final_integrate();
}

/* ---------------------------------------------------------------------- */

void FixNVE::reset_dt()
{
  dtv = update->dt;
  dtf = 0.5 * update->dt * force->ftm2v;
}
