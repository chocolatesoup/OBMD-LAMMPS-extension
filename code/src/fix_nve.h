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

#ifdef FIX_CLASS
// clang-format off
FixStyle(nve,FixNVE);
// clang-format on
#else

#ifndef LMP_FIX_NVE_H
#define LMP_FIX_NVE_H

#include "fix.h"
#include <fstream>

namespace LAMMPS_NS {

class FixNVE : public Fix {
 public:
  FixNVE(class LAMMPS *, int, char **);
  ~FixNVE();

  int setmask() override;
  void init() override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void initial_integrate_respa(int, int, int) override;
  void final_integrate_respa(int, int) override;
  void reset_dt() override;

 private:
   double linear_momentum[3];

   std::ofstream px_file;
   std::ofstream py_file;
   std::ofstream pz_file;
   
 protected:
  double dtv, dtf;
  double *step_respa;
  int mass_require;
};

}    // namespace LAMMPS_NS

#endif
#endif
