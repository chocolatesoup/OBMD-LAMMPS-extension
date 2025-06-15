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
FixStyle(adress,FixAdResS);
// clang-format on
#else

#ifndef LMP_FIX_ADRESS_H
#define LMP_FIX_ADRESS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdResS : public Fix {
 public:
  FixAdResS(class LAMMPS *, int, char **);
  ~FixAdResS() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_integrate() override;
    
 protected:
  double *step_respa;
  int mass_require;
    
  double len_AT;
  double len_HY;
  int flag_inverse;
  int topo;
  tagint idlo,idhi;
  int nmolecules;
  int *molmap_tmp;
  
  char *id_c;
  
  void options(int, char **);

  int find_mols(tagint &, tagint &);

  /* tagint maxtag_all, maxmol_all;
  void find_maxid(); */ // PP
  tagint maxmol_all;

  double x0lo,x0hi,x1lo,x1hi,x2lo,x2hi;
  double len_CG;
  double center_box[3];
};

}    // namespace LAMMPS_NS

#endif
#endif
