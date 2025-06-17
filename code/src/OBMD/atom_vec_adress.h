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

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(adress,AtomVecAdResS);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_ADRESS_H
#define LMP_ATOM_VEC_ADRESS_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecAdResS : public AtomVec {
 public:
  AtomVecAdResS(class LAMMPS *);
  ~AtomVecAdResS() override;

  void grow_pointers() override;
  void pack_restart_pre(int) override;
  void pack_restart_post(int) override;
  void unpack_restart_init(int) override;
  void data_atom_post(int) override;

 private:
  int *num_bond, *num_angle, *num_dihedral, *num_improper;
  int **bond_type, **angle_type, **dihedral_type, **improper_type;
  int **nspecial;

  int any_bond_negative, any_angle_negative, any_dihedral_negative, any_improper_negative;
  int bond_per_atom, angle_per_atom, dihedral_per_atom, improper_per_atom;
  int *bond_negative, *angle_negative, *dihedral_negative, *improper_negative;
  
  double *lambdaF;
  double **cms_mol;
  int *rep_atom;

  double **vcms_mol; // PPapez COMMENT: added
};

}    // namespace LAMMPS_NS

#endif
#endif
