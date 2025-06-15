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

#include "fix_tdforce.h"

#include "domain.h"
#include "error.h"
#include "input.h"
#include "memory.h" 
#include "lattice.h"
#include "modify.h"
#include "respa.h"
#include "update.h"
#include "variable.h"
#include "atom.h"
#include "tokenizer.h" // split
#include "table_file_reader.h"
#include "fix_adress.h"

#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

enum {LINEAR, SPLINE };
#define BIG MAXTAGINT

/* ---------------------------------------------------------------------- */

FixTDForce::FixTDForce(LAMMPS *lmp, int narg, char **arg) : 
    Fix(lmp, narg, arg), ntables(0), tables(nullptr) 
{ 

  x_lo = domain->boxlo[0];
  x_hi = domain->boxhi[0];
  // not used
  y_lo = domain->boxlo[1];
  y_hi = domain->boxhi[1];
  z_lo = domain->boxlo[2];
  z_hi = domain->boxhi[2];

  center_box[0] = (x_lo + x_hi) / 2.0;
  //not used
  center_box[1] = (y_lo + y_hi) / 2.0;
  center_box[2] = (z_lo + z_hi) / 2.0;

  if (narg != 7) error->all(FLERR, "Illegal fix tdforce command");
     
// tabstyle
  if (strcmp(arg[3], "linear") == 0) {
    tabstyle = LINEAR;
  }
  else if (strcmp(arg[3], "spline") == 0) {
    tabstyle = SPLINE;
  }
  else {
    error->all(FLERR, "Unknown table style {} in fix {}", arg[3], style);
  }

  // number of entries
  tablength = utils::inumeric(FLERR, arg[4], false, lmp);
  if (tablength < 2) error->all(FLERR, "Illegal number of fix {} table entries", style);

  int me;
  tables = nullptr;
  MPI_Comm_rank(world, &me);
  tables = (Table *)
  memory->srealloc(tables,(ntables+1)*sizeof(Table), "tdforce:tables"); // srealloc ! ! ! free ! ! !
  Table *tb = &tables[ntables]; 
  null_table(tb); 
  if (me == 0) read_table(tb, arg[5], arg[6]); // name of a table name of force to read
  bcast_table(tb);

  if (tb->ninput <= 1) {
    error->all(FLERR, "Invalid fix {} table length: {}", style, tb->ninput);
  }
      
  rlo = tb->rfile[0]; 
  rhi = tb->rfile[tb->ninput - 1];

  spline_table(tb);
  compute_table(tb);
  ntables++;
}

/* ---------------------------------------------------------------------- */

FixTDForce::~FixTDForce()
{
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);
}

/* ---------------------------------------------------------------------- */

int FixTDForce::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTDForce::init()
{
  ;
}

/* ---------------------------------------------------------------------- */

void FixTDForce::setup(int vflag)
{
  post_force(vflag);
} 

/* ---------------------------------------------------------------------- */

/* void FixTDForce::min_setup(int vflag)
{
  post_force(vflag);
} */

/* ---------------------------------------------------------------------- */

/* void FixTDForce::pre_exchange()
{
  // nothing in pre_exchange()  ... can be removed 
  ;
} */

/* ---------------------------------------------------------------------- */
// PP ... accounts external forces other than bonded and non-bonded forces and adds them to the force vectors of selected atoms
void FixTDForce::post_force(int vflag)
{
  // std::cout<<"fix tdforce  [post_force()]"<<"\n";
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  int *rep_atom = atom->rep_atom;
  double **cms = atom->cms_mol; // ? ? ? 
  double *rmass = atom->rmass;  
  int *mask = atom->mask;
  tagint *tag = atom->tag;

  int i;
  tagint imol;
  double massone;

  /* added */
  nmolecules = find_mols(idlo, idhi); 
  // instead find maxid ... back to find_mols where I'm setting maxmol_all to idhi
  // find_maxid();
  maxmol_all = idhi; 

  // mass
  double *massproc_tmp, *masstotal_tmp;
  memory->create(massproc_tmp,maxmol_all,"tdforce:massproc_tmp");
  memory->create(masstotal_tmp,maxmol_all,"tdforce:masstotal_tmp");
  // fill with zeros
  for (i = 0; i < maxmol_all; i++) {
        massproc_tmp[i] = 0.0;
  }

  // compute masses
  for (i = 0; i < nlocal; i++) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else {
        massone = mass[type[i]];
      }

      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;
      massproc_tmp[imol] += massone;

      /* if (imol == 3515) {
        std::cout<<"0"<<"\n";
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"massone: "<<massone<<"\n";
        std::cout<<"mass[type[i]]: "<<mass[type[i]]<<"\n";
        std::cout<<"massproc_tmp[imol]: "<< massproc_tmp[imol]<<"\n";
      }

      if (imol == 3516) {
        std::cout<<"0"<<"\n";
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"massone: "<<massone<<"\n";
        std::cout<<"mass[type[i]]: "<<mass[type[i]]<<"\n";
        std::cout<<"massproc_tmp[imol]: "<< massproc_tmp[imol]<<"\n";
      } */

      // PP ... just testing
      /* double **x = atom->x;
      tagint *tag = atom->tag;
      int G = tag[i];
      if (tag[i] == 10546) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      }

      if (tag[i] == 10547) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      }

      if (tag[i] == 10548) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      }

      if (tag[i] == 10549) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      }

      if (tag[i] == 10550) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      }

      if (tag[i] == 10551) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      }

      if (tag[i] == 7386) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      } */

    }
  }

  MPI_Allreduce(massproc_tmp, masstotal_tmp, maxmol_all, MPI_DOUBLE, MPI_SUM, world); // sum \& distribute

  /* std::cout<<"\n";
  std::cout<<"masstotal_tmp[1641]: "<<masstotal_tmp[1641]<<"\n";
  std::cout<<"masstotal_tmp[1642]: "<<masstotal_tmp[1642]<<"\n"; */

  // apply thermodynamic force
  double xtmp;
  double delx;
  double mol_mass, mass_frac;
  double dist, td_force;

  double **mol_f_td, **mol_f_all_td;
  memory->create(mol_f_td, maxmol_all, 3, "tdforce:mol_f_td");
  memory->create(mol_f_all_td, maxmol_all, 3, "tdforce:mol_f_all_td");
  // fill with zeros
  for (int i = 0; i < maxmol_all; i++) {
    mol_f_td[i][0] = mol_f_td[i][1] = mol_f_td[i][2] = 0.0; 
    mol_f_all_td[i][0] = mol_f_all_td[i][1] = mol_f_all_td[i][2] = 0.0;
  }

  // loop over particles \& find cms
  for (int i = 0; i < nlocal; i++) {
    if (rep_atom[i] == 0) {
      continue; // not rep atom ... using only rep atoms ! ! ! 
    }

    // PP ... just testing
    /* double **x = atom->x;
    tagint *tag = atom->tag;
    int G = tag[i];
    if (tag[i] == 7386) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
      std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
    } */

    /* if (tag[i] == 7601) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 7602) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    } 

    if (tag[i] == 7471) {
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 7472) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 7473) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    } 

    if (tag[i] == 2488) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 2491) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 2492) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    } 

    if (tag[i] == 124) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 125) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    }

    if (tag[i] == 126) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
    } */     

    /* --- */
    imol = molecule[i];
    if (molmap_tmp) imol = molmap_tmp[imol-idlo];
    else imol--;

    xtmp = cms[i][0];
    delx = xtmp - center_box[0];
    dist = std::abs(delx);
    // TD force as provided
    force_lookup(dist,td_force); // get force
    td_force /= dist;
    mol_f_td[imol][0] += delx * td_force;

    /* if (imol == 1641) {
        std::cout<<"1"<<"\n";
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"xtmp: "<<xtmp<<"\n";
        std::cout<<"mol_f_td[imol][0]: "<<mol_f_td[imol][0]<<"\n";
    } */
  }

  MPI_Allreduce(&mol_f_td[0][0], &mol_f_all_td[0][0], 3*maxmol_all, MPI_DOUBLE, MPI_SUM, world); // sum \& distribute

  // distribute among particles
  for (int i = 0; i < nlocal; i++) {
    imol = molecule[i];
    if (molmap_tmp) imol = molmap_tmp[imol-idlo]; // important ! ! !
    else imol--;

    mol_mass = masstotal_tmp[imol];
    mass_frac = mass[type[i]] / mol_mass;
    
    f[i][0] += mass_frac * mol_f_all_td[imol][0];

    /* if (imol == 3515) {
        std::cout<<"2"<<"\n";
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"mol_mass: "<<mol_mass<<"\n";
        std::cout<<"mol_f_all_td[imol][0]: "<<mol_f_all_td[imol][0]<<"\n";
    } 

    if (imol == 3516) {
        std::cout<<"2"<<"\n";
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"mol_mass: "<<mol_mass<<"\n";
        std::cout<<"mol_f_all_td[imol][0]: "<<mol_f_all_td[imol][0]<<"\n";
    } */

    // PP ... just testing
    /* double **x = atom->x;
    tagint *tag = atom->tag;
    int G = tag[i];
    if (tag[i] == 7386) {
      std::cout<<"tag[i]: "<<tag[i]<<"\n";
      std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
      std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
      std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
      std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
      std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
      std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
      std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
    } */
    
  }

  // free
  memory->destroy(massproc_tmp);
  memory->destroy(masstotal_tmp);
  memory->destroy(mol_f_td);
  memory->destroy(mol_f_all_td);

  // std::cout<<"fix tdforce  [post_force() ::: END]"<<"\n";
} 

/* ---------------------------------------------------------------------- */

// read table with TDForce data
void FixTDForce::read_table(Table *tb, char *file, char *keyword)
{

  TableFileReader reader(lmp, file, "fix tdforce");

  char *line = reader.find_section_start(keyword); 
  if (!line) {
    error->one(FLERR, "Dod not find keyword {} in table file", keyword);
  }

  line = reader.next_line();
  param_extract(tb, line);
  memory->create(tb->rfile, tb->ninput, "tdforce:rfile"); // PP ... distances
  memory->create(tb->ffile, tb->ninput, "tdforce:ffile"); // PP ... tdforce

  // read r and ftd from file
  reader.skip_line();
  for (int i = 0; i < tb->ninput; i++) {
    line = reader.next_line();

    if (!line) {
      error->one(FLERR, "Data missing when parsing table '{}' line {} of {}", keyword, i+1, tb->ninput);
    }
    try {
      ValueTokenizer values(line);
      values.next_int();
      tb->rfile[i] = values.next_double();
      tb->ffile[i] = values.next_double();
    }
    catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }
  }
}

/* ---------------------------------------------------------------------- */

// extract attributes from parameter line in table section (N value)
void FixTDForce::param_extract(Table *tb, char *line)
{
    // to set tb->ninput
    tb->ninput = 0;
    
    try {
      ValueTokenizer values(line);

      while(values.has_next()) {
        std::string word = values.next_string();

        if (word == "N") {
          tb->ninput = values.next_int();
        }
        else {
          error->one(FLERR, "Invalid keyword {} in fix {} table parameters", word, style);
        }
      }
    }

    catch (TokenizerException &e) {
      error->one(FLERR, e.what());
    }

    if (tb->ninput == 0) {
      error->one(FLERR, "Fix {} param_extract did not set N");
    }
}

/* ---------------------------------------------------------------------- */

// read table with TDForce data
void FixTDForce::null_table(Table *tb)
{
  tb->rfile = tb->ffile = nullptr;
  tb->f2file = nullptr;

  tb->r = tb->f = tb->df = tb->f2 = nullptr;
}

/* ---------------------------------------------------------------------- */

// broadcast read-in table info from proc 0 to other procs
void FixTDForce::bcast_table(Table *tb)
{
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);

  int me;
  MPI_Comm_rank(world, &me);
  if (me > 0) {
    memory->create(tb->rfile, tb->ninput, "tdforce:rfile");
    memory->create(tb->ffile, tb->ninput, "tdforce:ffile");
  }

  MPI_Bcast(tb->rfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

// read table with TDForce data
void FixTDForce::free_table(Table *tb)
{
  // memory management
  memory->destroy(tb->rfile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->f2file);

  memory->destroy(tb->f);
  memory->destroy(tb->df);
  memory->destroy(tb->f2);
}

/* ---------------------------------------------------------------------- */

// build spline representation
void FixTDForce::spline_table(Table *tb)
{
  memory->create(tb->f2file, tb->ninput, "tdforce:f2file");
  double fp0 = tb->fplo;
  double fpn = tb->fphi;
  spline(tb->rfile, tb->ffile, tb->ninput, fp0, fpn, tb->f2file);
}

/* ---------------------------------------------------------------------- */

// spline routine (interpolation)
void FixTDForce::spline(double *x, double *y, int n, double yp1, double ypn, double *y2)
{
  int i,k;
  double p,qn,sig,un;
  auto u = new double[n];

  if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
  else {
    y2[0] = -0.5;
    u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0]) / (x[1]-x[0]) - yp1);
  }
  for (i = 1; i < n-1; i++) {
    sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
    p = sig*y2[i-1] + 2.0;
    y2[i] = (sig-1.0) / p;
    u[i] = (y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i]-y[i-1]) / (x[i]-x[i-1]);
    u[i] = (6.0*u[i] / (x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }
  if (ypn > 0.99e30) qn = un = 0.0;
  else {
    qn = 0.5;
    un = (3.0/(x[n-1]-x[n-2])) * (ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]));
  }
  y2[n-1] = (un-qn*u[n-2]) / (qn*y2[n-2] + 1.0);
  for (k = n-2; k >= 0; k--) y2[k] = y2[k]*y2[k+1] + u[k];

  delete [] u;
}

/* ---------------------------------------------------------------------- */

// compute f vectors from splined values
void FixTDForce::compute_table(Table *tb)
{

  int tlm1 = tablength-1;
  tb->delta = (rhi - rlo) / tlm1; 
  tb->invdelta = 1.0 / tb->delta;
  tb->deltasq6 = tb->delta * tb->delta / 6.0;
  
  // N-1 evenly spaced bins in r from min to max
  // r, f values = value at lower edge of bin
  // df value = delta value of f
  // r, f are N in length so df can compute diff

  memory->create(tb->r, tablength, "tdforce:r");
  memory->create(tb->f, tablength, "tdforce:f");
  memory->create(tb->df, tablength, "tdforce:df");
  memory->create(tb->f2, tablength, "tdforce:f2");

  double a;
  for (int i = 0; i < tablength; i++) {
    a = i*tb->delta;
    tb->r[i] = a;
    tb->f[i] = splint(tb->rfile, tb->ffile, tb->f2file, tb->ninput, a);
  }

  // diff
  for (int i = 0; i < tlm1; i++) {
    tb->df[i] = tb->f[i+1] - tb->f[i];
  }

  tb->df[tlm1] = 2.0 * tb->df[tlm1 - 1] - tb->df[tlm1 - 2];

  spline(tb->r, tb->f, tablength, tb->fplo, tb->fphi, tb->f2);
}

/* ---------------------------------------------------------------------- */

// splint 
double FixTDForce::splint(double *xa, double *ya, double *y2a, int n, double x)
{
  int klo, khi, k;
  double h, b, a, y;

  klo = 0;
  khi = n - 1;
  while (khi - klo > 1) {
    k = (khi + klo) >> 1;
    if (xa[k] > x)
      khi = k;
    else
      klo = k;
  }
  h = xa[khi] - xa[klo];
  a = (xa[khi] - x) / h;
  b = (x - xa[klo]) / h;
  y = a * ya[klo] + b * ya[khi] +
      ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * (h * h) / 6.0;
  return y;

}

/* ---------------------------------------------------------------------- */

// calculate force f at distance x
void FixTDForce::force_lookup(double x, double &f) 
{
  double fraction, a, b;
  const Table *tb=&tables[0];
  int itable = static_cast<int> (x * tb->invdelta);

  if (tabstyle == LINEAR) {
    fraction = (x - tb->r[itable]) * tb->invdelta;
    f = tb->f[itable] + fraction * tb->df[itable];
  }
  else if (tabstyle == SPLINE) {
    b = (x - tb->r[itable]) * tb->invdelta;
    a = 1.0 -b;

    f = a * tb->f[itable] + b * tb->f[itable+1] + ((a * a * a - a) * tb->f2[itable] + (b * b * b - b) * tb->f2[itable+1]) * tb->deltasq6;
  } 
}

/* ---------------------------------------------------------------------- */

// find molecules
int FixTDForce::find_mols(tagint &idlo, tagint &idhi)
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

/* ---------------------------------------------------------------------- */

// void FixTDForce::min_post_force(int /*vflag*/)
// {
//   post_force(/*vflag*/);
// }
// PP ... not used in minimization procedure ! ! !

/* ----------------------------------------------------------------------
   maxmol_all = current max molecule ID for all atoms
------------------------------------------------------------------------- */

// PP using find_maxid()
/* void FixTDForce::find_maxid()
{ 
  
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  int max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,molecule[i]);
  MPI_Allreduce(&max,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
} */
