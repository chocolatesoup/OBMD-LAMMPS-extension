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

#include "pair_morse_adress_cg_ipc_org.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "domain.h"

#include <cmath>
#include <cstring>

#include <iostream>

using namespace LAMMPS_NS;
#define BIG MAXTAGINT

/* ---------------------------------------------------------------------- */

PairMorseAdResSCGIPCOrg::PairMorseAdResSCGIPCOrg(LAMMPS *lmp) : Pair(lmp)
{
  writedata = 1;

  massproc_tmp = nullptr;
  masstotal_tmp = nullptr;
  molmap_tmp = nullptr;
  mol_f_adr = nullptr;
  mol_f_all_adr = nullptr;

  nmolecules = find_mols(idlo, idhi);

  memory->create(massproc_tmp, nmolecules, "morse/adress/cg/ipc/org:massproc_tmp");
  memory->create(masstotal_tmp, nmolecules, "morse/adress/cg/ipc/org:masstotal_tmp");
  memory->create(mol_f_adr, nmolecules, 3, "morse/adress/cg/ipc/org:mol_f_adr"); // vector
  memory->create(mol_f_all_adr, nmolecules, 3, "morse/adress/cg/ipc/org:mol_f_all_adr"); // vector

  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;  
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  tagint *molecule = atom->molecule;
  int imol; // 
  double massone;

  for (int i = 0; i < nmolecules; i++) {
    massproc_tmp[i] = 0.0;
  }

  for (int i = 0; i < nlocal; i++) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];

      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;
      massproc_tmp[imol] += massone;
    }
  }

  MPI_Allreduce(massproc_tmp, masstotal_tmp, nmolecules, MPI_DOUBLE, MPI_SUM, world);

}

/* ---------------------------------------------------------------------- */

PairMorseAdResSCGIPCOrg::~PairMorseAdResSCGIPCOrg()
{
  if (copymode) return;

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);

    memory->destroy(cut);
    memory->destroy(d0);
    memory->destroy(alpha);
    memory->destroy(r0);
    memory->destroy(morse1);
    memory->destroy(offset);
  }
}

/* ---------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq, r, dr, dexp, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // AdResS
  tagint *molecule = atom->molecule;
  int imol, jmol;
  int *rep_atom = atom->rep_atom; // representative atom marked using 1
  double *lambda = atom->lambdaF;
  double ilambda, jlambda, lambda_factor;
  double **cms = atom->cms_mol; // vector
  double s_factor;

  for (int i = 0; i < nmolecules; i++) {
    mol_f_adr[i][0] = mol_f_adr[i][1] = mol_f_adr[i][2] = 0.0;
    mol_f_all_adr[i][0] = mol_f_all_adr[i][1] = mol_f_all_adr[i][2] = 0.0;
  }

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];

    if (rep_atom[i] == 0) {
      continue;
    }

    ilambda = lambda[i];

    xtmp = cms[i][0];  // x[i][0];
    ytmp = cms[i][1]; // x[i][1];
    ztmp = cms[i][2]; // x[i][2];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;
      jtype = type[j];

      if (rep_atom[j] == 0) {
        continue;
      }

      jlambda = lambda[j];

      delx = xtmp - cms[j][0]; // x[j][0];
      dely = ytmp - cms[j][1]; // x[j][1];
      delz = ztmp - cms[j][2]; // x[j][2];
      domain->minimum_image(delx, dely, delz); // ! ! !

      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq < cutsq[itype][jtype]) {
        // find corresponding molecules
        imol = molecule[i];
        jmol = molecule[j];
        if (molmap_tmp) {
          imol = molmap_tmp[imol-idlo];
          jmol = molmap_tmp[jmol-idlo];
        }
        else {
          imol--; // ! ! ! 
          jmol--; // ! ! ! 
        }

        r = sqrt(rsq);
        dr = r - r0[itype][jtype];
        dexp = exp(-alpha[itype][jtype] * dr);
        fpair = factor_lj * morse1[itype][jtype] * (dexp * dexp - dexp) / r;

        // AdResS 
        lambda_factor = 1.0 - ilambda * jlambda;
        s_factor = 4.0 * pow(sqrt(ilambda*jlambda) - 0.5,2.0);

        /* std::cout<<"STARTING CG HERE"<<"\n";
        std::cout<<"morse/adress/cg/test: cms[i][0]: "<<cms[i][0]<<" cms[j][0]: "<<cms[j][0]<<"\n";
        std::cout<<"morse/adress/cg/test: xtmp: "<<xtmp<<" x[j][0]: "<<x[j][0]<<" sqrt(rsq): "<<sqrt(rsq)<<" cutsq[itype][jtype]: "<<cutsq[itype][jtype]<<"\n";
        std::cout<<"morse/adress/cg/test: ilambda: "<<ilambda<<" jlambda: "<<jlambda<<"\n";
        std::cout<<"morse/adress/cg/test: lambda_factor: "<<lambda_factor<<"\n";
        std::cout<<"rep_atom[i]: "<<rep_atom[i]<<" rep_atom[j]: "<<rep_atom[j]<<"\n";
        tagint *tag = atom->tag;
        int Gi = tag[i];
        int Gj = tag[j];
        std::cout<<"morse/adress/cg/test: Gi: "<<Gi<<" Gj: "<<Gj<<"\n"; */
        
        // std::cout<<"morse/adress/cg/test: before: mol_f_adr[imol][0]: "<<mol_f_adr[imol][0]<<"\n";
        // std::cout<<"morse/adress/cg/test: delx: "<<delx<<" fpair: "<<fpair<<"\n",
        mol_f_adr[imol][0] += lambda_factor * s_factor * delx * fpair; // delx * fpair;
        mol_f_adr[imol][1] += lambda_factor * s_factor * dely * fpair; // dely * fpair;
        mol_f_adr[imol][2] += lambda_factor * s_factor * delz * fpair; // delz * fpair;
        // std::cout<<"morse/adress/cg/test: after: mol_f_adr[imol][0]: "<<mol_f_adr[imol][0]<<"\n";
        // std::cout<<"ENDING CG HERE"<<"\n";
        if (newton_pair || j < nlocal) {
          mol_f_adr[jmol][0] -= lambda_factor * s_factor * delx * fpair; // delx * fpair;
          mol_f_adr[jmol][1] -= lambda_factor * s_factor * dely * fpair; // dely * fpair;
          mol_f_adr[jmol][2] -= lambda_factor * s_factor * delz * fpair; // delz * fpair;
        }

        if (eflag) {
          evdwl = d0[itype][jtype] * (dexp * dexp - 2.0 * dexp) - offset[itype][jtype];
          evdwl *= factor_lj;
        }

        if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, 0.0, fpair, delx, dely, delz);
      }
    }
  }

  MPI_Allreduce(&mol_f_adr[0][0], &mol_f_all_adr[0][0], 3*nmolecules, MPI_DOUBLE, MPI_SUM, world);
  
  // AdResS
  // distribute force among the explicit atoms of molecule 
  double mol_mass, mass_frac;
  for (int i = 0; i < nlocal; i++) {
    imol = molecule[i];
    if (molmap_tmp) imol = molmap_tmp[imol-idlo]; 
    else imol--; // ! ! ! 

    // std::cout<<"morse/adress/cg/test: x[i][0]: "<<x[i][0]<<" cms[i][0]: "<<cms[i][0]<<"\n";

    mol_mass = masstotal_tmp[imol];
    mass_frac = mass[type[i]] / mol_mass;

    f[i][0] += mass_frac * mol_f_all_adr[imol][0];
    f[i][1] += mass_frac * mol_f_all_adr[imol][1];
    f[i][2] += mass_frac * mol_f_all_adr[imol][2];

    if (f[i][0] > 1000 || f[i][1] > 1000 || f[i][2] > 1000) { // check
      tagint *tag = atom->tag;
      int G = tag[i];
      std::cout<<"morse/adress/cg/test ::: JUST GO"<<"\n";
      std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      std::cout<<"index: "<<G<<"\n"; // index as in input.data
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  for (int i = 1; i < np1; i++)
    for (int j = i; j < np1; j++) setflag[i][j] = 0;

  memory->create(cutsq, np1, np1, "pair:cutsq");

  memory->create(cut, np1, np1, "pair:cut");
  memory->create(d0, np1, np1, "pair:d0");
  memory->create(alpha, np1, np1, "pair:alpha");
  memory->create(r0, np1, np1, "pair:r0");
  memory->create(morse1, np1, np1, "pair:morse1");
  memory->create(offset, np1, np1, "pair:offset");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Illegal pair_style command");

  cut_global = utils::numeric(FLERR, arg[0], false, lmp);

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::coeff(int narg, char **arg)
{
  if (narg < 5 || narg > 6) error->all(FLERR, "Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  double d0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double alpha_one = utils::numeric(FLERR, arg[3], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[4], false, lmp);

  double cut_one = cut_global;
  if (narg == 6) cut_one = utils::numeric(FLERR, arg[5], false, lmp);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo, i); j <= jhi; j++) {
      d0[i][j] = d0_one;
      alpha[i][j] = alpha_one;
      r0[i][j] = r0_one;
      cut[i][j] = cut_one;
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMorseAdResSCGIPCOrg::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  morse1[i][j] = 2.0 * d0[i][j] * alpha[i][j];

  if (offset_flag) {
    double alpha_dr = -alpha[i][j] * (cut[i][j] - r0[i][j]);
    offset[i][j] = d0[i][j] * (exp(2.0 * alpha_dr) - 2.0 * exp(alpha_dr));
  } else
    offset[i][j] = 0.0;

  d0[j][i] = d0[i][j];
  alpha[j][i] = alpha[i][j];
  r0[j][i] = r0[i][j];
  morse1[j][i] = morse1[i][j];
  offset[j][i] = offset[i][j];

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::write_restart(FILE *fp)
{
  /* write_restart_settings(fp);

  int i, j;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j], sizeof(int), 1, fp);
      if (setflag[i][j]) {
        fwrite(&d0[i][j], sizeof(double), 1, fp);
        fwrite(&alpha[i][j], sizeof(double), 1, fp);
        fwrite(&r0[i][j], sizeof(double), 1, fp);
        fwrite(&cut[i][j], sizeof(double), 1, fp);
      }
    }
  } */
 ;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::read_restart(FILE *fp)
{
  /* read_restart_settings(fp);

  allocate();

  int i, j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++) {
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR, &setflag[i][j], sizeof(int), 1, fp, nullptr, error);
      MPI_Bcast(&setflag[i][j], 1, MPI_INT, 0, world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR, &d0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &alpha[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &r0[i][j], sizeof(double), 1, fp, nullptr, error);
          utils::sfread(FLERR, &cut[i][j], sizeof(double), 1, fp, nullptr, error);
        }
        MPI_Bcast(&d0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&alpha[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&r0[i][j], 1, MPI_DOUBLE, 0, world);
        MPI_Bcast(&cut[i][j], 1, MPI_DOUBLE, 0, world);
      }
    }
  } */
 ;
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::write_restart_settings(FILE *fp)
{
  /* fwrite(&cut_global, sizeof(double), 1, fp);
  fwrite(&offset_flag, sizeof(int), 1, fp);
  fwrite(&mix_flag, sizeof(int), 1, fp); */
  ;
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::read_restart_settings(FILE *fp)
{
  /* if (comm->me == 0) {
    utils::sfread(FLERR, &cut_global, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &offset_flag, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &mix_flag, sizeof(int), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_global, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&offset_flag, 1, MPI_INT, 0, world);
  MPI_Bcast(&mix_flag, 1, MPI_INT, 0, world); */
  ;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::write_data(FILE *fp)
{
  /* for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp, "%d %g %g %g\n", i, d0[i][i], alpha[i][i], r0[i][i]); */
  ;
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMorseAdResSCGIPCOrg::write_data_all(FILE *fp)
{
  /* for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp, "%d %d %g %g %g %g\n", i, j, d0[i][j], alpha[i][j], r0[i][j], cut[i][j]); */
  ;
}

/* ---------------------------------------------------------------------- */

double PairMorseAdResSCGIPCOrg::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                         double /*factor_coul*/, double factor_lj, double &fforce)
{
  double r, dr, dexp, phi;

  r = sqrt(rsq);
  dr = r - r0[itype][jtype];
  dexp = exp(-alpha[itype][jtype] * dr);
  fforce = factor_lj * morse1[itype][jtype] * (dexp * dexp - dexp) / r;

  phi = d0[itype][jtype] * (dexp * dexp - 2.0 * dexp) - offset[itype][jtype];
  return factor_lj * phi;
}

/* ---------------------------------------------------------------------- */

void *PairMorseAdResSCGIPCOrg::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "d0") == 0) return (void *) d0;
  if (strcmp(str, "r0") == 0) return (void *) r0;
  if (strcmp(str, "alpha") == 0) return (void *) alpha;
  return nullptr;
}

/* ----------------------------------------------------------------------
  find molecules
------------------------------------------------------------------------- */

int PairMorseAdResSCGIPCOrg::find_mols(tagint &idlo, tagint &idhi) {

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
  
  memory->create(molmap_tmp,nlen,"morse/adress/cg:molmap");
  for (i = 0; i < nlen; i++) molmap_tmp[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i])
      molmap_tmp[molecule[i]-idlo] = 1;

  int *molmapall;
  memory->create(molmapall,nlen,"morse/adress/cg:molmapall");
  MPI_Allreduce(molmap_tmp,molmapall,nlen,MPI_INT,MPI_MAX,world);

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap_tmp[i] = nmolecules++;
    else molmap_tmp[i] = -1;
    
  memory->destroy(molmapall);  
  memory->destroy(molmap_tmp);  
  
  return nmolecules;

}