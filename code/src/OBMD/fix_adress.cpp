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

#include "fix_adress.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "respa.h"
#include "update.h"
#include "domain.h" 
#include "memory.h" 
#include "compute.h" 
#include "compute_chunk_atom.h" 
#include "modify.h"
#include "math_const.h" 

#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
#define BIG MAXTAGINT
#define EPS 1.0e-12

/* ---------------------------------------------------------------------- */

// Adapted from H-AdResS code

/* ---------------------------------------------------------------------- */

FixAdResS::FixAdResS(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // std::cout<<"fix_adress"<<"\n";
  len_AT = utils::numeric(FLERR,arg[3],false,lmp);
  len_HY = utils::numeric(FLERR,arg[4],false,lmp);
  
  x0lo = domain->boxlo[0];
  x0hi = domain->boxhi[0];
  x1lo = domain->boxlo[1];
  x1hi = domain->boxhi[1];
  x2lo = domain->boxlo[2];
  x2hi = domain->boxhi[2];

  center_box[0] = (x0hi + x0lo)/2.0;
  center_box[1] = (x1hi + x1lo)/2.0;
  center_box[2] = (x2hi + x2lo)/2.0;
    
  len_CG = 0.5*(x0hi-x0lo-len_AT-2.0*len_HY);
  
  options(narg-5,&arg[5]);
  
  dynamic_group_allow = 1;
  time_integrate = 1;
}

/* ---------------------------------------------------------------------- */

FixAdResS::~FixAdResS()
{
  
  ; // empty

}

/* ---------------------------------------------------------------------- */

int FixAdResS::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdResS::init()
{
  post_integrate();
}

/* ----------------------------------------------------------------------
   allow for both per-type and per-atom mass
------------------------------------------------------------------------- */

void FixAdResS::setup(int vflag)
{
	post_integrate();
}

/* ---------------------------------------------------------------------- */

void FixAdResS::post_integrate()
{ 
  // std::cout<<"fix_adress [post_integrate()]"<<"\n";
  double **x = atom->x;
  double **v = atom->v;
  double *mass = atom->mass;
  double *rmass=atom->rmass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  tagint *molecule = atom->molecule;
  imageint *image = atom->image;
  tagint imol;
  double massone;
  double unwrap[3];  
  double *lambdaF = atom->lambdaF;
  double **cms_mol = atom->cms_mol;
  double **vcms_mol = atom->vcms_mol; 
  int i, cnt;

  // PP 
  nmolecules = find_mols(idlo,idhi); 
  // find_maxid();
  maxmol_all = idhi;
  // std::cout<<"nmolecules: "<<nmolecules<<" maxmol_all: "<<maxmol_all<<"\n";
  
  // mass
  double *massproc_tmp, *masstotal_tmp;
  memory->create(massproc_tmp,maxmol_all,"adress:massproc_tmp");
  memory->create(masstotal_tmp,maxmol_all,"adress:masstotal_tmp");

  // center-of-mass and weight factor
  double **com_tmp, **comall_tmp;
  double **vcom_tmp, **vcomall_tmp;
  double *lambda_tmp;
  memory->create(com_tmp,maxmol_all,3,"adress:com_tmp");
  memory->create(comall_tmp,maxmol_all,3,"adress:comall_tmp");
  memory->create(vcom_tmp,maxmol_all,3,"adress:com_tmp");
  memory->create(vcomall_tmp,maxmol_all,3,"adress:vcomall_tmp");
  memory->create(lambda_tmp,maxmol_all,"adress:lambda_tmp");

  // init
  for (cnt = 0; cnt < maxmol_all; cnt++) {
    massproc_tmp[cnt] = 0.0;

    com_tmp[cnt][0] = 0.0;
    com_tmp[cnt][1] = 0.0;
    com_tmp[cnt][2] = 0.0;

    vcom_tmp[cnt][0] = 0.0;
    vcom_tmp[cnt][1] = 0.0;
    vcom_tmp[cnt][2] = 0.0;

    lambda_tmp[cnt] = 0.0;
  }

  // compute masses
  for (i = 0; i < nlocal; i++) {
    if (mask[i]) {
      if (rmass) massone = rmass[i];
      else {
        massone = mass[type[i]];
      }

      /* tagint *tag = atom->tag;
      int G = tag[i];
      double **f = atom->f;
      if (tag[i] == 10923) {
        std::cout<<"fix adress 1"<<"\n";
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
        std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
        std::cout<<"x[i][0]: "<<atom->x[i][0]<<" x[i][1]: "<<atom->x[i][1]<<" x[i][2]: "<<atom->x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
      } */

      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;
      massproc_tmp[imol] += massone;
    }
  }

  MPI_Allreduce(massproc_tmp, masstotal_tmp, maxmol_all, MPI_DOUBLE, MPI_SUM, world);

  // std::cout<<"masstotal_tmp[152]: "<<masstotal_tmp[152]<<"\n";
  /* std::cout<<"masstotal_tmp[1642]: "<<masstotal_tmp[1642]<<"\n"; */

  double **f = atom->f;
  // compute centre-of-mass
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      massone = mass[type[i]];
      
      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;

      /* if (x[i][0] < x0lo || x[i][0] > x0hi) {
        std::cout<<"x0lo: "<<x0lo<<" x0hi: "<<x0hi<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<"\n";
      } */

      domain->unmap(x[i],image[i],unwrap); // returns unwrap vector
      massone /= masstotal_tmp[imol]; // fraction
      com_tmp[imol][0] += unwrap[0] * massone;
      com_tmp[imol][1] += unwrap[1] * massone;
      com_tmp[imol][2] += unwrap[2] * massone;

      vcom_tmp[imol][0] += v[i][0] * massone;
      vcom_tmp[imol][1] += v[i][1] * massone;
      vcom_tmp[imol][2] += v[i][2] * massone;

      /* if (imol == 3515) {
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"i: "<<i<<"\n";
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"lambdaF[imol]: "<<lambdaF[imol]<<"\n";
        std::cout<<"massone: "<<massone<<" unwrap[0]: "<<unwrap[0]<<" unwrap[1]: "<<unwrap[1]<<" unwrap[2]: "<<unwrap[2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
        std::cout<<"com_tmp[imol][0]: "<<com_tmp[imol][0]<<" com_tmp[imol][1]: "<<com_tmp[imol][2]<<" com_tmp[imol][2]: "<<com_tmp[imol][2]<<"\n";
      } 

      if (imol == 3516) {
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"i: "<<i<<"\n";
        std::cout<<"tag[i]: "<<tag[i]<<"\n";
        std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
        std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
        std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
        std::cout<<"lambdaF[imol]: "<<lambdaF[imol]<<"\n";
        std::cout<<"massone: "<<massone<<" unwrap[0]: "<<unwrap[0]<<" unwrap[1]: "<<unwrap[1]<<" unwrap[2]: "<<unwrap[2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"f[i][0]: "<<f[i][0]<<" f[i][1]: "<<f[i][1]<<" f[i][2]: "<<f[i][2]<<"\n";
        std::cout<<"com_tmp[imol][0]: "<<com_tmp[imol][0]<<" com_tmp[imol][1]: "<<com_tmp[imol][2]<<" com_tmp[imol][2]: "<<com_tmp[imol][2]<<"\n";
      } */


      /* if (imol == 3640) {
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"lambdaF[imol]: "<<lambdaF[imol]<<"\n";
        std::cout<<"massone: "<<massone<<" unwrap[0]: "<<unwrap[0]<<" unwrap[1]: "<<unwrap[1]<<" unwrap[2]: "<<unwrap[2]<<"\n";
        std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";
        std::cout<<"com_tmp[imol][0]: "<<com_tmp[imol][0]<<" com_tmp[imol][1]: "<<com_tmp[imol][2]<<" com_tmp[imol][2]: "<<com_tmp[imol][2]<<"\n";
      } */
    }
  }

  MPI_Allreduce(&com_tmp[0][0],&comall_tmp[0][0],3*maxmol_all,MPI_DOUBLE,MPI_SUM,world);
  MPI_Allreduce(&vcom_tmp[0][0],&vcomall_tmp[0][0],3*maxmol_all,MPI_DOUBLE,MPI_SUM,world);
  
  // compute weights ... lambda
  double sinx, xtmp, ytmp, ztmp, rtmp, phase;

  for (i = 0; i < maxmol_all; i++) { // using molecules 
        
    // std::cout<<"i: "<<i<<" comall_tmp[i][0]: "<<comall_tmp[i][0]<<" comall_tmp[i][1]: "<<comall_tmp[i][1]<<" comall_tmp[i][2]: "<<comall_tmp[i][2]<<"\n";
    domain->remap(comall_tmp[i]); // fold coordintaes back in box ! ! !

    if(topo==1) { // using & testing slab geometry
        xtmp = comall_tmp[i][0];
        if(xtmp < x0lo + len_CG) {
        lambda_tmp[i] = 0.0;	
      } else if(xtmp < x0lo + len_CG + len_HY && xtmp > x0lo + len_CG) {
        phase = MY_PI*(xtmp - len_CG - x0lo)/2/len_HY; 
        sinx = sin(phase);
        lambda_tmp[i] = sinx*sinx;	
      } else if(xtmp < x0hi - len_CG - len_HY && xtmp > x0lo + len_CG + len_HY) {
        lambda_tmp[i] = 1.0;		  
      } else if(xtmp < x0hi - len_CG && xtmp > x0hi - len_CG - len_HY) {
        phase = MY_PI*(xtmp - len_CG - len_AT - x0lo)/2/len_HY; 
        sinx = sin(phase);
        lambda_tmp[i] = sinx*sinx;		  
      } else if(xtmp > x0hi - len_CG) {
        lambda_tmp[i] = 0.0;
      }
    }
    if(topo==2) {
      xtmp = comall_tmp[i][0]-center_box[0];
      ytmp = comall_tmp[i][1]-center_box[1];
      ztmp = comall_tmp[i][2]-center_box[2];
      rtmp = sqrt(xtmp*xtmp + ytmp*ytmp + ztmp*ztmp);
      if(rtmp < len_AT) {
        lambda_tmp[i] = 1.0;	
      } else if(rtmp < len_AT + len_HY && rtmp > len_AT) {
        phase = MY_PI*(len_HY+len_AT-rtmp)/2/len_HY;
        sinx = sin(phase);
        lambda_tmp[i] = sinx*sinx;	
      } else if(xtmp > len_CG) {
        lambda_tmp[i] = 0.0;
      }
    }
    if(topo==3) {
      xtmp = comall_tmp[i][0]-center_box[0];
      ztmp = comall_tmp[i][2]-center_box[2];
      rtmp = sqrt(xtmp*xtmp + ztmp*ztmp);
      if(rtmp < len_AT) {
        lambda_tmp[i] = 1.0;	
      } else if(rtmp < len_AT + len_HY && rtmp > len_AT) {
        phase = MY_PI*(len_HY+len_AT-rtmp)/2/len_HY;
        sinx = sin(phase);
        lambda_tmp[i] = sinx*sinx;	
      } else if(xtmp > len_CG) {
        lambda_tmp[i] = 0.0;
      }
    }
  }  
  
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      imol = molecule[i];
      if (molmap_tmp) imol = molmap_tmp[imol-idlo];
      else imol--;

      lambdaF[i] = lambda_tmp[imol]; 
      cms_mol[i][0] = comall_tmp[imol][0];
      cms_mol[i][1] = comall_tmp[imol][1];
      cms_mol[i][2] = comall_tmp[imol][2];

      vcms_mol[i][0] = vcomall_tmp[imol][0];
      vcms_mol[i][1] = vcomall_tmp[imol][1];
      vcms_mol[i][2] = vcomall_tmp[imol][2];

      /* if (imol == 3488) {
        std::cout<<"imol: "<<imol<<"\n";
        std::cout<<"lambdaF[imol]: "<<lambdaF[imol]<<"\n";
        std::cout<<"comall_tmp[imol][0]: "<<comall_tmp[imol][0]<<"\n";
      } */
    }
  }

  // PP ... just testing
  /* tagint *tag = atom->tag;
  int G = tag[i];
  if (tag[i] == 7386) {
    std::cout<<"tag[i]: "<<tag[i]<<"\n";
    std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
    std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
    std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
    std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
    std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
    std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
  } 

  if (tag[i] == 10546) {
    std::cout<<"tag[i]: "<<tag[i]<<"\n";
    std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
    std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
    std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
    std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
    std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
    std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
  } 

  if (tag[i] == 10548) {
    std::cout<<"tag[i]: "<<tag[i]<<"\n";
    std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
    std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
    std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
    std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
    std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
    std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
  } 

  if (tag[i] == 10549) {
    std::cout<<"tag[i]: "<<tag[i]<<"\n";
    std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
    std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
    std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
    std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
    std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
    std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
  } 

  if (tag[i] == 10551) {
    std::cout<<"tag[i]: "<<tag[i]<<"\n";
    std::cout<<"atom->molecule[i]: "<<atom->molecule[i]<<"\n";
    std::cout<<"atom->type[i]: "<<atom->type[i]<<"\n";
    std::cout<<"atom->rep_atom[i]: "<<atom->rep_atom[i]<<"\n";
    std::cout<<"atom->lambdaF[i]: "<<atom->lambdaF[i]<<"\n";
    std::cout<<"atom->cms_mol[i][0]: "<<atom->cms_mol[i][0]<<" atom->cms_mol[i][1]: "<<atom->cms_mol[i][1]<<" atom->cms_mol[i][2]: "<<atom->cms_mol[i][2]<<"\n";
    std::cout<<"x[i][0]: "<<x[i][0]<<" x[i][1]: "<<x[i][1]<<" x[i][2]: "<<x[i][2]<<"\n";;
  } */

  /* if (tag[i] == 7471) {
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

  memory->destroy(massproc_tmp);
  memory->destroy(masstotal_tmp);
  memory->destroy(com_tmp);
  memory->destroy(comall_tmp);
  memory->destroy(vcom_tmp);
  memory->destroy(vcomall_tmp);
  memory->destroy(lambda_tmp);

  // std::cout<<"fix_adress [post_integrate() ::: END]"<<"\n";
}

/* ---------------------------------------------------------------------- */
// find number of molecules in a system
int FixAdResS::find_mols(tagint &idlo, tagint &idhi)
{
  int i;
  tagint *molecule = atom->molecule;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
   
  tagint lo = BIG;
  tagint hi = -BIG;
  
  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      lo = MIN(lo,molecule[i]);
      hi = MAX(hi,molecule[i]); 
    }
  }

  MPI_Allreduce(&lo,&idlo,1,MPI_LMP_TAGINT,MPI_MIN,world);
  MPI_Allreduce(&hi,&idhi,1,MPI_LMP_TAGINT,MPI_MAX,world);
  
  if (idlo == BIG) return 0;
  
  int nlen = (int) idhi-idlo+1;
  
  memory->create(molmap_tmp,nlen,"adress:molmap");
  for (i = 0; i < nlen; i++) molmap_tmp[i] = 0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      molmap_tmp[molecule[i]-idlo] = 1;

  int *molmapall;
  memory->create(molmapall,nlen,"adress:molmapall");
  MPI_Allreduce(molmap_tmp,molmapall,nlen,MPI_INT,MPI_MAX,world);

  // nmolecules = # of non-zero IDs in molmap
  // molmap[i] = index of molecule, skipping molecules not in group with -1

  int nmolecules = 0;
  for (i = 0; i < nlen; i++)
    if (molmapall[i]) molmap_tmp[i] = nmolecules++;
    else molmap_tmp[i] = -1;
    
  memory->destroy(molmapall);  
  memory->destroy(molmap_tmp);  
  
  return nmolecules;
}

/* ----------------------------------------------------------------------
   maxtag_all = current max atom ID for all atoms
   maxmol_all = current max molecule ID for all atoms
------------------------------------------------------------------------- */

/* void FixAdResS::find_maxid()
{ 
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;

  int max = 0;
  for (int i = 0; i < nlocal; i++) max = MAX(max,molecule[i]);
  MPI_Allreduce(&max,&maxmol_all,1,MPI_LMP_TAGINT,MPI_MAX,world);
} */

/* ---------------------------------------------------------------------- */

void FixAdResS::options(int narg, char **arg)
{
  
  topo = 1;
  flag_inverse = 0;
  
  int iarg = 0;
  while (iarg < narg) {
    if(strcmp(arg[iarg],"slab") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix adress command");
      topo = 1;
      iarg += 1;
    } else if (strcmp(arg[iarg],"sphere") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix adress command");		
      topo = 2;
      iarg += 1;  
    } else if (strcmp(arg[iarg],"cylinder") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix adress command");
      topo = 3;
      iarg += 1;  
    } else if (strcmp(arg[iarg],"inverse") == 0) {
      if (iarg+1 > narg) error->all(FLERR,"Illegal fix adress command");
      flag_inverse = 1;
      iarg += 1;  
    }
  }
}
