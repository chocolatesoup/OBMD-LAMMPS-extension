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

#include "pair_lj_cut_rf_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "error.h"
#include "force.h"
#include "kokkos.h"
#include "memory_kokkos.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cmath>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutRFKokkos<DeviceType>::PairLJCutRFKokkos(LAMMPS *lmp):PairLJCutRF(lmp)
{
  respa_enable = 0;

  kokkosable = 1;
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = X_MASK | F_MASK | TYPE_MASK | Q_MASK | ENERGY_MASK | VIRIAL_MASK;
  datamask_modify = F_MASK | ENERGY_MASK | VIRIAL_MASK;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
PairLJCutRFKokkos<DeviceType>::~PairLJCutRFKokkos()
{
  if (copymode) return;

  if (allocated) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->destroy_kokkos(k_cutsq,cutsq);
    memoryKK->destroy_kokkos(k_cut_ljsq,cut_ljsq);
    memoryKK->destroy_kokkos(k_cut_coulsq,cut_coulsq);
    memoryKK->destroy_kokkos(k_epsilon_rf,epsilon_rf);
  }
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutRFKokkos<DeviceType>::compute(int eflag_in, int vflag_in)
{
  eflag = eflag_in;
  vflag = vflag_in;

  if (neighflag == FULL) no_virial_fdotr_compute = 1;

  ev_init(eflag,vflag,0);

  // reallocate per-atom arrays if necessary

  if (eflag_atom) {
    memoryKK->destroy_kokkos(k_eatom,eatom);
    memoryKK->create_kokkos(k_eatom,eatom,maxeatom,"pair:eatom");
    d_eatom = k_eatom.view<DeviceType>();
  }
  if (vflag_atom) {
    memoryKK->destroy_kokkos(k_vatom,vatom);
    memoryKK->create_kokkos(k_vatom,vatom,maxvatom,"pair:vatom");
    d_vatom = k_vatom.view<DeviceType>();
  }

  atomKK->sync(execution_space,datamask_read);
  k_cutsq.template sync<DeviceType>();
  k_cut_ljsq.template sync<DeviceType>();
  k_cut_coulsq.template sync<DeviceType>();
  k_epsilon_rf.template sync<DeviceType>();
  k_params.template sync<DeviceType>();
  if (eflag || vflag) atomKK->modified(execution_space,datamask_modify);
  else atomKK->modified(execution_space,F_MASK);

  x = atomKK->k_x.view<DeviceType>();
  c_x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  q = atomKK->k_q.view<DeviceType>();
  type = atomKK->k_type.view<DeviceType>();
  nlocal = atom->nlocal;
  nall = atom->nlocal + atom->nghost;
  special_lj[0] = force->special_lj[0];
  special_lj[1] = force->special_lj[1];
  special_lj[2] = force->special_lj[2];
  special_lj[3] = force->special_lj[3];
  special_coul[0] = force->special_coul[0];
  special_coul[1] = force->special_coul[1];
  special_coul[2] = force->special_coul[2];
  special_coul[3] = force->special_coul[3];
  qqrd2e = force->qqrd2e;
  newton_pair = force->newton_pair;

  // loop over neighbors of my atoms

  EV_FLOAT ev = pair_compute<PairLJCutRFKokkos<DeviceType>,void > 
    (this,(NeighListKokkos<DeviceType>*)list);

  if (eflag) {
    eng_vdwl += ev.evdwl;
    eng_coul += ev.ecoul;
  }
  if (vflag_global) {
    virial[0] += ev.v[0];
    virial[1] += ev.v[1];
    virial[2] += ev.v[2];
    virial[3] += ev.v[3];
    virial[4] += ev.v[4];
    virial[5] += ev.v[5];
  }

  if (eflag_atom) {
    k_eatom.template modify<DeviceType>();
    k_eatom.template sync<LMPHostType>();
  }

  if (vflag_atom) {
    k_vatom.template modify<DeviceType>();
    k_vatom.template sync<LMPHostType>();
  }

  if (vflag_fdotr) pair_virial_fdotr_compute(this);

}

/* ----------------------------------------------------------------------
   compute LJ 12-6 pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutRFKokkos<DeviceType>::
compute_fpair(const F_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;
  F_FLOAT forcelj;

  forcelj = r6inv *
    ((STACKPARAMS?m_params[itype][jtype].lj1:params(itype,jtype).lj1)*r6inv -
     (STACKPARAMS?m_params[itype][jtype].lj2:params(itype,jtype).lj2));

  return forcelj*r2inv;
}

/* ----------------------------------------------------------------------
   compute coulomb pair force between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutRFKokkos<DeviceType>::
compute_fcoul(const F_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& itype, const int& jtype,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT rf_fctr_0 = 1.0/pow(sqrt(rsq),3.0);
  const F_FLOAT epsilon_rf_ij = STACKPARAMS?m_params[itype][jtype].epsilon_rf:params(itype,jtype).epsilon_rf;
  const F_FLOAT rf_fctr_1 = (epsilon_rf_ij-1.0);
  const F_FLOAT rf_fctr_2 = 1.0+2.0*epsilon_rf_ij;
  const F_FLOAT cut_coulsq_ij = STACKPARAMS?m_params[itype][jtype].cut_coulsq:params(itype,jtype).cut_coulsq;
  const F_FLOAT cut_coul_ij = sqrt(cut_coulsq_ij);

  F_FLOAT forcecoul;

  forcecoul = (qqrd2e*qtmp*q[j])*((r2inv*sqrt(r2inv))-(1.0/pow(cut_coul_ij,3.0)*(2.0*rf_fctr_1/rf_fctr_2)));

  return factor_coul*forcecoul;
}

/* ----------------------------------------------------------------------
   compute LJ 12-6 pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutRFKokkos<DeviceType>::
compute_evdwl(const F_FLOAT& rsq, const int& /*i*/, const int& /*j*/,
              const int& itype, const int& jtype) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT r6inv = r2inv*r2inv*r2inv;

  return r6inv*
    ((STACKPARAMS?m_params[itype][jtype].lj3:params(itype,jtype).lj3)*r6inv
     - (STACKPARAMS?m_params[itype][jtype].lj4:params(itype,jtype).lj4))
    -  (STACKPARAMS?m_params[itype][jtype].offset:params(itype,jtype).offset);

}

/* ----------------------------------------------------------------------
   compute coulomb pair potential energy between atoms i and j
   ---------------------------------------------------------------------- */
template<class DeviceType>
template<bool STACKPARAMS, class Specialisation>
KOKKOS_INLINE_FUNCTION
F_FLOAT PairLJCutRFKokkos<DeviceType>::
compute_ecoul(const F_FLOAT& rsq, const int& /*i*/, const int& j,
              const int& itype, const int& jtype,
              const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const {
  const F_FLOAT r2inv = 1.0/rsq;
  const F_FLOAT rf_fctr_0 = 1.0/pow(sqrt(rsq),3.0);
  const F_FLOAT epsilon_rf_ij = STACKPARAMS?m_params[itype][jtype].epsilon_rf:params(itype,jtype).epsilon_rf;
  const F_FLOAT rf_fctr_1 = (epsilon_rf_ij-1.0);
  const F_FLOAT rf_fctr_2 = 1.0+2.0*epsilon_rf_ij;
  const F_FLOAT cut_coulsq_ij = STACKPARAMS?m_params[itype][jtype].cut_coulsq:params(itype,jtype).cut_coulsq;
  const F_FLOAT cut_coul_ij = sqrt(cut_coulsq_ij);

  F_FLOAT ecoul;

  ecoul = qqrd2e*qtmp*q[j]*sqrt(r2inv)*(1.0+(rf_fctr_1/rf_fctr_2)*(pow(sqrt(rsq)/cut_coul_ij,3.0)))-(qqrd2e*qtmp*q[j])*(1.0/cut_coul_ij)*(3.0*epsilon_rf_ij/rf_fctr_2);

  return factor_coul * ecoul;

}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutRFKokkos<DeviceType>::allocate()
{
  PairLJCutRF::allocate();

  int n = atom->ntypes;
  memory->destroy(cutsq);
  memoryKK->create_kokkos(k_cutsq,cutsq,n+1,n+1,"pair:cutsq");
  d_cutsq = k_cutsq.template view<DeviceType>();
  memory->destroy(cut_ljsq);
  memoryKK->create_kokkos(k_cut_ljsq,cut_ljsq,n+1,n+1,"pair:cut_ljsq");
  d_cut_ljsq = k_cut_ljsq.template view<DeviceType>();
  memory->destroy(cut_coulsq);
  memoryKK->create_kokkos(k_cut_coulsq,cut_coulsq,n+1,n+1,"pair:cut_coulsq");
  d_cut_coulsq = k_cut_coulsq.template view<DeviceType>();
  memory->destroy(epsilon_rf);
  memoryKK->create_kokkos(k_epsilon_rf,epsilon_rf,n+1,n+1,"pair:epsilon_rf");
  d_epsilon_rf = k_epsilon_rf.template view<DeviceType>();
  k_params = Kokkos::DualView<params_lj_rf**,Kokkos::LayoutRight,DeviceType>("PairLJCutRF::params",n+1,n+1);
  params = k_params.template view<DeviceType>();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutRFKokkos<DeviceType>::settings(int narg, char **arg)
{
  if (narg > 2) error->all(FLERR,"Illegal pair_style command");

  PairLJCutRF::settings(1,arg);
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

template<class DeviceType>
void PairLJCutRFKokkos<DeviceType>::init_style()
{
  PairLJCutRF::init_style();

  // error if rRESPA with inner levels

  /* if (update->whichflag == 1 && utils::strmatch(update->integrate_style,"^respa")) {
    int respa = 0;
    if (((Respa *) update->integrate)->level_inner >= 0) respa = 1;
    if (((Respa *) update->integrate)->level_middle >= 0) respa = 2;
    if (respa)
      error->all(FLERR,"Cannot use Kokkos pair style with rRESPA inner/middle");
  } */

  // adjust neighbor list request for KOKKOS

  neighflag = lmp->kokkos->neighflag;
  auto request = neighbor->find_request(this);
  request->set_kokkos_host(std::is_same_v<DeviceType,LMPHostType> &&
                           !std::is_same_v<DeviceType,LMPDeviceType>);
  request->set_kokkos_device(std::is_same_v<DeviceType,LMPDeviceType>);
  if (neighflag == FULL) request->enable_full();
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

template<class DeviceType>
double PairLJCutRFKokkos<DeviceType>::init_one(int i, int j)
{
  double cutone = PairLJCutRF::init_one(i,j);
  double cut_ljsqm = cut_ljsq[i][j];
  double cut_coulsqm = cut_coulsq[i][j];
  double epsilon_rfm = epsilon_rf[i][j];

  k_params.h_view(i,j).lj1 = lj1[i][j];
  k_params.h_view(i,j).lj2 = lj2[i][j];
  k_params.h_view(i,j).lj3 = lj3[i][j];
  k_params.h_view(i,j).lj4 = lj4[i][j];
  k_params.h_view(i,j).offset = offset[i][j];
  k_params.h_view(i,j).cut_ljsq = cut_ljsqm;
  k_params.h_view(i,j).cut_coulsq = cut_coulsqm;
  k_params.h_view(i,j).epsilon_rf = epsilon_rfm;

  k_params.h_view(j,i) = k_params.h_view(i,j);
  if (i<MAX_TYPES_STACKPARAMS+1 && j<MAX_TYPES_STACKPARAMS+1) {
    m_params[i][j] = m_params[j][i] = k_params.h_view(i,j);
    m_cutsq[j][i] = m_cutsq[i][j] = cutone*cutone;
    m_cut_ljsq[j][i] = m_cut_ljsq[i][j] = cut_ljsqm;
    m_cut_coulsq[j][i] = m_cut_coulsq[i][j] = cut_coulsqm;
    m_epsilon_rf[j][i] = m_epsilon_rf[i][j] = epsilon_rfm;
  }

  k_cutsq.h_view(i,j) = k_cutsq.h_view(j,i) = cutone*cutone;
  k_cutsq.template modify<LMPHostType>();
  k_cut_ljsq.h_view(i,j) = k_cut_ljsq.h_view(j,i) = cut_ljsqm;
  k_cut_ljsq.template modify<LMPHostType>();
  k_cut_coulsq.h_view(i,j) = k_cut_coulsq.h_view(j,i) = cut_coulsqm;
  k_cut_coulsq.template modify<LMPHostType>();
  k_epsilon_rf.h_view(i,j) = k_epsilon_rf.h_view(j,i) = epsilon_rfm;
  k_epsilon_rf.template modify<LMPHostType>();
  k_params.template modify<LMPHostType>();

  return cutone;
}

namespace LAMMPS_NS {
template class PairLJCutRFKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class PairLJCutRFKokkos<LMPHostType>;
#endif
}

