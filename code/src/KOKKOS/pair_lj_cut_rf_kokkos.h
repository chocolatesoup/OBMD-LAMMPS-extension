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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/rf/kk,PairLJCutRFKokkos<LMPDeviceType>);
PairStyle(lj/cut/rf/kk/device,PairLJCutRFKokkos<LMPDeviceType>);
PairStyle(lj/cut/rf/kk/host,PairLJCutRFKokkos<LMPHostType>);
// clang-format on
#else

// clang-format off
#ifndef LMP_PAIR_LJ_CUT_RF_KOKKOS_H
#define LMP_PAIR_LJ_CUT_RF_KOKKOS_H

#include "pair_kokkos.h"
#include "pair_lj_cut_rf.h"
#include "neigh_list_kokkos.h"

namespace LAMMPS_NS {

template<class DeviceType>
class PairLJCutRFKokkos : public PairLJCutRF { // inherits from PairLJCutRF (changing only compute())
 public:
  enum {EnabledNeighFlags=FULL|HALFTHREAD|HALF};
  enum {COUL_FLAG=1};
  typedef DeviceType device_type;
  typedef ArrayTypes<DeviceType> AT;
  PairLJCutRFKokkos(class LAMMPS *);
  ~PairLJCutRFKokkos() override;

  void compute(int, int) override;

  void settings(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;

  /* create new struct */
  // struct params_lj_rf{
  //  KOKKOS_INLINE_FUNCTION
  //  params_lj_rf() {cut_ljsq=0;cut_coulsq=0;epsilon_rf=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  //  KOKKOS_INLINE_FUNCTION
  //  params_lj_rf(int /*i*/) {cut_ljsq=0;cut_coulsq=0;epsilon_rf=0;lj1=0;lj2=0;lj3=0;lj4=0;offset=0;};
  //  F_FLOAT cut_ljsq,cut_coulsq,epsilon_rf,lj1,lj2,lj3,lj4,offset;
  // };

 protected:
  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fpair(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_fcoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_evdwl(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype) const;

  template<bool STACKPARAMS, class Specialisation>
  KOKKOS_INLINE_FUNCTION
  F_FLOAT compute_ecoul(const F_FLOAT& rsq, const int& i, const int&j,
                        const int& itype, const int& jtype, const F_FLOAT& factor_coul, const F_FLOAT& qtmp) const;

  Kokkos::DualView<params_lj_rf**,Kokkos::LayoutRight,DeviceType> k_params;
  typename Kokkos::DualView<params_lj_rf**,
    Kokkos::LayoutRight,DeviceType>::t_dev_const_um params;
  // hardwired to space for 12 atom types
  params_lj_rf m_params[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];

  F_FLOAT m_cutsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_ljsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_cut_coulsq[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  F_FLOAT m_epsilon_rf[MAX_TYPES_STACKPARAMS+1][MAX_TYPES_STACKPARAMS+1];
  typename AT::t_x_array_randomread x;
  typename AT::t_x_array c_x;
  typename AT::t_f_array f;
  typename AT::t_int_1d_randomread type;
  typename AT::t_float_1d_randomread q;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  typename AT::t_efloat_1d d_eatom;
  typename AT::t_virial_array d_vatom;

  int newton_pair;

  typename AT::tdual_ffloat_2d k_cutsq;
  typename AT::t_ffloat_2d d_cutsq;
  typename AT::tdual_ffloat_2d k_cut_ljsq;
  typename AT::t_ffloat_2d d_cut_ljsq;
  typename AT::tdual_ffloat_2d k_cut_coulsq;
  typename AT::t_ffloat_2d d_cut_coulsq;
  typename AT::tdual_ffloat_2d k_epsilon_rf;
  typename AT::t_ffloat_2d d_epsilon_rf;


  int neighflag;
  int nlocal,nall,eflag,vflag;

  double special_coul[4];
  double special_lj[4];
  double qqrd2e;

  void allocate() override;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,FULL,true,0>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,FULL,true,1>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,HALF,true>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,HALFTHREAD,true>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,FULL,false,0>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,FULL,false,1>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,HALF,false>;
  friend struct PairComputeFunctor<PairLJCutRFKokkos,HALFTHREAD,false>;
  friend EV_FLOAT pair_compute_neighlist<PairLJCutRFKokkos,FULL,0>(PairLJCutRFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutRFKokkos,FULL,1>(PairLJCutRFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutRFKokkos,HALF>(PairLJCutRFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute_neighlist<PairLJCutRFKokkos,HALFTHREAD>(PairLJCutRFKokkos*,NeighListKokkos<DeviceType>*);
  friend EV_FLOAT pair_compute<PairLJCutRFKokkos,void>(PairLJCutRFKokkos*,
                                                            NeighListKokkos<DeviceType>*);
  friend void pair_virial_fdotr_compute<PairLJCutRFKokkos>(PairLJCutRFKokkos*);

};

}

#endif
#endif

