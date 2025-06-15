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
FixStyle(tdforce,FixTDForce);
// clang-format on
#else

#ifndef LMP_FIX_TDFORCE_H
#define LMP_FIX_TDFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTDForce : public Fix {
    public:
        FixTDForce(class LAMMPS *, int, char **);
        ~FixTDForce() override;

        int setmask() override;
        void init() override;
        void post_force(int) override;
        void setup(int) override;
        // void min_setup(int) override;
        // void pre_exchange();  
    
    protected:
        int tabstyle, tablength;
        struct Table {
            int ninput;
            double fplo, fphi;
            double *rfile, *ffile; // distances forces
            double *f2file; // build spline representation and set values in f2file
            double delta, invdelta, deltasq6;
            double *r, *f, *df, *f2;
        };
        int ntables;
        Table *tables;

        double x_lo, x_hi, y_lo, y_hi, z_lo, z_hi; // last four are not used ... only one dir
        double center_box[3];
        double len_AT, len_HY, len_CG;
        double rlo, rhi;

        void null_table(Table *);
        void param_extract(Table *, char *);
        void read_table(Table *, char *, char *);
        void bcast_table(Table *);
        void spline_table(Table *);
        void spline(double *, double *, int, double, double, double *);
        void compute_table(Table *);
        double splint(double *, double *, double *, int, double);
        void free_table(Table *);   
        void force_lookup(double, double &);

        tagint idlo, idhi;
        int nmolecules;
        int find_mols(tagint &, tagint &); // PP using maxmol_all
        tagint *molmap_tmp;

        tagint maxmol_all; 
        // void find_maxid(); 
};
}

#endif
#endif