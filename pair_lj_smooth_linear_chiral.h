/* ----------------------------------------------------------------------
    Chiral Tetramer Molecular Model (CTMM), Copyright (2021) Nikolai D. Petsev. 
    This software is distributed under the GNU General Public License, 
    and is derived from LAMMPS - Large-Scale Atomic/Molecular Massively Parallel Simulator
    by Steve Plimpton, Copyright (2003) Sandia Corporation, https://www.lammps.org.

    See the README file in the top-level directory.
   
    This file is part of CTMM.

    CTMM is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    CTMM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with CTMM.  If not, see <https://www.gnu.org/licenses/>.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/smooth/linear/chiral,PairLJSmoothLinearChiral)

#else

#ifndef PAIR_LJ_SMOOTH_LINEAR_CHIRAL_H
#define PAIR_LJ_SMOOTH_LINEAR_CHIRAL_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSmoothLinearChiral : public Pair {
 public:
  PairLJSmoothLinearChiral(class LAMMPS *);
  virtual ~PairLJSmoothLinearChiral();
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void  computezeta();				
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 protected:
  int nmax;
  double *zeta, *dzetax, *dzetay, *dzetaz;
  double *i1_vec, *i2_vec, *i3_vec, *i4_vec;
  double cut_global;
  double **cut;
  int nstep, first;
  double **epsilon,**sigma;
  double **ljcut,**dljcut;
  double **lj1,**lj2,**lj3,**lj4;
  double bias_global;
  double **bias;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

*/
