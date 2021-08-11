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

#ifdef COMPUTE_CLASS

ComputeStyle(zeta/ave,ComputeZetaAve)

#else

#ifndef LMP_COMPUTE_ZETA_AVE_H
#define LMP_COMPUTE_ZETA_AVE_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeZetaAve : public Compute {
 public:
  ComputeZetaAve(class LAMMPS *, int, char **);
  ~ComputeZetaAve();
  void init();
  void updatezeta();
  double compute_scalar();
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);
  double memory_usage();

 private:
  int nmax;
  double *zeta;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
