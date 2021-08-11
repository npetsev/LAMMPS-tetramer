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

#include "mpi.h"
#include "compute_zeta_ee.h"
#include "atom.h"
#include "comm.h"
#include "update.h"
#include "force.h"
#include "domain.h"
#include "group.h"
#include "error.h"
#include "atom.h"
#include "neighbor.h"
#include "force.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeZetaEE::ComputeZetaEE(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all(FLERR,"Illegal compute ee command");

  scalar_flag = 1;
  extscalar = 1;

  // setup necessary array	
  zeta = NULL;

  // set comm size needed by this compute
  comm_forward = 1;
  comm_reverse = 1;
}

/* ---------------------------------------------------------------------- */

ComputeZetaEE::~ComputeZetaEE()
{	
  memory->destroy(zeta);
}

/* ---------------------------------------------------------------------- */

void ComputeZetaEE::init()
{
}

/* ---------------------------------------------------------------------- */

double ComputeZetaEE::compute_scalar()
{
  int nlocal = atom->nlocal;
  int ndihedrals = atom->ndihedrals;

  // update zeta values
  updatezeta();

  // sum over all processors
  double nl = 0.0;
  for (int i = 0; i < nlocal; i++) {
    if (zeta[i] < 0.0) nl += 1.0; 
  }

  MPI_Allreduce(&nl,&scalar,1,MPI_DOUBLE,MPI_SUM,world);
  scalar = 2.0 * scalar / double(4*ndihedrals) - 1.0;

  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeZetaEE::updatezeta()
{	
  int i1,i2,i3,i4;
  int i,m,n,type;
  double vb1x,vb1y,vb1z,vb2x,vb2y,vb2z,vb3x,vb3y,vb3z;
  double b1mag2,b1mag,b2mag2,b2mag,b3mag2,b3mag;
  double c1f,c2f,c3f,zetaf; 

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int newton_bond = force->newton_bond;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
 
  nmax = atom->nmax;
  memory->create(zeta,nmax,"zeta/ee:zeta");

  // initialize zeta/dzeta arrays to zero
  if (newton_bond) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) {
      zeta[i] = 0.0;
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      zeta[i] = 0.0;
    }
  }

  // loop through list of dihedrals to find zeta/dzeta
  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    // 1st bond
    vb1x = x[i2][0] - x[i1][0];
    vb1y = x[i2][1] - x[i1][1];
    vb1z = x[i2][2] - x[i1][2];

    // 2nd bond
    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    // 3rd bond
    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    // 1st and 2nd angle
    b1mag2 = vb1x*vb1x + vb1y*vb1y + vb1z*vb1z;
    b1mag = sqrt(b1mag2);
    b2mag2 = vb2x*vb2x + vb2y*vb2y + vb2z*vb2z;
    b2mag = sqrt(b2mag2);
    b3mag2 = vb3x*vb3x + vb3y*vb3y + vb3z*vb3z;
    b3mag = sqrt(b3mag2);

    c1f = (vb2y*vb3z - vb2z*vb3y);
    c2f = (vb2z*vb3x - vb2x*vb3z);
    c3f = (vb2x*vb3y - vb2y*vb3x);

    zetaf = vb1x*c1f + vb1y*c2f + vb1z*c3f;
    zetaf = zetaf/(b1mag*b2mag*b3mag);

    // store in zeta gradient arrays
    zeta[i1] = zetaf;
    zeta[i2] = zetaf;
    zeta[i3] = zetaf;
    zeta[i4] = zetaf;    
  }

  // communicate zeta values between procs
  if (newton_bond) comm->reverse_comm_compute(this);
  comm->forward_comm_compute(this);
}

/* ----------------------------------------------------------------------
   forward/reverse packing and unpacking routines for proc communication
------------------------------------------------------------------------- */

int ComputeZetaEE::pack_comm(int n, int *list, double *buf, int pbc_flag,
    int *pbc) {
  int i, j, m;
	    
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = zeta[j];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeZetaEE::unpack_comm(int n, int first, double *buf) {
  int i, m, last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    zeta[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int ComputeZetaEE::pack_reverse_comm(int n, int first, double *buf)
{
  int i, m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = zeta[i];
  }
  return 1;
}

/* ---------------------------------------------------------------------- */

void ComputeZetaEE::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    zeta[j] += buf[m++];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double ComputeZetaEE::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}