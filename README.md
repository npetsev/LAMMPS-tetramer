# LAMMPS-tetramer
Molecular chiral model force field for LAMMPS simulation package.

Note that files in repository are not standalone; the included files
must be added to LAMMPS source coded and then compiled.

----------INSTALLATION----------

This force field should work with any version of LAMMPS, though 
it has been tested and debugged using 2013 version. To install:

1) Add files in this repository to main LAMMPS folder (style headers 
are optional, see step 2).

2) For convenience, we include style_pair.h and style_compute.h files,
which can replace the existing ones in the LAMMPS directory. However,
this may not be compatible with all versions of LAMMPS (different versions
include different pair and compute styles). A safer route is to simply
open the version of style_pair.h in your version of LAMMPS, and add the
following line:

#include "pair_lj_smooth_linear_chiral.h"

Similarly, to use the average of zeta and ee computes, you can add the
following lines to the style_compute.h file in your existing LAMMPS
directory:

#include "compute_zeta_ave.h"
#include "compute_zeta_ee.h"

3) Compile LAMMPS.

--------------NOTES-------------

1) Version of the force field in this repository uses a cut and shifted
potential; the shifted term can be commented out, though without a 
shifted force, there is noticable energy drift unless a prohibitively large
cut-off is used. The cut and shifted potential has been verified to
conserve energy exactly.

2) The computes for finding the average of the chirality measure and 
enantiomeric excess perform an average over all the molecules. However,
some versions of LAMMPS automatically normalize quantities found from
computes, which can result in incorrect output from the zeta and 
ee computes. To avoid this issue, include in the input script:

thermo_modify norm no

---------IMPLEMENTATION---------

The source code is set up as follows [brackets denote which function includes 
the operation being discussed]:

1)	We create arrays with appropriate dimensions for zeta, and the x-, y-, 
and z-components of the gradient of zeta. We also create so-called mapping arrays 
(described in more detail below) [computezeta()].

2)	The newly created arrays are initialized to zero [computezeta()].

3)	We loop over the dihedral list to identify sets of 4 atoms that make up a 
single molecule (e.g. atom indeces i1, i2, i3, i4), and we use this to calculate 
zeta and the gradient of zeta for i1, i2, i3, and i4. 

4)	We need to access which atoms belong to a given dihedral when we do the force 
computation. However, due to how the dihedral is stored across multiple processors, 
this can lead to issues where an atom is owned on one processor, but the dihedral 
it belongs to is stored on a different processor. Therefore, we create our own 
“mapping array” that stores the same information as dihedrallist, but in a way that 
is transferrable across multiple processors by instead storing the global ids of 
i1, i2, i3, and i4 [computezeta()].

5)	During our loop over dihedrals, it is possible that we assigned zeta (or other 
array) values to ghost atoms. It is also possible that some locally owned atoms had 
their zeta value updated on a different processor. In light of this, we need to 
perform forward and reverse communications between processors. Therefore, after we 
complete the loop through all dihedrals to update zeta, the gradient of zeta, and 
the mapping arrays, we perform a reverse communication, and then a forward 
communication. Note that order is important here! After performing these 
communications, zeta and all other arrays are now updated across all processors 
[computezeta(), pack_comm(), unpack_comm(), pack_reverse_comm(), unpack_reverse_comm()].

6)	We now have access to all relevant information to perform the force loop, and hence 
we loop across i-j pairs. To compute this force for pair i and j, we identify which 
atoms are bonded to i and which are bonded to j. In order to do this, we use the mapping 
array that we constructed in computezeta(). The way we constructed the mapping arrays 
i1_vec, i2_vec, i3_vec, and i4_vec, we can find the global ids of the atoms that make up 
the molecule containing i using:

i1 = i1_vec[i];
i2 = i1_vec[i];
i3 = i1_vec[i];
i4 = i1_vec[i];

Note that i1, i2, i3, and i4 above are global ids. We can find the local ids using the 
LAMMPS mapping function atom-> map(). Small detail: the buffer is set up to handle values 
of type double, and hence we actually convert indeces to doubles when creating the mapping 
arrays. We can convert back to integer using the function int() [compute()].

7)	We can similarly identify which local atom ids correspond to atoms that are part of 
the same molecule as j. At this point it is possible to compute the relevant forces acting 
on all atoms. Note that we include a check to ensure that atom-> map(i1_vec[i]) > -1, 
since we do not want to try to assign forces to atoms that are not known to the local 
processor.



