This LAMMPSChiralTetramerModel.txt file was generated on 2021-08-01 by Nikolai D. Petsev

# GENERAL INFORMATION

1. Title of Research Code: Chiral Tetramer Molecular Model

2. Author Information
	A. Principal Investigator Contact Information
		Name: Nikolai D. Petsev 
		Institution: Princeton University Department of Chemical & Biological Engineering
		Address: Princeton, NJ 08540, USA
		Email: npetsev@princeton.edu

	B. Associate or Co-investigator Contact Information
		Names: Frank H. Stillinger, Pablo G. Debenedetti 
		Institution: Princeton University
		Address: Princeton, NJ 08540, USA
		Emails: fhs@princeton.edu, pdebene@princeton.edu

	Source code is based on the pair_lj_smooth_linear force-field by Jonathan Zimmerman and relies on LAMMPS - 
	Large-scale Atomic/Molecular Massively Parallel Simulator by Steve Plimpton. 
	Detailed breakdown of contributions for chiral model implementation:

	Original function by Jonathan Zimmerman, significant modifications by Nikolai Petsev:
	
    - PairLJSmoothLinearChiral::compute()
 
	Functions added by Nikolai Petsev:
	
	- PairLJSmoothLinearChiral::computezeta()
	- PairLJSmoothLinearChiral::pack_comm()  
	- PairLJSmoothLinearChiral::unpack_comm()  
	- PairLJSmoothLinearChiral::pack_reverse_comm()  
	- PairLJSmoothLinearChiral::unpack_reverse_comm()  
	- PairLJSmoothLinearChiral::memory_usage()  

	Original function by Jonathan Zimmerman, minor modifications by Nikolai Petsev:
	
	- PairLJSmoothLinearChiral::compute() 
	- PairLJSmoothLinearChiral::settings() 
	- PairLJSmoothLinear::coeff() 
	- PairLJSmoothLinearChiral::init_one()  
	- PairLJSmoothLinearChiral::write_restart()
	- PairLJSmoothLinearChiral::read_restart() 
	- PairLJSmoothLinearChiral::write_restart_settings()    
	- PairLJSmoothLinearChiral::read_restart_settings() 
	- PairLJSmoothLinearChiral::single()    

3. Date of research code production: 2020-08-01 to 2021-11-01 

4. Geographic location of code production: Henderson NV, USA 

5. Information about funding sources that supported the production of the code: National Science Foundation (grant CHE-1856704).


# SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the code:

    Chiral Tetramer Molecular Model, Copyright (2021) Nikolai D. Petsev. This software is distributed under the GNU
    General Public License, and is derived from LAMMPS - Large-Scale Atomic/Molecular Massively Parallel Simulator
    by Steve Plimpton, Copyright (2003) Sandia Corporation, https://www.lammps.org.

    Full license is included in the top-level directory.

2. Links to publications that cite or use the code: 

    Petsev et al. Effect of configuration-dependent multi-body forces on interconversion kinetics of a chiral tetramer model
 (JCP 2021) (https://doi.org/10.1063/5.0060266). 
    Wang et al. Vapor-liquid-liquid phase transition of a chiral tetramer model (in preparation 2021).
    Uralcan et al. Interconversion-controlled liquid-liquid phase separation in a molecular chiral model (in preparation 2021).

3. Links to other publicly accessible locations of the code: https://github.com/npetsev/LAMMPS-tetramer

4. Links/relationships to ancillary codebases: The chiral tetramer molcular model force field works with the LAMMPS simulation package (https://www.lammps.org/)

5. Was code derived from another source? 

    This software uses the LAMMPS - Large-Scale Atomic/Molecular Massively Parallel Simulator package by Steve Plimpton, Copyright (2003)
    Sandia Corporation, https://www.lammps.org. The modified/new code is derived from the pair_lj_smooth_linear force-field
    by Jonathan Zimmerman.

    For additional details on LAMMPS, please see the official website (https://www.lammps.org). A detailed description of LAMMPS can be found at:

    S. Plimpton, Fast Parallel Algorithms for Short-Range Molecular Dynamics, J Comp Phys, 117, 1-19 (1995). 

6. Recommended citation for this dataset: 

    Nikolai D. Petsev, Frank H. Stillinger, and Pablo G. Debenedetti, Chiral Tetramer Molecular Model. Princeton DataSpace. Available at http://arks.princeton.edu/ark:/88435/dsp01z890rx36m. Deposited June 2021.


# FILE OVERVIEW

1. File List: 

	- compute_zeta_ave.cpp -- compute for calculating the average chiral measure zeta of the system
	- compute_zeta_ave.h -- header file for zeta_ave compute
	- compute_zeta_ee.cpp -- compute for finding the system enantiomeric excess
	- compute_zeta_ee.h -- header file for ee compute
	- pair_lj_smooth_linear_chiral.cpp -- force calculation for the chiral tetramer model
	- pair_lj_smooth_linear_chiral.h -- header file for force calculation
	- style_compute -- header file that ensures we include added compute routines when compiling
	- style_pair -- header file that ensures we include chiral tetramer force interaction when compiling
	- in.NVE -- sample input file for performing simulation at constant NVE conditions
	- RestartMe.50000000 -- sample restart file
	- forcefield-base.dat -- sample input file with force field parameters

2. Relationship between files, if important: N/A 

3. Additional source code that was not included in the current package:

    The source code for the chiral tetramer molecular model must be compiled with the LAMMPS package, available at
    https://www.lammps.org.

4. Are there multiple versions of the code? Yes (LAMMPS 2014 compatible and LAMMPS 2020 compatible versions available)


# INSTALLATION

This force field has been tested and debugged using the 2014 and 2020 versions of LAMMPS. Note that the 2014 version is the stable and vetted release. To install:

1. Add files in this repository to main LAMMPS folder (style headers are optional, see step 2). Make sure you have downloaded source code compatible with your version of LAMMPS.

2. For convenience, we include style_pair.h and style_compute.h files, which can replace the existing ones in the LAMMPS directory. However, a safer route is to open the existing style_pair.h in your copy of LAMMPS, and add the following line:

		#include "pair_lj_smooth_linear_chiral.h"

Similarly, to use the average of zeta and ee computes, you can add the following lines to the style_compute.h file in your existing LAMMPS directory:

		#include "compute_zeta_ave.h" 
		#include "compute_zeta_ee.h"

3. Compile LAMMPS.


# VALIDATION

Detailed description of our validation methods is available in the accompanying publication. One important check to see if everything installed properly is to perform an NVE simulation and confirm that energy is being conserved. We provide a simple input script for such an NVE test. Note that a longer interaction cut-off (rc ~4.0) is required to give energy conservation.


# IMPLEMENTATION DETAILS

The source code is set up as follows (see accompanying publication for additional context) [brackets indicate function that performs indicated action]:

1. Create arrays with appropriate dimensions for zeta, and the x-, y-, and z-components of the gradient of zeta. We also create so-called mapping arrays (described in more detail below) [computezeta()].
2. The newly created arrays are initialized to zero [computezeta()].
3. Loop over the dihedral list to identify sets of 4 atoms that make up a single molecule (e.g. atom indeces i1, i2, i3, i4), and use this to calculate zeta and the gradient of zeta for i1, i2, i3, and i4.
4. Create our own “mapping array” that stores the same information as dihedrallist, but in a way that is transferrable across multiple processors by instead storing the global ids of i1, i2, i3, and i4 [computezeta()].
5. Perform forward and reverse communications of relevant quantities between processors. After we complete the loop through all dihedrals to update zeta, the gradient of zeta, and the mapping arrays, we perform a reverse communication, and then a forward communication. Note that order is important here. After performing these communications, zeta and all other arrays are now updated across all processors [computezeta(), pack_comm(), unpack_comm(), pack_reverse_comm(), unpack_reverse_comm()].
6. Execute force calculation loop. To compute this force for pair i and j, we identify which atoms are bonded to i and which are bonded to j using the mapping array that we constructed in computezeta(). The way we constructed the mapping arrays i1_vec, i2_vec, i3_vec, and i4_vec, we can find the global ids of the atoms that make up the molecule containing i using:

i1 = i1_vec[i]; i2 = i1_vec[i]; i3 = i1_vec[i]; i4 = i1_vec[i];

Note that i1, i2, i3, and i4 above are global ids. We can find the local ids using the LAMMPS mapping function atom-> map(). Note: the buffer is set up to handle values of type double, and hence we need to convert indeces to doubles when creating the mapping arrays. [compute()].

We include a check to ensure that atom-> map(i1_vec[i]) > -1, since we do not want to try to assign forces to atoms that are not known to the local processor.


# VARIABLE LIST

This list contains important variables there are not normally computed by LAMMPS or included in the original pair_lj_smooth_linear:

- zeta = array storing the chiral measure zeta value for each molecule (note that array length is same as the number of monomers in the system, and monomers belonging to the same - molecule have identical zeta values)
- dzetax = array storing the x-component of the gradient of zeta for each monomer in the system
- dzetay = array storing the y-component of the gradient of zeta for each monomer in the system
- dzetaz = array storing the z-component of the gradient of zeta for each monomer in the system
- i1_vec = mapping array storing i = 1 monomer ids
- i2_vec = mapping array storing i = 2 monomer ids
- i3_vec = mapping array storing i = 3 monomer ids
- i4_vec = mapping array storing i = 4 monomer ids
