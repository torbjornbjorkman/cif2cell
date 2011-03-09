CIF2CELL RELEASE INFORMATION

VERSION 0.2.5

* Setting up symmetry operations in matrix/translation vector format+
  output of these to CASTEP. NOTE: Not yet available for supercells.
* Bugfixes in cell generation without space group information 
  and EMTO interface.
* Some tidying up of the source code, splitting over more files,
  introducing some more convenient classes etc.

-------------------------------------------------------------- 
VERSION 0.2.3

* Bugfixes in the CASTEP interface (thanks to Keith Refson).

-------------------------------------------------------------- 
VERSION 0.2.2

* Added reciprocal lattice vectors method to CrystalData (useful 
for setting up a k-space mesh).
* Bugfixes.

-------------------------------------------------------------- 
VERSION 0.2.1

* Fixed problem with Python < 2.6

--------------------------------------------------------------
VERSION 0.2.0

* Translation of supercell possible.
* Possibility of adding vacuum added to the supercell
  generator.
* New input format for RSPt program.

--------------------------------------------------------------
VERSION 0.1.1

* A couple of minor bugfixes.
* Improved error handling.
* Added possibility to specify which CIF grammar to use.

--------------------------------------------------------------
VERSION 0.1 

This is the first official release. Program features:
* Generation of principal or primitive cell.
* Generation of supercells.
* Output for : ABINIT, Siesta, CPMD, CASTEP, Crystal09, 
  elk, EMTO, exciting, Fleur, ncol, RSPt, Siesta and VASP.
