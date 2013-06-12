CIF2CELL RELEASE INFORMATION


VERSION 1.0.12

* Fixed security issues related to vectors/matrices input from the 
  command line.
* New interface: CP2K
* New interface: .coo files
* New interface: FHI-AIMS
* Minor fixes and tweaks

--------------------------------------------------------------   
VERSION 1.0.10

* Critical bugfix in non-diagonal supercell generation.

--------------------------------------------------------------   
VERSION 1.0.9

* Fixed bug in the CIF output when exporting the primitive
 unit cell.
* Fixed bug in the VASP output when exporting unit cells with 
 left-handed set of lattice vectors.
(with thanks to Jens Kunstmann)

--------------------------------------------------------------  
VERSION 1.0.8

* The code now correctly destroys the symmetry of the system when
  using the --random-displacements flag. Use --random-displacements=0.0
  to remove all symmetries from the output  without distorting the structure.
* Output to screen has been made consistent with the above change.
* Included CIF's for the crystal structures of all elements in the 
  periodic table to the CIF collection.

-------------------------------------------------------------- 
VERSION 1.0.7

* Bug with supercell generation with python 2.7.
* Miscellaneous other minor fixes.

-------------------------------------------------------------- 
VERSION 1.0.6

* Fixed bug with rhombohedral/hexagonal settings in the CRYSTAL interface.

-------------------------------------------------------------- 
VERSION 1.0.5

* Implemented support for virtual crystal approximation setups of alloys in CASTEP.
* New functionality for CASTEP interface: --castep-cartesian and --castep-atomic-units.
* The CASTEP interface now supplies a full, commented out, pseudopotential block for
  easy editing.
* The program can now take arguments and options in any order.
* Fixed buggy behaviour when inconsistent space group symbols are given.
* Allow for '?' and '.' for unknown space group symbols.
Thanks again to Keith Refson for bug reports and help with CASTEP features.

-------------------------------------------------------------- 
VERSION 1.0.2

* Extended support for --force flag.
* Bugfix in CRYSTAL interface.
* Bugfix in EMTO interface.

-------------------------------------------------------------- 
VERSION 1.0.0

MAJOR NEW RELEASE. 
CIF2Cell leaves beta stage. Wohoo! 

New features include:
* Usable manual/tutorial (in pdf format).
* Major rework of the supercell generator to support general
  map matrices. This should make it possible to generate any
  possible supercell.
* Possibility to realign the unit cell vectors. In combination with
  the general map matrix, this makes it easy to generate surface
  supercells. For convenience, there are also two predefined rotations
  for aligning the cubic (111) direction with the z axis and for
  aligning the threefold rotation axis of a rhombohedral system 
  with the pseudocubic (111) direction.
* Possibility to add random displacements to all atoms.
* Support for xyz format output.
* An experimental, more general --force flag to try to enforce cell
  creation despite any problems encountered. Not completely 
  supported yet.

-------------------------------------------------------------- 
VERSION 0.4.5

* New interface to xband/SPRKKR.
* Fixed yet another bug with rhombohedral settings.

-------------------------------------------------------------- 
VERSION 0.4.4.4

* Fixed bug in the Siesta interface.
* Fixed bug with naming of output files for RSPt. 
* Added non-standard settings.
* Improved handling of species designated with '?' or '.'.

-------------------------------------------------------------- 
VERSION 0.4.4.2

Another small patch fixing some problems provided/omitted by the
previous patch related to the H-M symbols.

-------------------------------------------------------------- 
VERSION 0.4.4_1

A small patch fixing that a lot of Hermann-Mauguin symbols with the 
standard choice of origin and unique axes where not supported in their 
abbreviated form (such as C2/m for C12/m1). 

-------------------------------------------------------------- 
VERSION 0.4.4

Back on track again after the little mishap in version 0.4.3. 
Mostly, this release greatly extends the support for non-standard 
space group settings. The program now handles almost all of the 
ICSD database, and should detect any setting that is not handled 
and exit with an error message. 

Under the hood, the handling of the space group data has been
extensively modified to make it easy to handle non-standard settings. 
The program now uses Hall symbols internally, because all the time
having to parse the poorly standardized and/or non-unique
Hermann-Mauguin symbols almost drove me insane.

-------------------------------------------------------------- 
VERSION 0.4.3

Due to problems with covering all possibilites for writing
Hermann-Mauguin symbols, version 0.4.3 was unable to support
large portions of popular CIF databases, and was removed
from the downloading section here. This will be fixed shortly
in an upcoming release, meanwhile use 0.4.2.

* Reworked the space group information stored internally to
  more flexible and pythonic formats.
* The code should now handle all reasonably normal settings 
  (anything selectable in the Bilbao Crystallographic server) 
  and choke in a controlled way on anything else (things like
  face centered monoclinic settings, which are by convention
  represented in an equivalent base-centered monoclinic setting).

With thanks to Henning Glawe for reporting some irregular behaviour.

-------------------------------------------------------------- 
VERSION 0.4.2

Fixed serious bug in the cell generation. All users should immediately 
upgrade from version 0.4.1.

-------------------------------------------------------------- 
VERSION 0.4.1

Major new release.
* A lot of tidying up under the hood.
* Now comes with an unfinished manual.
* A --setup-all flag that attempts to do a more "complete" setup. 
  Presently only available for VASP.
* Sorting of atoms in a supercell is possible.
* Primitive/conventional cell treatment completely consistent
  for trigonal/rhombohedral systems. By default the minimal
  rhombohedral cell is chosen, and with --no-reduce you get 
  the hexagonal cell.
* Outputting charge state. 
* Exporting reference data in BibTeX format.
* A LOT of bugfixes.

-------------------------------------------------------------- 
VERSION 0.2.6

* Fixed problem in setup script that made some files not install properly.
* Bugfixes in supercell generation.
* Speeding up the program a bit (previous version could be very
  slow in some circumstances).

-------------------------------------------------------------- 
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
