[![Build Status](https://github.com/torbjornbjorkman/cif2cell/workflows/ci/badge.svg)](https://github.com/torbjornbjorkman/cif2cell/actions)

A tool to generate the geometrical setup for various electronic
structure codes from a CIF (Crystallographic Information
Framework) file. The code will generate the crystal structure for
the primitive cell or the conventional cell.

## CURRENTLY SUPPORTS

|code           | alloy support |   output files|
|---------------|---------|-----------------------------------|
|ASE            |   no    | positions.py|
|ATAT           |  yes    | [compoundname].in|
|VASP           |  VCA    | POSCAR|
|ABINIT         |   no    | [compoundname].in|
|Siesta         |   no    | [compoundname].fdf|
|CPMD           |   no    | [compoundname].inp|
|CASTEP         |  VCA    | [compoundname].cell|
|Crystal09      |   no    | [compoundname].d12|
|quantum espresso|  no    | [compoundname].in|
|FHI-aims       |   no    | geometry.in|
|RSPt           |   no    | symt.inp|
|Fleur          |   no    | inp_[compoundname]|
|hutsepot       |   no    | [compoundname].sys|
|cellgen        |   no    | cellgen.inp|
|elk            |   no    | GEOMETRY.OUT|
|exciting       |   no    | input.xml|
|spacegroup     |   no    | spacegroup.in|
|ncol           |   no    | [spacegroupname/compoundname].dat|
|               |         | for bstr.|
|emto           |   yes   | [spacegroupname/compoundname].dat|
|               |         | for kstr, bmdl, shape, kgrn and kfcd|
|               |         | in separate directories.|
|spr-kkr        |   yes   | [compoundname].sys|
|xyz            |   no    | [compoundname].xyz|


## CONTENTS

The repository includes:

* This README file.
* The file LICENSE with the GPLv3 license.
* The python files cif2cell, uctools.py and spacegroupdata.py
* Installation files, setup.py and MANIFEST.
* A manual.
* The directory cifs/ containing a set of example CIF files
  as well as the crystal structures of the full periodic table
  from COD, the Crystallography Open Database <http://www.crystallography.net>
  and also a few from ICSD (with permission).


## INSTALLATION INSTRUCTIONS

### Prerequisites

The program requires Python 2.4 or higher and the PyCIFRW python package (which
will be installed automatically if not present).
Note however that the output may be slightly different (but formally
equivalent) with Python 2.4 than with later Python versions.


```
pip install cif2cell
```

The installation will also create a directory $PREFIX/lib/cif2cell
that contains the manual and sample cif files.


## DOCUMENTATION

The setup will install the manual, cif2cell.pdf, into the
$PREFIX/lib/cif2cell/docs directory.


## RUNNING

Run `cif2cell -h` to get a listing of the different options.
Example:

```
cif2cell Ni20Mn3P6.cif -p vasp --vasp-cartesian-positions
```

will generate a POSCAR file for VASP with the positions in cartesian format.


## LICENSE INFORMATION

cif2cell is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cif2cell is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.

## HOW TO CITE

Please use the following citation information:

Torbjorn Bjorkman, "CIF2Cell: Generating geometries for electronic structure programs",
Computer Physics Communications 182, 1183-1186 (2011)
doi: [10.1016/j.cpc.2011.01.013](https://doi.org/10.1016/j.cpc.2011.01.013)

My name is rendered in ascii above, bonus points for getting umlauts over both of the o's.
See also below for a BibTeX entry for use with LaTeX, which should be readable
for most scientific reference handling software.

```
@article{cif2cell,
title = "CIF2Cell: Generating geometries for electronic structure programs",
journal = "Computer Physics Communications",
volume = "182",
number = "5",
pages = "1183 - 1186",
year = "2011",
issn = "0010-4655",
doi = "10.1016/j.cpc.2011.01.013",
url = "http://www.sciencedirect.com/science/article/pii/S0010465511000336",
author = "Torbj\"orn Bj\"orkman"
}
```




Happy computing!

Torbjorn Bjorkman
https://orcid.org/0000-0002-1154-9846
