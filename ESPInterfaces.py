# Copyright 2010 Torbjorn Bjorkman
# This file is part of cif2cell
#
# cif2cell is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cif2cell is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cif2cell.  If not, see <http://www.gnu.org/licenses/>.
#
#******************************************************************************************
#  Description: Interfaces for a number of electronic structure programs. Currently only
#               reads CIF and outputs to the ESP's. Supported programs are: ABINIT, CASTEP,
#               CPMD, Crystal09, Elk, EMTO, Exciting, Fleur, NCOL, RSPt, Siesta, VASP
#               
#  Author:      Torbjorn Bjorkman, torbjorn.bjorkman(at)aalto.fi
#  Affiliation: COMP, Aaalto University School of Science,
#               Department of Applied Physics, Espoo, Finland
#******************************************************************************************
import copy
import os
import math
from utils import *
from elementdata import *
from uctools import *
from re import search

ed = ElementData()

################################################################################################
class GeometryOutputFile:
    """
    Parent class for electronic struture code files generated from geometrical information.
    A CrystalStructure object and a documentation string are required input.
    """
    def __init__(self, crystalstructure, string):
        self.cell = crystalstructure
        self.docstring = string

################################################################################################
# ASE FILE
class ASEFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting data for ASE
    and the method __str__ that outputs the contents of the file as a string of
    python code.
    """
    def __init__(self,crystalstructure,docstring):
        GeometryOutputFile.__init__(self,crystalstructure,docstring)
        # Variables
        self.cartesian = True  # Cartesian coordinates?
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
    def __str__(self):
        filestring = "from ase import *\n\n"
        # Cartesian or lattice coordinates?
        if self.cartesian:
            transmtx = []
            for i in range(3):
                transmtx.append([])
                for j in range(3):
                    transmtx[i].append(self.cell.latticevectors[i][j] * self.cell.lengthscale)
        else:
            transmtx = [[1,0,0],[0,1,0],[0,0,1]]
        # positions and number of species
        nspcs = []
        positionstring = ""
        for sp in self.species:
            nsp = 0
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        nsp += 1
                        p = Vector(mvmult3(transmtx,b.position))
                        positionstring += "%11f, %11f, %11f),\n               ("%(p[0],p[1],p[2])
            nspcs.append(nsp)
        positionstring = positionstring.rstrip("\n (,")+"],\n"

        # Atoms object
        filestring += "atoms = Atoms("
        # Species
        for i in range(len(self.species)):
            filestring += "['"+self.species[i]+"' for i in range("+str(nspcs[i])+")]+"
        filestring = filestring.rstrip("+")+",\n"
        # Positions
        filestring += "              [("
        filestring += positionstring
        # Boundary conditions
        filestring += "              pbc = (True,True,True))\n"
        # Set lattice vectors
        filestring += "atoms.set_cell([["
        for i in range(3):
            for j in range(3):
                filestring += "%12f, "%(self.cell.latticevectors[i][j]*self.cell.lengthscale)
            filestring = filestring.rstrip(", ")+"],\n                ["
        filestring = filestring.rstrip(",[ ]\n")+"]],\n"
        filestring += "                scale_atoms = %s)\n"%str(not self.cartesian)
        
        return filestring
    
################################################################################################
# COO FILE
class COOFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting a .coo file
    and the method __str__ that outputs the contents of the .coo file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Document string on first line after '//'
        self.programdoc = string.rstrip("\n")
    def __str__(self):
        filestring = "//"+self.programdoc+"\n"
        a = self.cell.latticevectors[0].length()*self.cell.lengthscale
        b = self.cell.latticevectors[1].length()*self.cell.lengthscale
        c = self.cell.latticevectors[2].length()*self.cell.lengthscale
        alpha = abs(self.cell.latticevectors[1].angle(self.cell.latticevectors[2]))*180/pi
        beta = abs(self.cell.latticevectors[2].angle(self.cell.latticevectors[0]))*180/pi
        gamma = abs(self.cell.latticevectors[0].angle(self.cell.latticevectors[1]))*180/pi
        filestring += " %10.7f %10.7f %10.7f"%(a,b,c)
        filestring += " %10.7f %10.7f %10.7f %i\n"%(alpha,beta,gamma,len(self.cell.atomset))
        for a in self.cell.atomdata:
            for b in a:
                filestring += str(b.position)+" %3i  0.500 0.000 1.000\n"%(ed.elementnr[b.spcstring()])
        return filestring
    
################################################################################################
# XYZ FILE
class XYZFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed for outputting an .xyz file
    and the method __str__ that outputs the contents of the .xyz file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # To be put on the second line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        filestring += "%i \n"%sum([len(v) for v in self.cell.atomdata])
        filestring += self.docstring+"\n"
        lv = []
        for i in range(3):
            lv.append([])
            for j in range(3):
                lv[i].append(self.cell.lengthscale*self.cell.latticevectors[i][j])
        for a in self.cell.atomdata:
            for b in a:
                t = Vector(mvmult3(lv,b.position))
                filestring += str(b).split()[0]+"  "+str(t)+"\n"
        return filestring
    
################################################################################################
# NCOL FILES
class OldNCOLFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the ncol program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.bstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
        # Set atomic units for length scale
        self.cell.newunit("bohr")
    def __str__(self):
        # Element data
        ed = ElementData()
        # l quantum number setup (same as from bstr)
        l = { "s" : 2, "p" : 2, "d" : 3, "f" : 4 }
        filestring = ""
        tmpstring = "BULK      IDSYST=  7 SCRATCH=R"
        tmpstring = tmpstring.ljust(25)+"    "+deletenewline(self.programdoc,replace=" ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 BSAVE..=N COLD...=Y DOS...=N SPO...=N ISM...=G RCLCR...=Y\n"
        filestring += tmpstring
        filestring += "FOR001=./"+self.bstrjobnam+".tfm\n"
        filestring += "FOR002=\n"
        filestring += "FOR003=\n"
        filestring += "FOR004=\n"
        filestring += "FOR006=\n"
        filestring += "FOR010=\n"
        filestring += "Band: 4 lines, "+deletenewline(self.docstring,replace=" ")+"\n"
        filestring += "NITER.=200 NOB..=  2 NPRN.=  0 NFIX.=  0 MIXKEY=  2 NCOL.=Y  PMODE=K\n"
        filestring += "REP.....=B FIXD...=Y CRT....=S NB...= 16 CLSIZE= 32 NPROW= 0 NPCOL= 0\n"
        filestring += "NKX...=  1 NKY..=  1 NKZ..=  1 TFERMI..= 2000.0(K)\n"
        filestring += "AMIX.....=     0.100 TOLE....= 0.0000100 TOLEL...= 0.0000010\n"
        # average wigner-seitz radius
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        volume = abs(det3(self.cell.latticevectors))
        wsr = self.cell.lengthscale * (3*volume/(nosites * 4 * pi))**third
        filestring += "SWS......= %9f NP...=  1 SMIX.= 0.500 TMIX.= 0.0000\n" % wsr
        filestring += "Setup: 3 + NQ*NS*NP lines\n"
        filestring += "EFGS.....=    0.0000 EFGS....=   0.00000 FTMAG...=  0.000000\n"
        filestring += "DEO(l)...=     0.020     0.010     0.005     0.001      0.02\n"
        filestring += "Symb IQ IT NL IP NSP   SWP  QTRO  SPLT NFIX NDWF     Eny(spdf)\n"
        # set first species
        if self.cell.atomdata[0][0].alloy():
            prevspecies = "??"
        else:
            for v in self.cell.atomdata[0][0].species:
                prevspecies = v
        # type loop
        iq = 1
        it = 1
        nsp = 1
        for a in self.cell.atomdata:
            for b in a:
                if b.alloy():
                    species = "??"
                else:
                    species = b.spcstring()
                if species != prevspecies:
                    prevspecies = species
                    nsp += 1
                tmpstring = species.ljust(2)+"    "+str(iq).ljust(3)+str(it).ljust(3)
                try:
                    tmpstring += str(l[ed.elementblock[species]]).ljust(3)+str(1).ljust(3)
                except KeyError:
                    tmpstring += "  ?  1"
                tmpstring += str(nsp).ljust(3)
                tmpstring += "  1.000 .000 0.00 0000 1111   .0   .0   .0   .0"
                if b.alloy():
                    # print alloy components at the end of the line
                    tmpstring += "       "+b.spcstring()
                filestring += tmpstring+"\n"
                iq += 1
            it += 1
        for a in self.cell.atomdata:
            filestring += "Theta....=     90.00 Phia....=      0.00 FIXMOM..=         N moment..=      0.0\n"
        filestring += "PQX......=      0.00 PQY.....=      0.00 PQZ.....=   0.00000 COORD...=L\n"
        filestring += "Atom: 4 lines + NT*6 lines\n"
        filestring += "IEX...=  4  NP..=500 NES..= 15 NITER=250 IWAT.=  0\n"
        filestring += "VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n"
        filestring += "DPAS.....=  0.049000 DR1.....=  1.00E-08 TEST....=  1.00E-08\n"
        filestring += "TESTE....=  1.00E-07 TESTY...=  1.00E-08 TESTV...=  1.00E-07\n"
        for a in self.cell.atomdata:
            for comp in a[0].species:
                filestring += comp+"\n"
                try:
                    filestring += ed.emtoelements[comp]
                except KeyError:
                    filestring += "\n\n\n\n\n"
        return filestring

class BSTRFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the bstr program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.a = 1
        self.b = 1
        self.c = 1
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        ed = ElementData()
        filestring = ""
        tmpstring = "BSTR      IDSYST=  7"
        tmpstring = tmpstring.ljust(40)+deletenewline(self.programdoc,replace=" ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(9)+" MSGL.=   1 \n"
        filestring += tmpstring
        filestring += "MODE....=B STORE..=Y SCREEN.=B CMBC...=Y\n"
        filestring += "FOR001=\n"
        filestring += "FOR006=\n"
        filestring += deletenewline(self.docstring,replace=" ")+"\n"
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        # Setting the real space summation cutoff to 4.5*(wigner-seitz radius)
        volume = abs(det3(self.cell.latticevectors))
        wsr = (3*volume/(nosites * 4 * pi))**third
        tmpstring = "IALF...= 0 NPRN..= 1 DMAX....=%10.5f \n" % (wsr*4.5)
        filestring += tmpstring
        filestring += "ALF(spdf)= 0.3205350 0.0413320 0.0084290 0.0015370\nDKAPPA...= 0.00010\n"
        tmpstring = "NQ3....=%3i LAT...= 0 IPRIM.= 0" % nosites
        filestring += tmpstring
        # Set up basis functions. Just setting lmax = 2 for s-/p-, 3 for d- and 4 for f- blocks
        tmpstring = "\nNLX(IQ)..="
        for a in self.cell.atomdata:
            for b in a:
                for k in b.species:
                    l = 1
                    if ed.elementblock[k] == "s" or ed.elementblock[k] == "p":
                        l = max(l,2)
                    elif ed.elementblock[k] == "d":
                        l = max(l,3)
                    elif ed.elementblock[k] == "f":
                        l = max(l,4)
                tmpstring += " %1i" % l
                if len(tmpstring) % 69 == 0:
                    tmpstring += "\n          "
        # Need to strip newline character if the last line was 69 characters long...
        tmpstring = tmpstring.rstrip(string.whitespace)
        tmpstring = tmpstring+"\n"
        filestring += tmpstring
        # Print lattice vectors
        coa = self.c / self.a
        boa = self.b / self.a
        filestring += "A........=  1.00000000 B.......=  1.00000000 C.......=  1.00000000\n"
        tmpstring = ""
        lv = self.cell.latticevectors
        for i in range(3):
            tmpstring += "BSX......=%12.7f BSY.....=%12.7f BSZ.....=%12.7f\n" % (lv[i][0],lv[i][1],lv[i][2])
        filestring += tmpstring
        # All positions
        it = 1
        for a in self.cell.atomdata:
            for b in a:
                pos = mvmult3(lv,b.position)
                tmpstring = "QX.......=%12.7f QY......=%12.7f QZ......=%12.7f" % (pos[0],pos[1],pos[2])
                tmpstring += "      "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        filestring += "LAMDA....=    2.5000 AMAX....=    5.5000 BMAX....=    5.5000\n"
        return filestring
    
################################################################################################
# RSPT FILES
class CellgenFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a cellgen.inp file and the method
    __str__ that outputs the contents of an cellgen.inp file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.supercellmap = [[1,0,0],[0,1,0],[0,0,1]]
        self.referencevector = [0,0,0]
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.: "+str(self.cell.lengthscale)+"\n"
        # RSPt reads the lattice vectors as columns...
        filestring +="# Lattice vectors (columns)\n"
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%self.cell.latticevectors[j][i]
            tmpstring += "\n"
        filestring += tmpstring
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        filestring += "# Sites\n"
        filestring += str(nosites)+"\n"
        it = 1
        for a in self.cell.atomdata:
            for b in a:
                tmpstring = ""
                tmpstring += str(b.position)+" "
                if b.alloy():
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    tmpstring += "%3i"%ed.elementnr[b.spcstring()]
                tmpstring += " l "+chr(it+96)+"   # "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        filestring += "# Supercell map\n"
        tmpstring = ""
        for i in self.supercellmap:
            for j in i:
                tmpstring += str(j).rjust(4)
            tmpstring += "\n"
        filestring += tmpstring
        filestring += "# Reference vector\n"
        tmpstring = ""
        for i in self.referencevector:
            tmpstring += "%19.15f "%i
        tmpstring += "\n"
        filestring += tmpstring
        return filestring

################################################################################################
class SymtFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an old format symt.inp file and the method
    __str__ that outputs the contents of an symt.inp file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # Default spin axis is [0,0,0]
        self.spinaxis = [0.0, 0.0, 0.0]
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.: "+str(self.cell.lengthscale)+"\n"
        # RSPt reads the lattice vectors as columns...
        filestring +="# Lattice vectors (columns)\n"
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%self.cell.latticevectors[j][i]
            tmpstring += "\n"
        filestring += tmpstring
        filestring += "# Spin axis\n"
        filestring += "%19.15f %19.15f %19.15f  l\n"%(self.spinaxis[0],self.spinaxis[1],self.spinaxis[2])
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        filestring += "# Sites\n"
        filestring += str(nosites)+"\n"
        it = 1
        for a in self.cell.atomdata:
            for b in a:
                tmpstring = ""
                tmpstring += str(b.position)+" "
                if b.alloy():
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    tmpstring +=  "%3i"%ed.elementnr[b.spcstring()]
                tmpstring += " l "+chr(it+96)+"   # "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        return filestring

################################################################################################
class SymtFile2(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a new format symt.inp file and the method
    __str__ that outputs the contents of an symt.inp file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # Default spin axis is [0,0,0]
        self.spinaxis = [0.0, 0.0, 0.0]
        # parameters for spin polarization
        self.spinpol = False
        self.relativistic = False
        self.mtradii = 0
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.\n"
        filestring += "lengthscale\n"
        filestring += str(self.cell.lengthscale)+"\n"
        if self.spinpol:
            filestring += "# Spin polarized calculation\nspinpol\n"
            filestring += "# Spin polarize atomic densities\nspinpol_atomdens\n"
        if self.relativistic:
            filestring += "# Relativistic symmetries\nrelativistic\n"
        if self.mtradii != 0:
            filestring += "# Choice of MT radii\n"
            filestring += "mtradii\n"+str(self.mtradii)+"\n"
        # RSPt reads the lattice vectors as columns...
        filestring += "# Lattice vectors (columns)\n"
        filestring += "latticevectors\n"
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%self.cell.latticevectors[j][i]
            tmpstring += "\n"
        filestring += tmpstring
        filestring += "# Spin axis\n"
        filestring += "spinaxis\n"
        filestring += "%19.15f %19.15f %19.15f  l\n"%(self.spinaxis[0],self.spinaxis[1],self.spinaxis[2])
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        filestring += "# Sites\n"
        filestring += "atoms\n"
        filestring += str(nosites)+"\n"
        it = 1
        for a in self.cell.atomdata:
            for b in a:
                tmpstring = ""
                tmpstring += str(b.position)+" "
                if b.alloy():
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    tmpstring += "%3i"%ed.elementnr[b.spcstring()]
                tmpstring += " l "+chr(it+96)+"   # "+b.spcstring()+"\n"
                filestring += tmpstring
            it += 1
        return filestring

################################################################################################
class Crystal09File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed by Crystal09 and the method
    __str__ that outputs the contents of an Crystal09 input file as a string.
    Presently only handles standard settings (space group numbers, not H-M symbols)
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.HermannMauguin = ""
        self.spacegroupnr = 0
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.trigonalsetting = "H"
        # Set atomic units for length scale
        self.cell.newunit("angstrom")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        filestring += "CRYSTAL\n"
        system = crystal_system(self.spacegroupnr)
        # Space group setting and crystal parameters
        if system == "triclinic":
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        elif system == "monoclinic":
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c, self.beta)
        elif system == "orthorhombic":
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c)
        elif system == "tetragonal":
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f\n"%(self.a, self.c)
        elif system == "trigonal":
            if self.trigonalsetting == "H":
                filestring += "0 0 0\n"
                filestring += str(self.spacegroupnr)+"\n"
                filestring += "%13.8f %13.8f\n"%(self.a, self.c)
            elif self.trigonalsetting == "R":
                filestring += "0 1 0\n"
                filestring += str(self.spacegroupnr)+"\n"
                filestring += "%13.8f %13.8f\n"%(self.a, self.alpha)
            else:
                return "***Error: No such trigonal setting : "+self.trigonalsetting
        elif system == "hexagonal":
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f %13.8f\n"%(self.a, self.c)
        elif system == "cubic":
            filestring += "0 0 0\n"
            filestring += str(self.spacegroupnr)+"\n"
            filestring += "%13.8f\n"%(self.a)
        else:
            if self.force:
                sys.stderr.write("***Warning: Could not determine crystal system corresponding to space group "+str(self.spacegroupnr)+".")
                filestring += "0 0 0\n"
                filestring += str(self.spacegroupnr)+"\n"
                filestring += "%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
            else:
                return "***Error: Could not determine crystal system corresponding to space group "+str(self.spacegroupnr)+"."
        # Number of atoms
        filestring += str(len(self.cell.ineqsites))+"\n"
        # Atomic numbers and representative positions
        for a in self.cell.atomdata:
            if len(a[0].species) > 1:
                # don't know what to put for an alloy
                filestring += "??"
            else:
                for k in a[0].species:
                    filestring += str(ed.elementnr[k]).rjust(2)
            filestring += "  "+str(a[0].position)+"      ! "+a[0].spcstring()+"\n"
        ## occ = self.cell.occupations
        ## for i in range(len(self.cell.ineqsites)):
        ##     spcstring = ""
        ##     if len(occ[i]) > 1:
        ##         # don't know what to put for alloy
        ##         filestring += "??"
        ##         for k in occ[i]:
        ##             spcstring += str(ed.elementnr[k])+"/"
        ##         spcstring = spcstring.rstrip("/")+"  "
        ##         for k in occ[i]:
        ##             spcstring += k+"/"
        ##         spcstring = spcstring.rstrip("/")
        ##     else:
        ##         for k in occ[i]:
        ##             filestring += str(ed.elementnr[k]).rjust(2)
        ##             spcstring = k
        ##     filestring += str(self.cell.ineqsites[i])+"      ! "+spcstring+"\n"
        filestring += "END\n"
        return filestring

################################################################################################
class SpacegroupFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a spacegroup.in file and the method
    __str__ that outputs the contents of an spacegroup.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.HermannMauguin = ""
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.supercelldims = [1, 1, 1]
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = ""
        if (self.HermannMauguin[-1] == 'R' or self.HermannMauguin[-1] == 'H') and self.HermannMauguin[-2] != ':':
            self.HermannMauguin = self.HermannMauguin[0:len(self.HermannMauguin)-1]+':'+self.HermannMauguin[-1]
        tmpstring=" '"+self.HermannMauguin+"'"
        tmpstring = tmpstring.ljust(50)+": hrmg\n"
        filestring += tmpstring
        tmpstring = ""
        tmpstring += " %15.11f" % (self.a)
        tmpstring += " %15.11f" % (self.b)
        tmpstring += " %15.11f" % (self.c)
        tmpstring = tmpstring.ljust(50)+": a, b, c\n"
        filestring += tmpstring
        tmpstring = " %15.9f %15.9f %15.9f"%(self.gamma,self.beta,self.alpha)
        tmpstring = tmpstring.ljust(50)+": ab, ac, bc\n"
        filestring += tmpstring
        tmpstring = ""
        for i in self.supercelldims:
            tmpstring += str(i)+"  "
        tmpstring = tmpstring.ljust(50)
        tmpstring += ": ncell\n"
        filestring += tmpstring
        filestring += ".true.".ljust(50)+": primcell\n"
        # Get species info
        species = set([])
        for occ in self.cell.occupations:
            spcstring = ""
            for k in occ:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            species.add(spcstring)
        tmpstring = str(len(species)).ljust(50)+": nspecies\n"
        filestring += tmpstring
        for spcs in species:
            # find number of representative sites for this species
            spcsites = 0
            positionstring = ""
            i = 0
            for occ in self.cell.occupations:
                spcstring = ""
                for k in occ:
                    spcstring += k+"/"
                spcstring = spcstring.rstrip("/")
                if spcstring == spcs:
                    spcsites += 1
                    positionstring += str(self.cell.ineqsites[i])+"\n"
                i += 1
            # output species info
            if len(spcs) > 2:
                # alloy
                spcsheader = "'??'".ljust(50)+"! "+spcs+"\n"+str(spcsites)+"\n"
            else:
                spcsheader = "'"+spcs+"'\n"+str(spcsites)+"\n"
            filestring += spcsheader
            filestring += positionstring
        filestring += "\n"+self.docstring
        return filestring

################################################################################################
class ElkFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an elk.in file and the method
    __str__ that outputs the contents of an elk.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # Make sure the docstring has the form of a f90 comment
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("!")
            string = "!"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = self.docstring
        # Lattice vectors
        filestring += "avec\n"
        tmpstring = ""
        for pos in self.cell.latticevectors:
            for coord in pos:
                tmpstring += "  %13.10f"%coord
            tmpstring += "\n"
        tmpstring += "\n"
        filestring += tmpstring
        # Scale factor
        filestring += "scale\n"
        filestring += "  %13.10f\n\n"%self.cell.lengthscale
        # Atoms
        filestring += "atoms\n"
        # Get number of species
        species = set([])
        for a in self.cell.atomdata:
            for b in a:
                species.add(b.spcstring())
        tmpstring = ("  "+str(len(species))).ljust(37)+": nspecies\n"
        filestring += tmpstring
        # local B-field string
        bfcmtstring = "   0.00000000  0.00000000  0.00000000"
        # initialize some stuff
        natoms = 0
        spcstring = self.cell.atomdata[0][0].spcstring()
        positionstring = ""
        for a in self.cell.atomdata:
            for b in a:
                spcs = b.spcstring()
                # Accumulate for this species
                if spcs == spcstring:
                    ## natoms += len(a)
                    natoms += 1
                    positionstring += str(b.position)+bfcmtstring+"\n"
                else:
                    # Print species
                    if len(spcstring) > 2:
                        # alloy
                        filestring += "'??.in'".ljust(37)+": spfname = "+spcstring+"\n"
                    else:
                        filestring += ("'"+spcstring+".in'").ljust(37)+": spfname \n"
                    filestring += "  "+str(natoms)+"\n"
                    filestring += positionstring
                    # Initialize next species
                    spcstring = spcs
                    natoms = 1
                    positionstring = str(b.position)+bfcmtstring+"\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += "'??.in'".ljust(37)+": spfname = "+spcstring+"\n"
        else:
            filestring += ("'"+spcstring+".in'").ljust(37)+": spfname\n"
        filestring += "  "+str(natoms)+"\n"
        filestring += positionstring
        return filestring

################################################################################################
class ExcitingFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an input.xml file and the method
    __str__ that outputs the contents of an input.xml file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.title = ""
        self.docstring = self.docstring.rstrip("\n")+"\n"
    def __str__(self):
        filestring = "<input>\n"
        filestring += "  <title>\n"
        filestring += self.docstring
        filestring += "  </title>\n"
        # Add title if there is one
        if self.title != "":
            filestring += "  <title>"+self.title+"</title>\n"
        filestring += "  <structure>\n"
        # scale factor
        filestring += "    <crystal scale="+str(self.cell.lengthscale)+">\n"
        # Lattice vectors
        tmpstring = ""
        for pos in self.cell.latticevectors:
            tmpstring += "      <basevect>"
            for coord in pos:
                tmpstring += " %13.10f"%coord
            tmpstring += "</basevect>\n"
        filestring += tmpstring
        filestring += "    </crystal>\n"
        # Atoms
        # local B-field string
        bfcmtstring = "   0.00000000  0.00000000  0.00000000"
        # initialize some stuff
        spcstring = self.cell.atomdata[0][0].spcstring()
        positionstring = ""
        for a in self.cell.atomdata:
            for b in a:
                spcs = b.spcstring()
                # Accumulate for this species
                if spcs == spcstring:
                    positionstring += "      <atom coord=\""
                    positionstring += str(b.position)+"\"/>\n"
                else:
                    # Print species
                    if len(spcstring) > 2:
                        # alloy
                        filestring += "    <species speciesfile=\"??.xml\">"
                        filestring += "       <!-- "+spcstring+" -->\n"
                    else:
                        filestring += "    <species speciesfile=\""+spcstring+".xml\">\n"
                    filestring += positionstring+"    </species>\n"
                    # Initialize next species
                    spcstring = spcs
                    positionstring = "      <atom coord=\""+str(b.position)+"\"/>\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += "    <species speciesfile=\"??.xml\">"
            filestring += "       <!-- "+spcstring+" -->\n"
        else:
            filestring += "    <species speciesfile=\""+spcstring+".xml\">\n"
        filestring += positionstring
        filestring += "    </species>\n"
        filestring += "  </structure>\n"
        filestring += "</input>\n"
        return filestring

################################################################################################
class FleurFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a Fleur input generator input file (how about
    that, we generate input for the generator of the input...) and the method
    __str__ that outputs the contents as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # make sure the docstring goes on one line
        self.docstring = self.docstring.replace("\n"," ")
        if len(self.docstring) > 80:
            self.docstring = self.docstring[0:78]+"...\n"
    def __str__(self):
        ed = ElementData()
        filestring = self.docstring+"\n"
        filestring += "&input cartesian=f oldfleur=f\n\n"
        # Lattice vectors
        tmpstring = ""
        n = 1
        for pos in self.cell.latticevectors:
            tmpstring += str(pos)
            tmpstring += "    !  a%1i\n"%n
            n += 1
        filestring += tmpstring
        # Scale factor
        filestring += "%13.9f    ! aa\n"%self.cell.lengthscale
        filestring += "1.0  1.0  1.0 ! scale(1), scale(2), scale(3)\n"
        # Atoms
        natom = 0
        for a in self.cell.atomdata:
            natom += len(a)
        filestring += str(natom)+"\n"
        nspcs = 0
        spcs = ""
        coordstring = ""
        for a in self.cell.atomdata:
            for b in a:
                # Check for alloy
                if b.alloy():
                    prestring = "??"
                    poststring = "  ! "
                    for k in b.species:
                        poststring += str(ed.elementnr[k])+"/"
                    poststring = poststring.rstrip("/")+" "+str(b.spcstring())+"\n"
                else:
                    prestring = str(ed.elementnr[b.spcstring()]).ljust(2)
                    poststring = "  ! "+b.spcstring()+"\n"
                coordstring += prestring+str(b.position)+poststring
        # To filestring
        filestring += coordstring
        return filestring

################################################################################################
class CASTEPFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a CASTEP run and the method
    __str__ that outputs to a .cell file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Cartesian units?
        self.cartesian = False
        # What units?
        self.unit = "angstrom"
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        # VCA calculation?
        self.vca = False
    def __str__(self):
        # Set units
        self.cell.newunit(self.unit)
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring+"\n"
        filestring += "%BLOCK LATTICE_CART\n"
        # units
        if self.cell.unit == "angstrom":
            filestring += "ang    # angstrom units\n"
        elif self.cell.unit == "bohr":
            filestring += "bohr   # atomic units\n"
        # lattice
        for vec in lattice:
            for coord in vec:
                filestring += " %19.15f"%(coord*a)
            filestring += "\n"
        # Cutoff
        filestring += "%ENDBLOCK LATTICE_CART\n\n"
        # The atom position info
        if self.cartesian:
            # Correct block name and units
            filestring += "%BLOCK POSITIONS_ABS\n"
            if self.cell.unit == "angstrom":
                filestring += "ang    # angstrom units\n"
            elif self.cell.unit == "bohr":
                filestring += "bohr   # atomic units\n"
            # Set transformation matrix
            transmat = LatticeMatrix(self.cell.latticevectors)
            scalfac = self.cell.a
        else:
            filestring += "%BLOCK POSITIONS_FRAC\n"
            transmat = LatticeMatrix([[1,0,0],[0,1,0],[0,0,1]])
            scalfac = 1.0
        i = 0
        for a in self.cell.atomdata:
            for b in a:
                pos = Vector(mvmult3(transmat,b.position)).scalmult(scalfac)
                # Check for VCA calculation
                if self.cell.alloy and self.vca:
                    if len(b.species) > 1:
                        i = i + 1
                        for sp,conc in b.species.iteritems():
                            filestring += sp.ljust(2)+" "+str(pos)+"  MIXTURE:( %i %6.5f )\n"%(i,conc)
                    else:
                        filestring += b.spcstring().ljust(2)+" "+str(pos)+"\n"
                else:
                    filestring += b.spcstring().ljust(2)+" "+str(pos)+"\n"
        if self.cartesian:
            filestring += "%ENDBLOCK POSITIONS_ABS\n"
        else:
            filestring += "%ENDBLOCK POSITIONS_FRAC\n"
        # pseudo-potential block
        species = set([])
        for a in self.cell.atomdata:
            if self.vca:
                for sp,conc in a[0].species.iteritems():
                    species.add(sp)
            else:
                species.add(a[0].spcstring())
        filestring += "\n"
        filestring += "# Commented out pseudopotential block for easy editing\n"
        filestring += "#%BLOCK SPECIES_POT\n"
        for sp in species:
            filestring += "# "+sp.ljust(2)+"  "+sp+"_00.usp\n"
        filestring += "#%ENDBLOCK SPECIES_POT\n"
        # Put in the symmetry operations
        filestring += "\n%BLOCK SYMMETRY_OPS\n"
        latvect = self.cell.conventional_latticevectors()
        # make list and make sure that identity comes first
        symoplist = sorted(list(self.cell.symops))
        k = 1
        for op in symoplist:
            filestring += "# Symm. op. %i\n"%k
            filestring += str(op)
            k += 1
        filestring += "%ENDBLOCK SYMMETRY_OPS\n"        
        return filestring

################################################################################################
# PWSCF (Quantum Espresso)
class PWSCFFile(GeometryOutputFile):
    """
    Class for storing the geometrical data for a PWSCF run and the method
    __str__ that outputs to a .in file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.docstring += "\n"
    def __str__(self):
        filestring = self.docstring
        a = self.cell.latticevectors[0].length()*self.cell.lengthscale
        boa = self.cell.latticevectors[1].length()/self.cell.latticevectors[0].length()
        coa = self.cell.latticevectors[2].length()/self.cell.latticevectors[0].length()
        filestring += "&SYSTEM\n"
        filestring += "  ibrav = %i\n"%self.ibrav()
        filestring += "  celldm(1) = %10.5f, celldm(2) = %10.5f, celldm(3) = %10.5f\n"%(a,boa,coa)
        filestring += "  nat = %i, ntyp = %i\n"
        filestring += "/\n"
        filestring += "ATOMIC_SPECIES"
        
        return filestring
    # Return the PWscf internal bravais lattice number
    def ibrav(self):
        system = self.cell.crystal_system()
        setting = self.cell.spacegroupsetting
        if self.supercell:
            return 14
        if system == 'cubic':
            if self.primcell:
                if setting == 'P':
                    return 1
                elif setting == 'F':
                    return 2
                elif setting == 'I':
                    return 3
            else:
                return 1
        if system == 'hexagonal':
            if self.primcell:
                if setting == 'P':
                    return 4
                elif setting == 'R':
                    return 5
        
################################################################################################
# CP2K
class CP2KFile(GeometryOutputFile):
    """
    Class for storing the geometrical data for a CP2k run and the method
    __str__ that outputs to a .inp file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
        self.docstring += "\n"
    def __str__(self):
        filestring = self.docstring
        filestring += "&CELL\n"
        filestring += "  UNIT ANGSTROM\n"
        filestring += "  PERIODIC XYZ\n"
        filestring += "  A "+str(self.cell.latticevectors[0].scalmult(self.cell.lengthscale))+"\n"
        filestring += "  B "+str(self.cell.latticevectors[1].scalmult(self.cell.lengthscale))+"\n"
        filestring += "  C "+str(self.cell.latticevectors[2].scalmult(self.cell.lengthscale))+"\n"
        filestring += "&END CELL\n\n"
        filestring += "&COORD\n"
        for a in self.cell.atomdata:
            for b in a:
                filestring += b.spcstring()+str(Vector(mvmult3(self.cell.latticevectors,b.position)).scalmult(self.cell.lengthscale))+"\n"
        filestring += "&END COORD\n"
        return filestring 

################################################################################################
class CPMDFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a CPMD run and the method
    __str__ that outputs to a .inp file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("bohr")
        self.cutoff = 100.0
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # Transformation to cartesian coordinates
        transmtx = []
        for i in range(3):
            transmtx.append([])
            for j in range(3):
                transmtx[i].append(lattice[i][j] * a)
        # docstring
        filestring = self.docstring+"\n"
        filestring += "&SYSTEM\n"
        # lattice
        filestring += " CELL VECTORS\n"
        for vec in transmtx:
            for coord in vec:
                filestring += " %19.15f"%coord
            filestring += "\n"
        # Cutoff
        filestring += " CUTOFF\n"
        filestring += " "+str(self.cutoff)+"\n"
        filestring += "&END\n\n"
        # The atom position info
        filestring += "&ATOMS\n"
        # get all species
        species = set([])
        for a in self.cell.atomdata:
            for b in a:
                species.add(b.spcstring())
        for spc in species:
            filestring += "*[pseudopotential file for "+spc+" here]\n"
            # Find maximal angular momentum
            spcs = spc.split("/")
            l = "s"
            for s in spcs:
                 if ed.angularmomentum[ed.elementblock[s]] > ed.angularmomentum[l]:
                     l = ed.elementblock[s]
            natoms = 0
            posstring = ""
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == spc:
                        natoms +=1
                        posstring += str(Vector(mvmult3(transmtx,b.position)))+"\n"
            # Print
            filestring += str(natoms)+"\n"
            filestring += posstring
        filestring += "&END\n"
        return filestring

################################################################################################
class SiestaFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a Siesta run and the method
    __str__ that outputs to a .fdf file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Assign some local variables
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring
        filestring += "AtomicCoordinatesFormat".ljust(28)+"Fractional\n"
        species = set([])
        natom = 0
        for a in self.cell.atomdata:
            natom += len(a)
            for b in a:
                species.add(b.spcstring())
        species = list(species)
        nspcs = len(species)
        filestring += "LatticeConstant".ljust(28)+str(self.cell.lengthscale)+" Ang\n"
        filestring += "NumberOfAtoms".ljust(28)+str(natom)+"\n"
        filestring += "NumberOfSpecies".ljust(28)+str(nspcs)+"\n"
        # lattice
        filestring += "%block LatticeVectors\n"
        for vec in lattice:
            filestring += str(vec)+"\n"
        filestring += "%endblock LatticeVectors\n"
        # Atomic coordinates
        filestring += "%block AtomicCoordinatesAndAtomicSpecies\n"
        i = 1
        for sp in species:
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        filestring += str(b.position)
                        filestring += "   %i\n"%i
            i += 1
        filestring += "%endblock AtomicCoordinatesAndAtomicSpecies\n"
        # Chemical species
        filestring += "%block ChemicalSpeciesLabel\n"
        i = 1
        for sp in species:
            filestring += str(i).ljust(8)
            if len(sp) > 2:
                filestring += "??      ??    # "
                tsp = sp.split("/")
                for t in tsp:
                    filestring += str(ed.elementnr[t])+"/"
                filestring = filestring.rstrip("/")
                filestring += "      "+sp+"\n"
            else:
                filestring += str(ed.elementnr[sp]).ljust(8)+sp.ljust(8)+"\n"
            i += 1
        filestring += "%endblock ChemicalSpeciesLabel\n"
        return filestring

################################################################################################
class ABINITFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in an abinit run and the method
    __str__ that outputs the contents of a abinit input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("bohr")
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring
        # length scale and lattice
        filestring += "# Structural parameters\n"
        filestring += "acell   3*"+str(a)+"\n\n"
        filestring += "rprim   "
        for vec in lattice:
            filestring += str(vec)+"\n        "
        filestring += "\n"
        # The atom position info
        alloy = False
        spcs = ""
        typatstring = "typat    "
        natom = 0
        ntypat = 0
        znuclstring = "znucl    "
        alloystring = ""
        xredstring = "xred   "
        for a in self.cell.atomdata:
            for b in a:
                natom += 1
                if spcs != b.spcstring():
                    ntypat += 1
                    if b.alloy():
                        znuclstring += "?? "
                        alloystring += b.spcstring()+" "
                        alloy = True
                    else:
                        znuclstring += str(ed.elementnr[b.spcstring()])+" "
                typatstring += str(ntypat)+" "
                xredstring += str(b.position)+"\n       "
                spcs = b.spcstring()
        filestring += "natom    "+str(natom)+"\n"
        filestring += "ntypat   "+str(ntypat)+"\n"
        filestring += typatstring+"\n"
        filestring += znuclstring
        if alloy:
            filestring += "    # "+alloystring
        filestring += "\n"
        filestring += xredstring+"\n"
        return filestring

################################################################################################
class AIMSFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a FHI-AIMS run and the method
    __str__ that outputs the contents of a abinit input file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        self.cartesian = False
        # Make sure the docstring has comment form
        self.docstring = self.docstring.rstrip("\n")
        tmpstrings = self.docstring.split("\n")
        self.docstring = ""
        for string in tmpstrings:
            string = string.lstrip("#")
            string = "#"+string+"\n"
            self.docstring += string
    def __str__(self):
        filestring = self.docstring+"\n"
        latvecs = self.cell.latticevectors
        for i in range(3):
            for j in range(3):
                latvecs[i][j] = latvecs[i][j] * self.cell.lengthscale
        for vec in latvecs:
            filestring += "lattice_vector  "+str(vec)+"\n"
        for a in self.cell.atomdata:
            for b in a:
                if self.cartesian:
                    filestring += "atom  "+str(Vector(mvmult3(latvecs,b.position)))+" "+b.spcstring()+"\n"
                else:
                    filestring += "atom_frac  "+str(b.position)+" "+b.spcstring()+"\n"
        return filestring
        
################################################################################################
class POSCARFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a POSCAR file and the method
    __str__ that outputs the contents of a POSCAR file as a string.
    If you want POSCAR to be printed with the atomic positions in Cartesian form,
    then set
    POSCARFile.printcartpos = True
    If you want to put the overall length scale on the lattice vectors and print 1.0
    for the length scale, then set
    POSCARFile.printcartvecs = True
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
        self.printcartvecs = False
        self.printcartpos = False
        self.vasp5format = False
        self.selectivedyn = False
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
        # make sure the docstring goes on one line
        self.docstring = self.docstring.replace("\n"," ")
    def SpeciesOrder(self):
        """
        Return a string with the species in the order they appear in POSCAR.
        """
        returnstring = ""
        for sp in self.species:
            returnstring += sp+" "
        return returnstring
    def __str__(self):
        # Assign some local variables
        lattice = self.cell.latticevectors
        # VASP needs lattice vector matrix to have positive triple product
        if det3(lattice) < 0:
            if lattice[0].length() == lattice[1].length() == lattice[2].length():
                # Shift the first and last for cubic lattices
                transmtx = [[0, 0, 1],
                            [0, 1, 0],
                            [1, 0, 0]]
            else:
                # Else shift the two shortest
                if lattice[0].length() > lattice[1].length() and lattice[0].length() > lattice[2].length():
                    transmtx = [[1, 0, 0],
                                [0, 0, 1],
                                [0, 1, 0]]
                elif lattice[1].length() > lattice[2].length() and lattice[1].length() > lattice[0].length():
                    transmtx = [[0, 0, 1],
                                [0, 1, 0],
                                [1, 0, 0]]
                else:
                    transmtx = [[0, 1, 0],
                                [1, 0, 0],
                                [0, 0, 1]]
        else:
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        lattice = mmmult3(transmtx,lattice)
                
        # For output of atomic positions
        a = self.cell.lengthscale
        if self.selectivedyn:
            positionunits += "Selective dynamics\n"
        if self.printcartpos:
            positionunits = "Cartesian\n"
            coordmat = []
            for i in range(3):
                coordmat.append([])
                for j in range(3):
                    coordmat[i].append(lattice[i][j] * a)
        else:
            coordmat = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
            positionunits = "Direct\n"
        # The first line with info from input docstring
        filestring = self.docstring
        if not self.vasp5format:
            filestring += " Species order: "
            for sp in self.species:
                filestring += sp+" "
        filestring += "\n"
        # Lattice parameter and vectors
        if self.printcartvecs:
            latticestring = " 1.0\n"
            for i in range(3):
                latticestring += "%19.15f %19.15f %19.15f\n" % (a*lattice[i][0], a*lattice[i][1], a*lattice[i][2])
        else:
            latticestring = " %10f\n" % a
            for i in range(3):
                latticestring += "%19.15f %19.15f %19.15f\n" % (lattice[i][0], lattice[i][1], lattice[i][2])
        filestring += latticestring
        # print species here if vasp 5 format
        if self.vasp5format:
            for sp in self.species:
                filestring += (" "+sp).rjust(4)
            filestring += "\n"
        # positions and number of species
        nspstring = ""
        positionstring = ""
        for sp in self.species:
            nsp = 0
            for a in self.cell.atomdata:
                for b in a:
                    if b.spcstring() == sp:
                        nsp += 1
                        p = Vector(mvmult3(coordmat,mvmult3(transmtx,b.position)))
                        positionstring += str(p)
                        if self.selectivedyn:
                            positionstring += "   T  T  T"
                        positionstring += "\n"
            nspstring += (" "+str(nsp)).rjust(4)
        filestring += nspstring+"\n"
        filestring += positionunits
        filestring += positionstring
        return filestring

class POTCARFile:
    """
    Class for representing and outputting a POTCAR file for VASP.
    """
    def __init__(self, crystalstructure, directory=""):
        self.cell = crystalstructure
        if directory != "":
            self.dir = directory
        else:
            try:
                self.dir = os.environ['VASP_PSEUDOLIB']
            except:
                try:
                    self.dir = os.environ['VASP_PAWLIB']
                except:
                    self.dir = ""
        # check directory
        if self.dir == "":
            raise SetupError("No path to the VASP pseudopotential library specified.\n")
        if not os.path.exists(self.dir):
            raise SetupError("The specified path to the VASP pseudopotential library does not exist.\n"+self.dir)
        # set up species list
        tmp = set([])
        for a in self.cell.atomdata:
            for b in a:
                tmp.add(b.spcstring())
        self.species = list(tmp)
    def __str__(self):
        # Make good selection of potcars
        prioritylist = ["_d", "_pv", "_sv", "", "_h", "_s"]
        # get all files
        potcarlist = []
        for a in self.species:
            for version in prioritylist:
                potcarfile = self.dir+"/"+a+version+"/POTCAR"
                if os.path.exists(potcarfile):
                    potcarlist.append(potcarfile)
                    break
        # read potcar files and put in outstring
        outstring = ""
        for f in potcarlist:
            potcar = open(f,"r")
            outstring += potcar.read()
            potcar.close()
        return outstring

class KPOINTSFile:
    """
    Class for representing and outputting a KPOINTS file for VASP.
    """
    def __init__(self, crystalstructure, docstring="",kresolution=0.2):
        self.docstring = docstring
        self.kresolution = kresolution
        # set reciprocal lattice vectors in reciprocal angstroms
        reclatvect = crystalstructure.reciprocal_latticevectors()
        for j in range(3):
            for i in range(3):
                reclatvect[j][i] = reclatvect[j][i] / crystalstructure.a
        # Lengths of reciprocal lattice vectors
        reclatvectlen = [elem.length() for elem in reclatvect]
        self.kgrid = [max(1,int(round(elem/self.kresolution))) for elem in reclatvectlen]
    def __str__(self):
        tmp = self.docstring
        tmp += " k-space resolution ~"+str(self.kresolution)+"/A\n"
        tmp += " 0\n"
        tmp += "Gamma\n"
        tmp += str(self.kgrid[0])+" "+str(self.kgrid[1])+" "+str(self.kgrid[2])+"\n"
        tmp += "0 0 0\n"
        return tmp
        
class INCARFile:
    """
    Class for representing and outputting a INCAR file for VASP.
    """
    def __init__(self, crystalstructure, docstring="",potcardir=""):
        self.cell = crystalstructure
        self.docstring = "# "+docstring.lstrip("#").rstrip("\n")+"\n"
        # we need the potcar directory
        if potcardir != "":
            self.potcardir = potcardir
        else:
            try:
                self.potcardir = os.environ['VASP_PSEUDOLIB']
            except:
                try:
                    self.potcardir = os.environ['VASP_PAWLIB']
                except:
                    self.potcardir = ""
        # check directory
        if self.potcardir == "":
            raise SetupError("No path to the VASP pseudopotential library specified.\n")
        if not os.path.exists(self.potcardir):
            raise SetupError("The specified path to the VASP pseudopotential library does not exist.\n"+self.dir)
        # set up species list
        self.species = dict([])
        for a in self.cell.atomdata:
            for b in a:
                spcstr = b.spcstring()
                if spcstr in self.species:
                    t = self.species[spcstr] + 1
                    self.species[spcstr] = t
                else:
                    self.species[spcstr] = 1
        # get potcar list
        prioritylist = ["_d", "_pv", "_sv", "", "_h", "_s"]
        potcars = dict([])
        specieslist = []
        for a in self.species:
            for version in prioritylist:
                potcarfile = self.potcardir+"/"+a+version+"/POTCAR"
                if os.path.exists(potcarfile):
                    potcars[a] = potcarfile
                    specieslist.append(a)
                    break        
        # get maximal encut and number of electrons from potcars
        enmaxs = dict([])
        zvals = dict([])
        for a,f in potcars.iteritems():
            potcar = open(f,"r")
            for line in potcar:
                if search("ZVAL",line):
                    zvals[a] = float(line.split("ZVAL")[1].lstrip("= ").split()[0].strip(string.punctuation))
                if search("ENMAX",line):
                    enmaxs[a] = float(line.split("ENMAX")[1].lstrip("= ").split()[0].strip(string.punctuation))
                if search("END of PSCTR",line):
                    break
            potcar.close()
        self.maxencut = max([k for v,k in enmaxs.iteritems()])
        # ecut = max(encuts found in potcars)*encutfac
        self.encutfac = 1.5
        # do we suspect that this might be magnetic?
        self.magnetic = False
        self.magmomlist = []
        suspiciouslist = set(["Cr", "Mn", "Fe", "Co", "Ni",
                              "Ce","Pr","Nd","Pm","Sm","Eu",
                              "Gd","Tb","Dy","Ho","Er","Tm"])
        initialmoments = {"Cr" : 3, "Mn" : 3, "Fe" : 3, "Co" : 3, "Ni" : 1,
                          "Ce" : 1, "Pr" : 2, "Nd" : 3, "Pm" : 4, "Sm" : 5, "Eu" : 6,
                          "Gd" : 7, "Tb" : 8, "Dy" : 9, "Ho" : 10, "Er" : 11, "Tm" : 12}
        for s in self.species:
            if s in suspiciouslist:
                self.magnetic = True
                for i in range(self.species[s]):
                    self.magmomlist.append(str(initialmoments[s]))
            else:
                for i in range(self.species[s]):
                    self.magmomlist.append("0")
        # Determine NBANDS
        nmag = sum([int(i) for i in self.magmomlist])
        nelect = 0.0
        for sp,z in zvals.iteritems():
            for a in self.cell.atomdata:
                for b in a:
                    if sp == b.spcstring():
                        nelect += z
        if nmag > 0:
            nstates = int(math.ceil(nelect))
        else:
            nstates = int(math.ceil(nelect/2))
        # NBANDS is max of the default VASP definition and occupied bands+20
        natoms = sum([len(a) for a in self.cell.atomdata])
        self.nbands = max(max(max(int(math.ceil(nelect/2))+int(natoms/2),3), math.ceil(0.6*nelect))+nmag, nstates+20)
            
    def __str__(self):
        tmp = self.docstring
        tmp += "ENCUT = "+str(self.maxencut*self.encutfac)+"\n"
        tmp += "NBANDS = "+str(self.nbands)+"\n"
        ## tmp += "IBRION = 1\n"
        ## tmp += "POTIM = 0.4\n"
        ## tmp += "ISIF = 2\n"
        ## tmp += "NSW = 30\n"
        ## tmp += "NELMIN = 4\n"
        tmp += "PREC = Accurate\n"
        tmp += "LREAL = Auto\n"
        tmp += "ISMEAR = 0\n"
        tmp += "SIGMA = 0.1\n"
        if self.magnetic:
            tmp += "ISPIN = 2\n"
            tmp += "MAGMOM = "
            for species in self.magmomlist:
                tmp += species+" "
            tmp += "\n"
        return tmp
        

################################################################################################
# EMTO 
class KFCDFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the kfcd program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.kstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        tmpstring = "KFCD      MSGL..=  0"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam+"\n"
        filestring += tmpstring
        tmpstring = "STRNAM...="+self.kstrjobnam+"\n"
        filestring += tmpstring
        filestring += "DIR001=../kstr/smx/\n"
        filestring += "DIR002=../kgrn/chd/\n"
        filestring += "DIR003=../shape/shp/\n"
        filestring += "DIR004=../bmdl/mdl/\n"
        filestring += "DIR006=./\n"
        filestring += "Lmaxs.= 30 NTH..= 41 NFI..= 81 FPOT..= N\n"
        filestring += "OVCOR.=  Y UBG..=  N NPRN.=  0\n"
        return filestring

class KGRNFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the kgrn program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.kstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.latticenr = 14
    def __str__(self):
        ed = ElementData()
        filestring = ""
        tmpstring = "KGRN"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM="+self.jobnam+"\n"
        filestring += tmpstring
        filestring += "STRT..=  A MSGL.=  0 EXPAN.= S FCD..=  Y FUNC..= SCA\n"
        tmpstring = "FOR001=../kstr/smx/"+self.kstrjobnam+".tfh\n"
        tmpstring += "FOR001=../kstr/smx/"+self.kstrjobnam+"10.tfh\n"
        filestring += tmpstring
        filestring += "DIR002=pot/\n"
        filestring += "DIR003=pot/\n"
        tmpstring = "FOR004=../bmdl/mdl/"+self.kstrjobnam+".mdl\n"
        filestring += tmpstring
        filestring += "DIR006=\n"
        filestring += "DIR009=pot/\n"
        filestring += "DIR010=chd/\n"
        # Use environment variable TMPDIR if possible
        tmpstring = "DIR011="
        if "TMPDIR" in os.environ:
            tmpstring += os.environ["TMPDIR"]
            # Make sure the string will end with a single /
            tmpstring = tmpstring.rstrip("/")
        else:
            # ...else check for /tmp
            if os.path.isdir("/tmp"):
                tmpstring += "/tmp"
            else:
                # ...and last resort is ./
                tmpstring += "."
        tmpstring += "/\n"
        filestring += tmpstring
        filestring += self.docstring.replace("\n"," ")+"\n"
        filestring += "Band: 10 lines\n"
        tmpstring = "NITER.= 50 NLIN.= 31 NPRN.=  0 NCPA.= 20 NT...=%3i"%len(self.cell.atomdata)+" MNTA.="
        # Work out maximal number of species occupying a site
        mnta = 1
        for a in self.cell.atomdata:
            for b in a:
                mnta = max(mnta,len(b.species))
        tmpstring += "%3i"%mnta+"\n"
        filestring += tmpstring
        filestring += "MODE..= 3D FRC..=  N DOS..=  N OPS..=  N AFM..=  P CRT..=  M\n"
        filestring += "Lmaxh.=  8 Lmaxt=  4 NFI..= 31 FIXG.=  2 SHF..=  0 SOFC.=  N\n"
        # Choose brillouin zone by lattice type
        # Output the smallest allowed n for each direction in this lattice type
        if self.latticenr == 1:
            nkx = 0
            nky = 2
            nkz = 0
        elif self.latticenr == 2:
            nkx = 0
            nky = 5
            nkz = 0
        elif self.latticenr == 3:
            nkx = 0
            nky = 3
            nkz = 0
        elif self.latticenr == 4:
            nkx = 0
            nky = 3
            nkz = 2
        elif self.latticenr == 5:
            nkx = 0
            nky = 2
            nkz = 2
        elif self.latticenr == 6:
            nkx = 0
            nky = 3
            nkz = 2
        elif self.latticenr == 7:
            nkx = 0
            nky = 3
            nkz = 3
        elif self.latticenr == 8:
            nkx = 2
            nky = 2
            nkz = 1
        elif self.latticenr == 9:
            nkx = 2
            nky = 2
            nkz = 1
        elif self.latticenr == 10:
            nkx = 2
            nky = 2
            nkz = 1
        elif self.latticenr == 11:
            nkx = 1
            nky = 1
            nkz = 1
        elif self.latticenr == 12:
            nkx = 1
            nky = 2
            nkz = 1
        elif self.latticenr == 13:
            nkx = 3
            nky = 3
            nkz = 0
        else:
            nkx = 2
            nky = 2
            nkz = 2
        filestring += "KMSH...= G IBZ..= %2i NKX..= %2i NKY..= %2i NKZ..= %2i FBZ..=  N\n"%(self.latticenr,nkx,nky,nkz)
        filestring += "KMSH2..= G IBZ2.=  1 NKX2.=  4 NKY2.=  0 NKZ2.= 51\n"
        filestring += "ZMSH...= C NZ1..= 16 NZ2..= 16 NZ3..=  8 NRES.=  4 NZD.= 500\n"
        filestring += "DEPTH..=  1.500 IMAGZ.=  0.020 EPS...=  0.200 ELIM..= -1.000\n"
        filestring += "AMIX...=  0.100 EFMIX.=  1.000 VMTZ..=  0.000 MMOM..=  0.000\n"
        filestring += "TOLE...= 1.d-07 TOLEF.= 1.d-07 TOLCPA= 1.d-06 TFERMI=  500.0 (K)\n"
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        # average wigner-seitz radius
        volume = abs(det3(self.cell.latticevectors))
        wsr = self.cell.lengthscale * 3*volume/(nosites * 4 * pi)**third
        filestring += "SWS......=%8f   NSWS.=  1 DSWS..=   0.05 ALPCPA= 0.9020\n"%wsr
        filestring += "Setup: 2 + NQ*NS lines\n"
        filestring += "EFGS...=  0.000 HX....=  0.100 NX...= 11 NZ0..=  6 STMP..= Y\n"
        # atom info
        filestring += "Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT\n"
        iq = 1
        it = 1
        for a in self.cell.atomdata:
            ita = 1
            # THIS MAKES ASSUMPTIONS ABOUT THE ORDERING OF ATOMDATA
            # But we're OK for all orderings implemented so far
            for comp in a[0].spcstring().split("/"):
                for b in a:
                    if comp in b.species:
                        tmpstring = comp.ljust(4)+"  "+"%3i%3i%3i"%(iq,it,ita)
                        tmpstring += "%4i"%ed.elementnr[comp]
                        tmpstring += "%7.3f%7.3f%7.3f%7.3f"%(a[0].species[comp],1,1,1)
                        tmpstring += "%5.2f%5.2f\n"%(0,0)
                        filestring += tmpstring
                        iq += 1
                ita += 1
                iq -= len(a)
            iq += len(a)
            it += 1
        filestring += "Atom:  4 lines + NT*NTA*6 lines\n"
        filestring += "IEX...=  4 NP..= 251 NES..= 15 NITER=100 IWAT.=  0 NPRNA=  0\n"
        filestring += "VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n"
        filestring += "DX.......=  0.030000 DR1.....=  0.002000 TEST....=  1.00E-12\n"
        filestring += "TESTE....=  1.00E-12 TESTY...=  1.00E-12 TESTV...=  1.00E-12\n"
        for a in self.cell.atomdata:
            for comp in a[0].species:
                filestring += comp+"\n"
                try:
                    filestring += ed.emtoelements[comp]
                except KeyError:
                    filestring += "\n\n\n\n\n"
        return filestring

class ShapeFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the shape program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.jobnam = "default"
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        filestring = ""
        tmpstring = "SHAPE     HP......=N"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1\n"
        filestring += tmpstring
        filestring += "FOR001=../kstr/smx/"+self.jobnam+".tfh\n"
        filestring += "DIR002=shp/\n"
        filestring += "DIR006=./\n"
        filestring += "Lmax..= 30 NSR..=129 NFI..= 11\n"
        filestring += "NPRN..=  0 IVEF.=  3\n"
        return filestring

class BMDLFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the bmdl program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.latticenr = 14
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        lv = self.cell.latticevectors
        ed = ElementData()
        filestring = ""
        tmpstring = "BMDL      HP......=N"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 NPRN.=  0\n"
        filestring += tmpstring
        filestring += "DIR001=mdl/\n"
        filestring += "DIR006=./\n"
        filestring += self.docstring.replace("\n"," ")+"\n"
        filestring += "NL.....= 7\n"
        filestring += "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n"
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        if self.latticenr == 0:
            tmpstring = "NQ3...=%3i LAT...= 0 IPRIM.= 0 NGHBP.=13 NQR2..= 0\n" % (nosites,self.latticenr)
        else:
            tmpstring = "NQ3...=%3i LAT...=%2i IPRIM.= 1 NGHBP.=13 NQR2..= 0\n" % (nosites,self.latticenr)
        filestring += tmpstring
        boa = self.b/self.a
        coa = self.c/self.a
        filestring += "A........= 1.0000000 B.......=%10f C.......=%10f\n"%(boa,coa)
        tmpstring = ""
        if self.latticenr == 0:
            for i in range(3):
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (lv[i][0],lv[i][1],lv[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        for a in self.cell.atomdata:
            for b in a:
                v = mvmult3(lv,b.position)
                filestring += "QX(IQ)...=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                filestring += "      "+b.spcstring()+"\n"
        return filestring

class KSTRFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the kstr program
    and the method __str__ that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.latticenr = 14
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.hardsphere = 0.67
        self.iprim = 0
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        lv = self.cell.latticevectors
        ed = ElementData()
        filestring = ""
        tmpstring = "KSTR      HP......=N"
        tmpstring = tmpstring.ljust(25)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 MODE...=B STORE..=Y HIGH...=Y\n"
        filestring += tmpstring
        filestring += "DIR001=smx/\n"
        filestring += "DIR006=./\n"
        filestring += self.docstring.replace("\n"," ")+"\n"
        # NL = maximal l from element blocks
        maxl = 1
        for a in self.cell.atomdata:
            for b in a:
                for i in b.species:
                    if ed.elementblock[i] == "p":
                        maxl = max(maxl,2)
                    elif ed.elementblock[i] == "d":
                        maxl = max(maxl,3)
                    elif ed.elementblock[i] == "f":
                        maxl = max(maxl,4)
        tmpstring = "NL.....= %1i NLH...=11 NLW...= 9 NDER..= 6 ITRANS= 3 NPRN..= 0\n" % maxl
        filestring += tmpstring
        # Get number of sites
        nosites = 0
        for a in self.cell.atomdata:
            nosites += len(a)
        # Setting the real space summation cutoff to 4.5*(wigner-seitz radius)
        volume = abs(det3(lv))
        wsr = (3*volume/(nosites * 4 * pi))**third
        tmpstring = "(K*W)^2..=  0.000000 DMAX....=%10f RWATS...=      0.10\n" % (wsr*4.5)
        filestring += tmpstring
        tmpstring = "NQ3...=%3i LAT...=%2i IPRIM.=%2i NGHBP.=13 NQR2..= 0\n" % (nosites,self.latticenr,self.iprim)
        filestring += tmpstring
        boa = self.b/self.a
        coa = self.c/self.a
        filestring += "A........= 1.0000000 B.......=%10f C.......=%10f\n"%(boa,coa)
        tmpstring = ""
        if self.iprim == 0:
            for i in range(3):
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (lv[i][0],lv[i][1],lv[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        for a in self.cell.atomdata:
            for b in a:
                v = mvmult3(lv,b.position)
                filestring += "QX(IQ)...=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                filestring += "      "+b.spcstring()+"\n"
        for i in range(nosites):
            filestring += "a/w(IQ)..="
            for i in range(4):
                filestring += "%5.2f"%self.hardsphere
            filestring += "\n"
        filestring += "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n"
        return filestring

################################################################################################
# SPRKKR
class XBandSysFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].sys file for the xband program
    and the method __str__ that outputs the contents of the .sys file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.jobnam = "default"
        self.a = 1
        self.b = 1
        self.c = 1
        self.alpha = 90
        self.beta = 90
        self.gamma = 90
        self.minangmom = None
        self.filename = ""
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
        ed = ElementData()
        # First identify any site which is not filled (concentrations add to 1.0)
        # and fill up with vacuum sphere (Vc)
        self.cell.fill_out_empty(label='Vc')
        # docstring on a single line
        filestring = deletenewline(self.docstring,replace=" ")
        filestring += "\n"+self.filename+"\n"
        filestring += "xband-version\n"
        filestring += "5.0\n"
        # It would be really cool to support lower dimensions...one day.
        filestring += "dimension\n"
        filestring += "3D\n"
        filestring += "Bravais lattice\n"
        if self.cell.crystal_system() == 'triclinic':
            filestring += "1  triclinic   primitive      -1     C_i \n"
        elif self.cell.crystal_system() == 'monoclinic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "2  monoclinic  primitive      2/m    C_2h\n"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "3  monoclinic  primitive      2/m    C_2h\n"
            else:
                print "xband only knows primitive and base-centered monoclinic settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == 'orthorhombic':
            if self.cell.spacegroupsetting == 'P':
                filestring += "4  orthorombic primitive      mmm    D_2h\n"
            elif self.cell.spacegroupsetting == 'A' or self.cell.spacegroupsetting == 'B' or \
                     self.cell.spacegroupsetting == 'C':
                filestring += "5  orthorombic body-centered  mmm    D_2h\n"
            elif self.cell.spacegroupsetting == "I":
                filestring += "6  orthorombic body-centered  mmm    D_2h\n"
            elif self.cell.spacegroupsetting == "F":
                filestring += "7  orthorombic face-centered  mmm    D_2h\n"
            else:
                print "xband does not know %1s centering of an orthorhombic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        elif self.cell.crystal_system() == "tetragonal":
            if self.cell.spacegroupsetting == "P":
                filestring += "8  tetragonal  primitive      4/mmm  D_4h\n"
            elif self.cell.spacegroupsetting == "I":
                filestring += "9  tetragonal  body-centered  4/mmm  D_4h\n"
            else:
                print "xband only knows primitive and body-centered tetragonal settings!"
                sys.exit(43)
        elif self.cell.crystal_system() == "trigonal":
            filestring += "10 trigonal    primitive      -3m    D_3d\n"
        elif self.cell.crystal_system() == "hexagonal":
            filestring += "11 hexagonal   primitive      6/mmm  D_6h\n"
        elif self.cell.crystal_system() == "cubic":
            if self.cell.spacegroupsetting == "P":
                filestring += "12 cubic       primitive      m3m    O_h \n"
            elif self.cell.spacegroupsetting == "F":
                filestring += "13 cubic       face-centered  m3m    O_h \n"
            elif self.cell.spacegroupsetting == "I":
                filestring += "14 cubic       body-centered  m3m    O_h \n"
            else:
                print "xband does not know %1s centering of a cubic cell."%self.cell.spacegroupsetting
                sys.exit(43)
        filestring += "space group number (ITXC and AP)\n"
        filestring += "%5i%5i"%(self.cell.spacegroupnr,Number2AP[self.cell.spacegroupnr])+"\n"
        filestring += "structure type\n"
        filestring += "UNKNOWN\n"
        filestring += "lattice parameter A  [a.u.]\n"
        filestring += "%18.12f\n"%self.cell.lengthscale
        filestring += "ratio of lattice parameters  b/a  c/a\n"
        filestring += "%18.12f%18.12f\n"%(self.cell.boa,self.cell.coa)
        filestring += "lattice parameters  a b c  [a.u.]\n"
        a = self.cell.lengthscale
        b = self.cell.b * self.cell.lengthscale / self.cell.a
        c = self.cell.c * self.cell.lengthscale / self.cell.a
        filestring += "%18.12f%18.12f%18.12f\n"%(a,b,c)
        filestring += "lattice angles  alpha beta gamma  [deg]\n"
        filestring += "%18.12f%18.12f%18.12f\n"%(self.cell.alpha,self.cell.beta,self.cell.gamma)
        filestring += "primitive vectors     (cart. coord.) [A]\n"
        for vec in self.cell.latticevectors:
            for p in vec:
                filestring += "%18.12f"%p
            filestring += "\n"
        # Get number of sites and fill out with empty spheres if the sites are not fully filled
        filestring += "number of sites NQ\n"
        nq = 0
        for a in self.cell.atomdata:
            nq += len(a)
        self.cell.fill_out_empty(label="Vc")
        filestring += "%3i\n"%nq
        filestring += " IQ ICL     basis vectors     (cart. coord.) [A]                      RWS [a.u.]  NLQ  NOQ ITOQ\n"
        # Average Wigner-Seitz radius
        rws = pow(3*self.cell.volume()/(4*pi*len(self.cell.atomset)),1.0/3.0)*self.cell.lengthscale
        iq = 0
        icl = 0
        itoq = 0
        for a in self.cell.atomdata:
            icl += 1
            itoqs = []
            for sp in a[0].species:
                itoq += 1
                itoqs.append(itoq)
            for b in a:
                iq += 1
                if self.minangmom:
                    angmom = max(max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1,self.minangmom)
                else:
                    angmom = max([ed.angularmomentum[ed.elementblock[spcs]] for spcs in b.species])+1
                v = mvmult3(self.cell.latticevectors,b.position)
                filestring += "%3i%4i%18.12f%18.12f%18.12f  %18.12f%4i%5i "%(iq,icl,v[0],v[1],v[2],rws,angmom,len(a[0].species))
                for i in itoqs:
                    filestring += "%3i"%i
                filestring += "\n"
        filestring += "number of sites classes NCL\n"
        filestring += "%3i\n"%len(self.cell.atomdata)
        filestring += "ICL WYCK NQCL IQECL (equivalent sites)\n"
        iq = 0
        icl = 0
        for a in self.cell.atomdata:
            icl += 1
            filestring += "%3i   %1s%5i"%(icl,'-',len(a))
            for b in a:
                iq += 1
                filestring += "%3i"%iq
            filestring += "\n"
        filestring += "number of atom types NT\n"
        nt = 0
        for a in self.cell.atomdata:
            nt += len(a[0].species)
        filestring += "%3i\n"%nt
        filestring += " IT  ZT  TXTT  NAT  CONC  IQAT (sites occupied)\n"
        iq = 0
        it = 0
        for a in self.cell.atomdata:
            corr = 0
            for sp,conc in a[0].species.iteritems():
                it += 1
                filestring += " %2i%4i  %8s%5i%6.3f"%(it,ed.elementnr[sp],sp,len(a),conc)
                iq -= corr*len(a)
                for b in a:
                    iq += 1
                    filestring += "%3i"%iq
                corr = 1
                filestring += "\n"
        return filestring
