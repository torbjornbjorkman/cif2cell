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
#  Description: A set of tools to generate the geometrical
#               setup for various electronic structure codes.
#               Contains some container classes for structural
#               data and methods to extract these from a CIF
#               file.
#               Currently supports standard (conventional)
#               cell settings and from that reduction to the
#               primitive cell.
#  Author: Torbjorn Bjorkman, torbjorn(at)cc.hut.fi
#  Affiliation: COMP, Aaalto University School of
#               Science and Technology, Department of
#               Applied Physics, Espoo, Finland
#******************************************************************************************
import copy
from elementdata import *
from uctools import *

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
        tmpstring = tmpstring.ljust(25)+"    "+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(10)+" MSGL.=  1 BSAVE..=N COLD...=Y DOS...=N SPO...=N ISM...=G RCLCR...=Y\n"
        filestring += tmpstring
        filestring += "FOR001=./"+self.bstrjobnam+".tfm\n"
        filestring += "FOR002=\n"
        filestring += "FOR003=\n"
        filestring += "FOR004=\n"
        filestring += "FOR006=\n"
        filestring += "FOR010=\n"
        filestring += "Band: 4 lines, "+self.docstring.replace("\n"," ")+"\n"
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
        if len(self.cell.atomdata[0][0].species) > 1:
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
                if len(b.species) > 1:
                    species = "??"
                else:
                    for v in b.species:
                        species = v
                if species != prevspecies:
                    prevspecies = species
                    nsp += 1
                tmpstring = species.ljust(2)+"  "+"%3i%3i"%(iq,it)
                try:
                    tmpstring += "%3i%3i"%(l[ed.elementblock[species]],1)
                except KeyError:
                    tmpstring += "  ?  1"
                tmpstring += "%3i"%nsp
                tmpstring += "    1.000 .000 0.00 0000 1111   .0   .0   .0   .0"
                if len(b.species) > 1:
                    # print alloy components at the end of the line
                    tmpstring += "       "
                    for comp in b.species:
                        tmpstring += comp+"/"
                    tmpstring = tmpstring.rstrip("/")
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
        tmpstring = tmpstring.ljust(40)+self.programdoc.replace("\n"," ")+"\n"
        filestring += tmpstring
        tmpstring = "JOBNAM...="+self.jobnam.ljust(9)+" MSGL.=   1 \n"
        filestring += tmpstring
        filestring += "MODE....=B STORE..=Y SCREEN.=B CMBC...=Y\n"
        filestring += "FOR001=\n"
        filestring += "FOR006=\n"
        filestring += self.docstring.replace("\n"," ")+"\n"
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
                    if ed.elementblock[k] == 's' or ed.elementblock[k] == 'p':
                        l = max(l,2)
                    elif ed.elementblock[k] == 'd':
                        l = max(l,3)
                    elif ed.elementblock[k] == 'f':
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
                tmpstring += "      "
                for k in b.species:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
                filestring += tmpstring
            it += 1
        filestring += "LAMDA....=    2.5000 AMAX....=    5.5000 BMAX....=    5.5000\n"
        return filestring
    
################################################################################################
class CrystalFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed by Crystal and the method
    __str__ that outputs the contents of the input file as a string.
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
        # Set unit for length scale
        self.cell.newunit("angstrom")
    def __str__(self):
        filestring = ""
        
        
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
                if len(b.species) > 1:
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    for k in b.species:
                        tmpstring += "%3i"%ed.elementnr[k]
                tmpstring += " l "+chr(it+96)+"   # "
                for k in b.species:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
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
                if len(b.species) > 1:
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    for k in b.species:
                        tmpstring += "%3i"%ed.elementnr[k]
                tmpstring += " l "+chr(it+96)+"   # "
                for k in b.species:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
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
                if len(b.species) > 1:
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    for k in b.species:
                        tmpstring += "%3i"%ed.elementnr[k]
                tmpstring += " l "+chr(it+96)+"   # "
                for k in b.species:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
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
        # Space group 
        filestring += str(self.spacegroupnr)+"\n"
        system = crystal_system(self.spacegroupnr)
        # Space group setting and crystal parameters
        if system == "triclinic":
            filestring += "0 0 0\n"
            filestring += "%13.8f %13.8f %13.8f %13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c, self.alpha, self.beta, self.gamma)
        elif system == "monoclinic":
            filestring += "0 0 0\n"
            filestring += "%13.8f %13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c, self.beta)
        elif system == "orthorhombic":
            filestring += "0 0 0\n"
            filestring += "%13.8f %13.8f %13.8f\n"%(self.a, self.b, self.c)
        elif system == "tetragonal":
            filestring += "0 0 0\n"
            filestring += "%13.8f %13.8f\n"%(self.a, self.c)
        elif system == "trigonal":
            if self.trigonalsetting == "H":
                filestring += "0 0 0\n"
                filestring += "%13.8f %13.8f\n"%(self.a, self.c)
            elif self.trigonalsetting == "R":
                filestring += "0 1 0\n"
                filestring += "%13.8f %13.8f\n"%(self.a, self.alpha)
            else:
                return "ERROR: No such trigonal setting : "+self.trigonalsetting
        elif system == "hexagonal":
            filestring += "0 0 0\n"
            filestring += "%13.8f %13.8f\n"%(self.a, self.c)
        elif system == "cubic":
            filestring += "0 0 0\n"
            filestring += "%13.8f\n"%(self.a)
        else:
            return "ERROR: Could not determine crystal system corresponding to space group "+str(self.spacegroupnr)+"."
        # Number of atoms
        filestring += str(len(self.cell.sitedata))+"\n"
        # Atomic numbers and positions
        for site in self.cell.sitedata:
            spcstring = ""
            if len(site[1]) > 1:
                # Don't know what to put for alloy
                filestring += "??"
                for k in site[1]:
                    spcstring += str(ed.elementnr[k])+"/"
                spcstring = spcstring.rstrip("/")+"  "
                for k in site[1]:
                    spcstring += k+"/"
                spcstring = spcstring.rstrip("/")
            else:
                for k in site[1]:
                    filestring += str(ed.elementnr[k]).rjust(2)
                    spcstring = k
            for coord in site[0]:
                filestring += " %13.10f"%coord
            filestring += "      ! "+spcstring+"\n"
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
        tmpstring=' \''+self.HermannMauguin+'\''
        tmpstring = tmpstring.ljust(50)+': hrmg\n'
        filestring += tmpstring
        tmpstring = ""
        tmpstring += ' %15.11f' % (self.a)
        tmpstring += ' %15.11f' % (self.b)
        tmpstring += ' %15.11f' % (self.c)
        tmpstring = tmpstring.ljust(50)+': a, b, c\n'
        filestring += tmpstring
        tmpstring = ' %15.9f %15.9f %15.9f'%(self.gamma,self.beta,self.alpha)
        tmpstring = tmpstring.ljust(50)+': ab, ac, bc\n'
        filestring += tmpstring
        tmpstring = ""
        for i in self.supercelldims:
            tmpstring += str(i)+"  "
        tmpstring = tmpstring.ljust(50)
        tmpstring += ": ncell\n"
        filestring += tmpstring
        filestring += '.true.'.ljust(50)+': primcell\n'
        # Get species info
        species = set([])
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            species.add(spcstring)
        tmpstring = str(len(species)).ljust(50)+': nspecies\n'
        filestring += tmpstring
        for spcs in species:
            # find number of representative sites for this species
            spcsites = 0
            positionstring = ""
            for site in self.cell.sitedata:
                spcstring = ""
                for k in site[1]:
                    spcstring += k+"/"
                spcstring = spcstring.rstrip("/")
                if spcstring == spcs:
                    spcsites += 1
                    for coord in site[0]:
                        positionstring += " %15.11f" % coord
                    positionstring += "\n"
            # output species info
            if len(spcs) > 2:
                # alloy
                spcsheader = '\'??\''.ljust(50)+'! '+spcs+'\n'+str(spcsites)+'\n'
            else:
                spcsheader = '\''+spcs+'\'\n'+str(spcsites)+'\n'
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
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            species.add(spcstring)
        tmpstring = "  "+str(len(species)).ljust(38)+': nspecies\n'
        filestring += tmpstring
        # local B-field string
        bfcmtstring = "   0.00000000  0.00000000  0.00000000"
        # initialize some stuff
        natoms = 0
        spcstring = ""
        for k in self.cell.sitedata[0][1]:
            spcstring += k+"/"
        spcstring = spcstring.rstrip("/")
        positionstring = ""
        tofilestring = ""
        for site in self.cell.sitedata:
            spcs = ""
            for k in site[1]:
                spcs += k+"/"
            spcs = spcs.rstrip("/")
            # Accumulate for this species
            if spcs == spcstring:
                natoms += len(site[2])
                for pos in site[2]:
                    for coord in pos:
                        positionstring += " %15.11f" % coord
                    positionstring += bfcmtstring+"\n"
            else:
                # Print species
                if len(spcstring) > 2:
                    # alloy
                    filestring += '\'??.in\''.ljust(37)+': spfname = '+spcstring+'\n'
                else:
                    filestring += '\''+spcstring+'.in\''.ljust(37)+': spfname\n'
                filestring += "  "+str(natoms)+"\n"
                filestring += positionstring
                # Initialize next species
                spcstring = spcs
                natoms = len(site[2])
                positionstring = ""
                for pos in site[2]:
                    for coord in pos:
                        positionstring += " %15.11f" % coord
                    positionstring += bfcmtstring+"\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += '\'??.in\''.ljust(37)+': spfname = '+spcstring+'\n'
        else:
            filestring += '\''+spcstring+'.in\''.ljust(37)+': spfname\n'
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
        spcstring = ""
        for k in self.cell.sitedata[0][1]:
            spcstring += k+"/"
        spcstring = spcstring.rstrip("/")
        positionstring = ""
        tofilestring = ""
        for site in self.cell.sitedata:
            spcs = ""
            for k in site[1]:
                spcs += k+"/"
            spcs = spcs.rstrip("/")
            # Accumulate for this species
            if spcs == spcstring:
                for pos in site[2]:
                    positionstring += "      <atom coord=\""
                    for coord in pos:
                        positionstring += " %15.11f" % coord
                    positionstring += "\" bfcmt=\""+bfcmtstring+"\"></atom>\n"
            else:
                # Print species
                if len(spcstring) > 2:
                    # alloy
                    filestring += '\'??.in\''.ljust(37)+': spfname = '+spcstring+'\n'
                else:
                    filestring += "    <species speciesfile=\""+spcstring+".xml\">\n"
                filestring += positionstring
                filestring += "    </species>\n"
                # Initialize next species
                spcstring = spcs
                positionstring = ""
                for pos in site[2]:
                    positionstring += "      <atom coord=\""
                    for coord in pos:
                        positionstring += " %15.11f" % coord
                    positionstring += "\" bfcmt=\""+bfcmtstring+"\"></atom>\n"
        # Print last species
        if len(spcstring) > 2:
            # alloy
            filestring += '\'??.in\''.ljust(37)+': spfname = '+spcstring+'\n'
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
    __str__ that outputs the contents of an elk/exciting.in file as a string.
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
            for coord in pos:
                tmpstring += "  %13.10f"%coord
            tmpstring += "    !  a%1i\n"%n
            n += 1
        filestring += tmpstring
        # Scale factor
        filestring += "%13.9f    ! aa\n"%self.cell.lengthscale
        filestring += "1.0  1.0  1.0 ! scale(1), scale(2), scale(3)\n"
        # Atoms
        natom = 0
        nspcs = 0
        spcs = ""
        coordstring = ""
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            natom += len(site[2])
            # Check for alloy
            if len(site[1]) > 1:
                prestring = "??"
                poststring = "  ! "
                for k in site[1]:
                    poststring += str(ed.elementnr[k])+"/"
                poststring = poststring.rstrip("/")+" "+spcstring+"\n"
            else:
                prestring = str(ed.elementnr[spcstring]).rjust(2)
                poststring = "  ! "+spcstring+"\n"
            for pos in site[2]:
                coordstring += prestring
                for coord in pos:
                    coordstring += " %13.10f"%coord
                coordstring += poststring
        # To filestring
        filestring += str(natom)+"\n"
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
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring+"\n"
        filestring += "%BLOCK LATTICE_CART\n"
        # units
        filestring += "ang    # angstrom units\n"
        # lattice
        for vec in lattice:
            for coord in vec:
                filestring += " %19.15f"%(coord*a)
            filestring += "\n"
        # Cutoff
        filestring += "%ENDBLOCK LATTICE_CART\n\n"
        # The atom position info
        filestring += "%BLOCK POSITIONS_FRAC\n"
        for site in self.cell.sitedata:
            spcstring = ""
            l = ""
            for k in site[1]:
                spcstring += k+"/"
                if l == "":
                    l = ed.elementblock[k]
                else:
                    if ed.angularmomentum[ed.elementblock[k]] > ed.angularmomentum[l]:
                        l = ed.elementblock[k]
            spcstring = spcstring.rstrip("/")
            for pos in site[2]:
                filestring += spcstring.rjust(max(2,len(spcstring)))+" "
                for coord in pos:
                    filestring += " %19.15f"%coord
                filestring += "\n"
        filestring += "%ENDBLOCK POSITIONS_FRAC\n"
        # pseudo-potential block
        filestring += "\n"
        filestring += "%BLOCK SPECIES_POT\n"
        filestring += "%ENDBLOCK SPECIES_POT\n"
        # Put the symmetry operations
        filestring += "\n%BLOCK SYMMETRY_OPS\n"
        latvect = self.cell.conventional_latticevectors()
        # make sure that identity comes first
        identity = SymmetryOperation(['x','y','z'])
        if self.cell.symops[0] != identity:
            symops = copy.deepcopy(self.cell.symops)
            symops.remove(identity)
            symops.insert(0,identity)
        else:
            symops = self.cell.symops
        k = 1
        for op in symops:
            filestring += "# Symm. op. %i\n"%k
            filestring += str(op)
            k += 1
        filestring += "%ENDBLOCK SYMMETRY_OPS\n"        
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
        for site in self.cell.sitedata:
            spcstring = ""
            l = ""
            for k in site[1]:
                spcstring += k+"/"
                if l == "":
                    l = ed.elementblock[k]
                else:
                    if ed.angularmomentum[ed.elementblock[k]] > ed.angularmomentum[l]:
                        l = ed.elementblock[k]
            spcstring = spcstring.rstrip("/")
            # pseudo-potential for each type
            filestring += "*[pseudopotential file for "+spcstring+" here]\n"
            filestring += "  LMAX="+l.upper()+"\n"
            filestring += "  "+str(len(site[2]))+"\n"
            for pos in site[2]:
                v = mvmult3(transmtx,pos)
                for coord in v:
                    filestring += " %19.15f"%coord
                filestring += "\n"
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
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring
        # The atom position info
        filestring += "AtomicCoordinatesFormat".ljust(28)+"Fractional\n"
        Alloy = False
        natom = 0
        nspcs = 0
        spcs = ""
        coordstring = ""
        specieslabels = ""
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            natom += len(site[2])
            if spcs != spcstring:
                nspcs += 1
                specieslabels += str(nspcs).ljust(4)+"   "
                # Check for alloy
                if len(site[1]) > 1:
                    species = "    ??    # "+spcstring
                    specieslabels += "??".rjust(4)+"??".rjust(7)+"     # "
                    for k in site[1]:
                        specieslabels += str(ed.elementnr[k])+"/"
                    specieslabels = specieslabels.rstrip("/")
                    specieslabels += "   "+spcstring+"\n"
                else:
                    for k in site[1]:
                        species = "  "+k.rjust(4)
                    specieslabels += str(ed.elementnr[species.lstrip()]).rjust(4)+" "+species+"\n"
            for pos in site[2]:
                for coord in pos:
                    coordstring += "%19.15f "%coord
                coordstring += species.rjust(2)+"\n"
            spcs = spcstring
        filestring += "LatticeConstant".ljust(28)+str(a)+" Ang\n"
        filestring += "NumberOfAtoms".ljust(28)+str(natom)+"\n"
        filestring += "NumberOfSpecies".ljust(28)+str(nspcs)+"\n"
        # lattice
        filestring += "%block LatticeVectors\n"
        for vec in lattice:
            for coord in vec:
                filestring += "%19.15f "%coord
            filestring += "\n"
        filestring += "%endblock LatticeVectors\n"
        # Atomic coordinates
        filestring += "%block AtomicCoordinatesAndAtomicSpecies\n"
        filestring += coordstring
        filestring += "%endblock AtomicCoordinatesAndAtomicSpecies\n"
        # Chemical species
        filestring += "%block ChemicalSpeciesLabel\n"
        filestring += specieslabels
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
            for coord in vec:
                filestring += "%19.15f "%coord
            filestring += "\n        "
        filestring += "\n"
        # The atom position info
        alloy = False
        spcs = ""
        typatstring = "typat   "
        natom = 0
        ntypat = 0
        znuclstring = "znucl   "
        alloystring = ""
        xredstring = "xred   "
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            if spcs == spcstring:
                natom += len(site[2])
            else:
                natom += len(site[2])
                ntypat += 1
                # Check for alloy
                if len(site[1]) > 1:
                    alloy = True
                    znuclstring += "?? "
                    alloystring += spcstring+" "
                else:
                    for k in site[1]:
                        znuclstring += str(ed.elementnr[k])+" "
            for pos in site[2]:
                typatstring += str(ntypat)+" "
                for coord in pos:
                    xredstring += "%19.15f "%coord
                xredstring += "\n       "
            spcs = spcstring
        filestring += "natom   "+str(natom)+"\n"
        filestring += "ntypat  "+str(ntypat)+"\n"
        filestring += typatstring+"\n"
        filestring += znuclstring
        if alloy:
            filestring += " # "+alloystring
        filestring += "\n"
        filestring += xredstring+"\n"
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
        # make sure the docstring goes on one line
        self.docstring = self.docstring.replace("\n"," ")
    def SpeciesOrder(self):
        """
        Return a string with the species in the cell in the order they come in the
        input CrystalStructure.sitedata.
        """
        returnstring = ""
        spcs = ""
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            if spcs == spcstring:
                pass
            else:
                returnstring += " "+spcstring
                spcs = spcstring
        return returnstring
    def __str__(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        # For output of atomic positions
        if self.printcartpos:
            positionunits = "Cartesian\n"
            transmtx = []
            for i in range(3):
                transmtx.append([])
                for j in range(3):
                    transmtx[i].append(lattice[i][j] * a)
        else:
            positionunits = "Direct\n"
            transmtx = [[1, 0, 0],
                        [0, 1, 0],
                        [0, 0, 1]]
        # The first line with info from input docstring
        firstline = self.docstring+" Species order: "
        # loop over sites incrementing the position string, number of sites for each species and
        # collect information about the species at the end of the firstline string
        spcs = ""
        nsitestring = ""
        for site in self.cell.sitedata:
            spcstring = ""
            for k in site[1]:
                spcstring += k+"/"
            spcstring = spcstring.rstrip("/")
            if spcs == spcstring:
                nsites += len(site[2])
            else:
                if spcs != "":
                    nsitestring += " "+str(nsites)
                nsites = len(site[2])
                firstline += " "+spcstring
                spcs = spcstring
            for pos in site[2]:
                v = mvmult3(transmtx,pos)
                for coord in v:
                    positionunits += "%19.15f "%coord
                positionunits += "\n"
        nsitestring += " "+str(nsites)+"\n"
        firstline += "\n"
        # Write first string
        filestring = firstline
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
        # Species info
        filestring += nsitestring
        # All the sites
        filestring += positionunits
        return filestring

################################################################################################
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
        self.latticenr = "14"
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
        if 'TMPDIR' in os.environ:
            tmpstring += os.environ['TMPDIR']
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
        tmpstring = "NITER.= 50 NLIN.= 31 NPRN.=  0 NCPA.= 20 NT...=%3i"%len(self.cell.sitedata)+" MNTA.="
        # Work out maximal number of species occupying a site
        mnta = 1
        for site in self.cell.sitedata:
            mnta = max(mnta,len(site[1]))
        tmpstring += "%2i"%mnta+"\n"
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
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
        # average wigner-seitz radius
        volume = abs(det3(self.cell.latticevectors))
        wsr = self.cell.lengthscale * 3*volume/(nosites * 4 * pi)**third
        filestring += "SWS......=%8f   NSWS.=  1 DSWS..=   0.05 ALPCPA= 0.9020\n"%wsr
        filestring += "Setup: 2 + NQ*NS lines\n"
        filestring += "EFGS...=  0.000 HX....=  0.100 NX...= 11 NZ0..=  6 STMP..= Y\n"
        filestring += "Symb   IQ IT ITA NZ  CONC   Sm(s)  S(ws) WS(wst) QTR SPLT\n"
        iq = 1
        it = 1
        # type loop
        for site in self.cell.sitedata:
            # alloy component loop
            ita = 1
            for comp in site[1]:
                # site loop
                for pos in site[2]:
                    tmpstring = comp.ljust(4)+"  "+"%3i%3i%3i"%(iq,it,ita)
                    tmpstring += "%4i"%ed.elementnr[comp]
                    tmpstring += "%7.3f%7.3f%7.3f%7.3f"%(site[1][comp],1,1,1)
                    tmpstring += "%5.2f%5.2f"%(0,0)
                    tmpstring += "\n"
                    filestring += tmpstring
                    iq += 1
                ita += 1
                iq -= len(site[2])
            iq += len(site[2])
            it += 1
        filestring += "Atom:  4 lines + NT*NTA*6 lines\n"
        filestring += "IEX...=  4 NP..= 251 NES..= 15 NITER=100 IWAT.=  0 NPRNA=  0\n"
        filestring += "VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n"
        filestring += "DX.......=  0.030000 DR1.....=  0.002000 TEST....=  1.00E-12\n"
        filestring += "TESTE....=  1.00E-12 TESTY...=  1.00E-12 TESTV...=  1.00E-12\n"
        for site in self.cell.sitedata:
            for comp in site[1]:
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
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
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
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (self.cell.latticevectors[i][0],self.cell.latticevectors[i][1],self.cell.latticevectors[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        v = [0.0,0.0,0.0]
        for site in self.cell.sitedata:
            for pos in site[2]:
                v = mvmult3(self.cell.latticevectors,pos)
                tmpstring = "QX(IQ)...=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                tmpstring +=  "      "
                for k in site[1]:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+" "+str(site[3])+"\n"
                filestring += tmpstring
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
        self.latticenr = 14
        # To be put on the first line
        self.programdoc = ""
    def __str__(self):
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
        for site in self.cell.sitedata:
            for i in site[1]:
                if ed.elementblock[i] == 'p':
                    maxl = max(maxl,2)
                elif ed.elementblock[i] == 'd':
                    maxl = max(maxl,3)
                elif ed.elementblock[i] == 'f':
                    maxl = max(maxl,4)
        tmpstring = "NL.....= %1i NLH...=11 NLW...= 9 NDER..= 6 ITRANS= 3 NPRN..= 0\n" % maxl
        filestring += tmpstring
        # Get number of sites
        nosites = 0
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
        # Setting the real space summation cutoff to 4.5*(wigner-seitz radius)
        volume = abs(det3(self.cell.latticevectors))
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
                tmpstring += "BSX......=%10f BSY.....=%10f BSZ.....=%10f\n" % (self.cell.latticevectors[i][0],self.cell.latticevectors[i][1],self.cell.latticevectors[i][2])
        else:
            tmpstring +=  "ALPHA....=%10f BETA....=%10f GAMMA...=%10f\n" % (self.alpha, self.beta, self.gamma)
        filestring += tmpstring
        v = [0.0,0.0,0.0]
        for site in self.cell.sitedata:
            for pos in site[2]:
                v = mvmult3(self.cell.latticevectors,pos)
                tmpstring = "QX(IQ)...=%10f QY......=%10f QZ......=%10f" % (v[0],v[1],v[2])
                tmpstring +=  "      "
                for k in site[1]:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+" "+str(site[3])+"\n"
                filestring += tmpstring
        for i in range(nosites):
            filestring += "a/w(IQ)..="
            for i in range(4):
                filestring += "%5.2f"%self.hardsphere
            filestring += "\n"
        filestring += "LAMDA....=    2.5000 AMAX....=    4.5000 BMAX....=    4.5000\n"
        return filestring
