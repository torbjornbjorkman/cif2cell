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
from __future__ import division
import os
import sys
import string
import copy
import CifFile
from types import *
from math import sin,cos,pi,sqrt
from spacegroupdata import *
from elementdata import *

################################################################################################
# Miscellaneous
zero = 0.0
one = 1.0
two = 2.0
three = 3.0
four = 4.0
six = 6.0
third = 1/3
half = one/two
fourth = one/four
sixth = one/six
occepsilon = 0.000001
floatlist = [third, 2*third, half, fourth, one, zero, sqrt(2.0),sixth,5*sixth]

################################################################################################
# Exception classes
class SymmetryError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class PositionError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class CellError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class GeometryObjectError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
class SetupError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)
    
################################################################################################
class GeometryObject:
    """
    Parent class for anything geometrical contains:
    compeps  : epsilon for determining when two floats are equal
    """
    def __init__(self):
        self.compeps = 0.0002

class Charge(float):
    """
    Class for representing the charge state/oxidation number of an atom (ion).
    It is just an integer, except for a modified routine for turning into a string,
    viz. a two plus ion gets the string representation '2+'
    """
    def __init__(self,i):
        float.__init__(self)
    def __str__(self):
        if abs(self-int(self)) < 0.0001:
            if int(self) == 0:
                return '0'
            elif self > 0:
                return str(abs(int(self)))+'+'
            elif self < 0:
                return str(abs(int(self)))+'-'
        else:
            if self > 0:
                return str(abs(self))+"+"
            else:
                return str(abs(self))+"-"

class Vector(list,GeometryObject):
    """
    Floating point vector describing a lattice point. Assumed to have length three.
    * Supports testing for equality using the self.compeps parameter inherited
    from the parent class GeometryObject
    * Supports testing for < using the euclidean norm (length) of the vector
    * Supports addition with another Vector, in which case the
      vectors are added component-wise.
    More methods:
    length            : returns the euclidean norm of the vector
    transform(matrix) : returns matrix*vector
    improveprecision  : identify some conspicuous numbers and improve precision
    """
    def __init__(self, vec):
        GeometryObject.__init__(self)
        list.__init__(self, vec)
    def __eq__(self,other):
        for i in range(3):
            if abs(self[i]-other[i]) > self.compeps:
                return False
        return True
    def length(self):
        return sqrt(self[0]**2+self[1]**2+self[2]**2)
    def __lt__(self, other):
        sl = self[0]**2+self[1]**2+self[2]**2
        ol = other[0]**2+other[1]**2+other[2]**2
        return sl < ol
    # Addition of two vectors
    def __add__(self, other):
        t = Vector([self[i]+other[i] for i in range(3)])
        return t
    def __str__(self):
        return "%19.15f %19.15f %19.15f"%(self[0],self[1],self[2])
    # multiplication by scalar
    def scalmult(self, a):
        for i in range(3):
            self[i] *= a
        return self
    # dot product
    def dot(self,a):
        t = 0.0
        for i in range(3):
            t += self[i]*a[i]
        return t
    # coordinate transformation
    def transform(self, matrix):
        t = Vector(mvmult3(matrix, self))
        return t
    def improveprecision(self):
        for i in range(3):
            for f in floatlist:
                if abs(self[i]-f) <= self.compeps:
                    # 0
                    self[i] = f
                    break

class LatticeVector(Vector):
    """
    Vector of length three that maps back things into the cell
    """
    def __init__(self, vec):
        Vector.__init__(self, vec)
        # Interval we wish to use for the coordinates.
        # In practice either [0,1] or [-.5, 0.5]
        self.interval = [0.0, 1.0]
        self.intocell()
        self.improveprecision()
    # Addition of two vectors, putting the result back
    # into the cell if necessary
    def __add__(self, other):
        if self.interval[0] != other.interval[0] or self.interval[1] != other.interval[1]:
            raise GeometryObjectError("LatticeVectors must have the same definition intervals to be added.")
        t = LatticeVector([self[i]+other[i] for i in range(3)])
        t.intocell()
        return t
    # Change interval
    def change_interval(self, interval):
        self.interval = interval
        t = LatticeVector([0,0,0])
        t.interval = interval
        self = self + t
    # coordinate transformation
    def transform(self, matrix):
        t = LatticeVector(mvmult3(matrix, self))
        t.intocell()
        return t
    # Put the vector components into the cell interval defined by self.interval
    def intocell(self):
        for i in range(3):
            while not (self.interval[0] <= self[i] < self.interval[1]):
                self[i] -= copysign(1,self[i])

class LatticeMatrix(GeometryObject, list):
    """
    Three by three matrix
    """
    def __init__(self,mat):
        GeometryObject.__init__(self)
        t = []
        for vec in mat:
            t.append(Vector(vec))
        list.__init__(self, t)
    def __str__(self):
        matstr = ""
        for l in self:
            matstr += str(l)+"\n"
        return matstr
    def __eq__(self,other):
        for i in range(3):
            for j in range(3):
                if abs(self[i][j]-other[i][j]) > self.compeps:
                    return False
        return True
    # coordinate transformation
    def transform(self, matrix):
        return LatticeMatrix(mmmult3(matrix, self))

class AtomSite(GeometryObject):
    """
    Class for describing an atomic site.

    Contains data:
    position  : a vector that gives the position
    species   : a dictionary with element-occupancy pairs (e.g. {Fe : 0.2, Co : 0.8})
    label     : any label
    charge    : the charge state (oxidation number) of the species (e.g. '2-' for oxygen in most compounds)
    index     : any integer

    Functions:
        __eq__    : compare equality
        __str__   : one line with species and position info
        spcstring : species string ('Mn', 'La/Sr' ...)
        alloy     : true if there are more than one species occupying the site
    """
    def __init__(self,position=None,species=None,label="",charge=0,index=None):
        GeometryObject.__init__(self)
        if position != None:
            self.position = LatticeVector(position)
        else:
            self.position = None
        if species != None:
            self.species = species
        else:
            self.species = {}
        self.label = label
        self.charge = Charge(charge)
        self.index = index
    def __eq__(self,other):
        return self.position == other.position and self.species == other.species
    # Species string
    def spcstring(self):
        tmp = ""
        for k in self.species:
            tmp += k+"/"
        tmp = tmp.rstrip("/")
        return tmp
    # Is there more than one species on this site?
    def alloy(self):
        return len(self.species) > 1
    # print site data in some informative way
    def __str__(self):
        # Element symbol
        tmp = self.spcstring().ljust(8)
        # Position
        tmp += " %19.15f %19.15f %19.15f   "%(self.position[0],self.position[1],self.position[2])
        # occupancy
        for k,v in self.species.iteritems():
            tmp += str(v)+"/"
        tmp = tmp.rstrip("/")
        return tmp
        
class SymmetryOperation(GeometryObject):
    """
    Class describing a symmetry operation, with a rotation matrix and a translation.
    """
    def __init__(self, eqsite=None):
        GeometryObject.__init__(self)
        self.eqsite = eqsite
        if self.eqsite != None:
            self.rotation = self.rotmat()
            self.translation = self.transvec()
        else:
            self.rotation = None
            self.translation = None
    # This way of printing was useful for outputting to CASTEP.
    def __str__(self):
        return str(self.rotation)+str(self.translation)+"\n"
    # Two symmetry operations are equal if rotation matrices and translation vector
    # differ by at most compeps
    def __eq__(self, other):
        eq = True
        for i in range(3):
            for j in range(3):
                eq = eq and abs(self.rotation[i][j] - other.rotation[i][j]) < self.compeps
            eq = eq and self.translation == other.translation
        return eq
    # comparison between operations made by comparing lengths of translation vectors
    def __lt__(self, other):
        return self.translation < other.translation
    # Return a rotation matrix from "x,y,z" representation of a symmetry operation
    def rotmat(self):
        mat = [[0,0,0],[0,0,0],[0,0,0]]
        for j in range(len(self.eqsite)):
            xyz = self.eqsite[j].replace('+',' +').replace('-',' -').split()
            for i in xyz:
                if i.strip("+-") == 'x':
                    mat[0][j] = float(i.strip('x')+"1")
                elif i.strip("+-") == 'y':
                    mat[1][j] = float(i.strip('y')+"1")
                elif i.strip("+-") == 'z':
                    mat[2][j] = float(i.strip('z')+"1")            
        return LatticeMatrix(mat)
    # Return a translation vector from "x,y,z" representation of a symmetry operation
    def transvec(self):
        vec = []
        for i in range(3): vec.append(0.0)
        for j in range(len(self.eqsite)):
            xyz = self.eqsite[j].replace('+',' +').replace('-',' -').split()
            for i in xyz:
                if i.strip("+-xyz") != "":
                    vec[j] = eval(i)
        return LatticeVector(vec)
    # True if the operation is diagonal
    def diagonal(self):
        if abs(self.rotation[0][1]) < self.compeps and \
           abs(self.rotation[0][2]) < self.compeps and \
           abs(self.rotation[1][0]) < self.compeps and \
           abs(self.rotation[1][2]) < self.compeps and \
           abs(self.rotation[2][0]) < self.compeps and \
           abs(self.rotation[2][1]) < self.compeps:
            return True
        else:
            return False
        
    
################################################################################################
class CellData(GeometryObject):
    """
    Class for a lot of stuff specifying a unit cell.
    latticevectors : The Bravais lattice vectors
    lengthscale    : an overall length scale that multiplies the lattice vectors.
    unit           : the unit of the lengthscale
    alloy          : True if the compound is an alloy
    atomdata       : An array of arrays of AtomSite objects, in other words, a collection of 
                     collections of atoms. The getCrystalStructure method sets up a collection of
                     all the inequivalent sites, and each such collection contains all the
                     sites generated by the representative site. So: atomdata[0][2] is the
                     third atom generated from the first wyckoff position.
    Methods:
        getFromCIF          : obtain data for setting up a cell from a CIF block
        crystal_system      : return a string with the name of the crystal system
        latticevectors      : Return the Bravais lattice vectors as a 3x3 matrix
        reciprocal_latticevectors: Return the reciprocal lattice vectors as a 3x3 matrix.
        volume              : Return the unit cell volume
        primitive           : Returns a CrystalStructure object for the primitive cell
        conventional        : Returns a CrystalStructure object for the conventional cell.
        getCrystalStructure : The 'primitive' and 'conventional' methods are just
                              wrappers around this method. It requires the following
                              to be set beforehand:
                                  a, b, c, alpha, beta, gamma : the lattice parameters
                                  spacegroupnr : The space group number
                                  ineqsites : The inequivalent sites (wyckoff positions)
                                  occupations : The occupations of the different inequivalent sites
                                                in the form of a list of dictionaries. Each
                                                inequivalent site is supposed to have a dictionary
                                                with the species occupying it and its occupancy.
                                                Example: Two inequivalent sites, one with iron
                                                and the other with 90% oxygen and 10% flourine:
                                                    occupations = [{'Fe':1.0}, {'O':0.9, 'F':0.1}]
        getSuperCell        : Return a supercell from input supercell dimensions [i,j,k]. The crystal
                              structure must first have been initialized by getCrystalStructure (or
                              'primitive' or 'conventional').
    """
    def __init__(self):
        GeometryObject.__init__(self)
        self.initialized = False
        self.quiet = False
        self.coordepsilon = 0.0002
        self.spacegroupnr = 0
        self.spacegroupsymbol = ""
        self.spacegroupsymboltype = ""
        self.spacegroupsetting = ""
        self.lattrans = LatticeMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.transvecs = [LatticeVector([zero, zero, zero])]
        self.symops = []
        self.ineqsites = []
        self.occupations = []
        self.atomdata = []
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 0
        self.beta = 0
        self.gamma = 0
        self.coa = 1
        self.boa = 1
        self.latticevectors = LatticeMatrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        self.lengthscale = 1
        self.unit = "angstrom"
        self.alloy = False
        self.numofineqsites = 0
        self.primcell=False

    def newunit(self,newunit="angstrom"):
        """ Set new unit for the length scale. Valid choices are:
            * angstrom
            * bohr  (bohr radii, or a.u. (atomic unit))
            * nm    (nanometer)
        """
        if self.unit == newunit:
            return
        if self.unit == "angstrom" and newunit == "bohr":
            fact = 1.8897261
        elif self.unit == "bohr" and newunit == "angstrom":
            fact = 0.52917721
        elif self.unit == "angstrom" and newunit == "nm":
            fact = 0.1
        elif self.unit == "nm" and newunit == "angstrom":
            fact = 10
        elif self.unit == "bohr" and newunit == "nm":
            fact = 0.052917721
        elif self.unit == "nm" and newunit == "bohr":
            fact = 18.897261
        self.lengthscale *= fact
        self.unit = newunit

    def crystal_system(self):
        return crystal_system(self.spacegroupnr)

    def conventional_latticevectors(self):
        # Set up Bravais lattice vectors of the conventional cell
        self.coa = self.c / self.a
        self.boa = self.b / self.a
        alphar = self.alpha*pi/180
        betar  = self.beta*pi/180
        gammar = self.gamma*pi/180
        if self.crystal_system() == 'cubic':
            latticevectors = LatticeMatrix([[one, zero, zero],
                                            [zero, one, zero],
                                            [zero, zero, one]])
        elif self.crystal_system() == 'hexagonal':
            latticevectors = LatticeMatrix([[sin(gammar), cos(gammar), zero],
                                            [zero, one, zero],
                                            [zero, zero, self.coa]])
        elif self.crystal_system() == 'tetragonal':
            latticevectors = LatticeMatrix([[one, zero, zero],
                                            [zero, one, zero],
                                            [zero, zero, self.coa]])
        elif self.crystal_system() == 'orthorhombic':
            latticevectors = LatticeMatrix([[one, zero, zero], 
                                            [zero, self.boa, zero], 
                                            [zero, zero, self.coa]])
        elif self.crystal_system() == 'monoclinic':
            latticevectors = LatticeMatrix([[one, zero, zero], 
                                            [zero, self.boa, zero], 
                                            [self.coa*cos(betar), zero, self.coa*sin(betar)]])
        elif self.crystal_system() == 'trigonal':
            # Hexagonal cell normally used
            if abs(self.gamma-120) < self.coordepsilon:
                latticevectors = LatticeMatrix([[sin(gammar), cos(gammar), zero],
                                                [zero, one, zero],
                                                [zero, zero, self.coa]])
            else:
                # Symmetric rhombohedral cell (stretching a cube along the main diagonal)
                c = sqrt((1 + 2*cos(alphar))/(1 - cos(alphar)))
                a = pow(1/c, 1/3)
                t = a * (c + 2) / 3
                u = a * (c - 1) / 3
                latticevectors = LatticeMatrix([[t, u, u],
                                                [u, t, u],
                                                [u, u, t]])
        elif self.crystal_system() == 'triclinic' or self.crystal_system() == 'unknown':
            angfac1 = (cos(alphar) - cos(betar)*cos(gammar))/sin(gammar)
            angfac2 = sqrt(sin(gammar)**2 - cos(betar)**2 - cos(alphar)**2 
                       + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)
            latticevectors = LatticeMatrix([[one, zero, zero], 
                                            [self.boa*cos(gammar), self.boa*sin(gammar), zero], 
                                            [self.coa*cos(betar), self.coa*angfac1, self.coa*angfac2]])
        else:
            raise SymmetryError("No support for "+self.crystal_system()+" crystal systems.")
        return latticevectors

    # The reciprocal lattice vectors corresponding to the lattice vectors of the structure
    # (b1, b2, b3)^T = 2 * pi * (a1, a2, a3)^{-1}
    def reciprocal_latticevectors(self):
        t = minv3(self.latticevectors)
        reclatvect = []
        for j in range(3):
            reclatvect.append([])
            for i in range(3):
                reclatvect[j].append(t[i][j]*2*pi)
        return LatticeMatrix(reclatvect)

    # Define comparison functions
    def poscomp(self, pos1, pos2):
        # Return True if two positions are the same
        if abs(pos1[0] - pos2[0]) < self.coordepsilon and \
               abs(pos1[1] - pos2[1]) < self.coordepsilon and \
               abs(pos1[2] - pos2[2]) < self.coordepsilon:
            return True
        else:
            return False
        
    def transveccomp(self, pos1,pos2):
        # Return True if two positions only differ by one
        # of the induced lattice translations
        match = False
        for tv in self.transvecs:
            if (abs(pos1[0]-(pos2[0]-tv[0]))<self.coordepsilon and \
                abs(pos1[1]-(pos2[1]-tv[1]))<self.coordepsilon and \
                abs(pos1[2]-(pos2[2]-tv[2]))<self.coordepsilon):
                match = True
        return match
    
    def duplicates(self, poslist, compfunc = transveccomp):
        # Return list of indices of duplicates in a list,
        # sorted in reverse order to be easy to use for removing the duplicates.
        # Optionally supply a comparison function for when two coordinates
        # are the same, else use 'samecoords' function from above.
        removeindices = set([])
        for i in range(len(poslist)):
            for j in range(len(poslist)-1,i,-1):
                if compfunc(self,poslist[i],poslist[j]):
                    removeindices.add(j)
        removeindices = list(removeindices)
        removeindices.sort(reverse=True)
        return removeindices

    def volume(self):
        return det3(self.latticevectors)
    
    def primitive(self):
        """ Return a CrystalStructure object for the primitive cell."""
        w = self.getCrystalStructure(reduce=True)
        return w
    
    def conventional(self):
        """ Return a CrystalStructure object for the conventional cell."""
        w = self.getCrystalStructure(reduce=False)
        return w
    
    def getCrystalStructure(self, reduce=False):
        """
        Return a CrystalStructure object, either as it is or reduced to the
        primitive cell.
        """
        # reduce to primitive cell?
        self.primcell = reduce
        
        ###################################
        #   INITIALIZE SPACE GROUP DATA   #
        ###################################
        if self.spacegroupsymbol == "":
            if 0 < self.spacegroupnr < 231:
                self.spacegroupsymbol = SpaceGroupData().HMSymbol[str(self.spacegroupnr)]
                self.spacegroupsymboltype = "(H-M)"
        else:
            self.spacegroupsetting = self.spacegroupsymbol[0]
        # Check if we know enough:
        # To produce the conventional cell (self.primcell=False) we don't need the space group
        # symbol or number as long as we have the symmetry operations (equivalent sites).
        if self.spacegroupsetting == "":
            try:
                self.spacegroupsetting = SpaceGroupData().HMSymbol[str(self.spacegroupnr)][0]
            except:
                pass
        if self.a!=0 and self.b!=0 and self.c!=0 and self.alpha!=0 and self.beta!=0 and self.gamma!=0:
            if self.spacegroupnr == -1:
                if len(self.symops) >= 1:
                    if self.primcell == True:
                        if self.spacegroupsetting == 'P':
                            self.spacegroupnr = 0
                        else:
                            raise SymmetryError("Insufficient symmetry information to reduce to primitive cell"+\
                                                " (need space group number or Hermann-Mauguin symbol).\n"+\
                                                "  Run with --no-reduce to generate cell in the conventional setting.")
                    else:
                        self.spacegroupnr = 0
        else:
            raise CellError("No crystallographic parameter may be zero.")
        self.latticevectors = self.conventional_latticevectors()
        self.lengthscale = self.a
        self.alloy = False
        
        ###############################
        #     LATTICE TRANSLATIONS    #
        ###############################
        #
        # Choices of lattice vectors made to largely coincide with the choices
        # made at http://cst-www.nrl.navy.mil/lattice/
        #
        # The induced translation vectors are from Zachariasen, "Theory of x-ray
        # diffraction in crystals". 
        #
        # Relations between rhombohedral and hexagonal settings of trigonal
        # space groups from Tilley, "Crystals and crystal structures"
        # Bravais lattice vectors:
        #
        # a_r =  2/3 a_h + 1/3 b_h + 1/3 c_h
        # b_r = -1/3 a_h + 1/3 b_h + 1/3 c_h
        # c_r = -1/3 a_h - 2/3 b_h + 1/3 c_h
        #
        # a_h = a_r - b_r
        # b_h = b_r - c_r
        # c_h = a_r + b_r + c_r
        #
        # a, c and rhombohedral angle alpha:
        #
        # a_h = 2 * a_r * sin(alpha/2)
        # c_h = a_r * sqrt(3 + 6*cos(alpha))
        # 
        # a_r = sqrt(3*a_h^2 + c_h^2) / 3
        # sin(alpha/2) = 3*a_h / (2 * sqrt(3*a_h^2 + c_h^2))
        #
        if self.primcell:
            if self.spacegroupsetting == 'I':
                # Body centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,half])]
                if self.crystal_system() == 'cubic':
                    self.lattrans = LatticeMatrix([[-half, half, half],
                                                   [half, -half, half],
                                                   [half, half, -half]])
                else:
                    self.lattrans = LatticeMatrix([[one, zero, zero],
                                                   [zero, one, zero],
                                                   [half, half, half]])
            elif self.spacegroupsetting == 'F':
                # face centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,zero]),
                                  LatticeVector([half,zero,half]),
                                  LatticeVector([zero,half,half])]
                self.lattrans = LatticeMatrix([[half, half, zero],
                                               [half, zero, half],
                                               [zero, half, half]])
            elif self.spacegroupsetting == 'A':
                # A-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([zero,half,half])]
                self.lattrans = LatticeMatrix([[one, zero, zero],
                                               [zero, half, -half],
                                               [zero, half, half]])
            elif self.spacegroupsetting == 'B':
                # B-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,zero,half])]
                self.lattrans = LatticeMatrix([[half, zero, -half],
                                               [zero, one, zero],
                                               [half, zero, half]])
            elif self.spacegroupsetting == 'C':
                # C-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,zero])]
                self.lattrans = LatticeMatrix([[half, -half, zero],
                                               [half, half, zero],
                                               [zero, zero, one]])
            elif self.spacegroupsetting == 'R':
                if abs(self.gamma - 120) < self.coordepsilon:
                    # rhombohedral from hexagonal setting
                    self.transvecs = [LatticeVector([zero,zero,zero]),
                                      LatticeVector([third, 2*third, 2*third]),
                                      LatticeVector([2*third, third, third])]
                    self.lattrans = LatticeMatrix([[2*third, third, third],
                                                   [-third, third, third],
                                                   [-third, -2*third, third]])
                else:
                    self.transvecs = [LatticeVector([zero,zero,zero])]
                    self.lattrans = LatticeMatrix([[1, 0, 0],
                                                   [0, 1, 0],
                                                   [0, 0, 1]])
            else:
                self.transvecs = [LatticeVector([zero,zero,zero])]
                self.lattrans = LatticeMatrix([[1, 0, 0],
                                               [0, 1, 0],
                                               [0, 0, 1]])
            # Transform to primitive cell
            tmp = []
            for i in range(3):
                tmp.append(mvmult3(self.latticevectors,self.lattrans[i]))
            self.latticevectors = LatticeMatrix(tmp)
            # Improve precision again...
            for i in range(3):
                for j in range(3):
                    self.latticevectors[i][j] = improveprecision(self.latticevectors[i][j],self.coordepsilon)
        else:
            # If no reduction is to be done
            self.transvecs = [LatticeVector([zero,zero,zero])]
            self.lattrans = LatticeMatrix([[1, 0, 0],
                                          [0, 1, 0],
                                          [0, 0, 1]])

        # Find inverse lattice transformation matrix
        invlattrans = LatticeMatrix(minv3(self.lattrans))
        
        #################################
        #       SYMMETRY OPERATIONS     #
        #################################
        # Get equivalent sites from tables if not present
        if self.symops == []:
            for op in SpaceGroupData().EquivalentPositions[self.spacegroupnr]:
                self.symops.append(SymmetryOperation(op))

        # sort the symmetry operations
        # 1. by length of translation vectors
        self.symops.sort()
        # 2. by determinant (1 first then -1)
        self.symops.sort(key = lambda op : det3(op.rotation), reverse=True)
        # 3. put identity first
        identity = SymmetryOperation(['x','y','z'])
        self.symops.remove(identity)
        self.symops.insert(0,identity)
        
        # If we reduce the cell, remove the symmetry operations that are the
        # same up to an induced lattice translation
        if self.primcell:
            if len(self.transvecs) > 1:
                removelist = []
                for j in range(len(self.symops)):
                    for i in range(j+1,len(self.symops)):
                        for vec in self.transvecs:
                            if self.symops[i].translation+vec == self.symops[j].translation:
                                if self.symops[i].rotation == self.symops[j].rotation:
                                    removelist.append(i)
                removelist = list(set(removelist))
                removelist.sort(reverse=True)
                for i in removelist:
                    self.symops.pop(i)

        # To cartesian representation
        lv = self.conventional_latticevectors()
        for op in self.symops:
            op.rotation = lv.transform(op.rotation)
            op.rotation = op.rotation.transform(minv3(lv))
            # transform translations
            op.translation = op.translation.transform(minv3(self.lattrans))
                                
        #########################
        #    CELL GENERATION    # 
        #########################
        # Atomic species and the number of each species. Site occupancies.
        for i in range(len(self.ineqsites)):
            # Set up atomdata
            self.atomdata.append([])
            self.atomdata[i].append(AtomSite(position=self.ineqsites[i]))
            # Add species and occupations to atomdata
            for k,v in self.occupations[i].iteritems():
                self.atomdata[i][0].species[k] = v
            # Add charge state
            self.atomdata[i][0].charge = self.charges[i]
            # Determine if we have an alloy
            for element in self.occupations[i]:
                v = self.occupations[i][element]
                if abs(1-v) > occepsilon:
                    self.alloy = True

        # Make sites unique with array of occupation info in case of an alloy
        if self.alloy:
            removeindices = []
            # set up atomdata
            for i in range(len(self.atomdata)):
                for j in range(len(self.atomdata)-1,i,-1):
                    if self.atomdata[i][0].position == self.atomdata[j][0].position:
                        # Add the dictionary of site j to that of site i and schedule index j for
                        # removal. If there is already an instance of the species on this site,
                        # then add the occupancies (this happens when different valencies has been
                        # recorded)
                        for k in self.atomdata[j][0].species:
                            if k in self.atomdata[i][0].species:
                                v = self.atomdata[j][0].species[k] + self.atomdata[i][0].species[k]
                                self.atomdata[i][0].species[k] = v
                            else:
                                self.atomdata[i][0].species[k] = self.atomdata[j][0].species[k]
                                removeindices.append(j)
                        # ...also fix self.occupations
                        self.occupations[i][self.occupations[j].keys()[0]] = self.occupations[j].values()[0]
            # Remove duplicate elements
            removeindices = list(set(removeindices))
            removeindices.sort(reverse=True)
            for i in removeindices:
                self.atomdata.pop(i)
                self.ineqsites.pop(i)
                self.occupations.pop(i)

        # Work out all sites in the cell for atomdata
        for a in self.atomdata:
            for op in self.symops:
                # position expression string
                posexpr = [s for s in op.eqsite]
                for k in range(3):
                    # position expression string, replacing x,y,z with numbers
                    posexpr[k] = posexpr[k].replace('x',str(a[0].position[0]))
                    posexpr[k] = posexpr[k].replace('y',str(a[0].position[1]))
                    posexpr[k] = posexpr[k].replace('z',str(a[0].position[2]))
                position = LatticeVector([eval(pos) for pos in posexpr])
                b = AtomSite(position=position,species=a[0].species,charge=a[0].charge)
                append = True
                for site in a:
                    for vec in self.transvecs:
                        t = vec + b.position
                        if site.position == t:
                            append=False
                            break
                    if not append:
                        break
                if append:
                    a.append(b)
        for a in self.atomdata:
            for b in a:
                b.position = LatticeVector(mvmult3(invlattrans,b.position))
                    
        ######################
        #    MISCELLANEOUS   #
        ######################
        # Set flag and return the CrystalStructure in the conventional setting
        self.initialized = True
        return self

    def getSuperCell(self, supercellmap, vacuum, transvec, sort=""):
        """
        Returns a supercell based on the input supercell map, vacuum layers and translation vector.
        The cell must have been initialized prior to calling getSuperCell.
        
        The cell will be padded with some number of original unit cells of vacuum
        by simple rescaling of the lattice vectors and positions. This is controlled by 'vacuum'.
        
        All coordinates will be translated by 'transvec', which is given in units
        of the original lattice vectors.
        """
        
        # Sanity checks
        if not self.initialized:
            raise CellError("The unit cell must be initialized before a supercell can be generated.")
        if len(supercellmap) != 3 or type(supercellmap[0]) != IntType or type(supercellmap[1]) != IntType \
               or type(supercellmap[2]) != IntType or supercellmap[0] <= 0 or supercellmap[1] <= 0 \
               or supercellmap[2] <= 0:
            raise CellError("The supercell map must be an array of three positive integers\n"+\
                            "or three arrays of three integers.")
        if len(vacuum) != 3 or vacuum[0] < 0 or vacuum[1] < 0 or vacuum[2] < 0:
            raise CellError("The vacuum padding must be an array of three numbers >= 0.")
        vectors = self.latticevectors

        # map matrix
        mapmatrix = LatticeMatrix([[supercellmap[0],0,0],[0,supercellmap[1],0],[0,0,supercellmap[2]]])
        invmapmatrix = LatticeMatrix(minv3(mapmatrix))
        volratio = det3(mapmatrix)
        
        # previous length of lattice vectors
        oldlatlen = [vec.length() for vec in self.latticevectors]
            
        # New latticevectors from supercell map
        t = mmmult3(mapmatrix,self.latticevectors)
        self.latticevectors = LatticeMatrix(t)
        invlatvects = minv3(self.latticevectors)
                
        # Set up new translation group (in new latice vector coordinates)
        offsets = []
        multfac = 0
        for k in range(supercellmap[2]):
            for j in range(supercellmap[1]):
                for i in range(supercellmap[0]):
                    offsets.append([i,j,k])
                    multfac += 1
        # Transform to new coordinates
        newtranslations = []
        for trans in offsets:
            t = LatticeVector(mvmult3(minv3(mapmatrix),trans))
            newtranslations.append(t)

        # Transform coordinates to new basis
        for a in self.atomdata:
            for b in a:
                b.position = LatticeVector(mvmult3(invmapmatrix,b.position))
        for op in self.symops:
            op.translation = LatticeVector(mvmult3(invmapmatrix,op.translation))

        # Operate with new translation group on coordinates to generate all positions
        newsites = []
        i = 0
        for a in self.atomdata:
            newsites.append([])
            for b in a:
                for translation in newtranslations[1:]:
                    position = b.position + translation
                    newsites[i].append(AtomSite(position=position,species=b.species))
            i += 1
        i = 0
        for a in newsites:
            for b in a:
                self.atomdata[i].append(b)
            i += 1

        # previous length of lattice vectors
        oldlatlen = [vec.length() for vec in self.latticevectors]
            
        # New latticevectors after vacuum padding
        vacuummapmatrix = LatticeMatrix([[1,0,0],[0,1,0],[0,0,1]])
        if reduce(lambda x,y: x+y, vacuum) > 0:
            # add the given number of unit cell units along the lattice vectors
            for j in range(len(vacuum)):
                for i in range(len(self.latticevectors[j])):
                    self.latticevectors[j][i] = self.latticevectors[j][i] + vectors[j][i]*vacuum[j]
                vacuummapmatrix[j][j] = vacuummapmatrix[j][j] + vacuum[j]
        invvacmat = LatticeMatrix(minv3(vacuummapmatrix))
        # Rescale coordinates after padding
        for a in self.atomdata:
            for b in a:
                b.position = LatticeVector(mvmult3(invvacmat, b.position))

        # Multiply group by new translation group
        newops = []
        i = 0
        for vec in newtranslations[1:]:
            for op in self.symops:
                newops.append(SymmetryOperation())
                newops[i].rotation = op.rotation
                newops[i].translation = vec + op.translation
                i += 1 
        for op in newops:
            self.symops.append(op)

        # Weed out rotations that are broken by supercell map
        lv = self.latticevectors
        removelist = []
        # Ugly set of if's and special cases
        # THIS WILL NOT WORK WITH GENERAL SUPERCELL MAP MATRIX 
        if self.crystal_system()=="hexagonal" or (self.crystal_system()=="trigonal" and not self.primcell):
            if abs(lv[0].length()-lv[1].length()) < lv.compeps:
                # if a and b are still the same, no rotation symmetry is broken
                pass
            else:
                # if a and b are different, all rotation symmetries except inversion are broken
                e = SymmetryOperation(["x","y","z"])
                i = SymmetryOperation(["-x","-y","-z"])
                j = 0
                for op in self.symops:
                    if op.rotation == e.rotation or op.rotation == i.rotation:
                        pass
                    else:
                        removelist.append(j)
                    j += 1
        elif self.crystal_system()=="trigonal" and self.primcell:
            # if any latticevector has different length, all symmetries except inversion are broken
            if lv[0].length() == lv[1].length() == lv[2].length():
                pass
            else:
                e = SymmetryOperation(["x","y","z"])
                i = SymmetryOperation(["-x","-y","-z"])
                j = 0
                for op in self.symops:
                    if op.rotation == e.rotation or op.rotation == i.rotation:
                        pass
                    else:
                        removelist.append(j)
                    j += 1
        else:
            i = 0
            for op in self.symops:
                # diagonal operations always ok (since we only have diagonal map for now)
                if op.diagonal():
                    pass
                else:
                    for vec in lv:
                        t = Vector(mvmult3(op.rotation,vec))
                        r = Vector([-u for u in t])
                        # Symmetry operation OK if it maps a lattice vector into one of the other lattice vectors
                        equivalent = t == lv[0] or t == lv[1] or t == lv[2] or \
                                     r == lv[0] or r == lv[1] or r == lv[2]
                        if not equivalent:
                            removelist.append(i)
                i += 1
        # Weed out translations broken by the vacuum
        j = 0
        for op in self.symops:
            # check if vacuum padding destroys this translation
            t = op.translation
            for i in range(3):
                if abs(t[i]) > self.compeps and abs(vacuum[i]) > self.compeps:
                    removelist.append(j)
            j += 1
        # Remove broken symmetries
        removelist = sorted(list(set(removelist)), reverse=True)
        for i in removelist:
            self.symops.pop(i)

        # Move all atoms by transvec 
        if reduce(lambda x,y: x+y, transvec) != 0:
            for a in self.atomdata:
                for b in a:
                    for k in range(3):
                        fac = oldlatlen[k]/self.latticevectors[k].length()
                        b.position[k] = b.position[k] + fac*transvec[k]
                        
        # Put stuff back in ]-1,1[ interval
        # THIS SHOULD BE DONE BETTER USING THE INTERVAL OF LATTICEVECTOR CLASS
        for a in self.atomdata:
            for b in a:
                for k in range(3):
                    while abs(b.position[k]) >= 1:
                        b.position[k] = b.position[k] - copysign(1,b.position[k])

        # Sort the atomic positions in some way
        if sort != "":
            # remove previous ordering
            tempdata = []
            i = 0
            for a in self.atomdata:
                for b in a:
                    tempdata.append([])
                    tempdata[i].append(b)
                    i += 1
            self.atomdata = tempdata
            # Identify and sort by z-layers
            if sort == "zlayer":
                # sort atoms by z-coordinate
                self.atomdata.sort(key = lambda a: a[0].position.transform(self.latticevectors)[2])
                #
            elif len(sort) == 3:
                # check if the string contains only 'x', 'y' or 'z'
                allbutxyz = string.printable.replace("x","").replace("y","").replace("z","")
                # check if the string contains only '1', '2' or '3'
                allbut123 = string.printable.replace("1","").replace("2","").replace("3","")
                # dummy table to work in python < 2.6
                table = string.maketrans("a","a")
                if len(string.translate(sort,table,allbutxyz)) == 3:
                    # Sort atoms by cartesian coordinate
                    sortnum = string.translate(sort,string.maketrans("xyz","012"))
                    for i in sortnum:
                        self.atomdata.sort(key = lambda a: a[0].position.transform(self.latticevectors)[int(i)])
                elif len(string.translate(sort,table,allbut123)) == 3:
                    # sort atoms by lattice coordinates
                    for i in sort:
                        self.atomdata.sort(key = lambda a: a[0].position[int(i)-1])
        return self

    # Get lattice information from CIF block
    def getFromCIF(self, cifblock=None):
        # Get space group number
        try:
            self.spacegroupnr = int(cifblock['_space_group_IT_number'])
        except KeyError:
            try:
                self.spacegroupnr = int(cifblock['_symmetry_Int_Tables_number'])
            except KeyError:
                self.spacegroupnr = -1
        # Get space group symbol
        try:
            self.spacegroupsymbol=cifblock['_space_group_name_H-M_alt'].translate(string.maketrans("",""),string.whitespace)
            self.spacegroupsymboltype='(H-M)'
        except:
            try:
                self.spacegroupsymbol=cifblock['_space_group_name_Hall'].translate(string.maketrans("",""),string.whitespace)
                self.spacegroupsymboltype='(Hall)'
            except:
                try:
                    self.spacegroupsymbol=cifblock['_symmetry_space_group_name_H-M'].translate(string.maketrans("", ""),string.whitespace)
                    self.spacegroupsymboltype='(H-M)'
                except:
                    try:
                        self.spacegroupsymbol=cifblock['_symmetry_space_group_name_Hall'].translate(string.maketrans("",""),string.whitespace)
                        self.spacegroupsymboltype='(Hall)'
                    except:
                        self.spacegroupsymbol = ""
                        self.spacegroupsymboltype = ""
        # Found symbol but not number?
        if self.spacegroupnr == -1 and self.spacegroupsymbol != "":
            if self.spacegroupsymboltype == "(H-M)":
                try:
                    self.spacegroupnr = int(SpaceGroupData().HMtoSGnr[self.spacegroupsymbol])
                except:
                    pass
        # Save spacegroup setting separately
        if self.spacegroupsymbol == "" and self.spacegroupnr != -1:
            self.spacegroupsymbol = SpaceGroupData().HMSymbol[str(self.spacegroupnr)]
            self.spacegroupsymboltype = "(H-M)"
        else:
            pass
        # Get symmetry equivalent positions (space group operations)
        try:
            eqsitedata = cifblock.GetLoop('_symmetry_equiv_pos_as_xyz')
            try:
                eqsitestrs = eqsitedata.get('_symmetry_equiv_pos_as_xyz')
                # funny exception which will occur for the P1 space group
                if type(eqsitestrs) == StringType:
                    eqsitestrs = [eqsitestrs]
                self.symops = []
                for i in range(len(eqsitestrs)):
                    tmp = eqsitestrs[i].split(',')
                    for j in range(len(tmp)):
                        tmp[j] = tmp[j].strip().lower()
                    self.symops.append(SymmetryOperation(tmp))
            except KeyError:
                self.symops = []
        except KeyError:
            self.symops = []
        # Get cell
        try:
            self.a      = float(removeerror(cifblock['_cell_length_a']))
            self.b      = float(removeerror(cifblock['_cell_length_b']))
            self.c      = float(removeerror(cifblock['_cell_length_c']))
            self.alpha  = float(removeerror(cifblock['_cell_angle_alpha']))
            self.beta   = float(removeerror(cifblock['_cell_angle_beta']))
            self.gamma  = float(removeerror(cifblock['_cell_angle_gamma']))
            self.coa = self.c / self.a
            self.boa = self.b / self.a
        except:
            self.a = 0
            self.b = 0
            self.c = 0
            self.alpha = 0
            self.beta = 0
            self.gamma = 0
            raise CellError("Unable to read crystallographic parameters")
        # Get info on atom positions
        try:
            tmpdata = cifblock.GetLoop('_atom_site_fract_x')
        except:
            raise PositionError("Unable to find irreducible positions.")
        # Positions
        try:
            sitexer = tmpdata.get('_atom_site_fract_x')
            siteyer = tmpdata.get('_atom_site_fract_y')
            sitezer = tmpdata.get('_atom_site_fract_z')
            if type(sitexer) == NoneType or type(siteyer) == NoneType or type(sitezer) == NoneType:
                raise PositionError("Irreducible positions not found.")
        except KeyError:
            raise PositionError("Irreducible positions not found.")
        # Element names
        elementsymbs = tmpdata.get('_atom_site_type_symbol')
        if type(elementsymbs) == NoneType:
            elementsymbs = tmpdata.get('_atom_site_label')
            if type(elementsymbs) == NoneType:
                # Fill up with question marks if not found
                print "***Warning: Could not find element names."
                elementsymbs = ["??" for site in sitexer]
        # Find charge state
        try:
            # This is usually encoded in _atom_site_type_symbol, but we go by _atom_type_oxidation_number first.
            tmpdata2 = cifblock.GetLoop('_atom_type_oxidation_number')
            symbs = tmpdata2.get('_atom_type_symbol')
            oxnums = tmpdata2.get('_atom_type_oxidation_number')
            self.chargedict = dict([])
            self.charges = []
            for element in elementsymbs:
                i = 0
                for symb in symbs:
                    if symb == element:
                        self.chargedict[element] = Charge(oxnums[i])
                        self.charges.append(Charge(oxnums[i]))
                    i+=1
        except:
            # Try _atom_site_type_symbol
            try:
                oxnums = [elem.strip(string.letters) for elem in elementsymbs]
                self.charges = []
                for i in range(len(oxnums)):
                    if oxnums[i].strip(string.digits) == '+':
                        self.charges.append(Charge(float(oxnums[i].strip(string.punctuation))))
                    elif oxnums[i].strip(string.digits) == '-':
                        self.charges.append(Charge(-float(oxnums[i].strip(string.punctuation))))
            except:
                self.charges = [Charge(0) for element in elementsymbs]
        # Remove stuff (usually charge state specification) from element symbol strings
        elements = []
        i = 0
        for elem in elementsymbs:
            elements.append(elem.strip(string.punctuation+string.digits))
            # Make it ?? if there was nothing left after removing junk
            if elements[i] == "":
                elements[i] = "??"
            i += 1
        # Make element name start with capital and then have lower case letters
        elements[:] = [element[0].upper()+element[1:].lower() for element in elements]
        for element in elements:
            if not element in ElementData().elementnr:
                print "***Warning: "+element+" is not a chemical element."
        # Find occupancies
        try:
            siteoccer = tmpdata.get('_atom_site_occupancy')
        except KeyError:
            raise PositionError("Error reading site occupancies.")
        if siteoccer == None:
            if not self.quiet:
                print "***Warning : Site occupancies not found, assuming all occupancies = 1."
            siteoccer = []
            for site in elements:
                siteoccer.append("1.0")
        #
        self.ineqsites = []
        self.occupations = []
        for i in range(len(elements)):
            self.ineqsites.append([])
            # Remove error estimation from coordinates
            for j in sitexer[i], siteyer[i], sitezer[i]:
                try:
                    self.ineqsites[i].append(float(removeerror(j)))
                except:
                    raise PositionError("Invalid atomic position value : "+j)
            # Improve precision
            self.ineqsites[i] = LatticeVector(self.ineqsites[i])
            # dictionary of elements and occupancies
            k = elements[i]
            try:
                v = float(removeerror(siteoccer[i]))
            except ValueError:
                v = 1.0
            self.occupations.append({ k : v })

################################################################################################
class ReferenceData:
    """
    Container class for a set of reference data strings.
    Also contains the getFromCIF method to obtain the data from a CIF block.

    database        : database from which the data was obtained
    databasecode    : some identifier used by the database
    databasestring  : a ready formatted string describing the database info.
    compound        : long compound name
    cpd             : short compound name
    authors         : list of strings with the author names
    authorstring    : string with author names formatted as:
                       single author  -  authors name
                       two authors    -  first author and second author
                       > two authors  -  first author et al.
    journal         :\
    volume          : |
    firstpage       :  - self-explanatory
    lastpage        : |
    year            :/
    journalstring   : journal, volume firstpage-lastpage (year)
    referencestring : authorstring, journalstring
    
    Various attempts at identifying the data are made and if they all fail
    an empty string is returned.
    
    """
    def __init__(self):
        self.database = ""
        self.databasecode = ""
        self.databasestring = ""
        self.compound = ""
        self.cpd = ""
        self.authors = []
        self.authorstring = ""
        self.journal = ""
        self.volume = ""
        self.firstpage = ""
        self.lastpage = ""
        self.year = ""
        # Some known databases
        self.databasenames = {
            'CAS'  : 'Chemical Abstracts',
            'CSD'  : 'Cambridge Structural Database',
            'ICSD' : 'Inorganic Crystal Structure Database',
            'MDF'  : 'Metals Data File (metal structures)',
            'NBS'  : '(NIST) Crystal Data Database',
            'PDB'  : 'Protein Data Bank',
            'PDF'  : 'Powder Diffraction File (JCPDS/ICDD)',
            'COD'  : 'Crystallography Open Database'
            }
        self.databaseabbr = {
             'Chemical Abstracts'                   : 'CAS',
             'Cambridge Structural Database'        : 'CSD',
             'Inorganic Crystal Structure Database' : 'ICSD',
             'Metals Data File (metal structures)'  : 'MDF',
             '(NIST) Crystal Data Database'         : 'NBS',
             'Protein Data Bank'                    : 'PDB',
             'Powder Diffraction File (JCPDS/ICDD)' : 'PDF',
             'Crystallography Open Database'        : 'COD'
            }

    def journalstring(self):
        try:
            if not (self.journal == "" and self.volume == "" and self.firstpage=="" and self.year==""):
                journalstring = self.journal
                journalstring += " "+self.volume
                journalstring += ", "+self.firstpage+"-"+self.lastpage
                journalstring += " ("+self.year+")"
            else:
                journalstring = "No journal information"
        except:
            journalstring = "Failed to create journal string"
        return journalstring

    def referencestring(self):
        try:
            if self.authorstring == "":
                referencestring = "No author information. "
            else:
                referencestring = self.authorstring+", "
            referencestring += self.journalstring()
        except:
            referencestring = "Failed to create reference string"
        return referencestring

    def getFromCIF(self, cifblock=None):
        # Get long compound name
        self.compound = cifblock.get('_chemical_name_systematic')
        if type(self.compound) == NoneType:
            self.compound = cifblock.get('_chemical_name_mineral')
            if type(self.compound) == NoneType:
                self.compound = ""
        # Get short compound name
        self.cpd = cifblock.get('_chemical_formula_structural')
        if type(self.cpd) == NoneType:
            self.cpd = cifblock.get('_chemical_formula_sum')
            if type(self.cpd) == NoneType:
                self.cpd = ""
        if type(self.compound) != StringType:
            self.compound = ""
        if type(self.cpd) != StringType:
            self.cpd = ""
        # Try to identify a source database for the CIF
        for db in self.databasenames:
            # First all the standard ones
            try:
                tmp = cifblock.get('_database_code_'+db)
                if type(tmp) != NoneType:
                    self.databasecode = tmp
                    self.database = self.databasenames[db]
                    self.databasestring = "CIF file exported from "+self.database+\
                                     ".\nDatabase reference code: "+self.databasecode+"."
                    break
                else:
                    pass
            except KeyError:
                pass
            try:
                self.databasestring
            except AttributeError:
                self.databasecode = ""
                self.database = ""
                self.databasestring = ""
        # Check if it is a COD file
        if self.databasecode == "":
            try:
                tmp = cifblock.get('_cod_database_code')
                if type(tmp) != NoneType:
                    self.databasecode = tmp
                    self.database = self.databasenames["COD"]
                    self.databasestring = "CIF file exported from "+self.database+\
                                          ".\nDatabase reference code: "+self.databasecode+"."
            except:
                pass
        #
        # Get bibliographic information
        #
        try:
        # Authors
            authorsloop = cifblock.GetLoop('_publ_author_name')
            self.authors = authorsloop.get('_publ_author_name')
            if type(self.authors) == StringType:
                self.authorstring = self.authors
            else:                
                if len(self.authors) == 1:
                    self.authorstring = self.authors[0]
                elif len(self.authors) == 2:
                    self.authorstring = self.authors[0]+" and "+self.authors[1]
                elif len(self.authors) > 2:
                    self.authorstring = self.authors[0]+" et al."
        except KeyError:
            self.authors = []
            self.authorstring = "Failed to get author information"
        # get rid of newline characters
        self.authorstring = self.authorstring.replace("\n","")
        # Journal details
        failed = False
        try:
            # Look for citation block
            references = cifblock.GetLoop('_citation_id')
            # Pick primary reference (or the first in the list, if not found)
            i = 0
            try:
                while references.get('_citation_id')[i] != "primary":
                    i = i + 1
            except IndexError:
                # No primary reference found, using the first one.
                i = 0
            # journal/book title
            if type(references.get('_citation_journal_full')) != NoneType:
                self.journal = references.get('_citation_journal_full')[i]
            else:
                if type(references.get('_citation_journal_abbrev')) != NoneType:
                    self.journal = references.get('_citation_journal_abbrev')[i]
                else:
                    if type(references.get('_citation_book_title')) != NoneType:
                        self.journal = references.get('_citation_book_title')[i]
                    else:
                        self.journal = ""
            # volume
            if type(references.get('_citation_journal_volume')) != NoneType:
                self.volume = references.get('_citation_journal_volume')[i]
            else:
                self.volume = ""
            if type(self.volume) == NoneType:
                self.volume = ""
            # first page
            if type(references.get('_citation_page_first')) != NoneType:
                self.firstpage = references.get('_citation_page_first')[i]
            else:
                self.firstpage = ""
            if type(self.firstpage) == NoneType:
                self.firstpage = ""
            # last page
            if type(references.get('_citation_page_last')) != NoneType:
                self.lastpage = references.get('_citation_page_last')[i]
            else:
                self.lastpage = ""
            if type(self.lastpage) == NoneType:
                self.lastpage = ""
            # year
            if type(references.get('_citation_year')) != NoneType:
                self.year = references.get('_citation_year')[i]
            else:
                self.year = ""
            if type(self.year) == NoneType:
                self.year = ""
        except KeyError:
            try:
                # journal
                self.journal = cifblock.get('_journal_name_full')
                # volume
                self.volume = cifblock.get('_journal_volume')
                # pages
                self.firstpage = cifblock.get('_journal_page_first')
                self.lastpage = cifblock.get('_journal_page_last')
                # year
                self.year = cifblock.get('_journal_year')
            except KeyError:
                failed = True
        if self.journal == None:
            failed = True
            self.journal = ""
        if self.volume == None:
            failed = True
            self.volume = ""
        if self.firstpage == None:
            failed = True
            self.firstpage = ""
        if self.lastpage == None:
            failed = True
            self.lastpage = ""
        if self.year == None:
            failed = True
            self.year = ""
        if not failed:
            # get rid of newline characters
            try:
                self.journal = self.journal.replace("\n","")
            except:
                self.journal = "??????"
            try:
                self.volume = self.volume.replace("\n","")
            except:
                self.volume = "??"
            try:
                self.firstpage = self.firstpage.replace("\n","")
            except:
                self.firstpage = "??"
            try:
                self.lastpage = self.lastpage.replace("\n","")
            except:
                self.lastpage = "??"
            try:
                self.year = self.year.replace("\n","")
            except:
                self.year = "????"

################################################################################################
# Local functions
def removeerror(string):
    # Remove error estimates at the end of a number (as in 3.28(5))
    splitstr=string.split('(')
    return splitstr[0]

# Guess the "true" values of some conspicuous numbers
def improveprecision(x,eps):
    for f in floatlist:
        if abs(x-f) <= eps:                
            # 0
            return f
    # if no match found, return x
    return x

def latvectadd(a,b):
    t = []
    for i in range(3):
        t.append(a[i]+b[i])
        if abs(t[i]) >= 1-occepsilon:
            t[i] = t[i] - copysign(1,t[i])
        t[i] = improveprecision(t[i],occepsilon)
    return t

def putincell(coords,coordepsilon):
    # Put coordinates in the interval 0 <= x < 1
    for i in range(3):
        # first make the coordinate positive
        while coords[i] < 0:
            coords[i] = coords[i] + 1
        # then put it in the primitive cell
        while coords[i] > 1-coordepsilon:
            coords[i] = coords[i] - 1

# Determinant of 3x3 dimensional matrix
def det3(m):
    a = m[1][1]*m[2][2]-m[1][2]*m[2][1]
    b = m[1][2]*m[2][0]-m[1][0]*m[2][2]
    c = m[1][0]*m[2][1]-m[1][1]*m[2][0]
    return m[0][0]*a + m[0][1]*b + m[0][2]*c

# Inverse of 3x3 dimensional matrix
def minv3(m):
    det = det3(m)
    w = [[(m[1][1]*m[2][2]-m[1][2]*m[2][1])/det, (m[0][2]*m[2][1]-m[0][1]*m[2][2])/det, (m[0][1]*m[1][2]-m[0][2]*m[1][1])/det],
         [(m[1][2]*m[2][0]-m[1][0]*m[2][2])/det, (m[0][0]*m[2][2]-m[0][2]*m[2][0])/det, (m[0][2]*m[1][0]-m[0][0]*m[1][2])/det],
         [(m[1][0]*m[2][1]-m[1][1]*m[2][0])/det, (m[0][1]*m[2][0]-m[0][0]*m[2][1])/det, (m[0][0]*m[1][1]-m[0][1]*m[1][0])/det]]
    return w

# matrix-vector multiplication
def mvmult3(mat,vec):
    w = [0.,0.,0.]
    t = 0
    for i in range(3):
        r = mat[i]
        for j in range(3):
            t += r[j]*vec[j]
        w[i],t = t,0
    return w

# matrix-matrix multiplication
def mmmult3(m1,m2):
    w = []
    for i in range(3):
        w.append([])
        for j in range(3):
            t = 0
            for k in range(3):
                t += m1[i][k]*m2[k][j]
            w[i].append(t)
    return w

def crystal_system(spacegroupnr):
    # Determine crystal system
    if 0 < spacegroupnr <= 2:
        return "triclinic"
    elif 2 < spacegroupnr <=15:
        return "monoclinic"
    elif 15 < spacegroupnr <= 74:
        return "orthorhombic"
    elif 75 < spacegroupnr <= 142:
        return "tetragonal"
    elif 142 < spacegroupnr <= 167:
        return "trigonal"
    elif 167 < spacegroupnr <= 194:
        return "hexagonal"
    elif 194 < spacegroupnr <= 230:
        return "cubic"
    else:
        return "unknown"

# Return x with the sign of y
def copysign(x, y):
    if y >= 0:
        return x
    else:
        return -x
