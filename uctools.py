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
#
# TODO:
#   - Make a lot of things inherit from GeometryObject
#   - More output formats (always)
#   - Improve exception handling (probably also always...)
#   - Make the getSuperCell method to also set proper symmetries
#     of the new cell (probably requires using some more advanced
#     symmetry library like computational crystallographics toolbox).
#
# KNOWN BUGS:
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
    
################################################################################################
class GeometryObject:
    """
    Parent class for anything geometrical contains:
    compeps  : epsilon for determining when two floats are equal
    """
    def __init__(self):
        self.compeps = 0.0002

class LatticeVector(list,GeometryObject):
    """
    Floating point vector describing a lattice point. Assumed to have length three.
    * Supports testing for equality using the self.compeps parameter inherited
    from the parent class GeometryObject
    * Supports testing for < using the euclidean norm (length) of the vector
    * Supports addition with another LatticeVector, in which case the
      vectors are added component-wise and then translated back into
      the interval defined by self.interval (default=[0,1[)
    More methods:
    length            : returns the euclidean norm of the vector
    transform(matrix) : returns matrix*vector
    set_interval      : change the definition interval for the lattice vector
    """
    def __init__(self, vec):
        GeometryObject.__init__(self)
        list.__init__(self, vec)
        # Interval we wish to use for the coordinates.
        # In practice either [0,1] or [-.5, 0.5]
        self.interval = [0.0, 1.0]
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
    # Addition of two vectors, putting the result back
    # into the cell if necessary
    def __add__(self, other):
        if self.interval[0] != other.interval[0] or self.interval[1] != other.interval[1]:
            raise GeometryObjectError("LatticeVectors must have the same definition intervals to be added.")
        t = LatticeVector([])
        for i in range(3):
            t.append(self[i]+other[i])
        t.intocell()
        return t
    def __str__(self):
        return "%19.15f %19.15f %19.15f"%(self[0],self[1],self[2])
    # Put the vector components into the cell interval defined by self.interval
    def intocell(self):
        for i in range(3):
            if not (self.interval[0] <= self[i] < self.interval[1]):
                self[i] -= copysign(1,self[i])
    # multiplication by scalar
    def scalmult(self, a):
        for i in range(3):
            self[i] *= a
        return self
    # coordinate transformation
    def transform(self, matrix):
        t = LatticeVector(mvmult3(matrix, self))
        t.intocell()
        return t
    # Change interval
    def change_interval(self, interval):
        self.interval = interval
        t = LatticeVector([0,0,0])
        t.interval = interval
        self = self + t

class LatticeMatrix(list, GeometryObject):
    """
    Three by three matrix consisting of LatticeVector objects
    """
    def __init__(self,mat):
        GeometryObject.__init__(self)
        lmat = [LatticeVector(vec) for vec in mat]
        list.__init__(self, lmat)
    def __str__(self):
        matstr = ""
        for l in self:
            matstr += str(l)+"\n"
        return matstr
    def __eq__(self,other):
        for i in range(3):
            if self[i] != other[i]:
                return False
        return True
    # coordinate transformation
    def transform(self, matrix):
        return LatticeMatrix(mmmult3(matrix, self))

class SymmetryOperation(GeometryObject):
    """
    Class describing a symmetry operation, with a rotation matrix and a translation.
    """
    def __init__(self, eqsite):
        GeometryObject.__init__(self)
        self.eqsite = eqsite
        self.rotation_matrix = self.rotmat()
        self.translation = self.transvec()
    # This way of printing was useful for outputting to CASTEP.
    def __str__(self):
        return str(self.rotation_matrix)+str(self.translation)+"\n"
    # Two symmetry operations are equal if rotation matrices and translation vector
    # differ by at most compeps
    def __eq__(self, other):
        eq = True
        for i in range(3):
            for j in range(3):
                eq = eq and abs(self.rotation_matrix[i][j] - other.rotation_matrix[i][j]) < self.compeps
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
    
################################################################################################
class CellData(GeometryObject):
    """
    Class for a lot of stuff specifying a unit cell.
    latticevectors : The Bravais lattice vectors
    lengthscale    : an overall length scale that multiplies the lattice vectors.
    unit           : the unit of the lengthscale
    alloy          : True if the compound is an alloy
    sitedata       : master array for data about the sites
                      sitedata[site][0] : the positions of the sample sites (lattice coordinates)
                      sitedata[site][1] : dictionary with the atoms occupying this site and
                                          the occupancies.
                      sitedata[site][2] : positions of all atoms of site (lattice coordinates)
                      sitedata[site][3] : species index for when a species occupies more than
                                          one site (as the '2' in Pb2 for the second instance of
                                          a site occupied by a Pb atom)
                     Example: you get the second coordinate of the third lattice site of
                              the second type of site by:
                              sitedata[1][2][2][1]
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
        self.lattrans = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.transvecs = [LatticeVector([zero, zero, zero])]
        self.eqsites = []
        self.symops = []
        self.sitedata = []
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 0
        self.beta = 0
        self.gamma = 0
        self.coa = 1
        self.boa = 1
        self.latticevectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.lengthscale = 1
        self.unit = "angstrom"
        self.sitedata = []
        self.alloy = False
        self.numofineqsites = 0
        self.ineqsites = []
        self.occupations = []

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
        return reclatvect

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
        # Initialize 
        if self.spacegroupsymbol == "":
            if 0 < self.spacegroupnr < 231:
                self.spacegroupsymbol = SpaceGroupData().HMSymbol[str(self.spacegroupnr)]
                self.spacegroupsymboltype = "(H-M)"
        # Check if we know enough:
        # To produce the conventional cell (reduce=False) we don't need the space group
        # symbol or number as long as we have the symmetry operations (equivalent sites).
        if self.spacegroupsetting == "":
            try:
                self.spacegroupsetting = SpaceGroupData().HMSymbol[str(self.spacegroupnr)][0]
            except:
                pass
        if self.a!=0 and self.b!=0 and self.c!=0 and self.alpha!=0 and self.beta!=0 and self.gamma!=0:
            if self.spacegroupnr == -1:
                if len(self.eqsites) >= 1:
                    if reduce == True:
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
        self.sitedata = []
        self.alloy = False
        # REDUCTION TO PRIMITIVE CELL
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
        if reduce:
            if self.spacegroupsetting == 'I':
                # Body centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,half])]
                if self.crystal_system() == 'cubic':
                    self.lattrans = [[-half, half, half],
                                     [half, -half, half],
                                     [half, half, -half]]
                else:
                    self.lattrans = [[one, zero, zero],
                                     [zero, one, zero],
                                     [half, half, half]]
            elif self.spacegroupsetting == 'F':
                # face centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,zero]),
                                  LatticeVector([half,zero,half]),
                                  LatticeVector([zero,half,half])]
                self.lattrans = [[half, half, zero],
                                 [half, zero, half],
                                 [zero, half, half]]
            elif self.spacegroupsetting == 'A':
                # A-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([zero,half,half])]
                self.lattrans = [[one, zero, zero],
                                 [zero, half, -half],
                                 [zero, half, half]]
            elif self.spacegroupsetting == 'B':
                # B-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,zero,half])]
                self.lattrans = [[half, zero, -half],
                                 [zero, one, zero],
                                 [half, zero, half]]
            elif self.spacegroupsetting == 'C':
                # C-centered
                self.transvecs = [LatticeVector([zero,zero,zero]),
                                  LatticeVector([half,half,zero])]
                self.lattrans = [[half, -half, zero],
                                 [half, half, zero],
                                 [zero, zero, one]]
            elif self.spacegroupsetting == 'R':
                if abs(self.gamma - 120) < self.coordepsilon:
                    # rhombohedral from hexagonal setting
                    self.transvecs = [LatticeVector([zero,zero,zero]),
                                      LatticeVector([third, 2*third, 2*third]),
                                      LatticeVector([2*third, third, third])]
                    self.lattrans = [[2*third, third, third],
                                     [-third, third, third],
                                     [-third, -2*third, third]]
                else:
                    self.transvecs = [LatticeVector([zero,zero,zero])]
                    self.lattrans = [[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]
            else:
                self.transvecs = [LatticeVector([zero,zero,zero])]
                self.lattrans = [[1, 0, 0],
                                 [0, 1, 0],
                                 [0, 0, 1]]
            # Find inverse lattice transformation matrix
            invlattrans = minv3(self.lattrans)
            # Transform to primitive cell
            tmp = []
            for i in range(3):
                tmp.append(mvmult3(self.latticevectors,self.lattrans[i]))
            self.latticevectors = tmp
            # Improve precision again...
            for i in range(3):
                for j in range(3):
                    self.latticevectors[i][j] = improveprecision(self.latticevectors[i][j],self.coordepsilon)
        else:
            # If no reduction is to be done
            self.transvecs = [LatticeVector([zero,zero,zero])]
            self.lattrans = [[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]]
                
        # Atomic species and the number of each species. Site occupancies.
        for i in range(len(self.ineqsites)):
            self.sitedata.append([])
            self.sitedata[i].append(self.ineqsites[i])
            # Add a copy of occupations[i] to sitedata[i][1]
            for k,v in self.occupations[i].iteritems():
                self.sitedata[i].append({k:v})
            # Determine if we have an alloy
            for element in self.occupations[i]:
                v = self.occupations[i][element]
                if abs(1-v) > occepsilon:
                    self.alloy = True  
        # Make sites unique with array of occupation info in case of an alloy
        if self.alloy:
            removeindices = []
            # 
            for i in range(len(self.sitedata)):
                for j in range(len(self.sitedata)-1,i,-1):
                    if self.poscomp(self.sitedata[i][0], self.sitedata[j][0]):
                        # Add the dictionary of site j to that of site i and schedule index j for
                        # removal. If there is already an instance of the species on this site,
                        # then add the occupancies (this happens when different valencies has been
                        # recorded)
                        for k in self.sitedata[j][1]:
                            if k in self.sitedata[i][1]:
                                v = self.sitedata[j][1][k] + self.sitedata[i][1][k]
                                self.sitedata[i][1][k] = v
                            else:
                                self.sitedata[i][1][k] = self.sitedata[j][1][k]
                                removeindices.append(j)
            # Remove duplicate elements
            removeindices = list(set(removeindices))
            removeindices.sort(reverse=True)
            for i in removeindices:
                self.sitedata.pop(i)

        # Work out all sites in the cell
        if self.eqsites == []:
            self.eqsites = SpaceGroupData().EquivalentPositions[self.spacegroupnr]
        posexpr = []
        for i in range(len(self.sitedata)):
            posexpr.append([])
            self.sitedata[i].append([])
            for j in range(len(self.eqsites)):
                posexpr[i].append([])
                self.sitedata[i][2].append([])
                # position expression string (x,y,z)
                posexpr[i][j].append(self.eqsites[j][0])
                posexpr[i][j].append(self.eqsites[j][1])
                posexpr[i][j].append(self.eqsites[j][2])
                for k in range(3):
                    # position expression string, replacing x,y,z with numbers
                    posexpr[i][j][k] = posexpr[i][j][k].replace('x',str(self.sitedata[i][0][0]))
                    posexpr[i][j][k] = posexpr[i][j][k].replace('y',str(self.sitedata[i][0][1]))
                    posexpr[i][j][k] = posexpr[i][j][k].replace('z',str(self.sitedata[i][0][2]))
                # positions as decimal point numbers (this only works with __future__ division)
                self.sitedata[i][2][j].append(eval(posexpr[i][j][0]))
                self.sitedata[i][2][j].append(eval(posexpr[i][j][1]))
                self.sitedata[i][2][j].append(eval(posexpr[i][j][2]))
            for j in range(len(self.sitedata[i][2])):
                if reduce:
                    # Transform to primitive cell lattice vectors
                    self.sitedata[i][2][j] = mvmult3(invlattrans, self.sitedata[i][2][j])
                # Make all coordinates go in the interval [0,1)
                putincell(self.sitedata[i][2][j], self.coordepsilon)
            # Remove sites that are the same (up to one of the induced lattice translations)
            removelist = self.duplicates(self.sitedata[i][2])
            for j in removelist:
                self.sitedata[i][2].pop(j)

        # Symmetry operations
        for eqsite in self.eqsites:
            self.symops.append(SymmetryOperation(eqsite))
        # sort the symmetry operations (by length of translation vectors)
        self.symops.sort()
        
        # If we reduce the cell, remove the symmetry operations that are the
        # same up to an induced lattice translation
        if reduce:
            if len(self.transvecs) > 1:
                removelist = []
                for j in range(len(self.symops)):
                    for i in range(j+1,len(self.symops)):
                        for vec in self.transvecs:
                            if self.symops[i].translation+vec == self.symops[j].translation:
                                if self.symops[i].rotation_matrix == self.symops[j].rotation_matrix:
                                    removelist.append(i)
                removelist = list(set(removelist))
                removelist.sort(reverse=True)
                for i in removelist:
                    self.symops.pop(i)
                    
        # work out the number of times that a species appears and store index in sitedata[site][3]
        i = 0
        for site1 in self.sitedata:
            site1.append(1)
            species1 = ""
            occ1 = ""
            # get species string
            for k in site1[1]:
                species1 += k
            if len(site1[1]) > 1:
                # get occupancy string
                for v in site1[1].itervalues():
                    occ1 += str(v)
            # loop over sites up to this one
            for site2 in self.sitedata[:i]:
                species2 = ""
                occ2 = ""
                # get species string
                for k in site2[1]:
                    species2 += k
                # add one if species and occupancies matches
                if species1 == species2:
                    if len(site2[1]) > 1:
                        # get occupancy string
                        for v in site2[1].itervalues():
                            occ2 += str(v)
                    if occ1 == occ2:
                        site1[3] += 1
            i += 1
        # Finally, set flag and return the CrystalStructure in the conventional setting
        self.initialized = True
        return self

    def getSuperCell(self, supercellmap, vacuum, transvec):
        """
        Returns a supercell based on the input supercell map. The cell must
        have been initialized.
        The cell will be padded with some number of original unit cells of vacuum
        by simple rescaling of the lattice vectors and positions. Controlled by 'vacuum'.
        All coordinates will be translated by 'transvec', which is given in units
        of the original lattice vectors.
        """
        if not self.initialized:
            raise CellError("The unit cell must be initialized before a supercell can be generated.")
        if len(supercellmap) != 3 or type(supercellmap[0]) != IntType or type(supercellmap[1]) != IntType \
               or type(supercellmap[2]) != IntType or supercellmap[0] <= 0 or supercellmap[1] <= 0 \
               or supercellmap[2] <= 0:
            raise CellError("The supercell map must be an array of three positive integers.")
        if len(supercellmap) != 3 or vacuum[0] < 0 or vacuum[1] < 0 or vacuum[2] < 0:
            raise CellError("The vacuum padding must be an array of three integers >= 0.")
        struct = copy.deepcopy(self.structure)
        vectors = self.structure.latticevectors
        # original length of lattice vectors
        orglatlen = []
        for vec in self.latticevectors:
            leng = sqrt(vec[0]**2+vec[1]**2+vec[2]**2)
            orglatlen.append(leng)
        # Set up offsets
        offsets = []
        multfac = 0
        for k in range(supercellmap[2]):
            for j in range(supercellmap[1]):
                for i in range(supercellmap[0]):
                    offsets.append([i,j,k])
                    multfac += 1
        # generate additional positions
        for i in range(1,len(offsets)):
            for j in range(len(self.structure.sitedata)):
                for l in range(len(self.structure.sitedata[j][2])):
                    position = [self.structure.sitedata[j][2][l][m] for m in range(3)]
                    for m in range(3):
                        position[m] += offsets[i][m]
                    self.sitedata[j][2].append(position)
        # transform lattice vectors
        for i in range(len(supercellmap)):
            factor = supercellmap[i]
            for j in range(len(vectors[i])):
                self.latticevectors[i][j] = vectors[i][j] * factor
        # Pad with vacuum
        if reduce(lambda x,y: x+y, vacuum) > 0:
            # add the given number of unit cell units along the lattice vectors
            for j in range(len(vacuum)):
                for i in range(len(self.latticevectors[j])):
                    self.latticevectors[j][i] = self.latticevectors[j][i] + vectors[j][i]*vacuum[j]
        # new length of lattice vectors
        newlatlen = []
        for vec in self.latticevectors:
            leng = sqrt((vec[0])**2+(vec[1])**2+(vec[2])**2)
            newlatlen.append(leng)
        # Rescale coordinates
        for j in range(len(self.sitedata)):
            for l in range(len(self.sitedata[j][2])):
                for k in range(3):
                    fac = orglatlen[k]/newlatlen[k]
                    self.sitedata[j][2][l][k] = self.sitedata[j][2][l][k] * fac
        # Move all atoms by transvec
        if reduce(lambda x,y: x+y, transvec) != 0:
            for j in range(len(self.sitedata)):
                for l in range(len(self.sitedata[j][2])):
                    for k in range(3):
                        fac = orglatlen[k]/newlatlen[k]
                        self.sitedata[j][2][l][k] = self.sitedata[j][2][l][k] + fac*transvec[k]
        # Put stuff back in ]-1,1[ interval
        for j in range(len(self.sitedata)):
            for l in range(len(self.sitedata[j][2])):
                for k in range(3):
                    while abs(self.sitedata[j][2][l][k]) >= 1:
                        self.sitedata[j][2][l][k] = self.sitedata[j][2][l][k] - copysign(1,self.sitedata[j][2][l][k])
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
                # funny exception which will occur if there is only one position and there is still a loop
                if type(eqsitestrs) == StringType:
                    eqsitestrs = [eqsitestrs]
                self.eqsites = []
                for i in range(len(eqsitestrs)):
                    self.eqsites.append(eqsitestrs[i].split(','))
                    for j in range(len(self.eqsites[i])):
                        self.eqsites[i][j] = self.eqsites[i][j].strip().lower()
            except KeyError:
                self.eqsites = []
        except KeyError:
            self.eqsites = []
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
        elements = tmpdata.get('_atom_site_type_symbol')
        if type(elements) == NoneType:
            elements = tmpdata.get('_atom_site_label')
            if type(elements) == NoneType:
                # Fill up with question marks if not found
                print "***Warning: Could not find element names."
                elements = []
                elements[:] = ["??" for site in sitexer]
        # Remove junk from strings
        for i in range(len(elements)):
            elements[i] = elements[i].strip(string.punctuation).strip(string.digits).strip(string.punctuation)
            # Make it ?? if there was nothing left after removing junk
            if elements[i] == "":
                elements[i] = "??"
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
            # Remove error estimation from coordinates and store in sitedata[site][0]
            for j in sitexer[i], siteyer[i], sitezer[i]:
                try:
                    self.ineqsites[i].append(float(removeerror(j)))
                except:
                    raise PositionError("Invalid atomic position value : "+j)
            # Improve precision 
            for k in range(3):
                self.ineqsites[i][k] = improveprecision(self.ineqsites[i][k],self.coordepsilon)
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
floatlist = [third, 2*third, half, fourth, one, zero, sqrt(2.0),sixth,5*sixth]
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
    w = []
    for i in range(3):
        t = 0
        for j in range(3):
            t += mat[j][i]*vec[j]
        w.append(t)
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
