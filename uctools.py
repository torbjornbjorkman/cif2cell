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
#  Revision history:
#      2010-07-16 Torbjorn Bjorkman
#        - First version
#      2010-07-18 Torbjorn Bjorkman
#        - Added output for kgrn, kfcd and ncol.
#        - Added some comments.
#      2010-08-20 Torbjorn Bjorkman
#        - Added datastructures for inequivalent positions and
#          space group symbols (Hermann-Mauguin) in the spacegroupdata file.
#          The cell generation thus no longer require these to be in
#          the CIF file.
#      2010-08-20 Torbjorn Bjorkman
#        - Added output for ABINIT.
#      2010-08-22 Torbjorn Bjorkman
#        - Added output for Siesta.
#      2010-08-23 Torbjorn Bjorkman
#        - Added output for CPMD.
#      2010-08-23 Torbjorn Bjorkman
#        - Added output for CASTEP.
#      2010-08-26 Torbjorn Bjorkman
#        - Added output for Fleur.
#      2010-08-27 Torbjorn Bjorkman
#        - Added output for Crystal09.
#      2010-12-11 Torbjorn Bjorkman
#        - Added getSuperCell method to CellData.
#      2011-01-14 Torbjorn Bjorkman
#        - Added vacuum padding to getSuperCell method.
#      2011-01-19 Torbjorn Bjorkman
#        - Added translation vector to getSuperCell method.
#      2011-02-04 Torbjorn Bjorkman
#        - Added new format SymtFile2 class
#
# TODO:
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
from math import sin,cos,pi,sqrt,copysign
from spacegroupdata import *

################################################################################################
# Miscellaneous
zero = 0.0
one = 1.0
two = 2.0
three = 3.0
four = 4.0
third = 1/3
half = one/two
fourth = one/four
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
    
################################################################################################
class CrystalStructure:
    """
    Class for the minimal set of data for outputting the geometrical information
    to an electronic structure program.

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
    """
    def __init__(self):
        self.latticevectors = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.lengthscale = 1
        self.unit = "angstrom"
        self.sitedata = []
        self.alloy = False
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
        
################################################################################################
class CellData:
    """
    Class for a lot of stuff specifying a unit cell.
    Methods:
        getFromCIF          : obtain data for setting up a cell from a CIF block
        crystal_system      : return a string with the name of the crystal system
        latticevectors      : Return the Bravais lattice vectors as a 3x3 matrix
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
        self.initialized = False
        self.quiet = False
        self.coordepsilon = 0.0002
        self.spacegroupnr = 0
        self.spacegroupsymbol = ""
        self.spacegroupsymboltype = ""
        self.spacegroupsetting = ""
        self.lattrans = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        self.transvecs = [[zero, zero, zero]]
        self.a = 0
        self.b = 0
        self.c = 0
        self.alpha = 0
        self.beta = 0
        self.gamma = 0
        self.coa = 1
        self.boa = 1
        self.eqsites = []
        self.alloy = False
        self.numofineqsites = 0
        self.ineqsites = []
        self.occupations = []
        self.structure = CrystalStructure()

    def crystal_system(self):
        return crystal_system(self.spacegroupnr)

    def latticevectors(self):
        # Set up Bravais lattice vectors of the conventional cell
        self.coa = self.c / self.a
        self.boa = self.b / self.a
        alphar = self.alpha*pi/180
        betar  = self.beta*pi/180
        gammar = self.gamma*pi/180
        if self.crystal_system() == 'cubic':
            latticevectors = [[one, zero, zero], 
                              [zero, one, zero], 
                              [zero, zero, one]]
        elif self.crystal_system() == 'hexagonal':
            latticevectors = [[sin(gammar), cos(gammar), zero],
                              [zero, one, zero],
                              [zero, zero, self.coa]]
        elif self.crystal_system() == 'tetragonal':
            latticevectors = [[one, zero, zero],
                              [zero, one, zero],
                              [zero, zero, self.coa]]
        elif self.crystal_system() == 'orthorhombic':
            latticevectors = [[one, zero, zero], 
                              [zero, self.boa, zero], 
                              [zero, zero, self.coa]]
        elif self.crystal_system() == 'monoclinic':
            latticevectors = [[one, zero, zero], 
                              [zero, self.boa, zero], 
                              [self.coa*cos(betar), zero, self.coa*sin(betar)]]
        elif self.crystal_system() == 'trigonal':
            # Hexagonal cell normally used
            if abs(self.gamma-120) < self.coordepsilon:
                latticevectors = [[sin(gammar), cos(gammar), zero],
                                  [zero, one, zero],
                                  [zero, zero, self.coa]]
            else:
                # Symmetric rhombohedral cell (stretching a cube along the main diagonal)
                c = sqrt((1 + 2*cos(alphar))/(1 - cos(alphar)))
                a = pow(1/c, 1/3)
                t = a * (c + 2) / 3
                u = a * (c - 1) / 3
                latticevectors = [[t, u, u],
                                  [u, t, u],
                                  [u, u, t]]
        elif self.crystal_system() == 'triclinic' or self.crystal_system() == 'default':
            angfac1 = (cos(alphar) - cos(betar)*cos(gammar))/sin(gammar)
            angfac2 = sqrt(sin(gammar)**2 - cos(betar)**2 - cos(alphar)**2 
                       + 2*cos(alphar)*cos(betar)*cos(gammar))/sin(gammar)
            latticevectors = [[one, zero, zero], 
                              [self.boa*cos(gammar), self.boa*sin(gammar), zero], 
                              [self.coa*cos(betar), self.coa*angfac1, self.coa*angfac2]]
        else:
            raise SymmetryError("No support for "+self.crystal_system()+" crystal systems.")
        return latticevectors

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
        return det3(self.latticevectors())
    
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
        struct = self.structure
        if self.spacegroupsymbol == "":
            if 0 < self.spacegroupnr < 231:
                self.spacegroupsymbol = SpaceGroupData().HMSymbol[str(self.spacegroupnr)]
                self.spacegroupsymboltype = "(H-M)"
        # Check if we know enough:
        # To produce the conventional cell (reduce=False) we don't need the space group
        # symbol or number as long as we have the symmetry operations (equivalent sites).
        if self.a!=0 and self.b!=0 and self.c!=0 and self.alpha!=0 and self.beta!=0 and self.gamma!=0:
            if self.spacegroupnr == -1:
                if len(self.eqsites) >= 1:
                    if reduce == True:
                        if self.spacegroupsetting == 'P':
                            self.spacegroupnr = 0
                        else:
                            raise SymmetryError("Insufficient symmetry information to reduce to primitive cell.")
                    else:
                        self.spacegroupnr = 0
        else:
            raise CellError("No crystallographic parameter may be zero.")
        struct.latticevectors = self.latticevectors()
        struct.lengthscale = self.a
        struct.sitedata = []
        struct.alloy = False
        # For reduction to primitive cell
        if self.spacegroupsetting == "":
            self.spacegroupsetting = SpaceGroupData().HMSymbol[str(self.spacegroupnr)][0]
        if reduce:
            if self.spacegroupsetting == 'I':
                # Body centered
                self.transvecs = [[zero,zero,zero],
                                  [half,half,half]]
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
                self.transvecs = [[zero,zero,zero],
                                  [half,half,zero],
                                  [half,zero,half],
                                  [zero,half,half]]
                self.lattrans = [[half, half, zero],
                                 [half, zero, half],
                                 [zero, half, half]]
            elif self.spacegroupsetting == 'A':
                # A-centered
                self.transvecs = [[zero,zero,zero],
                                  [zero,half,half]]
                self.lattrans = [[one, zero, zero],
                                 [zero, half, -half],
                                 [zero, half, half]]
            elif self.spacegroupsetting == 'B':
                # B-centered
                self.transvecs = [[zero,zero,zero],
                                  [half,zero,half]]
                self.lattrans = [[half, zero, -half],
                                 [zero, one, zero],
                                 [half, zero, half]]
            elif self.spacegroupsetting == 'C':
                # C-centered
                self.transvecs = [[zero,zero,zero],
                                  [half,half,zero]]
                self.lattrans = [[half, -half, zero],
                                 [half, half, zero],
                                 [zero, zero, one]]
            elif self.spacegroupsetting == 'R':
                if abs(self.gamma - 120) < self.coordepsilon:
                    # rhombohedral from hexagonal setting
                    self.transvecs = [[zero,zero,zero],
                                      [third, 2*third, 2*third],
                                      [2*third, third, third]]
                    self.lattrans = [[2*third, third, third],
                                     [-third, third, third],
                                     [-third, -2*third, third]]
                else:
                    self.transvecs = [[zero,zero,zero]]
                    self.lattrans = [[1, 0, 0],
                                     [0, 1, 0],
                                     [0, 0, 1]]
            else:
                self.transvecs = [[zero,zero,zero]]
                self.lattrans = [[1, 0, 0],
                                 [0, 1, 0],
                                 [0, 0, 1]]
            # Find inverse lattice transformation matrix
            invlattrans = minv3(self.lattrans)
            # Transform to primitive cell
            tmp = []
            for i in range(3):
                tmp.append(mvmult3(struct.latticevectors,self.lattrans[i]))
            struct.latticevectors = tmp
            # Improve precision again...
            for i in range(3):
                for j in range(3):
                    struct.latticevectors[i][j] = improveprecision(struct.latticevectors[i][j],self.coordepsilon)
        else:
            # If no reduction is to be done
            self.transvecs = [[zero,zero,zero]]
            self.lattrans = [[1, 0, 0],
                             [0, 1, 0],
                             [0, 0, 1]]
                
        # Atomic species and the number of each species. Site occupancies.
        for i in range(len(self.ineqsites)):
            struct.sitedata.append([])
            struct.sitedata[i].append(self.ineqsites[i])
            # Add a copy of occupations[i] to sitedata[i][1]
            for k,v in self.occupations[i].iteritems():
                struct.sitedata[i].append({k:v})
            # Determine if we have an alloy
            for element in self.occupations[i]:
                v = self.occupations[i][element]
                if abs(1-v) > occepsilon:
                    struct.alloy = True  
        # Make sites unique with array of occupation info in case of an alloy
        if struct.alloy:
            removeindices = []
            # 
            for i in range(len(struct.sitedata)):
                for j in range(len(struct.sitedata)-1,i,-1):
                    if self.poscomp(struct.sitedata[i][0], struct.sitedata[j][0]):
                        # Add the dictionary of site j to that of site i and schedule index j for
                        # removal. If there is already an instance of the species on this site,
                        # then add the occupancies (this happens when different valencies has been
                        # recorded)
                        for k in struct.sitedata[j][1]:
                            if k in struct.sitedata[i][1]:
                                v = struct.sitedata[j][1][k] + struct.sitedata[i][1][k]
                                struct.sitedata[i][1][k] = v
                            else:
                                struct.sitedata[i][1][k] = struct.sitedata[j][1][k]
                                removeindices.append(j)
            # Remove duplicate elements
            removeindices = list(set(removeindices))
            removeindices.sort(reverse=True)
            for i in removeindices:
                struct.sitedata.pop(i)

        # Work out all sites in the cell
        if self.eqsites == []:
            self.eqsites = SpaceGroupData().EquivalentPositions[self.spacegroupnr]
        posexpr = []
        for i in range(len(struct.sitedata)):
            posexpr.append([])
            struct.sitedata[i].append([])
            for j in range(len(self.eqsites)):
                posexpr[i].append([])
                struct.sitedata[i][2].append([])
                # position expression string (x,y,z)
                posexpr[i][j].append(self.eqsites[j][0])
                posexpr[i][j].append(self.eqsites[j][1])
                posexpr[i][j].append(self.eqsites[j][2])
                for k in range(3):
                    # position expression string, replacing x,y,z with numbers
                    posexpr[i][j][k] = posexpr[i][j][k].replace('x',str(struct.sitedata[i][0][0]))
                    posexpr[i][j][k] = posexpr[i][j][k].replace('y',str(struct.sitedata[i][0][1]))
                    posexpr[i][j][k] = posexpr[i][j][k].replace('z',str(struct.sitedata[i][0][2]))
                # positions as decimal point numbers (this only works with __future__ division)
                struct.sitedata[i][2][j].append(eval(posexpr[i][j][0]))
                struct.sitedata[i][2][j].append(eval(posexpr[i][j][1]))
                struct.sitedata[i][2][j].append(eval(posexpr[i][j][2]))
            for j in range(len(struct.sitedata[i][2])):
                if reduce:
                    # Transform to primitive cell lattice vectors
                    struct.sitedata[i][2][j] = mvmult3(invlattrans, struct.sitedata[i][2][j])
                # Make all coordinates go in the interval [0,1)
                putincell(struct.sitedata[i][2][j], self.coordepsilon)
            # Remove sites that are the same (up to one of the induced lattice translations)
            removelist = self.duplicates(struct.sitedata[i][2])
            for j in removelist:
                struct.sitedata[i][2].pop(j)
                
        # work out the number of times that a species appears and store index in sitedata[site][3]
        i = 0
        for site1 in struct.sitedata:
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
            for site2 in struct.sitedata[:i]:
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
        return struct

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
        for vec in struct.latticevectors:
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
                    struct.sitedata[j][2].append(position)
        # transform lattice vectors
        for i in range(len(supercellmap)):
            factor = supercellmap[i]
            for j in range(len(vectors[i])):
                struct.latticevectors[i][j] = vectors[i][j] * factor
        # Pad with vacuum
        if reduce(lambda x,y: x+y, vacuum) > 0:
            # add the given number of unit cell units along the lattice vectors
            for j in range(len(vacuum)):
                for i in range(len(struct.latticevectors[j])):
                    struct.latticevectors[j][i] = struct.latticevectors[j][i] + vectors[j][i]*vacuum[j]
        # new length of lattice vectors
        newlatlen = []
        for vec in struct.latticevectors:
            leng = sqrt((vec[0])**2+(vec[1])**2+(vec[2])**2)
            newlatlen.append(leng)
        # Rescale coordinates
        for j in range(len(struct.sitedata)):
            for l in range(len(struct.sitedata[j][2])):
                for k in range(3):
                    fac = orglatlen[k]/newlatlen[k]
                    struct.sitedata[j][2][l][k] = struct.sitedata[j][2][l][k] * fac
        # Move all atoms by transvec
        if reduce(lambda x,y: x+y, transvec) != 0:
            for j in range(len(struct.sitedata)):
                for l in range(len(struct.sitedata[j][2])):
                    for k in range(3):
                        fac = orglatlen[k]/newlatlen[k]
                        struct.sitedata[j][2][l][k] = struct.sitedata[j][2][l][k] + fac*transvec[k]
        # Put stuff back in ]-1,1[ interval
        for j in range(len(struct.sitedata)):
            for l in range(len(struct.sitedata[j][2])):
                for k in range(3):
                    while abs(struct.sitedata[j][2][l][k]) >= 1:
                        struct.sitedata[j][2][l][k] = struct.sitedata[j][2][l][k] - copysign(1,struct.sitedata[j][2][l][k])
        return struct

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
        # Define setting
        try:
            self.spacegroupsetting = self.spacegroupsymbol[0]
        except:
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
            if len(self.authors) > 2:
                self.authorstring = self.authors[0]+" et al."
            elif len(self.authors) == 2:
                self.authorstring = self.authors[0]+" and "+self.authors[1]
            elif len(self.authors) == 1:
                self.authorstring = self.authors[0]
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
    and the method FileString that outputs the contents of the .dat file as a string.
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
    def FileString(self):
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
        for site in self.cell.sitedata:
            nosites += len(site[2])
        volume = abs(det3(self.cell.latticevectors))
        wsr = self.cell.lengthscale * (3*volume/(nosites * 4 * pi))**third
        filestring += "SWS......= %9f NP...=  1 SMIX.= 0.500 TMIX.= 0.0000\n" % wsr
        filestring += "Setup: 3 + NQ*NS*NP lines\n"
        filestring += "EFGS.....=    0.0000 EFGS....=   0.00000 FTMAG...=  0.000000\n"
        filestring += "DEO(l)...=     0.020     0.010     0.005     0.001      0.02\n"
        filestring += "Symb IQ IT NL IP NSP   SWP  QTRO  SPLT NFIX NDWF     Eny(spdf)\n"
        iq = 1
        it = 1
        nsp = 1
        if len(self.cell.sitedata[0][1]) > 1:
            prevspecies = "??"
        else:
            for v in self.cell.sitedata[0][1]:
                prevspecies = v
        # type loop
        for site in self.cell.sitedata:
            # Check for alloy
            if len(site[1]) > 1:
                species = "??"
            else:
                for v in site[1]:
                    species = v
            if species != prevspecies:
                prevspecies = species
                nsp += 1
            # site loop
            for pos in site[2]:
                tmpstring = species.ljust(2)+"  "+"%3i%3i"%(iq,it)
                try:
                    tmpstring += "%3i%3i"%(l[ed.elementblock[species]],1)
                except KeyError:
                    tmpstring += "  ?  1"
                tmpstring += "%3i"%nsp
                tmpstring += "    1.000 .000 0.00 0000 1111   .0   .0   .0   .0"
                if len(site[1]) > 1:
                    # print alloy components at the end of the line
                    tmpstring += "       "
                    for comp in site[1]:
                        tmpstring += comp+"/"
                    tmpstring = tmpstring.rstrip("/")
                filestring += tmpstring+"\n"
                iq += 1
            it += 1
        for site in self.cell.sitedata:
            filestring += "Theta....=     90.00 Phia....=      0.00 FIXMOM..=         N moment..=      0.0\n"
        filestring += "PQX......=      0.00 PQY.....=      0.00 PQZ.....=   0.00000 COORD...=L\n"
        filestring += "Atom: 4 lines + NT*6 lines\n"
        filestring += "IEX...=  4  NP..=500 NES..= 15 NITER=250 IWAT.=  0\n"
        filestring += "VMIX.....=  0.300000 RWAT....=  3.500000 RMAX....= 20.000000\n"
        filestring += "DPAS.....=  0.049000 DR1.....=  1.00E-08 TEST....=  1.00E-08\n"
        filestring += "TESTE....=  1.00E-07 TESTY...=  1.00E-08 TESTV...=  1.00E-07\n"
        for site in self.cell.sitedata:
            for comp in site[1]:
                filestring += comp+"\n"
                try:
                    filestring += ed.emtoelements[comp]
                except KeyError:
                    filestring += "\n\n\n\n\n"
        return filestring

class BSTRFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a [filename].dat file for the bstr program
    and the method FileString that outputs the contents of the .dat file as a string.
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
    def FileString(self):
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
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
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
        for site in self.cell.sitedata:
            for k in site[1]:
                l = 1
                if ed.elementblock[k] == 's' or ed.elementblock[k] == 'p':
                    l = max(l,2)
                elif ed.elementblock[k] == 'd':
                    l = max(l,3)
                elif ed.elementblock[k] == 'f':
                    l = max(l,4)
                lstring = " %1i" % l
            for pos in site[2]:
                tmpstring += lstring
                if len(tmpstring) % 69 == 0:
                    tmpstring += "\n          "
        # Need to strip newline character if it the last line was 69 characters long...
        tmpstring = tmpstring.rstrip(string.whitespace)
        tmpstring = tmpstring+"\n"
        filestring += tmpstring
        coa = self.c / self.a
        boa = self.b / self.a
        filestring += "A........=  1.00000000 B.......=  1.00000000 C.......=  1.00000000\n"
        tmpstring = ""
        for i in range(3):
            tmpstring += "BSX......=%12.7f BSY.....=%12.7f BSZ.....=%12.7f\n" % (self.cell.latticevectors[i][0],self.cell.latticevectors[i][1],self.cell.latticevectors[i][2])
        filestring += tmpstring
        for site in self.cell.sitedata:
            for pos in site[2]:
                pos = mvmult3(self.cell.latticevectors,pos)
                tmpstring = "QX.......=%12.7f QY......=%12.7f QZ......=%12.7f" % (pos[0],pos[1],pos[2])
                tmpstring += "      "
                for k in site[1]:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+" "+str(site[3])+"\n"
                filestring += tmpstring
        filestring += "LAMDA....=    2.5000 AMAX....=    5.5000 BMAX....=    5.5000\n"
        return filestring
    
################################################################################################
class CrystalFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed by Crystal and the method
    FileString that outputs the contents of the input file as a string.
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
    def FileString(self):
        filestring = ""
        
        
################################################################################################
# RSPT FILES
class CellgenFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a cellgen.inp file and the method
    FileString that outputs the contents of an cellgen.inp file as a string.
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
    def FileString(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.: "+str(self.cell.lengthscale)+"\n"
        filestring +="# Lattice vectors (columns)\n"
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%self.cell.latticevectors[j][i]
            tmpstring += "\n"
        filestring += tmpstring
        # Get number of sites
        nosites = 0
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
        filestring += "# Sites\n"
        filestring += str(nosites)+"\n"
        for site in self.cell.sitedata:
            for pos in site[2]:
                tmpstring = ""
                for k in range(3):
                    x = improveprecision(pos[k],0.00001)
                    tmpstring += "%19.15f " % x
                if len(site[1]) > 1:
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    for k in site[1]:
                        tmpstring += "%3i" % ed.elementnr[k]
                tmpstring += " l "+chr(site[3]+96)+"   # "
                for k in site[1]:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
                filestring += tmpstring
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
    FileString that outputs the contents of an symt.inp file as a string.
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
    def FileString(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.: "+str(self.cell.lengthscale)+"\n"
        filestring +="# Lattice vectors (columns)\n"
        tmpstring = ""
        for i in range(3):
            for j in range(3):
                tmpstring += "%19.15f "%self.cell.latticevectors[j][i]
            tmpstring += "\n"
        filestring += tmpstring
        filestring += "# Spin axis\n"
        filestring += "%19.15f %19.15f %19.15f  l\n"%(0,0,0)
        # Get number of sites
        nosites = 0
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
        filestring += "# Sites\n"
        filestring += str(nosites)+"\n"
        for site in self.cell.sitedata:
            for pos in site[2]:
                tmpstring = ""
                for k in range(3):
                    x = improveprecision(pos[k],0.00001)
                    tmpstring += "%19.15f " % x
                if len(site[1]) > 1:
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    for k in site[1]:
                        tmpstring += "%3i" % ed.elementnr[k]
                tmpstring += " l "+chr(site[3]+96)+"   # "
                for k in site[1]:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
                filestring += tmpstring
        return filestring

################################################################################################
class SymtFile2(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a new format symt.inp file and the method
    FileString that outputs the contents of an symt.inp file as a string.
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
    def FileString(self):
        # Initialize element data
        ed = ElementData()
        # Add docstring
        filestring = self.docstring
        # Add lattice constant
        filestring += "# Lattice constant in a.u.\n"
        filestring += "lengthscale\n"
        filestring += str(self.cell.lengthscale)+"\n"
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
        filestring += "%19.15f %19.15f %19.15f  l\n"%(0,0,0)
        # Get number of sites
        nosites = 0
        for site in self.cell.sitedata:
            for pos in site[2]:
                nosites += 1
        filestring += "# Sites\n"
        filestring += "atoms\n"
        filestring += str(nosites)+"\n"
        for site in self.cell.sitedata:
            for pos in site[2]:
                tmpstring = ""
                for k in range(3):
                    x = improveprecision(pos[k],0.00001)
                    tmpstring += "%19.15f " % x
                if len(site[1]) > 1:
                    # don't know what to put for an alloy
                    tmpstring += "???"
                else:
                    for k in site[1]:
                        tmpstring += "%3i" % ed.elementnr[k]
                tmpstring += " l "+chr(site[3]+96)+"   # "
                for k in site[1]:
                    tmpstring += k+"/"
                tmpstring = tmpstring.rstrip("/")+"\n"
                filestring += tmpstring
        return filestring

################################################################################################
class Crystal09File(GeometryOutputFile):
    """
    Class for storing the geometrical data needed by Crystal09 and the method
    FileString that outputs the contents of an Crystal09 input file as a string.
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
    def FileString(self):
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
    FileString that outputs the contents of an spacegroup.in file as a string.
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
    def FileString(self):
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
    FileString that outputs the contents of an elk.in file as a string.
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
    def FileString(self):
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
    FileString that outputs the contents of an input.xml file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        self.title = ""
        self.docstring = self.docstring.rstrip("\n")+"\n"
    def FileString(self):
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
    FileString that outputs the contents of an elk/exciting.in file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.cell.newunit("bohr")
        # make sure the docstring goes on one line
        self.docstring = self.docstring.replace("\n"," ")
        if len(self.docstring) > 80:
            self.docstring = self.docstring[0:78]+"...\n"
    def FileString(self):
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
    FileString that outputs to a .cell file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("angstrom")
    def FileString(self):
        # Assign some local variables
        a = self.cell.lengthscale
        lattice = self.cell.latticevectors
        ed = ElementData()
        # docstring
        filestring = self.docstring+"\n"
        filestring += "%BLOCK LATTICE_CART\n"
        # lattice
        for vec in lattice:
            for coord in vec:
                filestring += " %19.15f"%(coord*a)
            filestring += "\n"
        # Cutoff
        filestring += "&ENDBLOCK LATTICE_CART\n\n"
        # The atom position info
        filestring += "&BLOCK POSITIONS_FRAC\n"
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
        filestring += "&ENDBLOCK POSITIONS_FRAC\n"
        return filestring

################################################################################################
class CPMDFile(GeometryOutputFile):
    """
    Class for storing the geometrical data needed in a CPMD run and the method
    FileString that outputs to a .inp file as a string.
    """
    def __init__(self, crystalstructure, string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.cell.newunit("bohr")
        self.cutoff = 100.0
    def FileString(self):
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
        filestring += "%SYSTEM\n"
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
    FileString that outputs to a .fdf file as a string.
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
    def FileString(self):
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
    FileString that outputs the contents of a abinit input file as a string.
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
    def FileString(self):
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
    FileString that outputs the contents of a POSCAR file as a string.
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
    def FileString(self):
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
    and the method FileString that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        # Set atomic units for length scale
        self.jobnam = "default"
        self.kstrjobnam = "default"
        # To be put on the first line
        self.programdoc = ""
    def FileString(self):
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
    and the method FileString that outputs the contents of the .dat file as a string.
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
    def FileString(self):
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
            tmptring += "/tmp"
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
    and the method FileString that outputs the contents of the .dat file as a string.
    """
    def __init__(self,crystalstructure,string):
        GeometryOutputFile.__init__(self,crystalstructure,string)
        self.jobnam = "default"
        # To be put on the first line
        self.programdoc = ""
    def FileString(self):
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
    and the method FileString that outputs the contents of the .dat file as a string.
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
    def FileString(self):
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
    and the method FileString that outputs the contents of the .dat file as a string.
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
    def FileString(self):
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

################################################################################################
class ElementData:
    """
    Class for storing some data about the chemical elements.

    elementnr       :  dictionary of the element numbers.
                       Example: elementnr['O'] is 8.
    elementblock    :  dictionary of which block an element belong to (spdf).
                       Example: elementblock['Fe'] is 'd'
    angularmomentum :  dictionary for the angular momentum quantum number of s,p,d,f states
    emtoelements    :  element setups for EMTO
    """
    def __init__(self):
        # Element numbers
        self.elementnr = {
            'H'  : 1  ,
            'D'  : 1  ,
            'He' : 2  ,
            'Li' : 3  ,
            'Be' : 4  ,
            'B'  : 5  ,
            'C'  : 6  ,
            'N'  : 7  ,
            'O'  : 8  ,
            'F'  : 9  ,
            'Ne' : 10 ,
            'Na' : 11 ,
            'Mg' : 12 ,
            'Al' : 13 ,
            'Si' : 14 ,
            'P'  : 15 ,
            'S'  : 16 ,
            'Cl' : 17 ,
            'Ar' : 18 ,
            'K'  : 19 ,
            'Ca' : 20 ,
            'Sc' : 21 ,
            'Ti' : 22 ,
            'V'  : 23 ,
            'Cr' : 24 ,
            'Mn' : 25 ,
            'Fe' : 26 ,
            'Co' : 27 ,
            'Ni' : 28 ,
            'Cu' : 29 ,
            'Zn' : 30 ,
            'Ga' : 31 ,
            'Ge' : 32 ,
            'As' : 33 ,
            'Se' : 34 ,
            'Br' : 35 ,
            'Kr' : 36 ,
            'Rb' : 37 ,
            'Sr' : 38 ,
            'Y'  : 39 ,
            'Zr' : 40 ,
            'Nb' : 41 ,
            'Mo' : 42 ,
            'Tc' : 43 ,
            'Ru' : 44 ,
            'Rh' : 45 ,
            'Pd' : 46 ,
            'Ag' : 47 ,
            'Cd' : 48 ,
            'In' : 49 ,
            'Sn' : 50 ,
            'Sb' : 51 ,
            'Te' : 52 ,
            'I'  : 53 ,
            'Xe' : 54 ,
            'Cs' : 55 ,
            'Ba' : 56 ,
            'La' : 57 ,
            'Ce' : 58 ,
            'Pr' : 59 ,
            'Nd' : 60 ,
            'Pm' : 61 ,
            'Sm' : 62 ,
            'Eu' : 63 ,
            'Gd' : 64 ,
            'Tb' : 65 ,
            'Dy' : 66 ,
            'Ho' : 67 ,
            'Er' : 68 ,
            'Tm' : 69 ,
            'Yb' : 70 ,
            'Lu' : 71 ,
            'Hf' : 72 ,
            'Ta' : 73 ,
            'W'  : 74 ,
            'Re' : 75 ,
            'Os' : 76 ,
            'Ir' : 77 ,
            'Pt' : 78 ,
            'Au' : 79 ,
            'Hg' : 80 ,
            'Tl' : 81 ,
            'Pb' : 82 ,
            'Bi' : 83 ,
            'Po' : 84 ,
            'At' : 85 ,
            'Rn' : 86 ,
            'Fr' : 87 ,
            'Ra' : 88 ,
            'Ac' : 89 ,
            'Th' : 90 ,
            'Pa' : 91 ,
            'U'  : 92 ,
            'Np' : 93 ,
            'Pu' : 94 ,
            'Am' : 95 ,
            'Cm' : 96 ,
            'Bk' : 97 ,
            'Cf' : 98 ,
            'Es' : 99 ,
            'Fm' : 100,
            'Me' : 101,
            'No' : 102,
            'Lr' : 103,
            'Rf' : 104,
            'Db' : 105,
            'Sg' : 106,
            'Bh' : 107,
            'Hs' : 108,
            'Mt' : 109,
            'Uun': 110,
            'Uuu': 111,
            'Uub': 112}
        
        # Element classification in s, p, d and f blocks
        self.elementblock = { 
            'H'  : 's'  ,
            'D'  : 's'  ,
            'He' : 's'  ,
            'Li' : 's'  ,
            'Be' : 's'  ,
            'B'  : 'p'  ,
            'C'  : 'p'  ,
            'N'  : 'p'  ,
            'O'  : 'p'  ,
            'F'  : 'p'  ,
            'Ne' : 'p' ,
            'Na' : 's' ,
            'Mg' : 's' ,
            'Al' : 'p' ,
            'Si' : 'p' ,
            'P'  : 'p' ,
            'S'  : 'p' ,
            'Cl' : 'p' ,
            'Ar' : 'p' ,
            'K'  : 's' ,
            'Ca' : 's' ,
            'Sc' : 'd' ,
            'Ti' : 'd' ,
            'V'  : 'd' ,
            'Cr' : 'd' ,
            'Mn' : 'd' ,
            'Fe' : 'd' ,
            'Co' : 'd' ,
            'Ni' : 'd' ,
            'Cu' : 'd' ,
            'Zn' : 'd' ,
            'Ga' : 'p' ,
            'Ge' : 'p' ,
            'As' : 'p' ,
            'Se' : 'p' ,
            'Br' : 'p' ,
            'Kr' : 'p' ,
            'Rb' : 's' ,
            'Sr' : 's' ,
            'Y'  : 'd' ,
            'Zr' : 'd' ,
            'Nb' : 'd' ,
            'Mo' : 'd' ,
            'Tc' : 'd' ,
            'Ru' : 'd' ,
            'Rh' : 'd' ,
            'Pd' : 'd' ,
            'Ag' : 'd' ,
            'Cd' : 'd' ,
            'In' : 'p' ,
            'Sn' : 'p' ,
            'Sb' : 'p' ,
            'Te' : 'p' ,
            'I'  : 'p' ,
            'Xe' : 'p' ,
            'Cs' : 's' ,
            'Ba' : 's' ,
            'La' : 'f' ,
            'Ce' : 'f' ,
            'Pr' : 'f' ,
            'Nd' : 'f' ,
            'Pm' : 'f' ,
            'Sm' : 'f' ,
            'Eu' : 'f' ,
            'Gd' : 'f' ,
            'Tb' : 'f' ,
            'Dy' : 'f' ,
            'Ho' : 'f' ,
            'Er' : 'f' ,
            'Tm' : 'f' ,
            'Yb' : 'f' ,
            'Lu' : 'f' ,
            'Hf' : 'd' ,
            'Ta' : 'd' ,
            'W'  : 'd' ,
            'Re' : 'd' ,
            'Os' : 'd' ,
            'Ir' : 'd' ,
            'Pt' : 'd' ,
            'Au' : 'd' ,
            'Hg' : 'd' ,
            'Tl' : 'p' ,
            'Pb' : 'p' ,
            'Bi' : 'p' ,
            'Po' : 'p' ,
            'At' : 'p' ,
            'Rn' : 'p' ,
            'Fr' : 's' ,
            'Ra' : 's' ,
            'Ac' : 'f' ,
            'Th' : 'f' ,
            'Pa' : 'f' ,
            'U'  : 'f' ,
            'Np' : 'f' ,
            'Pu' : 'f' ,
            'Am' : 'f' ,
            'Cm' : 'f' ,
            'Bk' : 'f' ,
            'Cf' : 'f' ,
            'Es' : 'f' ,
            'Fm' : 'f',
            'Me' : 'f',
            'No' : 'f',
            'Lr' : 'd',
            'Rf' : 'd',
            'Db' : 'd',
            'Sg' : 'd',
            'Bh' : 'd',
            'Hs' : 'd',
            'Mt' : 'd',
            'Uun': 'd',
            'Uuu': 'd',
            'Uub': 'd' }

        # Angular momentum quantum numbers
        self.angularmomentum = { 's' : 0, 'p' : 1, 'd' : 2, 'f' : 3 }
        
        # EMTO element configurations (Couldn't think of any way to generate these on the fly...)
        self.emtoelements ={
            "Va":
            "Iz=   0 Norb=  0 Ion=  0 Config= 1s0 \n\
n      1\n\
Kappa -1\n\
Occup  0\n\
Valen  1\n",
            "Em":
            "Iz=   0 Norb=  0 Ion=  0 Config= 1s0 \n\
n      1\n\
Kappa -1\n\
Occup  0\n\
Valen  1\n",
            "H":
            "Iz=   1 Norb=  1 Ion=  0 Config= 1s1\n\
n      1\n\
Kappa -1\n\
Occup  1\n\
Valen  1\n",
            "D":
            "Iz=   1 Norb=  1 Ion=  0 Config= 1s1\n\
n      1\n\
Kappa -1\n\
Occup  1\n\
Valen  1\n",
            "He":
            "Iz=   2 Norb=  1 Ion=  0 Config= 1s2\n\
n      1\n\
Kappa -1\n\
Occup  2\n\
Valen  1\n",
            "Li":
            "Iz=   3 Norb=  2 Ion=  0 Config= 2s1\n\
n      1  2\n\
Kappa -1 -1\n\
Occup  2  1\n\
Valen  0  1\n",
            "Be":
            "Iz=   4 Norb=  2 Ion=  0 Config= 2s2\n\
n      1  2\n\
Kappa -1 -1\n\
Occup  2  2\n\
Valen  0  1\n",
            "B":
            "Iz=   5 Norb=  3 Ion=  0 Config= 2s1_2p1\n\
n      1  2  2\n\
Kappa -1 -1  1\n\
Occup  2  2  1\n\
Valen  0  1  1\n",
            "C":
            "Iz=   6 Norb=  3 Ion=  0 Config= 2s2 2p2\n\
n      1  2  2\n\
Kappa -1 -1  1\n\
Occup  2  2  2\n\
Valen  0  1  1\n",
            "N":
            "Iz=   7 Norb=  4 Ion=  0 Config= 2p3\n\
n      1  2  2  2\n\
Kappa -1 -1  1 -2\n\
Occup  2  2  2  1\n\
Valen  0  0  1  1\n",
            "O":
            "Iz=   8 Norb=  4 Ion=  0 Config= 2s2_2p4\n\
n      1  2  2  2\n\
Kappa -1 -1  1 -2\n\
Occup  2  2  2  2\n\
Valen  0  1  1  1\n",
            "O-2":
            "Iz=   8 Norb=  4 Ion= -2 Config= 2p6\n\
n      1  2  2  2\n\
Kappa -1 -1  1 -2\n\
Occup  2  2  2  4\n\
Valen  0  0  1  1\n",
            "F":
            "Iz=   9 Norb=  4 Ion=  0 Config= 2s2_2p5\n\
n      1  2  2  2\n\
Kappa -1 -1  1 -2\n\
Occup  2  2  2  3\n\
Valen  0  1  1  1\n",
            "Ne":
            "Iz=  10 Norb=  4 Ion=  0 Config= 2s2_2p6\n\
n      1  2  2  2\n\
Kappa -1 -1  1 -2\n\
Occup  2  2  2  4\n\
Valen  0  1  1  1\n",
            "Na":
            "Iz=  11 Norb=  5 Ion=  0 Config= 3s1\n\
n      1  2  2  2  3\n\
Kappa -1 -1  1 -2 -1\n\
Occup  2  2  2  4  1\n\
Valen  0  0  0  0  1\n",
            "Mg":
            "Iz=  12 Norb=  5 Ion=  0 Config= 3s2\n\
n      1  2  2  2  3\n\
Kappa -1 -1  1 -2 -1\n\
Occup  2  2  2  4  2\n\
Valen  0  0  0  0  1\n",
            "Al":
            "Iz=  13 Norb=  6 Ion=  0 Config= 3s2_3p1\n\
n      1  2  2  2  3  3\n\
Kappa -1 -1  1 -2 -1  1\n\
Occup  2  2  2  4  2  1\n\
Valen  0  0  0  0  1  1\n",
            "Si":
            "Iz=  14 Norb=  6 Ion=  0 Config= 3s2_3p2\n\
n      1  2  2  2  3  3\n\
Kappa -1 -1  1 -2 -1  1\n\
Occup  2  2  2  4  2  2\n\
Valen  0  0  0  0  1  1\n",
            "P":
            "Iz=  15 Norb=  7 Ion=  0 Config= 3s2_3p3\n\
n      1  2  2  2  3  3  3\n\
Kappa -1 -1  1 -2 -1  1 -2\n\
Occup  2  2  2  4  2  2  1\n\
Valen  0  0  0  0  1  1  1\n",
            "S":
            "Iz=  16 Norb=  7 Ion=  0 Config= 3s2_3p4\n\
n      1  2  2  2  3  3  3\n\
Kappa -1 -1  1 -2 -1  1 -2\n\
Occup  2  2  2  4  2  2  2\n\
Valen  0  0  0  0  1  1  1\n",
            "Cl":
            "Iz=  17 Norb=  7 Ion=  0 Config= 3s2_3p5\n\
n      1  2  2  2  3  3  3\n\
Kappa -1 -1  1 -2 -1  1 -2\n\
Occup  2  2  2  4  2  2  3\n\
Valen  0  0  0  0  1  1  1\n",
            "Ar":
            "Iz=  18 Norb=  7 Ion=  0 Config= 3s2_3p6\n\
n      1  2  2  2  3  3  3\n\
Kappa -1 -1  1 -2 -1  1 -2\n\
Occup  2  2  2  4  2  2  4\n\
Valen  0  0  0  0  1  1  1\n",
            "K":
            "Iz=  19 Norb=  8 Ion=  0 Config 4s1\n\
n      1  2  2  2  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  1\n\
Valen  0  0  0  0  0  0  0  1\n",
            "Ca":
            "Iz=  20 Norb=  8 Ion=  0 Config= 4s2\n\
n      1  2  2  2  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  2\n\
Valen  0  0  0  0  0  0  0  1\n",
            "Sc":
            "Iz=  21 Norb=  9 Ion=  0 Config= 3d1_4s2\n\
n      1  2  2  2  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  1  1\n",
            "Ti":
            "Iz=  22 Norb=  9 Ion=  0 Config= 3d2_4s2\n\
n      1  2  2  2  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  2  2\n\
Valen  0  0  0  0  0  0  0  1  1\n",
            "V":
            "Iz=  23 Norb=  9 Ion=  0 Config= 3d3_4s2\n\
n      1  2  2  2  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  3  2\n\
Valen  0  0  0  0  0  0  0  1  1\n",
            "Cr":
            "Iz=  24 Norb=  9 Ion=  0 Config= 3d4_4s2\n\
n      1  2  2  2  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  2\n\
Valen  0  0  0  0  0  0  0  1  1\n",
            "Mn":
            "Iz=  25 Norb= 10 Ion=  0 Config= 3d5_4s2\n\
n      1  2  2  2  3  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  1  2\n\
Valen  0  0  0  0  0  0  0  1  1  1\n",
            "Fe":
            "Iz=  26 Norb= 10 Ion=  0 Config= 3d7_4s1\n\
n      1  2  2  2  3  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  3  1\n\
Valen  0  0  0  0  0  0  0  1  1  1\n",
            "Co":
            "Iz=  27 Norb= 10 Ion=  0 Config= 3d7_4s2\n\
n      1  2  2  2  3  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  3  2\n\
Valen  0  0  0  0  0  0  0  1  1  1\n",
            "Ni":
            "Iz=  28 Norb= 10 Ion=  0 Config= 3d8_4s2\n\
n      1  2  2  2  3  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  4  2\n\
Valen  0  0  0  0  0  0  0  1  1  1\n",
            "Cu":
            "Iz=  29 Norb= 10 Ion=  0 Config= 3d10_4s1\n\
n      1  2  2  2  3  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  1\n\
Valen  0  0  0  0  0  0  0  1  1  1\n",
            "Zn":
            "Iz=  30 Norb= 10 Ion=  0 Config= 3d10_4s2\n\
n      1  2  2  2  3  3  3  3  3  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2\n\
Valen  0  0  0  0  0  0  0  1  1  1\n",
            "Ga":
            "Iz=  31 Norb= 11 Ion=  0 Config= 3d10_4s2_4p1\n\
n      1  2  2  2  3  3  3  3  3  4  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1\n\
Occup  2  2  2  4  2  2  4  4  6  2  1\n\
Valen  0  0  0  0  0  0  0  1  1  1  1\n",
            "Ge":
            "Iz=  32 Norb= 11 Ion=  0 Config= 4s2_4p2\n\
n      1  2  2  2  3  3  3  3  3  4  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  1  1\n",
            "As":
            "Iz=  33 Norb= 12 Ion=  0 Config= 3d10_4s2_4p3\n\
n      1  2  2  2  3  3  3  3  3  4  4  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  1\n\
Valen  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Se":
            "Iz=  34 Norb= 12 Ion=  0 Config= 3d10_4s2_4p4\n\
n      1  2  2  2  3  3  3  3  3  4  4  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  2\n\
Valen  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Br":
            "Iz=  35 Norb= 12 Ion=  0 Config= 3d10_4s2_4p5\n\
n      1  2  2  2  3  3  3  3  3  4  4  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  3\n\
Valen  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Kr":
            "Iz=  36 Norb= 12 Ion=  0 Config= 3d10_4s2_4p6\n\
n      1  2  2  2  3  3  3  3  3  4  4  4\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4\n\
Valen  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Rb":
            "Iz=  37 Norb= 13 Ion=  0 Config= 5s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1\n",
            "Sr":
            "Iz=  38 Norb= 13 Ion=  0 Config= 5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1\n",
            "Y":
            "Iz=  39 Norb= 14 Ion=  0 Config= 4d1_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Zr":
            "Iz=  40 Norb= 14 Ion=  0 Config= 4d2_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Nb":
            "Iz=  41 Norb= 14 Ion=  0 Config= 4d3_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  3  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Mo":
            "Iz=  42 Norb= 14 Ion=  0 Config= 4d4_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Tc":
            "Iz=  43 Norb= 15 Ion=  0 Config= 4d5_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Ru":
            "Iz=  44 Norb= 15 Ion=  0 Config= 4d6_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Rh":
            "Iz=  45 Norb= 15 Ion=  0 Config= 4d7_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  3  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Pd":
            "Iz=  46 Norb= 15 Ion=  0 Config= 4d8_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Ag":
            "Iz=  47 Norb= 15 Ion=  0 Config= 4d10_5s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Cd":
            "Iz=  48 Norb= 15 Ion=  0 Config= 4d10_5s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "In":
            "Iz=  49 Norb= 16 Ion=  0 Config= 4d10_5s2_5p1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n",
            "Sn":
            "Iz=  50 Norb= 16 Ion=  0 Config= 4d10_5s2_5p2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n",
            "Sb":
            "Iz=  51 Norb= 17 Ion=  0 Config= 4d10_5s2_5p3\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Te":
            "Iz=  52 Norb= 17 Ion=  0 Config= 4d10_5s2_5p4\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "I":
            "Iz=  53 Norb= 17 Ion=  0 Config= 4d10_5s2_5p5\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  3\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Xe":
            "Iz=  54 Norb= 17 Ion=  0 Config= 4d10_5s2_5p6\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Cs":
            "Iz=  55 Norb= 18 Ion=  0 Config= 6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n",
            "Ba":
            "Iz=  56 Norb= 18 Ion=  0 Config= 6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n",
            "La":
            "Iz=  57 Norb= 19 Ion=  0 Config= 5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Ce":
            "Iz=  58 Norb= 20 Ion=  0 Config= 4f1_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  1  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n",
            "Pr":
            "Iz=  59 Norb= 20 Ion=  0 Config= 4f2_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  2  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Nd":
            "Iz=  60 Norb= 20 Ion=  0 Config= 4f3_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  3  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Pm":
            "Iz=  61 Norb= 20 Ion=  0 Config= 4f4_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  4  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Sm":
            "Iz=  62 Norb= 20 Ion=  0 Config= 4f5_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  5  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Eu":
            "Iz=  63 Norb= 21 Ion=  0 Config= 4f7_5d1_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  1  2  2  4  1  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Gd":
            "Iz=  64 Norb= 21 Ion=  0 Config= 4f7_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  1  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Tb":
            "Iz=  65 Norb= 21 Ion=  0 Config= 4f8_5d2_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  2  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Dy":
            "Iz=  66 Norb= 21 Ion=  0 Config= 4f9_5d1_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  3  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Ho":
            "Iz=  67 Norb= 21 Ion=  0 Config= 4f10_5d1_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  4  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Er":
            "Iz=  68 Norb= 21 Ion=  0 Config= 4f11_5d1_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  5  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Tm":
            "Iz=  69 Norb= 21 Ion=  0 Config= 4f12_5d1_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  6  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Yb":
            "Iz=  70 Norb= 21 Ion=  0 Config= 4f14_5d1_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  1  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Lu":
            "Iz=  71 Norb= 21 Ion=  0 Config= 5d1_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Hf":
            "Iz=  72 Norb= 21 Ion=  0 Config= 5d2_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Ta":
            "Iz=  73 Norb= 21 Ion=  0 Config= 5d3_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  3  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "W":
            "Iz=  74 Norb= 21 Ion=  0 Config= 5d4_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Re":
            "Iz=  75 Norb= 22 Ion=  0 Config= 5d5_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  1  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Os":
            "Iz=  76 Norb= 22 Ion=  0 Config= 5d6_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Ir":
            "Iz=  77 Norb= 22 Ion=  0 Config= 5d7_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  3  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Pt":
            "Iz=  78 Norb= 22 Ion=  0 Config= 5d8_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Au":
            "Iz=  79 Norb= 22 Ion=  0 Config= 5d10_6s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Hg":
            "Iz=  80 Norb= 22 Ion=  0 Config= 5d10_6s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1\n",
            "Tl":
            "Iz=  81 Norb= 23 Ion=  0 Config= 5d10_6s2_6p1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n",
            "Pb":
            "Iz=  82 Norb= 23 Ion=  0 Config= 6s2_6p2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1\n",
            "Bi":
            "Iz=  83 Norb= 24 Ion=  0 Config= 5d10_6s2_6p3\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Po":
            "Iz=  84 Norb= 24 Ion=  0 Config= 5d10_6s2_6p4\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "At":
            "Iz=  85 Norb= 24 Ion=  0 Config= 5d10_6s2_6p5\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  3\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Rn":
            "Iz=  86 Norb= 24 Ion=  0 Config= 5d10_6s2_6p6\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  1  1  1\n",
            "Fr":
            "Iz=  87 Norb= 25 Ion=  0 Config= 7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n",
            "Ra":
            "Iz=  88 Norb= 25 Ion=  0 Config= 7s2\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  2\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1\n",
            "Ac":
            "Iz=  89 Norb= 26 Ion=  0 Config= 6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Th":
            "Iz=  90 Norb= 26 Ion=  0 Config= 6d3_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  4  3  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1\n",
            "Pa":
            "Iz=  91 Norb= 27 Ion=  0 Config= 5f1_6d3_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  1  2  2  4  3  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n",
            "U":
            "Iz=  92 Norb= 27 Ion=  0 Config= 5f2_6d3_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  2  2  2  4  3  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n",
            "Np":
            "Iz=  93 Norb= 27 Ion=  0 Config= 5f4_6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  4  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n",
            "Pu":
            "Iz=  94 Norb= 27 Ion=  0 Config= 5f5_6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  5  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n",
            "Am":
            "Iz=  95 Norb= 27 Ion=  0 Config= 5f6_6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  1  1\n",
            "Cm":
            "Iz=  96 Norb= 28 Ion=  0 Config= 5f7_6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  1  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  1  1\n",
            "Bk":
            "Iz=  97 Norb= 28 Ion=  0 Config= 5f8_6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  2  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  1  1\n",
            "Cf":
            "Iz=  98 Norb= 28 Ion=  0 Config= 5f9_6d2_7s1\n\
n      1  2  2  2  3  3  3  3  3  4  4  4  4  4  4  4  5  5  5  5  5  5  5  6  6  6  6  7\n\
Kappa -1 -1  1 -2 -1  1 -2  2 -3 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -3  3 -4 -1  1 -2  2 -1\n\
Occup  2  2  2  4  2  2  4  4  6  2  2  4  4  6  6  8  2  2  4  4  6  6  3  2  2  4  2  1\n\
Valen  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  1  0  0  0  1  1\n",
            }

################################################################################################
# Local functions
def removeerror(string):
    # Remove error estimates at the end of a number (as in 3.28(5))
    splitstr=string.split('(')
    return splitstr[0]

# Guess the "true" values of some conspicuous numbers
floatlist = [third, 2*third, half, fourth, one, zero, sqrt(2.0)]
def improveprecision(x,eps):
    for f in floatlist:
        if abs(x-f) <= eps:
            # 0
            return f
    # if no match found, return x
    return x

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
            t = t + mat[j][i]*vec[j]
        w.append(t)
    return w

def crystal_system(spacegroupnr):
    # Determine crystal system
    if spacegroupnr == 0:
        return "default"
    elif 0 < spacegroupnr <= 2:
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

