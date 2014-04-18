#!/usr/bin/python
#####################################
# Written by Gavin Heverly-Coulson
# Email: gavin <at> quantumgeranium.com
#####################################
# Converts Cartesian coordinates into fractional coordinates or vice versa
# If input is in Cartesian coordinates, should have .cart extension
# If input is in fractional coordinates, should have .frac extension
# 
# Assumes that all input data is in angstroms
# Requires the use of the coordinates.py library file
####################################
# Input file should be formatted as:
# Number of atoms: <numAtoms>
# 
# <coordType> coordinates:
# X  x/a1  y/b1  z/c1
# X  x/a2  y/b2  z/c2
# ....
#
# Lattice vectors:
# a:   a_x  a_y  a_z
# b:   b_x  b_y  b_z
# c:   c_x  c_y  c_z
#####################################
# Output will be formatted like input with the transformed 
# coordinates substituted for the original
#####################################
#####################################
# This work is licensed under a Simplified BSD License
# Copyright (c) 2014, Gavin Heverly-Coulson
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met: 
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer. 
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution. 
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
# ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import sys

#####################################
# Calculates the determinant of a 3x3 matrix
def det3(mat):
  return ((mat[0][0]*mat[1][1]*mat[2][2]) + (mat[0][1]*mat[1][2]*mat[2][0]) + (mat[0][2]*mat[1][0]*mat[2][1]) - (mat[0][2]*mat[1][1]*mat[2][0]) - (mat[0][1]*mat[1][0]*mat[2][2]) - (mat[0][0]*mat[1][2]*mat[2][1]))

#####################################
# a function that takes the cell parameters, in angstrom, and a list of fractional coordinates
# and returns the structure in Cartesian coordinates
#
# Assumes the cellParam matrix is of the form:
# | ax  ay  az |
# | bx  by  bz |
# | cx  cy  cz |
#
def frac2cart(cellParam, fracCoords):
  cartCoords = []
  for i in fracCoords:
    xPos = i[1]*cellParam[0][0] + i[2]*cellParam[1][0] + i[3]*cellParam[2][0]
    yPos = i[1]*cellParam[0][1] + i[2]*cellParam[1][1] + i[3]*cellParam[2][1]
    zPos = i[1]*cellParam[0][2] + i[2]*cellParam[1][2] + i[3]*cellParam[2][2]
    cartCoords.append([i[0], xPos, yPos, zPos])
  return cartCoords

#####################################
# a function that takes the cell parameters, in angstrom, and a list of Cartesian coordinates
# and returns the structure in fractional coordinates
#
# Uses Cramer's Rule to solve for the fractional coordinates
#
# Assumes the cellParam matrix is of the form:
# | ax  ay  az |
# | bx  by  bz |
# | cx  cy  cz |
#
# Need to use the transpose of this matrix in calculation, so call transpose function first
#
# Assumes cartCoords are of the form:
# | X  x0  y0  z0 |
# | X  x1  y1  z1 |
# | X  x2  y2  z2 |
# | ............. |
# where X is the element symbol
def cart2frac(cellParam, cartCoords):
  latCnt = [x[:] for x in [[None]*3]*3]
  for a in range(3):
    for b in range(3):
	  latCnt[a][b] = cellParam[b][a]

  fracCoords = []
  detLatCnt = det3(latCnt)
  for i in cartCoords:
    aPos = (det3([[i[1], latCnt[0][1], latCnt[0][2]], [i[2], latCnt[1][1], latCnt[1][2]], [i[3], latCnt[2][1], latCnt[2][2]]])) / detLatCnt
    bPos = (det3([[latCnt[0][0], i[1], latCnt[0][2]], [latCnt[1][0], i[2], latCnt[1][2]], [latCnt[2][0], i[3], latCnt[2][2]]])) / detLatCnt
    cPos = (det3([[latCnt[0][0], latCnt[0][1], i[1]], [latCnt[1][0], latCnt[1][1], i[2]], [latCnt[2][0], latCnt[2][1], i[3]]])) / detLatCnt
    fracCoords.append([i[0], aPos, bPos, cPos])
  return fracCoords

#####################################
reader = open(sys.argv[1], 'r')
sourceFile = reader.readlines()
reader.close()

# Determine which direction the transformation needs to happen
# filename[0] will hold the base filename
# filename[1] will hold the file extension
#   cart means the input is in Cartesian coordinates
#   frac means the input is in fractional coordinates
filename = sys.argv[1].split('.')

# Find the starting positions for the three pertinent sections of file
counter = 0
atomNumPos = 0
coordPos = 0
latVectPos = 0
for line in sourceFile:
  if "Number of atoms" in line:
    atomNumPos = counter
  elif "coordinates" in line:
    coordPos = counter + 1
  elif "Lattice vectors" in line:
    latVectPos = counter + 1
  counter += 1

numAtoms = int(sourceFile[atomNumPos].split(":")[-1].strip())

coords = []
for i in sourceFile[coordPos:coordPos+numAtoms]:
  temp = i.split()
  temp[1] = float(temp[1])
  temp[2] = float(temp[2])
  temp[3] = float(temp[3])
  coords.append(temp)

latVectsRaw = sourceFile[latVectPos:latVectPos+3]
latVects = []
temp1 = latVectsRaw[0].split()
latVects.append([float(temp1[1]), float(temp1[2]), float(temp1[3])])
temp2 = latVectsRaw[1].split()
latVects.append([float(temp2[1]), float(temp2[2]), float(temp2[3])])
temp3 = latVectsRaw[2].split()
latVects.append([float(temp3[1]), float(temp3[2]), float(temp3[3])])

# Convert the input coordinates into the other type
if (filename[-1] == "cart"):
  print "Converting to fractional coordinates..."
  newCoords = cart2frac(latVects, coords)
elif (filename[-1] == "frac"):
  print "Converting to Cartesian coordinates..."
  newCoords = frac2cart(latVects, coords) 

# Create the output strings and write them to the file
outputData = []
outputData.append(sourceFile[atomNumPos])
outputData.append('\n')
outputData.append(sourceFile[latVectPos-1])
outputData.append(sourceFile[latVectPos])
outputData.append(sourceFile[latVectPos+1])
outputData.append(sourceFile[latVectPos+2])
outputData.append('\n')
if (filename[-1] == "cart"):
  outputData.append("Fractional coordinates:\n")
elif (filename[-1] == "frac"):
  outputData.append("Cartesian coordinates:\n")
for k in newCoords:
  temp = "{0[0]:<2}   {0[1]:> 10.6f}   {0[2]:> 10.6f}   {0[3]:> 10.6f}\n".format(k)
  outputData.append(temp)
outputData.append('\n')

if (filename[-1] == "cart"):
  outFilename = filename[0] + ".frac"
elif (filename[-1] == "frac"):
  outFilename = filename[0] + ".cart"

writer = open(outFilename, 'w')
writer.write(''.join(outputData))
writer.close()
