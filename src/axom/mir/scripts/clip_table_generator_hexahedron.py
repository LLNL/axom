# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#--------------------------------------------------------------------------------
#  Clip Table Generator - Hexahedron | Marko Sterbentz 7/10/2019
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                             Midpoint Lookup Table
#--------------------------------------------------------------------------------

hexMidpoint = {
    (0,1) : 8,
    (1,0) : 8,
    (1,2) : 9,
    (2,1) : 9,
    (2,3) : 10,
    (3,2) : 10,
    (3,0) : 11,
    (0,3) : 11,

    (0,4) : 16,
    (4,0) : 16,
    (5,1) : 17,
    (1,5) : 17,
    (2,6) : 18,
    (6,2) : 18,
    (3,7) : 19,
    (7,3) : 19,

    (4,5) : 12,
    (5,4) : 12,
    (5,6) : 13,
    (6,5) : 13,
    (6,7) : 14,
    (7,6) : 14,
    (7,4) : 15,
    (4,7) : 15
}

#--------------------------------------------------------------------------------
#                              Utility Functions
#--------------------------------------------------------------------------------

def isBitOn(binaryNum, position, maxNumBits):
    bitMask = 1 << (maxNumBits - 1 - position)
    return (binaryNum & bitMask) == bitMask

#--------------------------------------------------------------------------------

def getOnBitHexIndices( caseNum ):
    numHexBits = 8
    onBits = []
    offBits = []
    for position in range(0, numHexBits):
        if isBitOn(caseNum, position, numHexBits):
            onBits.append( position )
        else:
            offBits.append( position )
    return onBits, offBits

#--------------------------------------------------------------------------------

def getAdjacentVertices( vertex ):
    if vertex == 0:
        return [1,3,4]
    elif vertex == 1:
        return [0,2,5]
    elif vertex == 2:
        return [1,3,6]
    elif vertex == 3:
        return [0,2,7]
    elif vertex == 4:
        return [0,5,7]
    elif vertex == 5:
        return [1,4,6]
    elif vertex == 6:
        return [2,5,7]
    elif vertex == 7:
        return [3,4,6]

#--------------------------------------------------------------------------------

def areAdjacent( v0, v1 ):
    v1AdjacentVertices = getAdjacentVertices( v1 )
    if v0 in v1AdjacentVertices:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def areOnSameFace( v0, v1 ):
    
    v0Neighbors = getAdjacentVertices( v0 )
    v1Neighbors = getAdjacentVertices( v1 )

    # If vertices are adjacent, then they are on the same face
    if areAdjacent( v0, v1 ):
        return True

    # If they are not adjacent, but they share two neighbors, then they are on the same face
    numSharedNeighbors = 0
    for v0n in v0Neighbors:
        if v0n in v1Neighbors:
            numSharedNeighbors += 1
    
    if numSharedNeighbors == 2:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

# Checks if the given list of four vertices are all on the same face.
def fourVertsOnSameFace( vertices ):
    faceOne = [0, 1, 2, 3] # bottom face
    faceTwo = [4, 5, 6, 7] # top face
    faceThree = [1, 2, 5, 6] # right face
    faceFour = [0, 3, 4, 7] # left face
    faceFive = [0, 1, 4, 5] # front face
    faceSix = [2, 3, 6, 7] # back face

    if vertices[0] in faceOne and vertices[1] in faceOne and vertices[2] in faceOne and vertices[3] in faceOne:
        return True
    elif vertices[0] in faceTwo and vertices[1] in faceTwo and vertices[2] in faceTwo and vertices[3] in faceTwo:
        return True
    elif vertices[0] in faceThree and vertices[1] in faceThree and vertices[2] in faceThree and vertices[3] in faceThree:
        return True
    elif vertices[0] in faceFour and vertices[1] in faceFour and vertices[2] in faceFour and vertices[3] in faceFour:
        return True
    elif vertices[0] in faceFive and vertices[1] in faceFive and vertices[2] in faceFive and vertices[3] in faceFive:
        return True
    elif vertices[0] in faceSix and vertices[1] in faceSix and vertices[2] in faceSix and vertices[3] in faceSix:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

# Checks if the given list of three vertices are all on the same face.
def threeVertsOnSameFace( vertices ):
    faceOne = [0, 1, 2, 3] # bottom face
    faceTwo = [4, 5, 6, 7] # top face
    faceThree = [1, 2, 5, 6] # right face
    faceFour = [0, 3, 4, 7] # left face
    faceFive = [0, 1, 4, 5] # front face
    faceSix = [2, 3, 6, 7] # back face

    if vertices[0] in faceOne and vertices[1] in faceOne and vertices[2] in faceOne:
        return True
    elif vertices[0] in faceTwo and vertices[1] in faceTwo and vertices[2] in faceTwo:
        return True
    elif vertices[0] in faceThree and vertices[1] in faceThree and vertices[2] in faceThree:
        return True
    elif vertices[0] in faceFour and vertices[1] in faceFour and vertices[2] in faceFour:
        return True
    elif vertices[0] in faceFive and vertices[1] in faceFive and vertices[2] in faceFive:
        return True
    elif vertices[0] in faceSix and vertices[1] in faceSix and vertices[2] in faceSix:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

# Returns a list of vertex index analogs as if the given vertex is zero and the rest
# of the vertices are remapped accordingly.
# Note: There are 8 cases; one for each vertex of the hexahedron.
def getAnalogsFromVertex( vertexZero ):
    v0Analog = vertexZero
    if vertexZero == 0:
        v1Analog = 1
        v2Analog = 2
        v3Analog = 3
        v4Analog = 4
        v5Analog = 5
        v6Analog = 6
        v7Analog = 7
    elif vertexZero == 1:
        v1Analog = 2
        v2Analog = 3
        v3Analog = 0
        v4Analog = 5
        v5Analog = 6
        v6Analog = 7
        v7Analog = 4
    elif vertexZero == 2:
        v1Analog = 3
        v2Analog = 0
        v3Analog = 1
        v4Analog = 6
        v5Analog = 7
        v6Analog = 4
        v7Analog = 5
    elif vertexZero == 3:
        v1Analog = 0
        v2Analog = 1
        v3Analog = 2
        v4Analog = 7
        v5Analog = 4
        v6Analog = 5
        v7Analog = 6
    elif vertexZero == 4:
        v1Analog = 5
        v2Analog = 1
        v3Analog = 0
        v4Analog = 7
        v5Analog = 6
        v6Analog = 2
        v7Analog = 3
    elif vertexZero == 5:
        v1Analog = 1
        v2Analog = 0
        v3Analog = 4
        v4Analog = 6
        v5Analog = 2
        v6Analog = 3
        v7Analog = 7
    elif vertexZero == 6:
        v1Analog = 2
        v2Analog = 1
        v3Analog = 5
        v4Analog = 7
        v5Analog = 3
        v6Analog = 0
        v7Analog = 4
    elif vertexZero == 7:
        v1Analog = 6
        v2Analog = 5
        v3Analog = 4
        v4Analog = 3
        v5Analog = 2
        v6Analog = 1
        v7Analog = 0

    return [ v0Analog, v1Analog, v2Analog, v3Analog, v4Analog, v5Analog, v6Analog, v7Analog ]

#--------------------------------------------------------------------------------

# Returns a list of vertex index analogs as if the given two vertices are vertex zero
# and four. The rest of the vertices are remapped accordingly.
# Note: There are 12*2 cases; two for each edge of the hexahedron (depending on direction of the edge).
# Note: Reversing the order of the edge vertices simply reverses the analog order.
def getAnalogsFromEdge( vertexZero, vertexFour ):
    v0Analog = vertexZero
    v4Analog = vertexFour

    if (v0Analog, v4Analog) == (0, 4):
        v1Analog = 1; v2Analog = 2; v3Analog = 3; v5Analog = 5; v6Analog = 6; v7Analog = 7
    elif (v0Analog, v4Analog) == (4, 0):
        v1Analog = 7; v2Analog = 6; v3Analog = 5; v5Analog = 3; v6Analog = 2; v7Analog = 1

    elif (v0Analog, v4Analog) == (0, 1):
        v1Analog = 3; v2Analog = 7; v3Analog = 4; v5Analog = 2; v6Analog = 6; v7Analog = 5
    elif (v0Analog, v4Analog) == (1, 0):
        v1Analog = 5; v2Analog = 6; v3Analog = 2; v5Analog = 4; v6Analog = 7; v7Analog = 3

    elif (v0Analog, v4Analog) == (0, 3):
        v1Analog = 4; v2Analog = 5; v3Analog = 1; v5Analog = 7; v6Analog = 6; v7Analog = 2
    elif (v0Analog, v4Analog) == (3, 0):
        v1Analog = 2; v2Analog = 6; v3Analog = 7; v5Analog = 1; v6Analog = 5; v7Analog = 4

    elif (v0Analog, v4Analog) == (1, 2):
        v1Analog = 0; v2Analog = 4; v3Analog = 5; v5Analog = 3; v6Analog = 7; v7Analog = 6
    elif (v0Analog, v4Analog) == (2, 1):
        v1Analog = 6; v2Analog = 7; v3Analog = 3; v5Analog = 5; v6Analog = 4; v7Analog = 0

    elif (v0Analog, v4Analog) == (1, 5):
        v1Analog = 2; v2Analog = 3; v3Analog = 0; v5Analog = 6; v6Analog = 7; v7Analog = 4
    elif (v0Analog, v4Analog) == (5, 1):
        v1Analog = 4; v2Analog = 7; v3Analog = 6; v5Analog = 0; v6Analog = 3; v7Analog = 2

    elif (v0Analog, v4Analog) == (2, 3):
        v1Analog = 1; v2Analog = 5; v3Analog = 6; v5Analog = 0; v6Analog = 4; v7Analog = 7
    elif (v0Analog, v4Analog) == (3, 2):
        v1Analog = 7; v2Analog = 4; v3Analog = 0; v5Analog = 6; v6Analog = 5; v7Analog = 1

    elif (v0Analog, v4Analog) == (2, 6):
        v1Analog = 3; v2Analog = 0; v3Analog = 1; v5Analog = 7; v6Analog = 4; v7Analog = 5
    elif (v0Analog, v4Analog) == (6, 2):
        v1Analog = 5; v2Analog = 4; v3Analog = 7; v5Analog = 1; v6Analog = 0; v7Analog = 3

    elif (v0Analog, v4Analog) == (3, 7):
        v1Analog = 0; v2Analog = 1; v3Analog = 2; v5Analog = 4; v6Analog = 5; v7Analog = 6
    elif (v0Analog, v4Analog) == (7, 3):
        v1Analog = 6; v2Analog = 5; v3Analog = 4; v5Analog = 2; v6Analog = 1; v7Analog = 0

    elif (v0Analog, v4Analog) == (4, 5):
        v1Analog = 0; v2Analog = 3; v3Analog = 7; v5Analog = 1; v6Analog = 2; v7Analog = 6
    elif (v0Analog, v4Analog) == (5, 4):
        v1Analog = 6; v2Analog = 2; v3Analog = 1; v5Analog = 7; v6Analog = 3; v7Analog = 0

    elif (v0Analog, v4Analog) == (4, 7):
        v1Analog = 5; v2Analog = 1; v3Analog = 0; v5Analog = 6; v6Analog = 2; v7Analog = 3
    elif (v0Analog, v4Analog) == (7, 4):
        v1Analog = 3; v2Analog = 2; v3Analog = 6; v5Analog = 0; v6Analog = 1; v7Analog = 5

    elif (v0Analog, v4Analog) == (5, 6):
        v1Analog = 1; v2Analog = 0; v3Analog = 4; v5Analog = 2; v6Analog = 3; v7Analog = 7
    elif (v0Analog, v4Analog) == (6, 5):
        v1Analog = 7; v2Analog = 3; v3Analog = 2; v5Analog = 4; v6Analog = 0; v7Analog = 1

    elif (v0Analog, v4Analog) == (6, 7):
        v1Analog = 2; v2Analog = 1; v3Analog = 5; v5Analog = 3; v6Analog = 0; v7Analog = 4
    elif (v0Analog, v4Analog) == (7, 6):
        v1Analog = 4; v2Analog = 0; v3Analog = 3; v5Analog = 5; v6Analog = 1; v7Analog = 2


    return [ v0Analog, v1Analog, v2Analog, v3Analog, v4Analog, v5Analog, v6Analog, v7Analog ]

#--------------------------------------------------------------------------------

# Returns a list of vertex index analogs as if the given four vertices are the face defined
# by vertices zero, one, two, and three. The rest of the vertices are remapped accordingly.
# Note: There are 6*4 cases; one for each face of the hexahedron (plus 3 more for each rotation of the face).
# Note: The user is responsible for ensuring vertices are properly ordered.
def getAnalogsFromFace( vertexZero, vertexOne, vertexTwo, vertexThree ):
    v0Analog = vertexZero
    v1Analog = vertexOne
    v2Analog = vertexTwo
    v3Analog = vertexThree

    # print("in getAnalogsFromFace(): " + str((v0Analog, v1Analog, v2Analog, v3Analog)) )

    # Face 1 (bottom, with normal facing out)
    if (v0Analog, v1Analog, v2Analog, v3Analog) == ( 0, 1, 2, 3 ):
        v4Analog = 4; v5Analog = 5; v6Analog = 6; v7Analog = 7
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 1, 2, 3, 0 ):
        v4Analog = 5; v5Analog = 6; v6Analog = 7; v7Analog = 4
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 2, 3, 0, 1 ):
        v4Analog = 6; v5Analog = 7; v6Analog = 4; v7Analog = 5
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 3, 0, 1, 2 ):
        v4Analog = 7; v5Analog = 4; v6Analog = 5; v7Analog = 6

    # Face 2 (top, with normal facing out)
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 5, 4, 7, 6 ):
        v4Analog = 1; v5Analog = 0; v6Analog = 3; v7Analog = 2
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 4, 7, 6, 5 ):
        v4Analog = 0; v5Analog = 3; v6Analog = 2; v7Analog = 1
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 7, 6, 5, 4 ):
        v4Analog = 3; v5Analog = 2; v6Analog = 1; v7Analog = 0
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 6, 5, 4, 7 ):
        v4Analog = 2; v5Analog = 1; v6Analog = 0; v7Analog = 3

    # Face 3 (right, with normal facing out)
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 1, 5, 6, 2 ):
        v4Analog = 0; v5Analog = 4; v6Analog = 7; v7Analog = 3
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 5, 6, 2, 1 ):
        v4Analog = 4; v5Analog = 7; v6Analog = 3; v7Analog = 0
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 6, 2, 1, 5 ):
        v4Analog = 7; v5Analog = 3; v6Analog = 0; v7Analog = 4
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 2, 1, 5, 6 ):
        v4Analog = 3; v5Analog = 0; v6Analog = 4; v7Analog = 7

    # Face 4 (left, with normal facing out) 
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 4, 0, 3, 7 ):
        v4Analog = 5; v5Analog = 1; v6Analog = 2; v7Analog = 6
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 0, 3, 7, 4 ):
        v4Analog = 1; v5Analog = 2; v6Analog = 6; v7Analog = 5
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 3, 7, 4, 0 ):
        v4Analog = 2; v5Analog = 6; v6Analog = 5; v7Analog = 1
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 7, 4, 0, 3 ):
        v4Analog = 6; v5Analog = 5; v6Analog = 1; v7Analog = 2

    # Face 5 (front, with normal facing out)
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 4, 5, 1, 0):
        v4Analog = 7; v5Analog = 6; v6Analog = 2; v7Analog = 3
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 5, 1, 0, 4 ):
        v4Analog = 6; v5Analog = 2; v6Analog = 3; v7Analog = 7
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 1, 0, 4, 5) :
        v4Analog = 2; v5Analog = 3; v6Analog = 7; v7Analog = 6
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 0, 4, 5, 1 ):
        v4Analog = 3; v5Analog = 7; v6Analog = 6; v7Analog = 2

    # Face 6 (back, with normal facing out)
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 3, 2, 6, 7):
        v4Analog = 0; v5Analog = 1; v6Analog = 5; v7Analog = 4
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 2, 6, 7, 3 ):
        v4Analog = 1; v5Analog = 5; v6Analog = 4; v7Analog = 0
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 6, 7, 3, 2 ):
        v4Analog = 5; v5Analog = 4; v6Analog = 0; v7Analog = 1
    elif (v0Analog, v1Analog, v2Analog, v3Analog) == ( 7, 3, 2, 6 ):
        v4Analog = 4; v5Analog = 0; v6Analog = 1; v7Analog = 5

    return [ v0Analog, v1Analog, v2Analog, v3Analog, v4Analog, v5Analog, v6Analog, v7Analog ]

#--------------------------------------------------------------------------------

# Wrapper for getAnalogsFromFace() which requires the user to only specify 
# the two vertices that form the diagonal from vertex zero to vertex two of
# the canonical cube.
def getAnalogsFromFaceDiagonal( vertexZero, vertexTwo ):
    # Face 1 (bottom)
    if (vertexZero, vertexTwo) == ( 0, 2 ):
        return getAnalogsFromFace(vertexZero, 1, vertexTwo, 3)
    elif (vertexZero, vertexTwo) == ( 2, 0 ):
        return getAnalogsFromFace(vertexZero, 3, vertexTwo, 1)
    elif (vertexZero, vertexTwo) == ( 1, 3 ):
        return getAnalogsFromFace(vertexZero, 2, vertexTwo, 0)
    elif (vertexZero, vertexTwo) == ( 3, 1 ):
        return getAnalogsFromFace(vertexZero, 0, vertexTwo, 2) 

    # Face 2 (top)
    elif (vertexZero, vertexTwo) == ( 5, 7 ):
        return getAnalogsFromFace(vertexZero, 4, vertexTwo, 6)
    elif (vertexZero, vertexTwo) == ( 7, 5 ):
        return getAnalogsFromFace(vertexZero, 6, vertexTwo, 4)
    elif (vertexZero, vertexTwo) == ( 4, 6 ):
        return getAnalogsFromFace(vertexZero, 7, vertexTwo, 5) 
    elif (vertexZero, vertexTwo) == ( 6, 4 ):
        return getAnalogsFromFace(vertexZero, 5, vertexTwo, 7) 

    # Face 3 (right)
    elif (vertexZero, vertexTwo) == ( 1, 6 ):
        return getAnalogsFromFace(vertexZero, 5, vertexTwo, 2) 
    elif (vertexZero, vertexTwo) == ( 6, 1 ):
        return getAnalogsFromFace(vertexZero, 2, vertexTwo, 5) 
    elif (vertexZero, vertexTwo) == ( 2, 5 ):
        return getAnalogsFromFace(vertexZero, 1, vertexTwo, 6) 
    elif (vertexZero, vertexTwo) == ( 5, 2 ):
        return getAnalogsFromFace(vertexZero, 6, vertexTwo, 1) 

    # Face 4 (left)
    elif (vertexZero, vertexTwo) == ( 0, 7 ):
        return getAnalogsFromFace(vertexZero, 3, vertexTwo, 4) 
    elif (vertexZero, vertexTwo) == ( 7, 0 ):
        return getAnalogsFromFace(vertexZero, 4, vertexTwo, 3)
    elif (vertexZero, vertexTwo) == ( 3, 4 ):
        return getAnalogsFromFace(vertexZero, 7, vertexTwo, 0) 
    elif (vertexZero, vertexTwo) == ( 4, 3 ):
        return getAnalogsFromFace(vertexZero, 0, vertexTwo, 7) 

    # Face 5 (front)
    elif (vertexZero, vertexTwo) == ( 1, 4 ):
        return getAnalogsFromFace(vertexZero, 0, vertexTwo, 5) 
    elif (vertexZero, vertexTwo) == ( 4, 1 ):
        return getAnalogsFromFace(vertexZero, 5, vertexTwo, 0)
    elif (vertexZero, vertexTwo) == ( 0, 5 ):
        return getAnalogsFromFace(vertexZero, 4, vertexTwo, 1) 
    elif (vertexZero, vertexTwo) == ( 5, 0 ):
        return getAnalogsFromFace(vertexZero, 1, vertexTwo, 4)

    # Face 6 (back)
    elif (vertexZero, vertexTwo) == ( 3, 6 ):
        return getAnalogsFromFace(vertexZero, 2, vertexTwo, 7)
    elif (vertexZero, vertexTwo) == ( 6, 3 ):
        return getAnalogsFromFace(vertexZero, 7, vertexTwo, 2) 
    elif (vertexZero, vertexTwo) == ( 2, 7 ):
        return getAnalogsFromFace(vertexZero, 6, vertexTwo, 3) 
    elif (vertexZero, vertexTwo) == ( 7, 2 ):
        return getAnalogsFromFace(vertexZero, 3, vertexTwo, 6) 

#--------------------------------------------------------------------------------
#                              Unique Case One
#                          All bits "on" or "off".
#--------------------------------------------------------------------------------

def isHexCaseOne( onBits, offBits ):
    if len( onBits ) == 8 or len( offBits ) == 8:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseOne(onBits, offBits):
    return [8,0,1,2,3,4,5,6,7,-1]

#--------------------------------------------------------------------------------
#                              Unique Case Two
#                         One of the corners is "on".
#--------------------------------------------------------------------------------

def isHexCaseTwo( onBits, offBits ):
    if len( onBits ) == 1 or len( offBits ) == 1:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseTwo(onBits, offBits):
    return handleHexCaseTough(onBits, offBits)

#--------------------------------------------------------------------------------
#                              Unique Case Three
#                       Two adjacent corners are "on".               
#--------------------------------------------------------------------------------

def isHexCaseThree( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:
        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        if clipBits[1] in getAdjacentVertices( clipBits[0] ):
            return True
        else:
            return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseThree(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    analogs = getAnalogsFromEdge( clipBits[0], clipBits[1] )

    # Add the wedge adjacent to the two clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( analogs[4] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[7] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[5] ) ] )

    # Add the hex adjacent to the wedge
    tableEntry.append( 8 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[5] ) ] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[7] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[7] ) ] )

    # Add the wedge that is not adjacent to the two clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[6] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Four
#                      Two opposite vertices are "on". 
#                   (i.e. they share no adjacent vertices)            
#--------------------------------------------------------------------------------

def isHexCaseFour( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:
        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Ensure the clip bits don't share any adjacent vertices
        v0AdjacentVertices = getAdjacentVertices( clipBits[0] )
        v1AdjacentVertices = getAdjacentVertices( clipBits[1] )

        if clipBits[0] in v1AdjacentVertices:
            return False

        if clipBits[1] in v0AdjacentVertices:
            return False

        for v0a in v0AdjacentVertices:
            for v1a in v1AdjacentVertices:
                if v0a == v1a:
                    return False
        return True

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseFour(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    analogs = getAnalogsFromVertex( clipBits[0] )

    # Add tet adjacent to first clipped vertex
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( analogs[0] )

    # Add wedge adjacent to this first tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[3] )

    # Add the first of the two central pyramids
    tableEntry.append( 5 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[5] )

    # Add the second of the two central pyramids
    tableEntry.append( 5 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[3] )

    # Add the wedge adjacent to the second tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[2] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[7] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[5] ) ] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[5] )

    # Add the tet adjacent to the second clipped vertex
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[2] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[7] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[5] ) ] )
    tableEntry.append( analogs[6] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry
  
#--------------------------------------------------------------------------------
#                              Unique Case Five
#             Two non-adjacent vertices on the same face are "on".            
#--------------------------------------------------------------------------------

def isHexCaseFive( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:
        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Ensure that the two clip vertices are not adjacent
        v0AdjacentVertices = getAdjacentVertices( clipBits[0] )
        if clipBits[1] in v0AdjacentVertices:
            return False

        # Ensure the clip bits share some adjacent vertices
        v1AdjacentVertices = getAdjacentVertices( clipBits[1] )

        for v0a in v0AdjacentVertices:
            for v1a in v1AdjacentVertices:
                if v0a == v1a:
                    return True
        return False

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseFive(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    analogs = getAnalogsFromFaceDiagonal( clipBits[0], clipBits[1] )

    # Add the tet adjacent to the first clip vertex
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( analogs[0] )

    # Add the wedge adjacent to the first tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[3] )

    # Add the middle tet that is on the other side of this wedge
    tableEntry.append( 4 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[5] )

    # Add the middle tet that is on the other side of the other wedge
    tableEntry.append( 4 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[6] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[5] )

    # Add the middle pyramid
    tableEntry.append( 5 )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[6] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[3] )

    # Add the tet adjacent to the second clip vertex
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[6] ) ] )
    tableEntry.append( analogs[2] )

    # Add the wedge adjacent to the second tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[6] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[6] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Six
#            Three vertices are "on", and all three are on one face.        
#--------------------------------------------------------------------------------

def isHexCaseSix( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        if len( onBits ) == 3:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Check if the three "on" bits are on the same face.
        v0 = clipBits[0]
        v1 = clipBits[1]
        v2 = clipBits[2]
        for offVertex in nonClipBits:
            faceVertices = [v0, v1, v2]
            faceVertices.append( offVertex )
            if fourVertsOnSameFace( faceVertices ):
                return True

        return False

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseSix(onBits, offBits):
    return handleHexCaseTough(onBits, offBits)

#--------------------------------------------------------------------------------
#                              Unique Case Seven
#        Three vertices are "on", but only two are adjacent to each other.
#--------------------------------------------------------------------------------

def isHexCaseSeven( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        if len( onBits ) == 3:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        v0 = clipBits[0]
        v1 = clipBits[1]
        v2 = clipBits[2]

        if areAdjacent(v0, v1) and not areAdjacent(v0, v2) and not areAdjacent(v1, v2):
            return True
        elif areAdjacent(v1, v2) and not areAdjacent(v1, v0) and not areAdjacent(v2, v0):
            return True
        elif areAdjacent(v2, v0) and not areAdjacent(v2, v1) and not areAdjacent(v0, v1):
            return True
        else:
            return False

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseSeven(onBits, offBits):
     # Init the table entry list
    tableEntry = []

    if len( onBits ) == 3:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # vertexZero = one of the two vertices adjacent to each other AND on the same face as the non-adjacent vertex
    # vertexTwo = the non-adjacent vertex
    # Note: Doing this should ensure that the other vertex adjacent to vertexZero corresponds to vertexFour in the canonical hex

    v0 = clipBits[0]
    v1 = clipBits[1]
    v2 = clipBits[2]

    if areAdjacent(v0, v1) and not areAdjacent(v0, v2) and not areAdjacent(v1, v2):
        # Vertex two is always the vertex not adjacent to the others
        vertexTwo = v2

        # Determine whether v0 or v1 is on the same face as v2, and set that one as vertexZero
        if areOnSameFace( v0, vertexTwo ):
            vertexZero = v0
        else:
            vertexZero = v1
        
    elif areAdjacent(v1, v2) and not areAdjacent(v1, v0) and not areAdjacent(v2, v0):
        # Vertex two is always the vertex not adjacent to the others
        vertexTwo = v0

        # Determine whether v2 or v1 is on the same face as v0, and set that one as vertexZero
        if areOnSameFace( v1, vertexTwo ):
            vertexZero = v1
        else:
            vertexZero = v2

    elif areAdjacent(v2, v0) and not areAdjacent(v2, v1) and not areAdjacent(v0, v1):
        # Vertex two is always the vertex not adjacent to the others
        vertexTwo = v1

        # Determine whether v0 or v2 is on the same face as v1, and set that one as vertexZero
        if areOnSameFace( v2, vertexTwo):
            vertexZero = v2
        else:
            vertexZero = v0

    analogs = getAnalogsFromFaceDiagonal( vertexZero, vertexTwo )

    # Add the wedge with the two-adjacent clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( analogs[4] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[7] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[5] ) ] )

    # Add the hex adjacent to the wedge
    tableEntry.append( 8 )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[7] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[5] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[7] ) ] )

    # Add the pyramid in the middle
    tableEntry.append( 5 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[6] )

    # Add the wedge adjacent to the tet
    tableEntry.append( 6 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[6] )
    tableEntry.append( analogs[3] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[6] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[3] ) ] )
    

    # Add the tet with the non-adjacent clip bit
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[6] ) ] )
    tableEntry.append( analogs[2] )


    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Eight
#          Three vertices are "on", but none are adjacent to each other.                             
#--------------------------------------------------------------------------------

def isHexCaseEight( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        if len( onBits ) == 3:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        v0 = clipBits[0]
        v1 = clipBits[1]
        v2 = clipBits[2]

        if not areAdjacent(v0, v1) and not areAdjacent(v0, v2) and not areAdjacent(v1, v2):
            return True
        else:
            return False

        # Alternate Method: Find "off" vertex that is adjacent to all three "on" vertices.
        # for offVertex in offBits:
        #     offVertexNeighbors = getAdjacentVertices( offVertex )
        #     if set(offVertexNeighbors).issubset( onBits ):
        #         return True
        
        # return False

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseEight(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 3:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Find the non-clip vertex that is adjacent to all three of the clip bits, set this to vertexZero, and remap the indices using just the single vertexZero
    for vertIndex in range (0, 8):
        if vertIndex in getAdjacentVertices( clipBits[0] ) and vertIndex in getAdjacentVertices( clipBits[1] ) and vertIndex in getAdjacentVertices( clipBits[2] ) :
            vertexZero = vertIndex

    analogs = getAnalogsFromVertex( vertexZero )

    # Add the first tet associated with one of the clip bits
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[1] , analogs[0]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[1] , analogs[2]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[1] , analogs[5]) ] )
    tableEntry.append( analogs[1] )

    # Add the wedge adjacent to the first tet
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[2] )
    tableEntry.append( hexMidpoint[ ( analogs[1] , analogs[0] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[1] , analogs[5] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[1] , analogs[2] ) ] )

    # Add the second tet associated with the a clip bit
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[4] , analogs[0] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4] , analogs[5] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4] , analogs[7] ) ] )
    tableEntry.append( analogs[4] )

    # Add the wedge adjacent to the second tet
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[5] )
    tableEntry.append( hexMidpoint[ ( analogs[4] , analogs[0] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4] , analogs[7] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4] , analogs[5] ) ] )

    # Add the third tet associated with one of the clip bits
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[3] , analogs[0] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[3] , analogs[7] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[3] , analogs[2] ) ] )
    tableEntry.append( analogs[3] )

    # Add the wedge adjacent to the third tet
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[7] )
    tableEntry.append( hexMidpoint[ ( analogs[3] , analogs[0] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[3] , analogs[2] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[3] , analogs[7] ) ] )
    
    # Add the first of the central tets
    tableEntry.append( 4 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[5] )

    # Add the second of the central tets
    tableEntry.append( 4 )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[7] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[6] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Nine
#                  Four vertices of a single face are "on", 
#                  four vertices of the opposite face are "off".                    
#--------------------------------------------------------------------------------

def isHexCaseNine( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        if len( onBits ) == 4:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Check if all four of the clipBits are on the same face
        if fourVertsOnSameFace( clipBits ):
            return True
        else:
            return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseNine(onBits, offBits):
     # Init the table entry list
    tableEntry = []

    if len( onBits ) == 4:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Find the vertices that form the diagonal of the bottom face of the canonical hex
    vertexZero = clipBits[0]

    otherFaceVertices = clipBits
    otherFaceVertices.remove( vertexZero )

    for ofv in otherFaceVertices:
        if not areAdjacent(vertexZero, ofv):
            vertexTwo = ofv

    analogs = getAnalogsFromFaceDiagonal( vertexZero, vertexTwo )

    # Add the hex with the clip bits
    tableEntry.append( 8 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[3] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[4]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[1], analogs[5]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[6]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[3], analogs[7]) ] )

    # Add the hex without the clip bits
    tableEntry.append( 8 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[4]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[1], analogs[5]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[6]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[3], analogs[7]) ] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[6] )
    tableEntry.append( analogs[7] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Ten
#           Four vertices are "on", one of which has no adjacent 
#           vertices that are also "on", AND the other three "on" 
#           vertices are on the same face.
#--------------------------------------------------------------------------------

def isHexCaseTen( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        if len( onBits ) == 4:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # If there is any "on" vertex that only has "off" neighbors, then this is the loneOnVertex
        loneOnVertex = None
        for onVert in clipBits:
            onVertNeighbors = getAdjacentVertices( onVert )
            numOffNeighbors = 0
            for onVertNeighbor in onVertNeighbors:
                if onVertNeighbor in nonClipBits:
                    numOffNeighbors += 1
            if numOffNeighbors == 3:
                loneOnVertex = onVert
    
        # Check if the other three "on" vertices are on the same face
        if not loneOnVertex == None:
            clipBitsCopy = list(clipBits)
            clipBitsCopy.remove( loneOnVertex )
            v0 = clipBitsCopy[0]
            v1 = clipBitsCopy[1]
            v2 = clipBitsCopy[2]
            for offVertex in nonClipBits:
                faceVertices = [v0, v1, v2]
                faceVertices.append( offVertex )
                if fourVertsOnSameFace( faceVertices ):
                    return True
        else:
            return False

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseTen(onBits, offBits):
    return handleHexCaseTough(onBits, offBits)

#--------------------------------------------------------------------------------
#                        Unique Cases Eleven and Twelve
#          Four vertices are "on", and two faces have three "on" vertices.                
#--------------------------------------------------------------------------------

def isHexCaseElevenOrTwelve( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        if len( onBits ) == 4:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # It is case eleven or twelve if there are two faces with three "on" vertices each
        numFacesWithThreeOnVertices = 0
        if threeVertsOnSameFace( [ clipBits[0], clipBits[1], clipBits[2] ] ):
            numFacesWithThreeOnVertices += 1
        if threeVertsOnSameFace( [ clipBits[1], clipBits[2], clipBits[3] ] ):
            numFacesWithThreeOnVertices += 1
        if threeVertsOnSameFace( [ clipBits[2], clipBits[3], clipBits[0] ] ):
            numFacesWithThreeOnVertices += 1
        if threeVertsOnSameFace( [ clipBits[3], clipBits[0], clipBits[1] ] ):
            numFacesWithThreeOnVertices += 1

        if numFacesWithThreeOnVertices == 2:
            return True
        else:
            return False

    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseElevenOrTwelve(onBits, offBits):
    return handleHexCaseTough(onBits, offBits)

#--------------------------------------------------------------------------------
#                              Unique Case Thirteen
#                Four vertices are "on", four vertices are "off".          
#                Each vertex has one neighbor in the same on/off set, 
#                 and two neighbors in the opposite on/off set.
#--------------------------------------------------------------------------------

def isHexCaseThirteen( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        if len( onBits ) == 4:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        for onVert in clipBits:
            onVertNeighbors = getAdjacentVertices( onVert )
            numVertsInSameSet = 0
            numVertsInDifferentSet = 0
            for onVertNeighbor in onVertNeighbors:
                if onVertNeighbor in clipBits:
                    numVertsInSameSet += 1
                elif onVertNeighbor in nonClipBits:
                    numVertsInDifferentSet += 1
            if not numVertsInSameSet == 1 or not numVertsInDifferentSet == 2:
                return False

        for offVert in nonClipBits:
            offVertNeighbors = getAdjacentVertices( offVert )
            numVertsInSameSet = 0
            numVertsInDifferentSet = 0
            for offVertNeighbor in offVertNeighbors:
                if offVertNeighbor in nonClipBits:
                    numVertsInSameSet += 1
                elif offVertNeighbor in clipBits:
                    numVertsInDifferentSet += 1
            if not numVertsInSameSet == 1 or not numVertsInDifferentSet == 2:
                return False
        
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseThirteen(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 4:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Arbitrarily choose one of the clipBits to be vertexZero
    vertexZero = clipBits[0]

    # Find the adjacent vertex to vertexZero that is also in clipBits
    vertexZeroNeighbors = getAdjacentVertices(vertexZero)
    for vertexZeroNeighbor in vertexZeroNeighbors:
        if vertexZeroNeighbor in clipBits:
            vertexFour = vertexZeroNeighbor

    analogs = getAnalogsFromEdge( vertexZero, vertexFour )

    # Add the wedge associated with the first two clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( analogs[4] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[7]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[5]) ] )

    # Add the hex adjacent to this wedge
    tableEntry.append( 8 )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( hexMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[5]) ] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[7] )
    tableEntry.append( hexMidpoint[ ( analogs[4], analogs[7]) ] )

    # Add the hex adjacent to the other wedge
    tableEntry.append( 8 )
    tableEntry.append( analogs[1] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[1]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[3]) ] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[5] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[5]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[7]) ] )
    tableEntry.append( analogs[7] )

    # Add the wedge associated with the other two clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[2] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[1]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2], analogs[3]) ] )
    tableEntry.append( analogs[6] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[5]) ] )
    tableEntry.append( hexMidpoint[ ( analogs[6], analogs[7]) ] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Fourteen
#                There are four "on" vertices, and four "off" vertices.
#                Each "on" vertex's adjacent neighbors are "off", and 
#                each "off" vertex's adjacent neighbors are "on."             
#--------------------------------------------------------------------------------

def isHexCaseFourteen( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        if len( onBits ) == 4:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Not Case Fourteen under two conditions:
            # 1. If there is any clip vertex that has neighbors that are also clip bits, then return False.
            # 2. if there is any clip vertex whose neighbors are not all in nonClipBits, then return False
        for cv in clipBits:
            cvNeighbors = getAdjacentVertices( cv )
            for cvNeighbor in cvNeighbors:
                if cvNeighbor in clipBits or cvNeighbor not in nonClipBits:
                    return False

        return True
    else:
        return False

#--------------------------------------------------------------------------------
 
def handleHexCaseFourteen(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 4:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Arbitrarily choose one of the clipBits to align to vertexZero
    vertexZero = clipBits[0]
    analogs = getAnalogsFromVertex( vertexZero )

    # Add tet associated with first clipBit
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[0] , analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0] , analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0] , analogs[3] ) ] )
    tableEntry.append( analogs[0] )

    # Add wedge adjacent to this first tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[0] , analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0] , analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[0] , analogs[3] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[3] )

    # Add tet associated with second clipBit
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[5] , analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[5] , analogs[6] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[5] , analogs[4] ) ] )
    tableEntry.append( analogs[5] )

    # Add wedge adjacent to this second tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[5] , analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[5] , analogs[6] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[5] , analogs[4] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[6] )
    tableEntry.append( analogs[4] )

    # Add tet associated with third clipBit
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[7] , analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[7] , analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[7] , analogs[6] ) ] )
    tableEntry.append( analogs[7] )

    # Add wedge adjacent to this third tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[7] , analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[7] , analogs[4] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[7] , analogs[6] ) ] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[6] )

    # Add tet associated with fourth clipBit
    tableEntry.append( 4 )
    tableEntry.append( hexMidpoint[ ( analogs[2] , analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2] , analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2] , analogs[6] ) ] )
    tableEntry.append( analogs[2] )

    # Add wedge adjacent to this fourth tet
    tableEntry.append( 6 )
    tableEntry.append( hexMidpoint[ ( analogs[2] , analogs[1] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2] , analogs[3] ) ] )
    tableEntry.append( hexMidpoint[ ( analogs[2] , analogs[6] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[6] )

    # Add the central tet
    tableEntry.append( 4 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[6] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Fifteen
#   There are four "on" vertices, one of which is adjacent to the other three.              
#--------------------------------------------------------------------------------

def isHexCaseFifteen( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        if len( onBits ) == 4:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Is Case Fifteen
        for cv in clipBits:
            cvNeighbors = getAdjacentVertices( cv )
            # Check if all of the current clip vertex's neighbors are also in clipBits
            if set(cvNeighbors).issubset(clipBits):
                return True 
        return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleHexCaseFifteen(onBits, offBits):
    return handleHexCaseTough(onBits, offBits)

#--------------------------------------------------------------------------------
#                                  Tough Case
#         There is no way to decompose cleanly using 2D clipping analogs.
#--------------------------------------------------------------------------------

def handleHexCaseTough(onBits, offBits):
    return [5,0,4,5,1,20,
            5,1,5,6,2,20,
            5,2,6,7,3,20,
            5,3,7,4,0,20,
            5,4,7,6,5,20,
            5,1,2,3,0,20,-1]

#--------------------------------------------------------------------------------
#                              Wrapper Functions
#--------------------------------------------------------------------------------

def determineHexCase( caseNum ):
    onBits, offBits = getOnBitHexIndices( caseNum )

    if isHexCaseOne( onBits, offBits ):
        return 1
    if isHexCaseTwo( onBits, offBits ):
        return 2
    if isHexCaseThree( onBits, offBits ):
        return 3
    if isHexCaseFour( onBits, offBits ):
        return 4
    if isHexCaseFive( onBits, offBits ):
        return 5
    if isHexCaseSix( onBits, offBits ):
        return 6
    if isHexCaseSeven( onBits, offBits ):
        return 7
    if isHexCaseEight( onBits, offBits ):
        return 8
    if isHexCaseNine( onBits, offBits ):
        return 9
    if isHexCaseTen( onBits, offBits ):
        return 10
    if isHexCaseElevenOrTwelve( onBits, offBits ):
        return 11
    if isHexCaseThirteen( onBits, offBits ):
        return 13
    if isHexCaseFourteen( onBits, offBits ):
        return 14
    if isHexCaseFifteen( onBits, offBits ):
        return 15

    # No case found. This shouldn't happen.
    return 16

#--------------------------------------------------------------------------------

def handleHexClipping( clippingCase, onBits, offBits ):
    if clippingCase == 1:
        return handleHexCaseOne( onBits, offBits )
    elif clippingCase == 2:
        return handleHexCaseTwo( onBits, offBits )
    elif clippingCase == 3:
        return handleHexCaseThree( onBits, offBits )
    elif clippingCase == 4:
        return handleHexCaseFour( onBits, offBits )
    elif clippingCase == 5:
        return handleHexCaseFive( onBits, offBits )
    elif clippingCase == 6:
        return handleHexCaseSix( onBits, offBits )
    elif clippingCase == 7:
        return handleHexCaseSeven( onBits, offBits )
    elif clippingCase == 8:
        return handleHexCaseEight( onBits, offBits )
    elif clippingCase == 9:
        return handleHexCaseNine( onBits, offBits )
    elif clippingCase == 10:
        return handleHexCaseTen( onBits, offBits )
    elif clippingCase == 11:
        return handleHexCaseElevenOrTwelve( onBits, offBits )
    elif clippingCase == 13:
        return handleHexCaseThirteen( onBits, offBits )
    elif clippingCase == 14:
        return handleHexCaseFourteen( onBits, offBits )
    elif clippingCase == 15:
        return handleHexCaseFifteen( onBits, offBits )
    else:
        return handleHexCaseTough( onBits, offBits )

#--------------------------------------------------------------------------------
#                                     Main
#--------------------------------------------------------------------------------

def main():
    # Determine which vertices correspond to the ones in the clip case, re-assign the vertex IDs based on this, and then handle the clipping case normally
    clipTable = []
    numCases = 2**8

    for clipTableIndex in range(0, numCases):
        # Determine which clipping case it currently is
        clippingCase = determineHexCase( clipTableIndex )

        # Determine which bits are "on" or "off" for the current case
        onBits, offBits = getOnBitHexIndices( clipTableIndex )

        # Call the appropriate function to generate the clipping table entry
        clipTableEntry = handleHexClipping( clippingCase, onBits, offBits )

        # Append the return clipping table entry to the clipTable list 
        clipTable.append( clipTableEntry )

    # Print the table to the console
    for entry in clipTable:
        print( entry )

    # Write out the table to a file
    file = open("hexahedron_clipping_table.txt", "w")
    file.write("{\n")
    for entry in clipTable:
        file.write("{")
        file.write(",".join(str(e) for e in entry))
        file.write("},\n")
    file.write("};")
    file.close()
 
#--------------------------------------------------------------------------------
    
if __name__ == "__main__":
    main()

#--------------------------------------------------------------------------------