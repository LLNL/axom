# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#--------------------------------------------------------------------------------
#  Clip Table Generator - Triangular Prism | Marko Sterbentz 7/5/2019
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                             Midpoint Lookup Table
#--------------------------------------------------------------------------------

wedgeMidpoint = {
    (0,1) : 6,
    (1,0) : 6,
    (1,2) : 7,
    (2,1) : 7,
    (2,0) : 8,
    (0,2) : 8,
    (0,3) : 9,
    (3,0) : 9,
    (1,4) : 10,
    (4,1) : 10,
    (2,5) : 11,
    (5,2) : 11,
    (3,4) : 12,
    (4,3) : 12,
    (4,5) : 13,
    (5,4) : 13,
    (3,5) : 14,
    (5,3) : 14
}

#--------------------------------------------------------------------------------
#                              Utility Functions
#--------------------------------------------------------------------------------

# Returns a list of node indices adjacent to the given node index.
def getAdjacentWedgeNodeIndices(node):
    if node == 0:
        return [1,2,3]
    elif node == 1:
        return [0,2,4]
    elif node == 2:
        return [0,1,5]
    elif node == 3:
        return [0,4,5]
    elif node == 4:
        return [1,3,5]
    elif node == 5:
        return [2,3,4]
    else:
        return [-1]

#--------------------------------------------------------------------------------

# Returns a list of vertex index analogs as if the given vertex is zero
def getAnalogsFromVertex( vertexZero ):
    v0Analog = vertexZero

    if vertexZero == 0:
        v1Analog = 1
        v2Analog = 2
        v3Analog = 3
        v4Analog = 4
        v5Analog = 5
    elif vertexZero == 1:
        v1Analog = 2
        v2Analog = 0
        v3Analog = 4
        v4Analog = 5
        v5Analog = 3
    elif vertexZero == 2:
        v1Analog = 0
        v2Analog = 1
        v3Analog = 5
        v4Analog = 3
        v5Analog = 4
    elif vertexZero == 3:
        v1Analog = 5
        v2Analog = 4
        v3Analog = 0
        v4Analog = 2
        v5Analog = 1
    elif vertexZero == 4:
        v1Analog = 3
        v2Analog = 5
        v3Analog = 1
        v4Analog = 0
        v5Analog = 2
    elif vertexZero == 5:
        v1Analog = 4
        v2Analog = 3
        v3Analog = 2
        v4Analog = 1
        v5Analog = 0

    return [ v0Analog, v1Analog, v2Analog, v3Analog, v4Analog, v5Analog ]

#--------------------------------------------------------------------------------

def areOnSameTriangle( v0, v1 ):
    triangleOne = [0,1,2]
    triangleTwo = [3,4,5]
    
    if v0 in triangleOne and v1 in triangleOne:
        return True
    elif v0 in triangleTwo and v1 in triangleTwo:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def areAdjacent( v0, v1 ):
    adjacentNodes = getAdjacentWedgeNodeIndices( v0 )
    if v1 in adjacentNodes:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def isBitOn(binaryNum, position, maxNumBits):
    bitMask = 1 << (maxNumBits - 1 - position)
    return (binaryNum & bitMask) == bitMask

#--------------------------------------------------------------------------------

def getOnBitWedgeIndices( caseNum ):
    numWedgeBits = 6
    onBits = []
    offBits = []
    for position in range(0, numWedgeBits):
        if isBitOn(caseNum, position, numWedgeBits):
            onBits.append( position )
        else:
            offBits.append( position )
    return onBits, offBits

#--------------------------------------------------------------------------------
#                              Unique Case One
#                             All bits "on" or "off".
#--------------------------------------------------------------------------------

def isWedgeCaseOne( onBits, offBits ):
    if len( onBits ) == 6 or len( offBits ) == 6:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseOne(onBits, offBits):
    return [6,0,1,2,3,4,5,-1]

#--------------------------------------------------------------------------------
#                              Unique Case Two
#                        One of the corners is "on"
#--------------------------------------------------------------------------------

def isWedgeCaseTwo( onBits, offBits ):
    if len( onBits ) == 1 or len( offBits ) == 1:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseTwo(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 1:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    vertexZero = clipBits[0]

    analogs = getAnalogsFromVertex( vertexZero ) 

    # Add the tetrahedron associated with the on vertex
    tableEntry.append( 4 )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[2] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( analogs[0] )

    # Add the wedge adjacent to the tetrahedron
    tableEntry.append( 6 )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[2] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[3] )

    # Add the pyramid adjacent to the wedge
    tableEntry.append( 5 )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] ) 

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Three
#                 Two corners on the same triangle are "on".
#--------------------------------------------------------------------------------

def isWedgeCaseThree( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:
        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        if areOnSameTriangle( clipBits[0], clipBits[1] ):
            return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseThree(onBits, offBits):
    return handleWedgeCaseTough( onBits, offBits )

#--------------------------------------------------------------------------------
#                              Unique Case Four
#          Two adjacent corners on two different triangles are "on".
#--------------------------------------------------------------------------------

def isWedgeCaseFour( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:
        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits
        
        if areAdjacent( clipBits[0], clipBits[1] ) and not areOnSameTriangle( clipBits[0], clipBits[1] ):
            return True
        else:
            return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseFour(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    vertexZero = clipBits[0]

    analogs = getAnalogsFromVertex( vertexZero )

    # Add the wedge associated with the two clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[2] ) ] )
    tableEntry.append( analogs[3] )
    tableEntry.append( wedgeMidpoint[ ( analogs[3], analogs[4] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[3], analogs[5] ) ] )

    # Add the hexahedron associated with the non-clip bits
    tableEntry.append( 8 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[4] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[2] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[3], analogs[5] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[3], analogs[4] ) ] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Five
#         Two non-adjacent corners on two different triangles are "on".
#--------------------------------------------------------------------------------

def isWedgeCaseFive( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:
        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits
        
        if not areAdjacent( clipBits[0], clipBits[1] ) and not areOnSameTriangle( clipBits[0], clipBits[1] ):
            return True
        else:
            return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseFive(onBits, offBits):
    return handleWedgeCaseTough( onBits, offBits )

#--------------------------------------------------------------------------------
#                              Unique Case Six
#                 All three corners of one triangle are "on".
#--------------------------------------------------------------------------------

def isWedgeCaseSix( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        if len( onBits ) == 3:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits
        
        if areAdjacent( clipBits[0], clipBits[1] ) and areAdjacent( clipBits[1], clipBits[2] ) and areAdjacent( clipBits[2], clipBits[0] ):
            if areOnSameTriangle( clipBits[0], clipBits[1] ) and areOnSameTriangle( clipBits[1], clipBits[2] ) and areOnSameTriangle( clipBits[2], clipBits[0] ):
                return True
            else:
                return False
        else:
            return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseSix(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    # Since this clipping case is symmetric, always consider the set with vertex 0 to be the clip bits
    if len( onBits ) == 3 and 0 in onBits:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Arbitrarily choose one of the clip bits to be vertex zero
    vertexZero = clipBits[0]

    analogs = getAnalogsFromVertex( vertexZero )

    # Add the first wedge with the clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[1], analogs[4] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[2], analogs[5] ) ] )

    # Add the second wedge with the non-clip bits
    tableEntry.append( 6 )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[5] )
    tableEntry.append( analogs[4] )
    tableEntry.append( wedgeMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[2], analogs[5] ) ] )
    tableEntry.append( wedgeMidpoint[ ( analogs[1], analogs[4] ) ] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Seven
#               Two corners on one triangle and an adjacent 
#                  corner on the other triangle are "on".
#--------------------------------------------------------------------------------

def isWedgeCaseSeven( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        if len( onBits ) == 3:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Find the lone vertex that is on the other triangle from the other two "on" vertices
        if areOnSameTriangle( clipBits[1], clipBits[2] ) and not areOnSameTriangle( clipBits[0], clipBits[1] ):
            loneVertex = clipBits[0]
        elif areOnSameTriangle( clipBits[0], clipBits[2] ) and not areOnSameTriangle( clipBits[1], clipBits[2] ):
            loneVertex = clipBits[1]
        elif areOnSameTriangle( clipBits[0], clipBits[1] ) and not areOnSameTriangle( clipBits[2], clipBits[0] ):
            loneVertex = clipBits[2]
        else:
            # all 3 "on" vertices are on the same triangle, so it is not Case 7.
            return False

        # If the loneVertex is adjacent to either of the other "on" vertices, then return True.
        clipBitsCopy = list( clipBits )
        clipBitsCopy.remove( loneVertex )
        if areAdjacent( loneVertex, clipBitsCopy[0] ) or areAdjacent( loneVertex, clipBitsCopy[1] ):
            return True
        else:
            return False
    else:
        return False

#--------------------------------------------------------------------------------

def handleWedgeCaseSeven(onBits, offBits):
    return handleWedgeCaseTough( onBits, offBits )

#--------------------------------------------------------------------------------
#                              Unique Case Eight
#               Two corners on one triangle and a non-adjacent 
#                   corner on the other triangle are "on".
#--------------------------------------------------------------------------------

def isWedgeCaseEight( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        if len( onBits ) == 3:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        # Find the lone vertex that is on the other triangle from the other two "on" vertices
        if areOnSameTriangle( clipBits[1], clipBits[2] ) and not areOnSameTriangle( clipBits[0], clipBits[1] ):
            loneVertex = clipBits[0]
        elif areOnSameTriangle( clipBits[0], clipBits[2] ) and not areOnSameTriangle( clipBits[1], clipBits[2] ):
            loneVertex = clipBits[1]
        elif areOnSameTriangle( clipBits[0], clipBits[1] ) and not areOnSameTriangle( clipBits[2], clipBits[0] ):
            loneVertex = clipBits[2]
        else:
            # all 3 "on" vertices are on the same triangle, so it is not Case 7.
            return False

        # If the loneVertex is adjacent to either of the other "on" vertices, then return False.
        clipBitsCopy = list( clipBits )
        clipBitsCopy.remove( loneVertex )
        if areAdjacent( loneVertex, clipBitsCopy[0] ) or areAdjacent( loneVertex, clipBitsCopy[1] ):
            return False
        else:
            return True
    else:
        return False

#--------------------------------------------------------------------------------

# Case Eight: Two corners on one triangle and a non-adjacent corner on the other triangle are on.
def handleWedgeCaseEight(onBits, offBits):
    return handleWedgeCaseTough( onBits, offBits )

#--------------------------------------------------------------------------------
#                                  Tough Case
#         There is no way to decompose cleanly using 2D clipping analogs.
#--------------------------------------------------------------------------------

# Wedge Tough Case: There is no way to decompose cleanly using 2D clipping analogs. 
def handleWedgeCaseTough(onBits, offBits):
    return [4,0,2,1,15,
            4,3,4,5,15,
            5,0,3,5,2,15,
            5,1,2,5,4,15,
            5,0,1,4,3,15,-1]

#--------------------------------------------------------------------------------
#                              Wrapper Functions
#--------------------------------------------------------------------------------

def determineWedgeCase( caseNum ):
    onBits, offBits = getOnBitWedgeIndices( caseNum )
    
    if isWedgeCaseOne( onBits, offBits ):
        return 1
    elif isWedgeCaseTwo( onBits, offBits ):
        return 2
    elif isWedgeCaseThree( onBits, offBits ):
        return 3
    elif isWedgeCaseFour( onBits, offBits ):
        return 4
    elif isWedgeCaseFive( onBits, offBits ):
        return 5
    elif isWedgeCaseSix( onBits, offBits ):
        return 6
    elif isWedgeCaseSeven( onBits, offBits ):
        return 7
    elif isWedgeCaseEight( onBits, offBits ):
        return 8

#--------------------------------------------------------------------------------

def handleWedgeClipping( clippingCase, onBits, offBits ):
    if clippingCase == 1:
        return handleWedgeCaseOne( onBits, offBits )
    elif clippingCase == 2:
        return handleWedgeCaseTwo( onBits, offBits )
    elif clippingCase == 3:
        return handleWedgeCaseThree( onBits, offBits )
    elif clippingCase == 4:
        return handleWedgeCaseFour( onBits, offBits )
    elif clippingCase == 5:
        return handleWedgeCaseFive( onBits, offBits )
    elif clippingCase == 6:
        return handleWedgeCaseSix( onBits, offBits )
    elif clippingCase == 7:
        return handleWedgeCaseSeven( onBits, offBits )
    elif clippingCase == 8:
        return handleWedgeCaseEight( onBits, offBits )
    else:
        return handleWedgeCaseTough( onBits, offBits )

#--------------------------------------------------------------------------------
#                                     Main
#--------------------------------------------------------------------------------

def main():
    # Determine which vertices correspond to the ones in the clip case, re-assign the vertex IDs based on this, and then handle the clipping case normally
    clipTable = []
    numCases = 2**6
    for clipTableIndex in range(0, numCases):
        # Determine which clipping case it currently is
        clippingCase = determineWedgeCase( clipTableIndex )
        
        # Determine which bits are "on" or "off" for the current case
        onBits, offBits = getOnBitWedgeIndices( clipTableIndex )

        # Call the appropriate function to generate the clipping table entry
        clipTableEntry = handleWedgeClipping( clippingCase, onBits, offBits )

        # Append the return clipping table entry to the clipTable list 
        clipTable.append( clipTableEntry )

    # Print the table to the console
    for entry in clipTable:
        print( entry )

    # Write out the table to a file
    file = open("wedge_clipping_table.txt", "w")
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