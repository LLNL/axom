# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#--------------------------------------------------------------------------------
#  Clip Table Generator - Pyramid | Marko Sterbentz 7/5/2019
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                             Midpoint Lookup Table
#--------------------------------------------------------------------------------

pyramidMidpoint = {
    (0,1) : 5,
    (1,0) : 5,
    (1,2) : 6,
    (2,1) : 6,
    (2,3) : 7,
    (3,2) : 7,
    (0,3) : 8,
    (3,0) : 8,
    (0,4) : 9,
    (4,0) : 9,
    (1,4) : 10,
    (4,1) : 10,
    (2,4) : 11,
    (4,2) : 11,
    (3,4) : 12,
    (4,3) : 12
}

#--------------------------------------------------------------------------------
#                              Utility Functions
#--------------------------------------------------------------------------------

def isBitOn(binaryNum, position, maxNumBits):
    bitMask = 1 << (maxNumBits - 1 - position)
    return (binaryNum & bitMask) == bitMask

#--------------------------------------------------------------------------------

def getOnBitPyramidIndices( caseNum ):
    numPyramidBits = 5
    onBits = []
    offBits = []
    for position in range(0, 5):
        if isBitOn(caseNum, position, numPyramidBits):
            onBits.append( position )
        else:
            offBits.append( position )
    return onBits, offBits

#--------------------------------------------------------------------------------

# Used to get the vertex analogs in unique case 2.
def getAnalogsFromOneBaseVertex( baseVertex ):
    v0Analog = baseVertex
    v4Analog = 4
    if baseVertex == 0:
        v1Analog = 1
        v2Analog = 2
        v3Analog = 3
    if baseVertex == 1:
        v1Analog = 2
        v2Analog = 3
        v3Analog = 0
    if baseVertex == 2:
        v1Analog = 3
        v2Analog = 0
        v3Analog = 1
    if baseVertex == 3:
        v1Analog = 0
        v2Analog = 1
        v3Analog = 2

    return [ v0Analog, v1Analog, v2Analog, v3Analog, v4Analog ]

#--------------------------------------------------------------------------------

# Returns a list of node indices adjacent to the given node index.
def getAdjacentPyramidNodeIndices(node):
    if node == 0:
        return [1,3,4]
    elif node == 1:
        return [0,2,4]
    elif node == 2:
        return [1,3,4]
    elif node == 3:
        return [0,2,4]
    else:
        return [0,1,2,3]

#--------------------------------------------------------------------------------

def isSquareBaseNode(node):
    return node == 0 or node == 1 or node == 2 or node == 3

#--------------------------------------------------------------------------------
#                              Unique Case One
#                           All bits "on" or "off".
#--------------------------------------------------------------------------------

def isPyramidCaseOne( onBits ):
    if len( onBits ) == 5 or len( onBits ) == 0:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handlePyramidCaseOne(onBits, offBits):
    return [5,0,1,2,3,4,-1]

#--------------------------------------------------------------------------------
#                              Unique Case Two
#                  One of the square base corners is "on".
#--------------------------------------------------------------------------------

def isPyramidCaseTwo( onBits ):
    if len( onBits ) == 1 and isSquareBaseNode( onBits[0] ):
        return True
    elif len( onBits ) == 4 and 4 in onBits:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handlePyramidCaseTwo(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 1:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    loneVertex = clipBits[0]

    analogs = getAnalogsFromOneBaseVertex( loneVertex )

    # Add the tet including the lone vertex to the table entry
    tableEntry.append( 4 )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[4]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( analogs[0] )

    # Add the wedge to the table entry
    tableEntry.append( 6 )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[4]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[3] )

    # Add the opposite tet to the table entry
    tableEntry.append( 4 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[2] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Three
#                         The top, tip corner is "on".
#--------------------------------------------------------------------------------

def isPyramidCaseThree( onBits ):
    if len( onBits ) == 1 and onBits[0] == 4:
        return True
    elif len( onBits ) == 4 and 4 not in onBits:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handlePyramidCaseThree(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    # Add the top pyramid to the table entry
    tableEntry.append( 5 )
    tableEntry.append( 9 )
    tableEntry.append( 10 )
    tableEntry.append( 11 )
    tableEntry.append( 12 )
    tableEntry.append( 4 )

    # Add the bottom hexahedron to the table entry
    tableEntry.append( 8 )
    tableEntry.append( 0 )
    tableEntry.append( 1 )
    tableEntry.append( 2 )
    tableEntry.append( 3 )
    tableEntry.append( 9 )
    tableEntry.append( 10 )
    tableEntry.append( 11 )
    tableEntry.append( 12 )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Four
#                 Two adjacent square base corners are "on".
#--------------------------------------------------------------------------------

def isPyramidCaseFour( onBits, offBits ):
    if len( onBits ) == 2 and 4 not in onBits:
        adjacentNodes = getAdjacentPyramidNodeIndices( onBits[0] )
        if onBits[1] in adjacentNodes:
            return True
    elif len( offBits ) == 2 and 4 not in offBits:
        adjacentNodes = getAdjacentPyramidNodeIndices( offBits[0] )
        if offBits[1] in adjacentNodes:
            return True
    else:
        return False

#--------------------------------------------------------------------------------

def handlePyramidCaseFour(onBits, offBits):
    return handlePyramidCaseTough( onBits, offBits )
    
#--------------------------------------------------------------------------------
#                              Unique Case Five
#       The top, tip corner and one of the square base corners are "on".
#--------------------------------------------------------------------------------

def isPyramidCaseFive( onBits ):
    if len( onBits ) == 2 and 4 in onBits:
        return True
    elif len( onBits ) == 3 and 4 not in onBits:
        return True
    else:
        return False

#--------------------------------------------------------------------------------
 
def handlePyramidCaseFive(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits
    
    # Get the lone, base vertex and the tip vertex to be clipped
    tipVertex = 4
    clipBits.remove( tipVertex )    # Remove the pyramid tip index, leaving just the lone vertex in the list
    loneVertex = clipBits[0]

    analogs = getAnalogsFromOneBaseVertex( loneVertex )

    # Add the wedge containing the two vertices to be clipped
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[1] ) ] )
    tableEntry.append( analogs[4] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[3] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[1] ) ] )

    # Add the tetrahedron adjacent to the first wedge
    tableEntry.append( 4 )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[1] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[2] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[3] ) ] )
    tableEntry.append( analogs[4] )

    # Add the wedge adjacent to the first wedge
    tableEntry.append( 6 )
    tableEntry.append( analogs[1] )
    tableEntry.append( pyramidMidpoint[ ( analogs[1], analogs[4] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[1], analogs[0] ) ] )
    tableEntry.append( analogs[3] )
    tableEntry.append( pyramidMidpoint[ ( analogs[3], analogs[4] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[3], analogs[0] ) ] )

    # Add the wedge that is not adjacent to the first wedge
    tableEntry.append( 6 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( analogs[2] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[1] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[3] ) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[4], analogs[2] ) ] )    

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Six
#                Two non-adjacent square base corners are "on".
#--------------------------------------------------------------------------------

def isPyramidCaseSix( onBits, offBits ):
    if len ( onBits ) == 2 and 4 not in onBits:
        adjacentNodes = getAdjacentPyramidNodeIndices( onBits[0] )
        if onBits[1] not in adjacentNodes:
            return True
    elif len ( offBits ) == 2 and 4 not in offBits:
        adjacentNodes = getAdjacentPyramidNodeIndices( offBits[0] )
        if offBits[1] not in adjacentNodes:
            return True
    else:
        return False

#--------------------------------------------------------------------------------

def handlePyramidCaseSix(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Get the base vertices to be clipped and the tip vertex
    tipVertex = 4
    nonClipBits.remove( tipVertex )

    if 0 in clipBits:
        # The case where 0 and 2 are on/off
        analogs = getAnalogsFromOneBaseVertex( 0 )
    else:
        # The case where 1 and 3 are on/off
        analogs = getAnalogsFromOneBaseVertex( 1 ) 

    # Add one of the tets adjacent to one of the base verties to be clipped
    tableEntry.append( 4 )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[4]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( analogs[0] )

    # Add the wedge adjacent to the previous tet
    tableEntry.append( 6 )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[4]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[3] )

    # Add the wedge adjacent to the next tet
    tableEntry.append( 6 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[4] )
    tableEntry.append( analogs[3] )
    tableEntry.append( pyramidMidpoint[ ( analogs[2], analogs[1]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[2], analogs[4]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[2], analogs[3]) ] )

    # Add one the other tet adjacent to the other base vertices to be clipped
    tableEntry.append( 4 )
    tableEntry.append( pyramidMidpoint[ ( analogs[2], analogs[1]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[2], analogs[3]) ] )
    tableEntry.append( pyramidMidpoint[ ( analogs[2], analogs[4]) ] )
    tableEntry.append( analogs[2] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                                  Tough Case
#         There is no way to decompose cleanly using 2D clipping analogs.
#--------------------------------------------------------------------------------

def handlePyramidCaseTough( onBits, offBits ):
    return [4, 0, 4, 1, 13,
            4, 1, 4, 2, 13,
            4, 2, 4, 3, 13,
            4, 3, 4, 0, 13,
            4, 0, 1, 3, 13,
            4, 1, 2, 3, 13, -1]

#--------------------------------------------------------------------------------
#                              Wrapper Functions
#--------------------------------------------------------------------------------

def determinePyramidCase( caseNum ):
    onBits, offBits = getOnBitPyramidIndices( caseNum )

    if isPyramidCaseOne( onBits ):
        return 1
    elif isPyramidCaseTwo( onBits ):
        return 2
    elif isPyramidCaseThree( onBits ):
        return 3
    elif isPyramidCaseFour( onBits, offBits ):
        return 4
    elif isPyramidCaseFive( onBits ):
        return 5
    elif isPyramidCaseSix( onBits, offBits ):
        return 6

#--------------------------------------------------------------------------------

def handlePyramidClipping( clippingCase, onBits, offBits ):
    if clippingCase == 1:
        return handlePyramidCaseOne( onBits, offBits )
    elif clippingCase == 2:
        return handlePyramidCaseTwo( onBits, offBits )
    elif clippingCase == 3:
        return handlePyramidCaseThree( onBits, offBits )
    elif clippingCase == 4:
        return handlePyramidCaseFour( onBits, offBits )
    elif clippingCase == 5:
        return handlePyramidCaseFive( onBits, offBits )
    elif clippingCase == 6:
        return handlePyramidCaseSix( onBits, offBits )

#--------------------------------------------------------------------------------
#                                     Main
#--------------------------------------------------------------------------------

def main():
    # Determine which vertices correspond to the ones in the clip case, re-assign the vertex IDs based on this, and then handle the clipping case normally
    clipTable = []
    numCases = 2**5
    for clipTableIndex in range(0, numCases):
        # Determine which clipping case it currently is
        clippingCase = determinePyramidCase( clipTableIndex )
        
        # Determine which bits are "on" or "off" for the current case
        onBits, offBits = getOnBitPyramidIndices( clipTableIndex )

        # Call the appropriate function to generate the clipping table entry
        clipTableEntry = handlePyramidClipping( clippingCase, onBits, offBits )

        # Append the return clipping table entry to the clipTable list 
        clipTable.append( clipTableEntry )

    # Print the table to the console
    for entry in clipTable:
        print( entry )

    # Write out the table to a file
    file = open("pyramid_clipping_table.txt", "w")
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