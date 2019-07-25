# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#--------------------------------------------------------------------------------
#  Clip Table Generator - Quadrilateral | Marko Sterbentz 7/25/2019
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                             Midpoint Lookup Table
#--------------------------------------------------------------------------------

quadMidpoint = {
    (0,1) : 4,
    (1,0) : 4,
    (1,2) : 5,
    (2,1) : 5,
    (2,3) : 6,
    (3,2) : 6,
    (3,0) : 7,
    (0,3) : 7
}

#--------------------------------------------------------------------------------
#                              Utility Functions
#--------------------------------------------------------------------------------

def isBitOn(binaryNum, position, maxNumBits):
    bitMask = 1 << (maxNumBits - 1 - position)
    return (binaryNum & bitMask) == bitMask

#--------------------------------------------------------------------------------

def getOnBitQuadIndices( caseNum ):
    numBits = 4
    onBits = []
    offBits = []
    for position in range(0, numBits):
        if isBitOn(caseNum, position, numBits):
            onBits.append( position )
        else:
            offBits.append( position )
    return onBits, offBits

#--------------------------------------------------------------------------------

def getAdjacentVertices( vertex ):
    if vertex == 0:
        return [1,3]
    elif vertex == 1:
        return [0,2]
    elif vertex == 2:
        return [1,3]
    elif vertex == 3:
        return [0,2]

#--------------------------------------------------------------------------------

# Remap the vertices such that the given vertex is vertex zero.
def getAnalogsFromVertex( vertexZero ):
    v0Analog = vertexZero

    if vertexZero == 0:
        v1Analog = 1
        v2Analog = 2
        v3Analog = 3
    elif vertexZero == 1:
        v1Analog = 2
        v2Analog = 3
        v3Analog = 0
    elif vertexZero == 2:
        v1Analog = 3
        v2Analog = 0
        v3Analog = 1
    elif vertexZero == 3:
        v1Analog = 0
        v2Analog = 1
        v3Analog = 2

    return [ v0Analog, v1Analog, v2Analog, v3Analog ]

#--------------------------------------------------------------------------------

# Remap the vertices such that the given vertices correspond to vertex zero and vertex one.
def getAnalogsFromEdge( vertexZero, vertexOne ):

    if ( vertexZero, vertexOne ) == (0, 1) or ( vertexZero, vertexOne ) == (1, 0):
        v0Analog = 0
        v1Analog = 1
        v2Analog = 2
        v3Analog = 3
    if ( vertexZero, vertexOne ) == (1, 2) or ( vertexZero, vertexOne ) == (2, 1):
        v0Analog = 1
        v1Analog = 2
        v2Analog = 3
        v3Analog = 0
    elif ( vertexZero, vertexOne ) == (2, 3) or ( vertexZero, vertexOne ) == (3, 2):
        v0Analog = 2
        v1Analog = 3
        v2Analog = 0
        v3Analog = 1
    elif ( vertexZero, vertexOne ) == (0, 3) or ( vertexZero, vertexOne ) == (3, 0):
        v0Analog = 3
        v1Analog = 0
        v2Analog = 1
        v3Analog = 2

    return [ v0Analog, v1Analog, v2Analog, v3Analog ]


#--------------------------------------------------------------------------------
#                              Unique Case One
#                        All bits are "on" or "off".
#--------------------------------------------------------------------------------

def isQuadCaseOne( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleQuadCaseOne(onBits, offBits):
    return [4,0,1,2,3,-1]

#--------------------------------------------------------------------------------
#                              Unique Case Two
#                        One of the corners is "on".
#--------------------------------------------------------------------------------

def isQuadCaseTwo( onBits, offBits ):
    if len( onBits ) == 1 or len ( offBits ) == 1:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleQuadCaseTwo(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 1:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits
    
    analogs = getAnalogsFromVertex( clipBits[0] )

    # Add the triangle containing the "on" vertex 
    tableEntry.append( 3 )
    tableEntry.append( analogs[0] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[3]) ] )

    # Add the quad adjacent to this triangle
    tableEntry.append( 4 )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[3]) ] )

    # Add the triangle not adjacent to the clipped triangle
    tableEntry.append( 3 )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[3] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Three
#                        Two non-adjacent corners "on".
#--------------------------------------------------------------------------------

def isQuadCaseThree( onBits, offBits ):
    if len( onBits ) == 2 or len( offBits ) == 2:

        if len( onBits ) == 2:
            clipBits = onBits
            nonClipBits = offBits
        else:
            clipBits = offBits
            nonClipBits = onBits

        if not clipBits[1] in getAdjacentVertices( clipBits[0] ):
            return True
        else:
            return False

    return False

#--------------------------------------------------------------------------------

def handleQuadCaseThree( onBits, offBits ):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits
    
    analogs = getAnalogsFromVertex( clipBits[0] )

    # Add the triangle adjacent to the first "on" vertex
    tableEntry.append( 3 )
    tableEntry.append( analogs[0] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[3]) ] )

    # Add the quad adjacent to the first triangle
    tableEntry.append( 4 )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[3] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[3]) ] )

    # Add the quad adjacent to the second triangle
    tableEntry.append( 4 )
    tableEntry.append( analogs[1] )
    tableEntry.append( quadMidpoint[ ( analogs[2], analogs[1]) ] )
    tableEntry.append( quadMidpoint[ ( analogs[2], analogs[3]) ] )
    tableEntry.append( analogs[3] )

    # Add the triangle adjacent to the second "on" vertex
    tableEntry.append( 3 )
    tableEntry.append( quadMidpoint[ ( analogs[2], analogs[1]) ] )
    tableEntry.append( analogs[2] )
    tableEntry.append( quadMidpoint[ ( analogs[2], analogs[3]) ] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Four
#                        Two adjacent corners are "on".
#--------------------------------------------------------------------------------

def isQuadCaseFour( onBits, offBits ):
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

    return False

#--------------------------------------------------------------------------------

def handleQuadCaseFour( onBits, offBits ):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits
    
    analogs = getAnalogsFromEdge( clipBits[0], clipBits[1] )

    # Add the quad adjacent to the "on" vertices
    tableEntry.append( 4 )
    tableEntry.append( analogs[0] )
    tableEntry.append( analogs[1] )
    tableEntry.append( quadMidpoint[ ( analogs[1], analogs[2]) ] )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[3]) ] )

    # Add the quad adjacent to the "off" vertices
    tableEntry.append( 4 )
    tableEntry.append( quadMidpoint[ ( analogs[0], analogs[3]) ] )
    tableEntry.append( quadMidpoint[ ( analogs[1], analogs[2]) ] )
    tableEntry.append( analogs[2] )
    tableEntry.append( analogs[3] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Wrapper Functions
#--------------------------------------------------------------------------------

def determineQuadCase( caseNum ):
    onBits, offBits = getOnBitQuadIndices( caseNum )
    
    if isQuadCaseOne( onBits, offBits ):
        return 1
    elif isQuadCaseTwo( onBits, offBits ):
        return 2
    elif isQuadCaseThree( onBits, offBits ):
        return 3
    elif isQuadCaseFour( onBits, offBits ):
        return 4

#--------------------------------------------------------------------------------

def handleQuadClipping( clippingCase, onBits, offBits ):
    if clippingCase == 1:
        return handleQuadCaseOne( onBits, offBits )
    elif clippingCase == 2:
        return handleQuadCaseTwo( onBits, offBits )
    elif clippingCase == 3:
        return handleQuadCaseThree( onBits, offBits )
    elif clippingCase == 4:
        return handleQuadCaseFour( onBits, offBits )

#--------------------------------------------------------------------------------
#                                     Main
#--------------------------------------------------------------------------------

def main():
    # Determine which vertices correspond to the ones in the clip case, re-assign the vertex IDs based on this, and then handle the clipping case normally
    clipTable = []
    numCases = 2**4
    for clipTableIndex in range(0, numCases):
        # Determine which clipping case it currently is
        clippingCase = determineQuadCase( clipTableIndex )
        
        # Determine which bits are "on" or "off" for the current case
        onBits, offBits = getOnBitQuadIndices( clipTableIndex )

        # Call the appropriate function to generate the clipping table entry
        clipTableEntry = handleQuadClipping( clippingCase, onBits, offBits )

        # Append the return clipping table entry to the clipTable list 
        clipTable.append( clipTableEntry )

    # Print the table to the console
    for entry in clipTable:
        print( entry )

    # Write out the table to a file
    file = open("quad_clipping_table.txt", "w")
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


# [4, 0, 1, 2, 3, -1]
# [3, 3, 7, 6, 4, 7, 0, 2, 6, 3, 0, 1, 2, -1]
# [3, 2, 6, 5, 4, 6, 3, 1, 5, 3, 3, 0, 1, -1]
# [4, 2, 3, 7, 5, 4, 5, 7, 0, 1, -1]  #
# [3, 1, 5, 4, 4, 5, 2, 0, 4, 3, 2, 3, 0, -1]
# [3, 1, 5, 4, 4, 5, 2, 0, 4, 4, 2, 6, 7, 0, 3, 6, 3, 7, -1]  #
# [4, 1, 2, 6, 4, 4, 4, 6, 3, 0, -1]  #
# [3, 0, 4, 7, 4, 4, 1, 3, 7, 3, 1, 2, 3, -1]
# [3, 0, 4, 7, 4, 4, 1, 3, 7, 3, 1, 2, 3, -1]
# [4, 3, 0, 4, 6, 4, 6, 4, 1, 2, -1]  #
# [3, 0, 4, 7, 4, 4, 1, 3, 7, 4, 1, 5, 6, 3, 3, 5, 2, 6, -1]  #
# [3, 1, 5, 4, 4, 5, 2, 0, 4, 3, 2, 3, 0, -1]
# [4, 0, 1, 5, 7, 4, 7, 5, 2, 3, -1]  #
# [3, 2, 6, 5, 4, 6, 3, 1, 5, 3, 3, 0, 1, -1]
# [3, 3, 7, 6, 4, 7, 0, 2, 6, 3, 0, 1, 2, -1]
# [4, 0, 1, 2, 3, -1]


[4, 0, 1, 2, 3, -1]
[3, 3, 7, 6, 4, 7, 0, 2, 6, 3, 0, 1, 2, -1]
[3, 2, 6, 5, 4, 6, 3, 1, 5, 3, 3, 0, 1, -1]
[4, 0, 1, 5, 7, 4, 7, 5, 2, 3, -1]
[3, 1, 5, 4, 4, 5, 2, 0, 4, 3, 2, 3, 0, -1]
[3, 0, 4, 7, 4, 4, 1, 3, 7, 4, 1, 5, 6, 3, 3, 5, 2, 6, -1]
[4, 3, 0, 4, 6, 4, 6, 4, 1, 2, -1]
[3, 0, 4, 7, 4, 4, 1, 3, 7, 3, 1, 2, 3, -1]
[3, 0, 4, 7, 4, 4, 1, 3, 7, 3, 1, 2, 3, -1]
[4, 3, 0, 4, 6, 4, 6, 4, 1, 2, -1]
[3, 0, 4, 7, 4, 4, 1, 3, 7, 4, 1, 5, 6, 3, 3, 5, 2, 6, -1]
[3, 1, 5, 4, 4, 5, 2, 0, 4, 3, 2, 3, 0, -1]
[4, 0, 1, 5, 7, 4, 7, 5, 2, 3, -1]
[3, 2, 6, 5, 4, 6, 3, 1, 5, 3, 3, 0, 1, -1]
[3, 3, 7, 6, 4, 7, 0, 2, 6, 3, 0, 1, 2, -1]
[4, 0, 1, 2, 3, -1]