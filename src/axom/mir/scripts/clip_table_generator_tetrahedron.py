# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#--------------------------------------------------------------------------------
#  Clip Table Generator - Tetrahedron | Marko Sterbentz 7/5/2019
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                             Midpoint Lookup Table
#--------------------------------------------------------------------------------

tetMidpoint = {
    (0,1) : 4,
    (1,0) : 4,
    (1,2) : 5,
    (2,1) : 5,
    (0,2) : 6,
    (2,0) : 6,
    (0,3) : 7,
    (3,0) : 7,
    (1,3) : 8,
    (3,1) : 8,
    (2,3) : 9,
    (3,2) : 9
}

#--------------------------------------------------------------------------------
#                              Utility Functions
#--------------------------------------------------------------------------------

def isBitOn(binaryNum, position, maxNumBits):
    bitMask = 1 << (maxNumBits - 1 - position)
    return (binaryNum & bitMask) == bitMask

#--------------------------------------------------------------------------------

def getOnBitTetrahedronIndices( caseNum ):
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

# Get the neighboring vertices of the tet vertex such that by appending loneVertex to the end of 
# the list, you get a tet is proper .vtk order.
def getOrderedNeighborsCaseTwo( loneVertex ):
    if loneVertex == 0:
        return [1,3,2]
    elif loneVertex == 1:
        return [0,2,3]
    elif loneVertex == 2:
        return [0,3,1]
    elif loneVertex == 3:
        return [0,1,2]

#--------------------------------------------------------------------------------

# Given any two pairs of onBits to clip, assume they correspond to vertex 2 and 3, 
# and map the two remaining vertices to vertex 0 and 1.
def getVertexAnalogs( onBits ):

    v2Analog = onBits[0]
    v3Analog = onBits[1]

    if (v2Analog, v3Analog) == (2, 3):
        v0Analog = 0
        v1Analog = 1
    elif (v2Analog, v3Analog) == (3, 2):
        v0Analog = 1
        v1Analog = 0
    elif (v2Analog, v3Analog) == (0, 1):
        v0Analog = 2
        v1Analog = 3
    elif (v2Analog, v3Analog) == (1, 0):
        v0Analog = 3
        v2Analog = 2
    elif (v2Analog, v3Analog) == (0, 2):
        v0Analog = 3
        v1Analog = 1
    elif (v2Analog, v3Analog) == (2, 0):
        v0Analog = 1
        v1Analog = 3
    elif (v2Analog, v3Analog) == (0, 3):
        v0Analog = 1
        v1Analog = 2
    elif (v2Analog, v3Analog) == (3, 0):
        v0Analog = 2
        v1Analog = 1
    elif (v2Analog, v3Analog) == (1, 2):
        v0Analog = 0
        v1Analog = 3
    elif (v2Analog, v3Analog) == (2, 1):
        v0Analog = 3
        v1Analog = 0
    elif (v2Analog, v3Analog) == (1, 3):
        v0Analog = 2
        v1Analog = 0
    elif (v2Analog, v3Analog) == (3, 1):
        v0Analog = 0
        v1Analog = 2

    return [ v0Analog, v1Analog, v2Analog, v3Analog ]

#--------------------------------------------------------------------------------
#                              Unique Case One
#                        All bits are "on" or "off".
#--------------------------------------------------------------------------------

def isTetrahedronCaseOne( onBits, offBits ):
    if len( onBits ) == 4 or len( offBits ) == 4:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleTetrahedronCaseOne(onBits, offBits):
    return [4,0,1,2,3,-1]

#--------------------------------------------------------------------------------
#                              Unique Case Two
#                        One of the corners is "on".
#--------------------------------------------------------------------------------

def isTetrahedronCaseTwo( onBits, offBits ):
    if len( onBits ) == 1 or len ( offBits ) == 1:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleTetrahedronCaseTwo(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 1:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits
    
    # Get the index of the lone vertex to clip off
    loneVertex = clipBits[0]

    # Get the midpoint nodes between the lone vertex and the others
    midpoints = []
    midpoints.append( tetMidpoint[ ( loneVertex, nonClipBits[0] ) ] )
    midpoints.append( tetMidpoint[ ( loneVertex, nonClipBits[1] ) ] )
    midpoints.append( tetMidpoint[ ( loneVertex, nonClipBits[2] ) ] )

    # Get the neighbors of the lone vertex in an order that conforms to .vtk file standard and avoids inverted tets and other rendering issues
    orderedNeighbors = getOrderedNeighborsCaseTwo( loneVertex )

    # Add the tet containing the lone vertex 
    tableEntry.append( 4 )
    tableEntry.append( tetMidpoint[ ( loneVertex, orderedNeighbors[0]) ] )
    tableEntry.append( tetMidpoint[ ( loneVertex, orderedNeighbors[1]) ] )
    tableEntry.append( tetMidpoint[ ( loneVertex, orderedNeighbors[2]) ] )
    tableEntry.append( loneVertex )

    # Add the wedge
    tableEntry.append( 6 )
    tableEntry.append( tetMidpoint[ ( loneVertex, orderedNeighbors[0]) ] )
    tableEntry.append( tetMidpoint[ ( loneVertex, orderedNeighbors[1]) ] )
    tableEntry.append( tetMidpoint[ ( loneVertex, orderedNeighbors[2]) ] )
    tableEntry.append( orderedNeighbors[0] )
    tableEntry.append( orderedNeighbors[1] )
    tableEntry.append( orderedNeighbors[2] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Unique Case Three
#                        Two of the corners are "on".
#--------------------------------------------------------------------------------

def isTetrahedronCaseThree( onBits, offBits ):
    if len( onBits ) == 2 or len ( offBits ) == 2:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleTetrahedronCaseThree(onBits, offBits):
    # Init the table entry list
    tableEntry = []

    if len( onBits ) == 2:
        clipBits = onBits
        nonClipBits = offBits
    else:
        clipBits = offBits
        nonClipBits = onBits

    # Map the given vertices to clip to known vertex indices and clip
    analogs = getVertexAnalogs( clipBits )

    # Add the first wedge that contains the two on vertices
    tableEntry.append( 6 )
    tableEntry.append( analogs[2] )
    tableEntry.append( tetMidpoint[ ( analogs[2], analogs[1] ) ] )
    tableEntry.append( tetMidpoint[ ( analogs[2], analogs[0] ) ] )
    tableEntry.append( analogs[3] )
    tableEntry.append( tetMidpoint[ ( analogs[3], analogs[1] ) ] )
    tableEntry.append( tetMidpoint[ ( analogs[3], analogs[0] ) ] )

    # Add the second wedge that contains the two off vertices
    tableEntry.append( 6 )
    tableEntry.append( analogs[0] )
    tableEntry.append( tetMidpoint[ ( analogs[0], analogs[3] ) ] )
    tableEntry.append( tetMidpoint[ ( analogs[0], analogs[2] ) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( tetMidpoint[ ( analogs[1], analogs[3] ) ] )
    tableEntry.append( tetMidpoint[ ( analogs[1], analogs[2] ) ] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Wrapper Functions
#--------------------------------------------------------------------------------

def determineTetrahedronCase( caseNum ):
    onBits, offBits = getOnBitTetrahedronIndices( caseNum )
    
    if isTetrahedronCaseOne( onBits, offBits ):
        return 1
    elif isTetrahedronCaseTwo( onBits, offBits ):
        return 2
    elif isTetrahedronCaseThree( onBits, offBits ):
        return 3

#--------------------------------------------------------------------------------

def handleTetrahedronClipping( clippingCase, onBits, offBits ):
    if clippingCase == 1:
        return handleTetrahedronCaseOne( onBits, offBits )
    elif clippingCase == 2:
        return handleTetrahedronCaseTwo( onBits, offBits )
    elif clippingCase == 3:
        return handleTetrahedronCaseThree( onBits, offBits )

#--------------------------------------------------------------------------------
#                                     Main
#--------------------------------------------------------------------------------

def main():
    # Determine which vertices correspond to the ones in the clip case, re-assign the vertex IDs based on this, and then handle the clipping case normally
    clipTable = []
    numCases = 2**4
    for clipTableIndex in range(0, numCases):
        # Determine which clipping case it currently is
        clippingCase = determineTetrahedronCase( clipTableIndex )
        
        # Determine which bits are "on" or "off" for the current case
        onBits, offBits = getOnBitTetrahedronIndices( clipTableIndex )

        # Call the appropriate function to generate the clipping table entry
        clipTableEntry = handleTetrahedronClipping( clippingCase, onBits, offBits )

        # Append the return clipping table entry to the clipTable list 
        clipTable.append( clipTableEntry )

    # Print the table to the console
    for entry in clipTable:
        print( entry )

    # Write out the table to a file
    file = open("tetrahedron_clipping_table.txt", "w")
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