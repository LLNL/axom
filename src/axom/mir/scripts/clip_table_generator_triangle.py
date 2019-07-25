# Copyright (c) 2017-2019, Lawrence Livermore National Security, LLC and
# other Axom Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (BSD-3-Clause)

#--------------------------------------------------------------------------------
#  Clip Table Generator - Triangle | Marko Sterbentz 7/25/2019
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
#                             Midpoint Lookup Table
#--------------------------------------------------------------------------------

triangleMidpoint = {
    (0,1) : 3,
    (1,0) : 3,
    (1,2) : 4,
    (2,1) : 4,
    (0,2) : 5,
    (2,0) : 5,
}

#--------------------------------------------------------------------------------
#                              Utility Functions
#--------------------------------------------------------------------------------

def isBitOn(binaryNum, position, maxNumBits):
    bitMask = 1 << (maxNumBits - 1 - position)
    return (binaryNum & bitMask) == bitMask

#--------------------------------------------------------------------------------

def getOnBitTriangleIndices( caseNum ):
    numBits = 3
    onBits = []
    offBits = []
    for position in range(0, numBits):
        if isBitOn(caseNum, position, numBits):
            onBits.append( position )
        else:
            offBits.append( position )
    return onBits, offBits

#--------------------------------------------------------------------------------

# Remap the vertices such that the given vertex is vertex zero.
def getAnalogsFromVertex( vertexZero ):
    v0Analog = vertexZero

    if vertexZero == 0:
        v1Analog = 1
        v2Analog = 2
    elif vertexZero == 1:
        v1Analog = 2
        v2Analog = 0
    elif vertexZero == 2:
        v1Analog = 0
        v2Analog = 1

    return [ v0Analog, v1Analog, v2Analog ]

#--------------------------------------------------------------------------------
#                              Unique Case One
#                        All bits are "on" or "off".
#--------------------------------------------------------------------------------

def isTriangleCaseOne( onBits, offBits ):
    if len( onBits ) == 3 or len( offBits ) == 3:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleTriangleCaseOne(onBits, offBits):
    return [3,0,1,2,-1]

#--------------------------------------------------------------------------------
#                              Unique Case Two
#                        One of the corners is "on".
#--------------------------------------------------------------------------------

def isTriangleCaseTwo( onBits, offBits ):
    if len( onBits ) == 1 or len ( offBits ) == 1:
        return True
    else:
        return False

#--------------------------------------------------------------------------------

def handleTriangleCaseTwo(onBits, offBits):
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
    tableEntry.append( triangleMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( triangleMidpoint[ ( analogs[0], analogs[2]) ] )

    # Add the quad containing the "off" vertices
    tableEntry.append( 4 )
    tableEntry.append( triangleMidpoint[ ( analogs[0], analogs[1]) ] )
    tableEntry.append( analogs[1] )
    tableEntry.append( analogs[2] )
    tableEntry.append( triangleMidpoint[ ( analogs[0], analogs[2]) ] )

    # Add the table entry delimiter
    tableEntry.append( -1 )

    return tableEntry

#--------------------------------------------------------------------------------
#                              Wrapper Functions
#--------------------------------------------------------------------------------

def determineTriangleCase( caseNum ):
    onBits, offBits = getOnBitTriangleIndices( caseNum )
    
    if isTriangleCaseOne( onBits, offBits ):
        return 1
    elif isTriangleCaseTwo( onBits, offBits ):
        return 2

#--------------------------------------------------------------------------------

def handleTriangleClipping( clippingCase, onBits, offBits ):
    if clippingCase == 1:
        return handleTriangleCaseOne( onBits, offBits )
    elif clippingCase == 2:
        return handleTriangleCaseTwo( onBits, offBits )

#--------------------------------------------------------------------------------
#                                     Main
#--------------------------------------------------------------------------------

def main():
    # Determine which vertices correspond to the ones in the clip case, re-assign the vertex IDs based on this, and then handle the clipping case normally
    clipTable = []
    numCases = 2**3
    for clipTableIndex in range(0, numCases):
        # Determine which clipping case it currently is
        clippingCase = determineTriangleCase( clipTableIndex )
        
        # Determine which bits are "on" or "off" for the current case
        onBits, offBits = getOnBitTriangleIndices( clipTableIndex )

        # Call the appropriate function to generate the clipping table entry
        clipTableEntry = handleTriangleClipping( clippingCase, onBits, offBits )

        # Append the return clipping table entry to the clipTable list 
        clipTable.append( clipTableEntry )

    # Print the table to the console
    for entry in clipTable:
        print( entry )

    # Write out the table to a file
    file = open("triangle_clipping_table.txt", "w")
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