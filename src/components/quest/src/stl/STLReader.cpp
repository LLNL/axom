/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */


/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file STLReader.cpp
 *
 * \date Dec 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */

#include "STLReader.hpp"

// ATK includes
#include "common/ATKMacros.hpp"
#include "common/CommonTypes.hpp"
#include "slic/slic.hpp"


// C/C++ includes
#include <cmath>    // for STL IO
#include <cstddef>  // for NULL
#include <cstdlib>  // for STL IO
#include <ctime>    // for STL IO
#include <fstream>  // for STL IO
#include <iomanip>  // for STL IO
#include <iostream> // for STL IO
#include <string>   // for STL string

using std::string;   // for STL Routines
using std::ifstream; // for STL Routines
using std::ofstream; // for STL Routines
using std::cout;     // for STL Routines

//------------------------------------------------------------------------------
//      STL IO Routine Definitions
//------------------------------------------------------------------------------
namespace
{

char ch_cap ( char c );
bool ch_eqi ( char c1, char c2 );
int ch_to_digit ( char c );

bool s_eqi ( string s1, string s2 );
int s_len_trim ( string s );
double s_to_r8 ( string s, int *lchar, bool *error );

int stla_offset_get ( void );

bool stla_read ( string input_file_name, int node_num, int face_num,
  double node_xyz[], int face_node[], double face_normal[] );
void stla_size ( string input_file_name, int *solid_num, int *node_num,
  int *face_num, int *text_num );

string word_next_read ( string s, bool *done );

//
//  This variable determines whether indices are 0 or 1 based.
//
static int stla_offset_value = 0;

} /* end anonymous namespace */

//------------------------------------------------------------------------------
//      STLReader Implementation
//------------------------------------------------------------------------------
namespace quest
{

STLReader::STLReader() :
        m_fileName(""),
        m_num_nodes(0),
        m_num_faces(0),
        m_nodes(ATK_NULLPTR),
        m_face_normals(ATK_NULLPTR),
        m_face_connectivity(ATK_NULLPTR)

{

}

//------------------------------------------------------------------------------
STLReader::~STLReader()
{
  this->clear();
}

//------------------------------------------------------------------------------
void STLReader::clear()
{

  if ( m_nodes != ATK_NULLPTR ) {
      delete [] m_nodes;
        m_nodes = ATK_NULLPTR;
  }

  if ( m_face_normals != ATK_NULLPTR ) {
      delete [] m_face_normals;
      m_face_normals = ATK_NULLPTR;
  }

  if ( m_face_connectivity != ATK_NULLPTR ) {
      delete [] m_face_connectivity;
      m_face_connectivity = ATK_NULLPTR;
  }

}

//------------------------------------------------------------------------------
void STLReader::read()
{
  SLIC_ASSERT( m_fileName != "" );

  // STEP 0: clear internal data-structures
  this->clear();

  // STEP 1: Query STL file for sizes
  int numSolids = 0;
  int numNodes  = 0;
  int numFaces  = 0;
  int numText   = 0;
  stla_size( m_fileName, &numSolids, &numNodes, &numFaces, &numText );

  m_num_nodes = numNodes;
  m_num_faces = numFaces;

  // STEP 2: Allocate internal data-structures
  m_nodes             = new double[ 3*numNodes ];
  m_face_connectivity = new int[ 3*numFaces ];
  m_face_normals      = new double[ 3*numFaces ];

  // STEP 3: Read in geometry to internal data-structures
  stla_read( m_fileName, numNodes, numFaces,
          m_nodes, m_face_connectivity, m_face_normals );

}

//------------------------------------------------------------------------------
void STLReader::getMesh(
        meshtk::UnstructuredMesh< meshtk::LINEAR_TRIANGLE >* mesh )
{
  /* Sanity checks */
  SLIC_ASSERT( mesh != ATK_NULLPTR );
  SLIC_ASSERT( m_nodes != ATK_NULLPTR );
  SLIC_ASSERT( m_face_connectivity != ATK_NULLPTR );

  for ( int i=0; i < m_num_nodes; ++i ) {
      mesh->insertNode( m_nodes[i*3], m_nodes[i*3+1], m_nodes[i*3+2] );
  }

  for ( int i=0; i < m_num_faces; ++i ) {
      mesh->insertCell( &m_face_connectivity[ i*3],meshtk::LINEAR_TRIANGLE,3);
  }

}

} /* namespace quest */

//------------------------------------------------------------------------------
//      STL IO Routine Implementation
//------------------------------------------------------------------------------
namespace
{
char ch_cap ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_CAP capitalizes a single character.
//
//  Discussion:
//
//    This routine should be equivalent to the library "toupper" function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    19 July 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the character to capitalize.
//
//    Output, char CH_CAP, the capitalized character.
//
{
  if ( 97 <= c && c <= 122 )
  {
    c = c - 32;
  }

  return c;
}
//****************************************************************************80

bool ch_eqi ( char c1, char c2 )

//****************************************************************************80
//
//  Purpose:
//
//    CH_EQI is true if two characters are equal, disregarding case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C1, C2, the characters to compare.
//
//    Output, bool CH_EQI, is true if the two characters are equal,
//    disregarding case.
//
{
  if ( 97 <= c1 && c1 <= 122 )
  {
    c1 = c1 - 32;
  }
  if ( 97 <= c2 && c2 <= 122 )
  {
    c2 = c2 - 32;
  }

  return ( c1 == c2 );
}
//****************************************************************************80

int ch_to_digit ( char c )

//****************************************************************************80
//
//  Purpose:
//
//    CH_TO_DIGIT returns the integer value of a base 10 digit.
//
//  Example:
//
//     C   DIGIT
//    ---  -----
//    '0'    0
//    '1'    1
//    ...  ...
//    '9'    9
//    ' '    0
//    'X'   -1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 June 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, char C, the decimal digit, '0' through '9' or blank are legal.
//
//    Output, int CH_TO_DIGIT, the corresponding integer value.  If C was
//    'illegal', then DIGIT is -1.
//
{
  int digit;

  if ( '0' <= c && c <= '9' )
  {
    digit = c - '0';
  }
  else if ( c == ' ' )
  {
    digit = 0;
  }
  else
  {
    digit = -1;
  }

  return digit;
}
//****************************************************************************


//****************************************************************************80

bool s_eqi ( string s1, string s2 )

//****************************************************************************80
//
//  Purpose:
//
//    S_EQI reports whether two strings are equal, ignoring case.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S1, S2, two strings.
//
//    Output, bool S_EQI, is true if the strings are equal.
//
{
  int i;
  int nchar;
  int s1_length;
  int s2_length;

  s1_length = s1.length ( );
  s2_length = s2.length ( );

  if ( s1_length < s2_length )
  {
    nchar = s1_length;
  }
  else
  {
    nchar = s2_length;
  }
//
//  The strings are not equal if they differ over their common length.
//
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return false;
    }
  }
//
//  The strings are not equal if the longer one includes nonblanks
//  in the tail.
//
  if ( nchar < s1_length )
  {
    for ( i = nchar; i < s1_length; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return false;
      }
    }
  }
  else if ( nchar < s2_length )
  {
    for ( i = nchar; i < s2_length; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return false;
      }
    }
  }

  return true;
}
//****************************************************************************80

int s_len_trim ( string s )

//****************************************************************************80
//
//  Purpose:
//
//    S_LEN_TRIM returns the length of a string to the last nonblank.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string.
//
//    Output, int S_LEN_TRIM, the length of the string to the last nonblank.
//    If S_LEN_TRIM is 0, then the string is entirely blank.
//
{
  int n;

  n = s.length ( );

  while ( 0 < n )
  {
    if ( s[n-1] != ' ' )
    {
      return n;
    }
    n = n - 1;
  }

  return n;
}
//****************************************************************************80

double s_to_r8 ( string s, int *lchar, bool *error )

//****************************************************************************80
//
//  Purpose:
//
//    S_TO_R8 reads an R8 from a string.
//
//  Discussion:
//
//    This routine will read as many characters as possible until it reaches
//    the end of the string, or encounters a character which cannot be
//    part of the real number.
//
//    Legal input is:
//
//       1 blanks,
//       2 '+' or '-' sign,
//       2.5 spaces
//       3 integer part,
//       4 decimal point,
//       5 fraction part,
//       6 'E' or 'e' or 'D' or 'd', exponent marker,
//       7 exponent sign,
//       8 exponent integer part,
//       9 exponent decimal point,
//      10 exponent fraction part,
//      11 blanks,
//      12 final comma or semicolon.
//
//    with most quantities optional.
//
//  Example:
//
//    S                 R
//
//    '1'               1.0
//    '     1   '       1.0
//    '1A'              1.0
//    '12,34,56'        12.0
//    '  34 7'          34.0
//    '-1E2ABCD'        -100.0
//    '-1X2ABCD'        -1.0
//    ' 2E-1'           0.2
//    '23.45'           23.45
//    '-4.2E+2'         -420.0
//    '17d2'            1700.0
//    '-14e-2'         -0.14
//    'e2'              100.0
//    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    05 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, the string containing the
//    data to be read.  Reading will begin at position 1 and
//    terminate at the end of the string, or when no more
//    characters can be read to form a legal real.  Blanks,
//    commas, or other nonnumeric data will, in particular,
//    cause the conversion to halt.
//
//    Output, int *LCHAR, the number of characters read from
//    the string to form the number, including any terminating
//    characters such as a trailing comma or blanks.
//
//    Output, bool *ERROR, is true if an error occurred.
//
//    Output, double S_TO_R8, the real value that was read from the string.
//
{
  char c;
  int ihave;
  int isgn;
  int iterm;
  int jbot;
  int jsgn;
  int jtop;
  int nchar;
  int ndig;
  double r;
  double rbot;
  double rexp;
  double rtop;
  char TAB = 9;

  nchar = s_len_trim ( s );
  *error = false;
  r = 0.0;
  *lchar = -1;
  isgn = 1;
  rtop = 0.0;
  rbot = 1.0;
  jsgn = 1;
  jtop = 0;
  jbot = 1;
  ihave = 1;
  iterm = 0;

  for ( ; ; )
  {
    c = s[*lchar+1];
    *lchar = *lchar + 1;
//
//  Blank or TAB character.
//
    if ( c == ' ' || c == TAB )
    {
      if ( ihave == 2 )
      {
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        iterm = 1;
      }
      else if ( 1 < ihave )
      {
        ihave = 11;
      }
    }
//
//  Comma.
//
    else if ( c == ',' || c == ';' )
    {
      if ( ihave != 1 )
      {
        iterm = 1;
        ihave = 12;
        *lchar = *lchar + 1;
      }
    }
//
//  Minus sign.
//
    else if ( c == '-' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
        isgn = -1;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
        jsgn = -1;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Plus sign.
//
    else if ( c == '+' )
    {
      if ( ihave == 1 )
      {
        ihave = 2;
      }
      else if ( ihave == 6 )
      {
        ihave = 7;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Decimal point.
//
    else if ( c == '.' )
    {
      if ( ihave < 4 )
      {
        ihave = 4;
      }
      else if ( 6 <= ihave && ihave <= 8 )
      {
        ihave = 9;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Exponent marker.
//
    else if ( ch_eqi ( c, 'E' ) || ch_eqi ( c, 'D' ) )
    {
      if ( ihave < 6 )
      {
        ihave = 6;
      }
      else
      {
        iterm = 1;
      }
    }
//
//  Digit.
//
    else if ( ihave < 11 && '0' <= c && c <= '9' )
    {
      if ( ihave <= 2 )
      {
        ihave = 3;
      }
      else if ( ihave == 4 )
      {
        ihave = 5;
      }
      else if ( ihave == 6 || ihave == 7 )
      {
        ihave = 8;
      }
      else if ( ihave == 9 )
      {
        ihave = 10;
      }

      ndig = ch_to_digit ( c );

      if ( ihave == 3 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
      }
      else if ( ihave == 5 )
      {
        rtop = 10.0 * rtop + ( double ) ndig;
        rbot = 10.0 * rbot;
      }
      else if ( ihave == 8 )
      {
        jtop = 10 * jtop + ndig;
      }
      else if ( ihave == 10 )
      {
        jtop = 10 * jtop + ndig;
        jbot = 10 * jbot;
      }

    }
//
//  Anything else is regarded as a terminator.
//
    else
    {
      iterm = 1;
    }
//
//  If we haven't seen a terminator, and we haven't examined the
//  entire string, go get the next character.
//
    if ( iterm == 1 || nchar <= *lchar + 1 )
    {
      break;
    }

  }
//
//  If we haven't seen a terminator, and we have examined the
//  entire string, then we're done, and LCHAR is equal to NCHAR.
//
  if ( iterm != 1 && (*lchar) + 1 == nchar )
  {
    *lchar = nchar;
  }
//
//  Number seems to have terminated.  Have we got a legal number?
//  Not if we terminated in states 1, 2, 6 or 7!
//
  if ( ihave == 1 || ihave == 2 || ihave == 6 || ihave == 7 )
  {
    *error = true;
    return r;
  }
//
//  Number seems OK.  Form it.
//
  if ( jtop == 0 )
  {
    rexp = 1.0;
  }
  else
  {
    if ( jbot == 1 )
    {
      rexp = pow ( 10.0, jsgn * jtop );
    }
    else
    {
      rexp = jsgn * jtop;
      rexp = rexp / jbot;
      rexp = pow ( 10.0, rexp );
    }

  }

  r = isgn * rexp * rtop / rbot;

  return r;
}
//****************************************************************************80


//****************************************************************************80

int stla_offset_get ( void )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_OFFSET_GET gets the STLA offset.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    22 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    3D Systems, Inc,
//    Stereolithography Interface Specification,
//    October 1989.
//
//  Parameters:
//
//    Output, int STLA_OFFSET_GET, the current value of the STLA offset.
//    This should only be 0 or 1.
//
{
  return stla_offset_value;
}

//****************************************************************************80

bool stla_read ( string input_file_name, int node_num, int face_num,
  double node_xy[], int face_node[], double face_normal[] )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_READ reads graphics information from an ASCII StereoLithography file.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    23 September 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    3D Systems, Inc,
//    Stereolithography Interface Specification,
//    October 1989.
//
//  Parameters:
//
//    Input, string INPUT_FILE_NAME, the name of the input file.
//
//    Input, int NODE_NUM, the number of vertices defined.
//
//    Input, int FACE_NUM, the number of faces defined.
//
//    Output, double NODE_XY[3*NODE_NUM], the coordinates of points.
//
//    Output, int FACE_NODE[3*FACE_NUM], the nodes that make up each face.
//
//    Output, double FACE_NORMAL[3*FACE_NUM], the normal vector
//    at each face.
//
//    Output, bool STLA_READ, is TRUE if an error occurred.
//
{
  bool done;
  double dval;
  bool error;
  int face;
  int i;
  int ierror;
  ifstream input;
  int lchar;
  int node;
  int offset;
  int state;
  double temp[3];
  string text;
  int text_num;
  int vertex;
  string word1;
  string word2;

  error = false;
  state = 0;
  offset = stla_offset_get ( );
  text_num = 0;

  face = 0;
  node = 0;
//
//  Open the file.
//
  input.open ( input_file_name.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "STLA_READ - Fatal error!\n";
    cout << "  Could not open the file \"" << input_file_name << "\".\n";
    error = true;
    return error;
  }
//
//  Read the next line of text.
//
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      if ( state != 0 && state != 1 )
      {
        cout << "\n";
        cout << "STLA_READ - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  End-of-file, but model not finished.\n";
        error = true;
        return error;
      }
      break;
    }

    text_num = text_num + 1;

    done = true;
//
//  Read the first word in the line.
//
    word1 = word_next_read ( text, &done );

    if ( done )
    {
      cout << "\n";
      cout << "STLA_READ - Fatal error!\n";
      cout << "  File line number = " << text_num << "\n";
      cout << "  No information on line.\n";
      error = true;
      return error;
    }
//
//  "Doctor" the text, changing a beginning occurrence of:
//
//      END FACET to ENDFACET
//      END LOOP to ENDLOOP
//      END SOLID to ENDSOLID
//      FACET NORMAL to FACETNORMAL
//      OUTER LOOP to OUTERLOOP
//
    if ( s_eqi ( word1, "END" ) )
    {
      word2 = word_next_read ( text, &done );

      if ( !s_eqi ( word2, "FACET" ) &&
           !s_eqi ( word2, "LOOP" ) &&
           !s_eqi ( word2, "SOLID" ) )
      {
        cout << "\n";
        cout << "STLA_READ - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  The tag END was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"FACET\", \"LOOP\", or \"SOLID\".\n";
        error = true;
        return error;
      }

      word1 = word1 + word2;
    }
    else if ( s_eqi ( word1, "FACET" ) )
    {
      word2 = word_next_read ( text, &done );

      if ( !s_eqi ( word2, "NORMAL" ) )
      {
        cout << "\n";
        cout << "STLA_READ - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  The tag FACET was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"NORMAL\".\n";
        error = true;
        return error;
      }
      word1 = word1 + word2;
    }
    else if ( s_eqi ( word1, "OUTER" ) )
    {
      word2 = word_next_read ( text, &done );

      if ( !s_eqi ( word2, "LOOP" ) )
      {
        cout << "\n";
        cout << "STLA_READ - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  The tag OUTER was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"LOOP\".\n";
        error = true;
        return error;
      }
      word1 = word1 + word2;
    }
//
//  This first word tells us what to do.
//
//  SOLID - begin a new solid.
//    Valid in state 0, moves to state 1.
//  ENDSOLID - end current solid.
//    Valid in state 1, moves to state 0.
//
//  FACETNORMAL - begin a new facet.
//    Valid in state 0 or 1, moves to state 2.
//  ENDFACET - end current facet.
//    Valid in state 2, moves to state 1.
//
//  OUTERLOOP - begin a list of vertices.
//    Valid in state 2, moves to state 3.
//  ENDLOOP - end vertex list.
//    Valid in state 3, moves to state 2.
//
//  VERTEX - give coordinates of next vertex.
//    Valid in state 3 if current vertex count is 0, 1 or 2.
//
//  End of file -
//    Valid in state 0 or 1.
//
    if ( s_eqi ( word1, "SOLID" ) )
    {
      if ( state != 0 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for SOLID.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }
      state = 1;
    }
    else if ( s_eqi ( word1, "ENDSOLID" ) )
    {
      if ( state != 1 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for ENDSOLID.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }
      state = 0;
    }
    else if ( s_eqi ( word1, "FACETNORMAL" ) )
    {
      if ( state != 0 && state != 1 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for FACET.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }

      state = 2;

      if ( face_num <= face )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  More faces being read than expected.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }

      for ( i = 0; i < 3; i++ )
      {
        face_normal[i+face*3] = 0.0;
        word2 = word_next_read ( text, &done );
        if ( !done )
        {
          dval = s_to_r8 ( word2, &lchar, &error );
          if ( error )
          {
            return error;
          }
          face_normal[i+face*3] = dval;
        }
      }
    }
    else if ( s_eqi ( word1, "ENDFACET" ) )
    {
      if ( state != 2 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for ENDFACET.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }
      face = face + 1;
      state = 1;
    }
    else if ( s_eqi ( word1, "OUTERLOOP" ) )
    {
      if ( state != 2 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for OUTERLOOP.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }

      state = 3;
      vertex = 0;
    }
    else if ( s_eqi ( word1, "ENDLOOP" ) )
    {
      if ( state != 3 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for ENDLOOP.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }
      state = 2;
    }
    else if ( s_eqi ( word1, "VERTEX" ) )
    {
      if ( state != 3 )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Model not in right state for VERTEX.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }

      if ( 3 <= vertex )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  Too many vertices for face.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }

      for ( i = 0; i < 3; i++ )
      {
        word2 = word_next_read ( text, &done );
        if ( done )
        {
          error = true;
          return error;
        }
        dval = s_to_r8 ( word2, &lchar, &error );
        if ( error )
        {
          return error;
        }
        temp[i] = dval;
      }

      if ( node_num <= node )
      {
        cout << "\n";
        cout << "STLA_READ - Warning!\n";
        cout << "  More nodes being read than expected.\n";
        cout << "  File line number = " << text_num << "\n";
        error = true;
        return error;
      }

      for ( i = 0; i < 3; i++ )
      {
        node_xy[i+node*3] = temp[i];
      }
      face_node[vertex+face*3] = node + offset;

      node = node + 1;
      vertex = vertex + 1;
    }
    else
    {
      cout << "\n";
      cout << "STLA_READ - Warning!\n";
      cout << "  Unrecognized line in file.\n";
      cout << "  File line number = " << text_num << "\n";
      error = true;
      return error;
    }
  }
//
//  Close the file.
//
  input.close ( );

  /* silence warnings about unused variables */
  static_cast< void >( ierror );

  return error;
}
//****************************************************************************80

void stla_size ( string input_file_name, int *solid_num, int *node_num,
  int *face_num, int *text_num )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_SIZE determines sizes associated with an STLA file.
//
//  Discussion:
//
//    This routine assumes that the file is a legal STLA file.
//
//    To perform checks on the file, call STLA_CHECK first.
//
//    Note that the counts for the number of nodes and edges are
//    overestimates, since presumably, most nodes will be defined several
//    times, once for each face they are part of, and most edges will
//    be defined twice.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 February 2007
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    3D Systems, Inc,
//    Stereolithography Interface Specification,
//    October 1989.
//
//  Parameters:
//
//    Input, string INPUT_FILE_NAME, the name of the input file.
//
//    Output, int *SOLID_NUM, the number of solids defined.
//    Presumably, this is 1.
//
//    Output, int *NODE_NUM, the number of vertices defined.
//
//    Output, int *FACE_NUM, the number of faces defined.
//
//    Output, int *TEXT_NUM, the number of lines of text.
//
{
  bool done;
  double dval;
  bool error;
  int i;
  int ierror;
  ifstream input;
  int lchar;
  int state;
  string text;
  int vertex;
  string word1;
  string word2;

  ierror = 0;

  state = 0;

  *text_num = 0;
  *solid_num = 0;
  *node_num = 0;
  *face_num = 0;
//
//  Open the file.
//
  input.open ( input_file_name.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "STLA_SIZE - Fatal error!\n";
    cout << "  Could not open the file \"" << input_file_name << "\".\n";
    return;
  }
//
//  Read the next line of text.
//
  for ( ; ; )
  {
    getline ( input, text );

    if ( input.eof ( ) )
    {
      if ( state != 0 &&
           state != 1 )
      {
        cout << "\n";
        cout << "STLA_SIZE - Fatal error!\n";
        cout << "  File line number = " << *text_num << "\n";
        cout << "  End-of-file, but model not finished.\n";
        return;
      }
      break;
    }

    *text_num = *text_num + 1;

    done = true;
//
//  Read the first word in the line.
//
    word1 = word_next_read ( text, &done );

    if ( done )
    {
      cout << "\n";
      cout << "STLA_CHECK - Fatal error!\n";
      cout << "  File line number = " << *text_num << "\n";
      cout << "  No information on line.\n";
      return;
    }
//
//  "Doctor" the text, changing a beginning occurrence of:
//
//      END FACET to ENDFACET
//      END LOOP to ENDLOOP
//      END SOLID to ENDSOLID
//      FACET NORMAL to FACETNORMAL
//      OUTER LOOP to OUTERLOOP
//
    if ( s_eqi ( word1, "END" ) )
    {
      word2 = word_next_read ( text, &done );

      if ( !s_eqi ( word2, "FACET" ) &&
           !s_eqi ( word2, "LOOP" ) &&
           !s_eqi ( word2, "SOLID" ) )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << *text_num << "\n";
        cout << "  The tag END was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"FACET\", \"LOOP\", or \"SOLID\".\n";
        return;
      }
      word1 = word1 + word2;
    }
    else if ( s_eqi ( word1, "FACET" ) )
    {
      word2 = word_next_read ( text, &done );

      if ( !s_eqi ( word2, "NORMAL" ) )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << *text_num << "\n";
        cout << "  The tag FACET was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"NORMAL\".\n";
        return;
      }
      word1 = word1 + word2;
    }
    else if ( s_eqi ( word1, "OUTER" ) )
    {
      word2 = word_next_read ( text, &done );

      if ( !s_eqi ( word2, "LOOP" ) )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << *text_num << "\n";
        cout << "  The tag OUTER was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"LOOP\".\n";
        return;
      }
      word1 = word1 + word2;
    }
//
//  This first word tells us what to do.
//
//  SOLID - begin a new solid.
//    Valid in state 0, moves to state 1.
//  ENDSOLID - end current solid.
//    Valid in state 1, moves to state 0.
//
//  FACETNORMAL - begin a new facet.
//    Valid in state 0 or 1, moves to state 2.
//  ENDFACET - end current facet.
//    Valid in state 2, moves to state 1.
//
//  OUTERLOOP - begin a list of vertices.
//    Valid in state 2, moves to state 3.
//  ENDLOOP - end vertex list.
//    Valid in state 3, moves to state 2.
//
//  VERTEX - give coordinates of next vertex.
//    Valid in state 3 if current vertex count is 0, 1 or 2.
//
//  End of file -
//    Valid in state 0 or 1.
//
    if ( s_eqi ( word1, "SOLID" ) )
    {
      if ( state != 0 )
      {
        return;
      }
      state = 1;
    }
    else if ( s_eqi ( word1, "ENDSOLID" ) )
    {
      if ( state != 1 )
      {
        return;
      }
      state = 0;

      *solid_num = *solid_num + 1;
    }
    else if ( s_eqi ( word1, "FACETNORMAL" ) )
    {
      if ( state != 0 && state != 1 )
      {
        return;
      }
      state = 2;

      for ( i = 1; i <= 3; i++ )
      {
        word2 = word_next_read ( text, &done );

        if ( done )
        {
          return;
        }

        dval = s_to_r8 ( word2, &lchar, &error );

        if ( error )
        {
          return;
        }
      }
    }
    else if ( s_eqi ( word1, "ENDFACET" ) )
    {
      if ( state != 2 )
      {
        return;
      }
      state = 1;
      *face_num = *face_num + 1;
    }
    else if ( s_eqi ( word1, "OUTERLOOP" ) )
    {
      if ( state != 2 )
      {
        return;
      }
      state = 3;
      vertex = 0;
    }
    else if ( s_eqi ( word1, "ENDLOOP" ) )
    {
      if ( state != 3 )
      {
        return;
      }
      state = 2;
    }
    else if ( s_eqi ( word1, "VERTEX" ) )
    {
      if ( state != 3 )
      {
        return;
      }

      if ( 3 <= vertex )
      {
        return;
      }

      for ( i = 1; i <= 3; i++ )
      {
        word2 = word_next_read ( text, &done );

        if ( done )
        {
          return;
        }

        dval = s_to_r8 ( word2, &lchar, &error );

        if ( error )
        {
          return;
        }

      }
      vertex = vertex + 1;
      *node_num = *node_num + 1;
    }
    else
    {
      return;
    }

  }
//
//  Close the file.
//
  input.close ( );

  /* silence warnings about unused variables */
  static_cast< void >( dval );
  static_cast< void >( ierror );

  return;
}
//****************************************************************************80

string word_next_read ( string s, bool *done )

//****************************************************************************80
//
//  Purpose:
//
//    WORD_NEXT_READ "reads" words from a string, one at a time.
//
//  Discussion:
//
//    This routine was written to process tokens in a file.
//    A token is considered to be an alphanumeric string delimited
//    by whitespace, or any of various "brackets".
//
//    The following characters are considered to be a single word,
//    whether surrounded by spaces or not:
//
//      " ( ) { } [ ]
//
//    Also, if there is a trailing comma on the word, it is stripped off.
//    This is to facilitate the reading of lists.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 October 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, string S, a string, presumably containing words
//    separated by spaces.
//
//    Input/output, bool *DONE.
//    On input with a fresh string, set DONE to TRUE.
//    On output, the routine sets DONE:
//      FALSE if another word was read,
//      TRUE if no more words could be read.
//
//    Output, string WORD_NEXT_READ.
//    If DONE is FALSE, then WORD contains the "next" word read.
//    If DONE is TRUE, then WORD is NULL, because there was no more to read.
//
{
  int i;
  int ilo;
  int j;
  static int lenc = 0;
  static int next = 0;
  char TAB = 9;
  string word;
  char *word_chstar;
//
//  We "remember" LENC and NEXT from the previous call.
//
//  An input value of DONE = TRUE signals a new line of text to examine.
//
  if ( *done )
  {
    next = 0;
    *done = false;
    lenc = s.length ( );
    if ( lenc <= 0 )
    {
      *done = true;
      word = "\n";;
      return word;
    }
  }
//
//  Beginning at index NEXT, search the string for the next nonblank,
//  which signals the beginning of a word.
//
  ilo = next;
//
//  ...S(NEXT:) is blank.  Return with WORD = ' ' and DONE = TRUE.
//
  for ( ; ; )
  {
    if ( lenc < ilo )
    {
      word = "\n";
      *done = true;
      next = lenc + 1;
      return word;
    }
//
//  If the current character is blank, skip to the next one.
//
    if ( s[ilo] != ' ' && s[ilo] != TAB )
    {
      break;
    }
    ilo = ilo + 1;
  }
//
//  ILO is the index of the next nonblank character in the string.
//
//  If this initial nonblank is a special character,
//  then that's the whole word as far as we're concerned,
//  so return immediately.
//
  if ( s[ilo] == '"' )
  {
    word = """";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '(' )
  {
    word = "(";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == ')' )
  {
    word = ")";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '{' )
  {
    word = "{";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '}' )
  {
    word = "}";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == '[' )
  {
    word = "[";
    next = ilo + 1;
    return word;
  }
  else if ( s[ilo] == ']' )
  {
    word = "]";
    next = ilo + 1;
    return word;
  }
//
//  Now search for the last contiguous character that is not a
//  blank, TAB, or special character.
//
  next = ilo + 1;

  while ( next <= lenc )
  {
    if ( s[next] == ' ' )
    {
      break;
    }
    else if ( s[next] == TAB )
    {
      break;
    }
    else if ( s[next] == '"' )
    {
      break;
    }
    else if ( s[next] == '(' )
    {
      break;
    }
    else if ( s[next] == ')' )
    {
      break;
    }
    else if ( s[next] == '{' )
    {
      break;
    }
    else if ( s[next] == '}' )
    {
      break;
    }
    else if ( s[next] == '[' )
    {
      break;
    }
    else if ( s[next] == ']' )
    {
      break;
    }

    next = next + 1;
  }
//
//  Allocate WORD, copy characters, and return.
//
  if ( s[next-1] == ',' )
  {
    word_chstar = new char[next-ilo];
    i = 0;
    for ( j = ilo; j <= next - 2; j++ )
    {
      word_chstar[i] = s[j];
      i = i + 1;
    }
    word_chstar[i] = '\0';
    word = string ( word_chstar );
    delete [] word_chstar;
  }
  else
  {
    word_chstar = new char[next+1-ilo];
    i = 0;
    for ( j = ilo; j <= next-1; j++ )
    {
      word_chstar[i] = s[j];
      i = i + 1;
    }
    word_chstar[i] = '\0';
    word = string ( word_chstar );
    delete [] word_chstar;
  }

  return word;
}

}

