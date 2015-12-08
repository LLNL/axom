# include <cstdlib>
# include <cmath>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>

using namespace std;

# include "stla_io.hpp"

//****************************************************************************80

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

int i4_min ( int i1, int i2 )

//****************************************************************************
//
//  Purpose:
//
//    I4_MIN returns the smaller of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two I4's to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//
//
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
//*********************************************************************

double r8_max ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//
{
  if ( y < x )
  {
    return x;
  } 
  else
  {
    return y;
  }
}
//*********************************************************************

double r8_min ( double x, double y )

//*********************************************************************
//
//  Purpose:
//
//    R8_MIN returns the minimum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MIN, the minimum of X and Y.
//
{
  if ( y < x )
  {
    return y;
  } 
  else
  {
    return x;
  }
}
//********************************************************************

double *r8vec_cross_3d ( double v1[3], double v2[3] )

//********************************************************************
//
//  Purpose:
//
//    R8VEC_CROSS_3D computes the cross product of two R8VEC's in 3D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    07 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double V1[3], V2[3], the coordinates of the vectors.
//
//    Output, double R8VEC_CROSS_3D[3], the cross product vector.
//
{
  double *v3;

  v3 = new double[3];

  v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
  v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
  v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

  return v3;
}
//********************************************************************

double r8vec_dot ( int n, double a1[], double a2[] )

//********************************************************************
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of two R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }

  return value;
}
//****************************************************************************80

double r8vec_length ( int dim_num, double x[] )

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_LENGTH returns the Euclidean length of the vector X
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the spatial dimension.
//
//    Input, double X[DIM_NUM], the vector.
//
//    Output, double R8VEC_LENGTH, the Euclidean length of the vector.
//
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < dim_num; i++ )
  {
    value = value + pow ( x[i], 2 );
  }
  value = sqrt ( value );

  return value;
}
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

bool stla_check ( string input_file_name )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_CHECK checks an ASCII StereoLithography file.
//
//  Example:
//
//    solid MYSOLID
//      facet normal 0.4 0.4 0.2
//        outerloop
//          vertex  1.0 2.1 3.2
//          vertex  2.1 3.7 4.5
//          vertex  3.1 4.5 6.7
//        end loop
//      end facet
//      ...
//      facet normal 0.2 0.2 0.4
//        outerloop
//          vertex  2.0 2.3 3.4
//          vertex  3.1 3.2 6.5
//          vertex  4.1 5.5 9.0
//        end loop
//      end facet
//    end solid MYSOLID
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
//    Input, string INPUT_FILE_NAME, the name of the input file.
//
//    Output, bool STLA_CHECK, is TRUE if the file is legal.
//
{
  bool check;
  bool done;
  double dval;
  bool error;
  int i;
  int ierror;
  ifstream input;
  int lchar;
  int state;
  string text;
  int text_num;
  int vertex;
  string word1;
  string word2;

  state = 0;
  text_num = 0;
//
//  Open the file.
//
  input.open ( input_file_name.c_str ( ) );

  if ( !input )
  {
    cout << "\n";
    cout << "STLA_CHECK - Fatal error!\n";
    cout << "  Could not open the file \"" << input_file_name << "\".\n";
    check = false;
    return check;
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
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  End-of-file, but model not finished.\n";
        check = false;
        return check;
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
      cout << "STLA_CHECK - Fatal error!\n";
      cout << "  File line number = " << text_num << "\n";
      cout << "  No information on line.\n";
      check = false;
      return check;
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
        cout << "  File line number = " << text_num << "\n";
        cout << "  The tag END was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"FACET\", \"LOOP\", or \"SOLID\".\n";
        check = false;
        return check;
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
        cout << "  File line number = " << text_num << "\n";
        cout << "  The tag FACET was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"NORMAL\".\n";
        check = false;
        return check;
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
        cout << "  File line number = " << text_num << "\n";
        cout << "  The tag OUTER was followed by an illegal word:\n";
        cout << "  \"" << word2 << "\"\n";
        cout << "  when expecting \"LOOP\".\n";
        check = false;
        return check;
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
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  A new SOLID statement was encountered, but we\n";
        cout << "  have not finished processing the current solid.\n";
        check = false;
        return check;
      }
      state = 1;
    }
    else if ( s_eqi ( word1, "ENDSOLID" ) )
    {
      if ( state != 1 )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  An END SOLID statement was encountered, but\n";
        cout << "  either we have not begun a solid at all, or we\n";
        cout << "  are not at an appropriate point to finish the\n";
        cout << "  current solid.\n";
        check = false;
        return check;
      }

      state = 0;
    }
    else if ( s_eqi ( word1, "FACETNORMAL" ) )
    {
      if ( state != 0 && state != 1 )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  Model not in right state for FACET.\n";
        check = false;
        return check;
      }

      state = 2;

      for ( i = 1; i <= 3; i++ )
      {
        word2 = word_next_read ( text, &done );

        if ( done )
        {
          cout << "\n";
          cout << "STLA_CHECK - Fatal error!\n";
          cout << "  File line number = " << text_num << "\n";
          cout << "  End of information while reading a component\n";
          cout << "  of the normal vector.\n";
          check = false;
          return check;
        }

        dval = s_to_r8 ( word2, &lchar, &error );

        if ( error )
        {
          cout << "\n";
          cout << "STLA_CHECK - Fatal error!\n";
          cout << "  File line number = " << text_num << "\n";
          cout << "  Error while reading a component of the normal vector.\n";
          check = false;
          return check;
        }
      }
    }
    else if ( s_eqi ( word1, "ENDFACET" ) )
    {
      if ( state != 2 )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  Model not in right state for ENDFACET.\n";
        check = false;
        return check;
      }
      state = 1;
    }
    else if ( s_eqi ( word1, "OUTERLOOP" ) )
    {
      if ( state != 2 )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  Model not in right state for OUTERLOOP.\n";
        check = false;
        return check;
      }

      state = 3;
      vertex = 0;
    }
    else if ( s_eqi ( word1, "ENDLOOP" ) )
    {
      if ( state != 3 )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  Model not in right state for ENDLOOP.\n";
        check = false;
        return check;
      }
      state = 2;
    }
    else if ( s_eqi ( word1, "VERTEX" ) )
    {
      if ( state != 3 )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  Model not in right state for VERTEX.\n";
        check = false;
        return check;
      }

      if ( 3 <= vertex )
      {
        cout << "\n";
        cout << "STLA_CHECK - Fatal error!\n";
        cout << "  File line number = " << text_num << "\n";
        cout << "  More than 3 vertices specified for a face.\n";
        check = false;
        return check;
      }

      for ( i = 1; i <= 3; i++ )
      {
        word2 = word_next_read ( text, &done );

        if ( done )
        {
          cout << "\n";
          cout << "STLA_CHECK - Fatal error!\n";
          cout << "  File line number = " << text_num << "\n";
          cout << "  The value of a vertex coordinate is missing.\n";
          check = false;
          return check;
        }

        dval = s_to_r8 ( word2, &lchar, &error );

        if ( error )
        {
          cout << "\n";
          cout << "STLA_CHECK - Fatal error!\n";
          cout << "  File line number = " << text_num << "\n";
          cout << "  The value of a vertex coordinate makes no sense.\n";
          check = false;
          return check;
        }
      }
      vertex = vertex + 1;
    }
    else
    {
      cout << "\n";
      cout << "STLA_CHECK - Fatal error!\n";
      cout << "  File line number = " << text_num << "\n";
      cout << "  Unrecognized line in file.\n";
      check = false;
      return check;
    }
  }
//
//  Close the file.
//
  input.close ( );

  check = true;

  return check;
}
//****************************************************************************80

void stla_face_node_print ( int face_num, int face_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_FACE_NODE_PRINT prints the node indices for each face.
//
//  Discussion:
//
//    If the global variable OFFSET is set to 1, then it is assumed that 
//    all indices are 1-based.  In that case, this routine will print
//    face numbers from 1 to FACE_NUM.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    21 September 2005
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
//    Input, int FACE_NUM, the number of faces.
//
//    Input, int FACE_NODE[3*FACE_NUM], the nodes that make up each face.
//
{
  int face;
  int offset;
  int vertex;

  offset = stla_offset_get ( );

  cout << "\n";
  cout << "    Face         Nodes\n";
  cout << "\n";

  for ( face = 0; face < face_num; face++ )
  {
    cout << "  " << setw(6) << ( face + offset );
    for ( vertex = 0; vertex < 3; vertex++ )
    {
      cout << "  " << setw(6) << face_node[vertex+face*3];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

double *stla_face_normal_compute ( int node_num, int face_num, double node_xyz[], 
  int face_node[] )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_FACE_NORMAL_COMPUTE computes normal vectors for an ASCII StereoLithography file.
//
//  Discussion:
//
//    This routine computes the normal vector to each triangular face
//    in the STLA solid.  If the nodes of each triangular face are
//    listed in counterclockwise order (as seen from outside the solid),
//    then the normal vectors will be properly outward facing.
//
//    The normal vectors will have unit Euclidean norm.
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
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the node coordinates.
//
//    Input, int FACE_NODE[3*FACE_NUM], the nodes making faces.
//
//    Input, int FACE_MAX, the maximum number of faces.
//
//    Output, double STLA_FACE_NORMAL_COMPUTE[3*FACE_NUM], the normal 
//    vector at each face.
//
{
  int face;
  double *face_normal;
  int i;
  int n1;
  int n2;
  int n3;
  double norm;
  int offset;
  double v1[3];
  double v2[3];
  double *v3;

  offset = stla_offset_get ( );

  face_normal = new double[3*face_num];
  
  for ( face = 0; face < face_num; face++ )
  {
    n1 = face_node[0+face*3] - offset;
    n2 = face_node[1+face*3] - offset;
    n3 = face_node[2+face*3] - offset;

    for ( i = 0; i < 3; i++ )
    {
      v1[i] = node_xyz[i+n2*3] - node_xyz[i+n1*3];
    }
    for ( i = 0; i < 3; i++ )
    {
      v2[i] = node_xyz[i+n3*3] - node_xyz[i+n1*3];
    }

    v3 = r8vec_cross_3d ( v1, v2 );

    norm = r8vec_length ( 3, v3 );

    if ( norm != 0.0 )
    {
      for ( i = 0; i < 3; i++ )
      {
        face_normal[i+face*3] = v3[i] / norm;
      }
    }
    else
    {
      for ( i = 0; i < 3; i++ )
      {
        face_normal[i+face*3] = v3[i];
      }
    }
    delete [] v3;
  }

  return face_normal;
}
//****************************************************************************80

void stla_face_normal_print ( int face_num, double face_normal[] )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_FACE_NORMAL_PRINT prints the normal vectors.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 September 2005
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
//    Input, int FACE_NUM, the number of faces.
//
//    Input, double FACE_NORMAL[3*FACE_NUM], the normal vector at each face.
//
{
  int face;
  int i;

  cout << "\n";
  cout << "    Face         Normal Vectors\n";
  cout << "\n";

  for ( face = 0; face < face_num; face++ )
  {
    cout << "  " << setw(6) << face;
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << setw(14) << face_normal[i+face*3];
    }
    cout << "\n";
  }

  return;
}
//****************************************************************************80

void stla_node_xyz_print ( int node_num, double node_xyz[] )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_NODE_XYZ_PRINT prints the node coordinates.
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
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, double NODE_XYZ[3*FACE_NUM], the normal vector at each face.
//
{
  int i;
  int node;

  cout << "\n";
  cout << "    Node         Coordinates\n";
  cout << "\n";

  for ( node = 0; node < node_num; node++ )
  {
    cout << "  " << setw(6) << node;
    for ( i = 0; i < 3; i++ )
    {
      cout << "  " << setw(14) << node_xyz[i+node*3];
    }
    cout << "\n";
  }

  return;
}
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

void stla_offset_set ( int offset )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_OFFSET_SET sets the STLA offset.
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
//    Input, int OFFSET, the new value for the STLA offset.
//    This should only be 0 or 1.
//
{
  if ( offset != 0 && offset != 1 )
  {
    cout << "\n";
    cout << "STLA_OFFSET_SET - Fatal error!\n";
    cout << "  Input values of OFFSET must be 0 or 1.\n";
    cout << "  Illegal input value was " << offset << "\n";
    exit ( 1 );
  }

  stla_offset_value = offset;

  return;
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

  return;
}
//****************************************************************************80

void stla_size_print ( string input_file_name, int solid_num, int node_num,
  int face_num, int text_num )

//****************************************************************************80
//
//  Purpose:
//
//    STLA_SIZE_PRINT prints sizes associated with an STLA file.
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
//    Input, string INPUT_FILE_NAME, the name of the input file, or the
//    name of the object.
//
//    Input, int SOLID_NUM, the number of solids defined.
//
//    Input, int NODE_NUM, the number of vertices defined.
//
//    Input, int FACE_NUM, the number of faces defined.
//
//    Input, int TEXT_NUM, the number of lines of text in the file.
//
{
  cout << "\n";
  cout << "  Sizes for STLA object \"" << input_file_name << "\".\n";
  cout << "\n";
  cout << "  Solids =                   " << solid_num << "\n";
  cout << "  Nodes (may be repeated) =  " << node_num << "\n";
  cout << "  Faces (triangular only) =  " << face_num << "\n";
  cout << "\n";
  cout << "  The index offset value =   " << stla_offset_get ( ) << "\n";
  cout << "  Number of lines of text =  " << text_num << "\n";

  return;
}
//****************************************************************************

void stla_write ( string output_file_name, int node_num, int face_num, 
  double node_xyz[], int face_node[], double face_normal[] )

//****************************************************************************
//
//  Purpose:
//   
//    STLA_WRITE writes an ASCII STL (stereolithography) file.
//
//  Example:
//
//    solid MYSOLID
//      facet normal 0.4 0.4 0.2
//        outerloop
//          vertex  1.0 2.1 3.2
//          vertex  2.1 3.7 4.5
//          vertex  3.1 4.5 6.7
//        end loop
//      end facet
//      ...
//      facet normal 0.2 0.2 0.4
//        outerloop
//          vertex  2.0 2.3 3.4
//          vertex  3.1 3.2 6.5
//          vertex  4.1 5.5 9.0
//        end loop
//      end facet
//    end solid MYSOLID
//
//  Discussion:
//
//    The polygons in an STL file should only be triangular.  This routine 
//    will try to automatically decompose higher-order polygonal faces into 
//    suitable triangles, without actually modifying the internal graphics 
//    data.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    22 September 1998
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
//    Input, string OUTPUT_FILE_NAME, the name of the output file.
//
//    Input, int NODE_NUM, the number of nodes.
//
//    Input, int FACE_NUM, the number of faces.
//
//    Input, double NODE_XYZ[3*NODE_NUM], the node coordinates.
//
//    Input, int FACE_NODE[3*FACE_NUM], the nodes making faces.
//
//    Input, int FACE_MAX, the maximum number of faces.
//
//    Input, double FACE_NORMAL[3*FACE_NUM], the normal vector at each face.
//
{
  int face;
  int i;
  int node;
  int offset;
  ofstream output_unit;
  int text_num;
  int vertex;

  offset = stla_offset_get ( );
//
//  Open the file.
//
  output_unit.open ( output_file_name.c_str ( ) );

  if ( !output_unit )
  {
    cout << "\n";
    cout << "STLA_WRITE - Fatal error!\n";
    cout << "  Could not open the file \"" << output_file_name << "\".\n";
    return;
  }
//
//  Initialize.
//
  output_unit << "solid MYSOLID\n";

  for ( face = 0; face < face_num; face++ )
  {
    output_unit << "  facet normal";
    for ( i = 0; i < 3; i++ )
    {
      output_unit << "  " << setw(10) << face_normal[i+face*3];
    }
    output_unit << "\n";
    output_unit << "    outer loop\n";
    for ( vertex = 0; vertex < 3; vertex++ )
    {
      node = face_node[vertex+face*3] - offset;
      output_unit << "      vertex  ";
      for ( i = 0; i < 3; i++ )
      {
        output_unit << "  " << setw(10) << node_xyz[i+node*3];
      }      
      output_unit << "\n";
    }
    output_unit << "    end loop\n";
    output_unit << "  end facet\n";
  }

  output_unit << "end solid MYSOLID\n";
//
//  Close the file.
//
  output_unit.close ( );

  return;
}
//****************************************************************************80

void timestamp ( void )

//****************************************************************************80
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
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
