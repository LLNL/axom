/*
 * $Id$
 */

/*!
 *******************************************************************************
 * \file
 *
 * \date Apr 8, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *
 *******************************************************************************
 */

#ifndef DETERMINANTS_H_
#define DETERMINANTS_H_

namespace quest {

namespace math {

/*!
 *******************************************************************************
 * \brief Computes a 2X2 determinant of the given matrix.
 * \tparam [in] a00 matrix element
 * \tparam [in] a01 matrix element
 * \tparam [in] a10 matrix element
 * \tparam [in] a11 matrix element
 * \return det the determinant of the 2X2 matrix
 *******************************************************************************
 */
template < class real >
inline real determinant( const real& a00, const real& a01,
                         const real& a10, const real& a11 )
{
   const real det = a00*a11 - a10*a01;
   return det;
}

/*!
 *******************************************************************************
 * \brief Computes the 3x3 determinant of the given matrix.
 * \param [in] a00 matrix element
 * \param [in] a01 matrix element
 * \param [in] a02 matrix element
 * \param [in] a10 matrix element
 * \param [in] a11 matrix element
 * \param [in] a12 matrix element
 * \param [in] a20 matrix element
 * \param [in] a21 matrix element
 * \param [in] a22 matrix element
 * \return det the determinant of the 3x3 matrix.
 *******************************************************************************
 */
template < class real >
inline real determinant( const real& a00,  const real& a01,  const real& a02,
                         const real& a10,  const real& a11,  const real& a12,
                         const real& a20,  const real& a21,  const real& a22 )
{
   const real m01 = a00*a11 - a10*a01;
   const real m02 = a00*a21 - a20*a01;
   const real m12 = a10*a21 - a20*a11;

   const real det = m01*a22 - m02*a12 + m12*a02;
   return det;
}

/*!
 *******************************************************************************
 * \brief Computes the 4x4 determinant of the given matrix.
 * \param [in] a00 matrix element
 * \param [in] a01 matrix element
 * \param [in] a02 matrix element
 * \param [in] a03 matrix element
 * \param [in] a10 matrix element
 * \param [in] a11 matrix element
 * \param [in] a12 matrix element
 * \param [in] a13 matrix element
 * \param [in] a20 matrix element
 * \param [in] a21 matrix element
 * \param [in] a22 matrix element
 * \param [in] a23 matrix element
 * \param [in] a30 matrix element
 * \param [in] a31 matrix element
 * \param [in] a32 matrix element
 * \param [in] a33 matrix element
 * \return det the determinant of the 4x4 matrix
 *******************************************************************************
 */
template < class real >
inline real determinant(
     const real& a00, const real& a01, const real& a02, const real& a03,
     const real& a10, const real& a11, const real& a12, const real& a13,
     const real& a20, const real& a21, const real& a22, const real& a23,
     const real& a30, const real& a31, const real& a32, const real& a33  )
{
   const real m01 = a10*a01 - a00*a11;
   const real m02 = a20*a01 - a00*a21;
   const real m03 = a30*a01 - a00*a31;
   const real m12 = a20*a11 - a10*a21;
   const real m13 = a30*a11 - a10*a31;
   const real m23 = a30*a21 - a20*a31;

   const real m012 = m12*a02 - m02*a12 + m01*a22;
   const real m013 = m13*a02 - m03*a12 + m01*a32;
   const real m023 = m23*a02 - m03*a22 + m02*a32;
   const real m123 = m23*a12 - m13*a22 + m12*a32;

   const real det = m123*a03 - m023*a13 + m013*a23 - m012*a33;
   return det;
}


} /* namespace math */

} /* namespace quest */

#endif /* DETERMINANTS_H_ */
