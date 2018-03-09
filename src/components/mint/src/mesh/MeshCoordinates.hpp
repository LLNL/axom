/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

#ifndef MESHCOORDINATES_HXX_
#define MESHCOORDINATES_HXX_

#include "axom/Macros.hpp"

// C/C++ includes
#include <vector> // for STL vector

#define X_COORDINATE 0
#define Y_COORDINATE 1
#define Z_COORDINATE 2

namespace axom
{
namespace mint
{

/*!
 * \class MeshCoordinates
 *
 * \brief Provides functionality for storing mesh coordinates.
 */
class MeshCoordinates
{
public:

  /*!
   * \brief Default constructor.
   */
  MeshCoordinates();

  /*!
   * \brief Custom constructor.
   * \param [in] dimension the dimension of the ambient space.
   * \pre dimension >= 1 && dimension <= 3
   */
  explicit MeshCoordinates(int dimension);

  /*!
   * \brief Creates a MeshCoordinates instance with the given number of points.
   * \param [in] dimension the dimension of the ambient space.
   * \param [in] npoints the total number of points to allocate space for.
   * \note The points are not initialized but, can be set using setPoint
   */
  MeshCoordinates( int dimension, int npoints );

  /*!
   * \brief Creates a MeshCoordinates instance where the coordinates are given
   *  by the cartesian product of the coordinates in each dimension.
   * \param [in] dimension the dimension of the ambient space.
   * \param [in] ndims number of coordinate values along each dimension.
   * \pre dimension >= 1 && dimension <= 3
   * \pre ndims != AXOM_NULLPTR.
   */
  MeshCoordinates( int dimension, int ndims[3] );

  /*!
   * \brief Destructor.
   */
  virtual ~MeshCoordinates();

  /*!
   * \brief Inserts a new point into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \pre m_ndims == 1.
   */
  void insertPoint( double x );

  /*!
   * \brief Inserts a new point into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] y the y--coordinate.
   * \pre m_ndims == 2.
   */
  void insertPoint( double x, double y );

  /*!
   * \brief Inserts a new point into this MeshCoordinates instance.
   * \param [in] x the x--coordinate.
   * \param [in] y the y--coordinate.
   * \param [in] z the z--coordinate.
   * \pre m_ndims == 3.
   */
  void insertPoint( double x, double y, double z );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 1.
   */
  void setPoint( int pntIdx, double x );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \param [in] y the y--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 2.
   */
  void setPoint( int pntIdx, double x, double y );

  /*!
   * \brief Sets the point at the supplied index to the given coordinates.
   * \param [in] pntIdx the index of the point to set
   * \param [in] x the x--coordinate to set
   * \param [in] y the y--coordinate to set
   * \param [in] z the z--coordinate to set
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   * \pre m_ndims == 3.
   */
  void setPoint( int pntIdx, double x, double y, double z );

  /*!
   * \brief Returns the coordinate of a point at the given dimension.
   * \param [in] pntIdx the index of the point in query.
   * \param [in] dim the dimension in query.
   * \return coord the coordinate of the point.
   * \pre dim < m_ndims
   * \pre (pntIdx >= 0) && (pntIdx < this->getNumberOfPoints())
   */
  double getCoordinate( int pntIdx, int dim );

  /*!
   * \brief Returns a const pointer to the coordinate array.
   * \param [in] dim the requested dimension.
   * \return coord_array const pointer to the coordinate array
   * \pre dim < m_ndims
   * \post coord_array != AXOM_NULLPTR.
   */
  double* getCoordinateArray( int dim );

  /*!
   * \brief Number of coordinates along the given dimension.
   * \param [in] idim the dimension in query.
   * \return N the number of coordinates along the given dimension.
   * \note Used for RectilinearMesh
   */
  int getCoordinateArraySize( int idim ) const
  { return static_cast< int >( m_coordinates[ idim ].size() ); };

  /*!
   * \brief Returns the number of points in this MeshCoordinates instance.
   * \return npoint the number points in this MeshCoordinates instance.
   */
  int getNumberOfPoints() const;

private:
  // TODO: support different memory layouts...
  // TODO: replace with flat array, just using vector for now...
  typedef std::vector< std::vector< double > > Coordinates;

  int m_ndims;
  Coordinates m_coordinates;

  /*!
   * \brief Helper method used to initialize the MeshCoordinates data-structure.
   * \note Called only from the constructor.
   */
  void initialize();

  /*!
   * \brief Helper method used to initialize the MeshCoordinate data-structure.
   *  It creates a MeshCoordinates instance with the given number of points.
   * \param [in] npoints the number of points.
   *
   * \overload
   */
  void initialize(int npoints);

  DISABLE_COPY_AND_ASSIGNMENT(MeshCoordinates);
  DISABLE_MOVE_AND_ASSIGNMENT(MeshCoordinates);
};

} /* namespace mint */
} /* namespace axom */

#endif /* MESHCOORDINATES_HXX_ */
