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
 * \file StructuredMesh.hxx
 *
 * \date Sep 20, 2015
 * \author George Zagaris (zagaris2@llnl.gov)
 *******************************************************************************
 */


#include "mint/StructuredMesh.hpp"

#include "mint/MeshType.hpp"


namespace mint
{


StructuredMesh::StructuredMesh() :
    Mesh(-1,MINT_UNDEFINED_MESH,-1,-1),
    m_extent( ATK_NULLPTR )
{
// TODO Auto-generated constructor stub

}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh(int meshType, int ndims, int ext[6]) :
   Mesh( ndims, meshType, 0, 0 ),
   m_extent( new Extent< int >(ndims, ext) )
{

}

//------------------------------------------------------------------------------
StructuredMesh::StructuredMesh( int meshType, int ndims, int ext[6],
                                int blockId,int partId) :
   Mesh( ndims, meshType, blockId, partId ),
   m_extent( new Extent< int >(ndims, ext) )

{

}

//------------------------------------------------------------------------------
StructuredMesh::~StructuredMesh()
{
  delete m_extent;
  m_extent = ATK_NULLPTR;
}



} /* namespace mint */
