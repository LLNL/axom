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

#ifndef MINT_CellTypes_HXX_
#define MINT_CellTypes_HXX_

#include "mint/config.hpp"

namespace axom
{
namespace mint
{

constexpr IndexType MAX_NUM_NODES = 27;
constexpr IndexType NEED_INDIRECTION = -1;

enum class CellTypes : signed char
{
  UNDEFINED_ELEMENT = -1, ///< UNDEFINED

  VERTEX,                 ///< VERTEX
  SEGMENT,                ///< LINE_SEGMENT

  TRIANGLE,               ///< LINEAR_TRIANGLE
  QUAD,                   ///< LINEAR_QUAD
  TET,                    ///< LINEAR_TET
  HEX,                    ///< LINEAR_HEX
  PRISM,                  ///< LINEAR_PRISM
  PYRAMID,                ///< LINEAR_PYRAMID

  QUAD9,                  ///< QUADRATIC QUAD
  HEX27,                  ///< QUADRATIC HEX

  MIXED,                  ///< MIXED
};

template< CellTypes TYPE > struct cell_info 
{
  static constexpr bool valid = false;  
};

template<> struct cell_info< CellTypes::VERTEX >
{
  static constexpr const char* name = "VERTEX";
  static constexpr const char* blueprint_name = "point";
  static constexpr int vtk_type = 1;
  static constexpr int num_nodes = 1;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::SEGMENT >
{
  static constexpr const char* name = "SEGMENT";
  static constexpr const char* blueprint_name = "LINE";
  static constexpr int vtk_type = 3;
  static constexpr int num_nodes = 2;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::TRIANGLE >
{
  static constexpr const char* name = "TRIANGLE";
  static constexpr const char* blueprint_name = "tri";
  static constexpr int vtk_type = 5;
  static constexpr int num_nodes = 3;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::QUAD >
{
  static constexpr const char* name = "QUAD";
  static constexpr const char* blueprint_name = "quad";
  static constexpr int vtk_type = 9;
  static constexpr int num_nodes = 4;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::TET >
{
  static constexpr const char* name = "TET";
  static constexpr const char* blueprint_name = "tet";
  static constexpr int vtk_type = 10;
  static constexpr int num_nodes = 4;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::HEX >
{
  static constexpr const char* name = "HEX";
  static constexpr const char* blueprint_name = "hex";
  static constexpr int vtk_type = 12;
  static constexpr int num_nodes = 8;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::PRISM >
{
  static constexpr const char* name = "PRISM";
  static constexpr const char* blueprint_name = "prism_non_bp";
  static constexpr int vtk_type = 13;
  static constexpr int num_nodes = 6;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::PYRAMID >
{
  static constexpr const char* name = "PYRAMID";
  static constexpr const char* blueprint_name = "pyramid_non_bp";
  static constexpr int vtk_type = 14;
  static constexpr int num_nodes = 5;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::QUAD9 >
{
  static constexpr const char* name = "QUAD9";
  static constexpr const char* blueprint_name = "quad9_non_bp";
  static constexpr int vtk_type = 28;
  static constexpr int num_nodes = 9;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::HEX27 >
{
  static constexpr const char* name = "HEX27";
  static constexpr const char* blueprint_name = "HEX27_non_bp";
  static constexpr int vtk_type = 29;
  static constexpr int num_nodes = 27;
  static constexpr bool valid = true;  
};

template<> struct cell_info< CellTypes::MIXED >
{
  static constexpr const char* name = "MIXED";
  static constexpr const char* blueprint_name = "mixed_non_bp";
  static constexpr int vtk_type = 0;
  static constexpr int num_nodes = NEED_INDIRECTION;
  static constexpr bool valid = true;  
};


} /* namespace mint */
} /* namespace axom */

#endif /* MINT_CellTypes_HXX_ */
