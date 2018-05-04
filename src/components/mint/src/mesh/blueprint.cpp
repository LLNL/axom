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
#include "mint/blueprint.hpp"

// Axom includes
#include "axom/Types.hpp"      // for AXOM_NULLPTR

// Mint includes
#include "mint/config.hpp"     // for MINT_USE_SIDRE compile-time definition
#include "mint/Extent.hpp"     // for mint::Extent
#include "mint/MeshTypes.hpp"  // for mesh types

// Slic includes
#include "slic/slic.hpp"    // for SLIC macros

#ifdef MINT_USE_SIDRE
#include "sidre/Group.hpp"  // for sidre::Group
#include "sidre/View.hpp"   // for sidre::View
#else

namespace axom
{
namespace sidre
{
struct Group { };
}
}

#endif

namespace axom
{
namespace mint
{
namespace blueprint
{

#ifdef MINT_USE_SIDRE

//------------------------------------------------------------------------------
bool validRootGroup( const sidre::Group* group )
{
  SLIC_ERROR_IF( group==AXOM_NULLPTR, "supplied group is NULL!" );

  const bool hasCoordsets  = group->hasChildGroup( "coordsets" );
  const bool hasTopologies = group->hasChildGroup( "topologies" );
  const bool hasFields     = group->hasChildGroup( "fields" );

  SLIC_WARNING_IF( !hasCoordsets, "sidre::Group " <<  group->getPathName() <<
                   " is missing coordsets group!" );
  SLIC_WARNING_IF( !hasTopologies, "sidre::Group " << group->getPathName() <<
                   " is missing topologies group!" );
  SLIC_WARNING_IF( !hasFields, "sidre::Group " << group->getPathName() <<
                   " is missing fields group!" );

  return ( ( hasCoordsets && hasTopologies && hasFields ) );
}

//------------------------------------------------------------------------------
bool validTopologyGroup( const sidre::Group* topo )
{
  SLIC_ERROR_IF( topo==AXOM_NULLPTR, "supplied topology group is NULL!" );

  const std::string path = topo->getPathName();

  const bool hasTypeView = topo->hasChildView( "type" );
  SLIC_WARNING_IF( !hasTypeView, "[" << path << "] is missing 'type' view!" );

  const bool isTypeAString =
      (hasTypeView) ? topo->getView( "type" )->isString() : false;
  SLIC_WARNING_IF( !isTypeAString,
                   "'type' view in [" << path << "] is not a string" );

  const bool hasCoordset = topo->hasChildView( "coordset" );
  SLIC_WARNING_IF( !hasCoordset,
                   "[" << path << "] is missing 'coordset' view!" );

  const bool isCoordsetAString =
      (hasCoordset) ? topo->getView( "coordset" )->isString() : false;
  SLIC_WARNING_IF( !isCoordsetAString,
                   "'coordset' view in [" << path << "] is not a string" );

  const bool status = hasTypeView   && hasCoordset &&
                      isTypeAString && isCoordsetAString;

  return ( status );
}

//------------------------------------------------------------------------------
bool validCoordsetGroup( const sidre::Group* coordset )
{
  SLIC_ERROR_IF( coordset==AXOM_NULLPTR, "supplied coordset group is NULL!" );

  const std::string path   = coordset->getPathName();

  const bool hasTypeView   = coordset->hasChildView( "type" );
  SLIC_WARNING_IF( !hasTypeView, "[" << path << "] is missing 'type' view!" );

  const bool isTypeAString =
      (hasTypeView) ? coordset->getView("type")->isString() : false;
  SLIC_WARNING_IF( !isTypeAString,
                   "'type' view in [" << path << "] is not a string" );

  return ( hasTypeView && isTypeAString );
}

//------------------------------------------------------------------------------
const sidre::Group* getCoordsetGroup( const sidre::Group* group,
                                      const sidre::Group* topology )
{
  SLIC_ERROR_IF( !blueprint::validRootGroup( group ),
                 "supplied group does not conform to the blueprint!" );
  SLIC_ERROR_IF( topology==AXOM_NULLPTR, "supplied topology group is null!" );
  SLIC_ERROR_IF( !blueprint::validTopologyGroup( topology ),
                 "supplied topology group does not conform to the blueprint!" );

  const sidre::Group* coordsets  = group->getGroup( "coordsets" );

  const char* coordset_name = topology->getView( "coordset" )->getString();
  SLIC_WARNING_IF( ! coordsets->hasChildGroup( coordset_name ),
                   "cannot find coordset [" << coordset_name << "] in " <<
                    coordsets->getPathName() );

  const sidre::Group* coords = coordsets->getGroup( coordset_name );
  SLIC_WARNING_IF( coords==AXOM_NULLPTR,
    "null coordset [" << coordset_name << "] in " <<coordsets->getPathName() );

  return coords;
}

//------------------------------------------------------------------------------
const sidre::Group* getCoordsetGroup( const sidre::Group* group,
                                      const std::string& coords )
{
  SLIC_ERROR_IF( !blueprint::validRootGroup( group ),
                 "supplied group does not conform to the blueprint!" );

  const sidre::Group* coordsets  = group->getGroup( "coordsets" );
  const std::string path         = coordsets->getPathName();

  // get coordset group
  const sidre::Group* coordset = AXOM_NULLPTR;
  if ( coords.empty() )
  {
    SLIC_ERROR_IF( coordsets->getNumGroups()==0,
                  "[" << coordsets->getPathName() << "] is empty!" );
    SLIC_WARNING_IF( coordsets->getNumGroups() > 1,
                     "multiple coordsets found!  " );
    coordset = coordsets->getGroup( 0 );
  }
  else
  {
    SLIC_ERROR_IF( !coordsets->hasChildGroup( coords ),
      "[" << path << "] is missing requested coordset group [" << coords << "]" );

    coordset = coordsets->getGroup( coords );
  }

  return( coordset );
}

//------------------------------------------------------------------------------
const sidre::Group* getTopologyGroup( const sidre::Group* group,
                                      const std::string& topo )
{
  SLIC_ERROR_IF( !blueprint::validRootGroup( group ),
                 "supplied group does not conform to the blueprint!" );

  const sidre::Group* topologies = group->getGroup( "topologies" );
  const std::string path         = topologies->getPathName();

  // get topology group
  const sidre::Group* topology = AXOM_NULLPTR;
  if ( topo.empty() )
  {
    SLIC_ERROR_IF( topologies->getNumGroups()==0,
                  "[" << topologies->getPathName() << "] is empty!" );
    SLIC_WARNING_IF( topologies->getNumGroups() > 1,
                     "multiple topologies found!  " );
    topology = topologies->getGroup( 0 );
  }
  else
  {
    SLIC_ERROR_IF( !topologies->hasChildGroup( topo ),
      "[" << path << "] is missing requested topology group [" << topo << "]" );

    topology = topologies->getGroup( topo );
  }

  return( topology );
}

//------------------------------------------------------------------------------
void initializeTopologyGroup( sidre::Group* group,
                              const std::string& topo,
                              const std::string& coordset,
                              const std::string& type )
{
  SLIC_ASSERT( group != AXOM_NULLPTR );
  sidre::Group* topo_group = group->getGroup( "topologies" );
  SLIC_ASSERT( topo_group != AXOM_NULLPTR );
  sidre::Group* cur_topo = topo_group->getGroup( topo );
  SLIC_ASSERT( cur_topo != AXOM_NULLPTR );

  cur_topo->createView( "type" )->setString( type );
  cur_topo->createView( "coordset" )->setString( coordset );
}

//------------------------------------------------------------------------------
void getMeshTypeAndDimension( int& mesh_type, int& dimension,
                               const sidre::Group* group,
                               const std::string& topo )
{
  SLIC_ERROR_IF( !blueprint::validRootGroup( group ),
                 "supplied group does not conform to the blueprint!" );

  const sidre::Group* topology = blueprint::getTopologyGroup( group, topo );
  SLIC_ERROR_IF( !blueprint::validTopologyGroup( topology ),
                 "mesh topology does not conform to the blueprint!" );

  const sidre::Group* coords = blueprint::getCoordsetGroup( group, topology );
  SLIC_ERROR_IF( !blueprint::validCoordsetGroup( coords ),
                 "mesh coordset does not conform to the blueprint!" );

  // get topology type
  const char* topo_type = topology->getView( "type" )->getString();
  SLIC_ASSERT( topo_type != AXOM_NULLPTR );

  // detect mesh type based on the topology type
  if ( strcmp( topo_type, "uniform" )==0 )
  {

    SLIC_ERROR_IF( !coords->hasChildGroup("origin"),
            "missing [origin] group from [" << coords->getPathName() <<
            "], required for a uniform mesh" );

    mesh_type = STRUCTURED_UNIFORM_MESH;
    dimension = coords->getGroup("origin")->getNumViews();

  } // END if UNIFORM MESH
  else if ( strcmp( topo_type, "rectilinear" )== 0 )
  {
    SLIC_ERROR_IF( !coords->hasChildGroup("values"),
                "missing [values] group from [" << coords->getPathName() <<
                "], required for a rectilinear mesh" );

    mesh_type = STRUCTURED_RECTILINEAR_MESH;
    dimension = coords->getGroup( "values" )->getNumViews();

  } // END if RECTILINEAR_MESH
  else if ( strcmp( topo_type, "structured" )==0 )
  {

    SLIC_ERROR_IF( !coords->hasChildGroup("values"),
                   "missing [values] group from [" << coords->getPathName() <<
                   "], required for a structured mesh" );

    mesh_type = STRUCTURED_CURVILINEAR_MESH;
    dimension = coords->getGroup( "values" )->getNumViews();

  } // END if STRUCTURED_MESH
  else if ( strcmp( topo_type, "particle" )==0 )
  {
    // NOTE: currently the blueprint doesn't provide a topology type
    // that indicates particles
    SLIC_ERROR_IF( !coords->hasChildGroup("values"),
                   "missing [values] group from [" << coords->getPathName() <<
                   "], required for a particle mesh" );

    mesh_type = PARTICLE_MESH;
    dimension = coords->getGroup( "values" )->getNumViews();

  } // END if PARTICLE_MESH
  else if ( strcmp( topo_type, "unstructured" )==0 )
  {
    SLIC_ERROR_IF( !coords->hasChildGroup("values"),
                   "missing [values] group from [" << coords->getPathName() <<
                   "], required for a unstructured mesh" );

    // check if this is a particle mesh stored as an unstructured mesh
    const char* shape = topology->getView( "elements/shape" )->getString();
    mesh_type = (strcmp(shape, "point")==0) ? PARTICLE_MESH : UNSTRUCTURED_MESH;
    dimension = coords->getGroup( "values" )->getNumViews();

  } // END if UNSTRUCTURED_MESH
  else
  {
    mesh_type = UNDEFINED_MESH;
    dimension = -1;
    SLIC_ERROR( "invalid mesh topology_type=[" << topo_type << "] " );
  }

}

//------------------------------------------------------------------------------
void getMeshTypeAndDimension( int& mesh_type, int& dimension,
                              const sidre::Group* group )
{
  blueprint::getMeshTypeAndDimension( mesh_type, dimension, group, "" );
}

//------------------------------------------------------------------------------
void getUniformMesh( int dimension,
                     const sidre::Group* coordset,
                     const sidre::Group* topology,
                     double* origin,
                     double* spacing,
                     int64* extent )
{
  SLIC_ERROR_IF( (dimension < 1) && ( dimension > 3), "invalid dimension!" );
  SLIC_ERROR_IF( !blueprint::validCoordsetGroup( coordset ),
                 "invalid coordset group!" );
  SLIC_ERROR_IF( !blueprint::validTopologyGroup( topology ),
                 "invalid topology group!" );
  SLIC_ERROR_IF( origin==AXOM_NULLPTR, "supplied null pointer for origin!");
  SLIC_ERROR_IF( spacing==AXOM_NULLPTR, "supplied null pointer for spacing!" );
  SLIC_ERROR_IF( extent==AXOM_NULLPTR, "supplied null pointer for extent!" );

  sidre::Group* c = const_cast< sidre::Group* >( coordset );
  sidre::Group* t = const_cast< sidre::Group* >( topology );

  const char* dim_names[]     = { "dims/i", "dims/j", "dims/k" };
  const char* origin_names[]  = { "origin/x" , "origin/y", "origin/z" };
  const char* spacing_names[] = { "spacing/dx", "spacing/dy", "spacing/dz" };
  const char* topo_names[]    = { "elements/origin/i0",
                                  "elements/origin/j0",
                                  "elements/origin/k0"   };

  for ( int i=0; i < dimension; ++i )
  {
    origin [ i ] = c->getView( origin_names[ i ] )->getScalar();
    spacing[ i ] = c->getView( spacing_names[ i ] )->getScalar();

    const int idx    = i*2;
    const int N      = c->getView( dim_names[ i ] )->getScalar();
    const int offset = t->getView( topo_names[ i ] )->getScalar();
    extent[ idx ]    = offset;
    extent[ idx+1 ]  = offset + N - 1;
  }

}

//------------------------------------------------------------------------------
void setUniformMesh( int dim,
                     const double* origin,
                     const double* spacing,
                     const mint::Extent* extent,
                     sidre::Group* coordset,
                     sidre::Group* topology )
{
  SLIC_ERROR_IF( (dim < 1) && (dim > 3), "invalid dimension!" );
  SLIC_ERROR_IF( origin==AXOM_NULLPTR, "supplied null pointer for origin!" );
  SLIC_ERROR_IF( spacing==AXOM_NULLPTR, "supplied null pointer for spacing!" );
  SLIC_ERROR_IF( extent==AXOM_NULLPTR, "supplied extent is null!" );
  SLIC_ERROR_IF( coordset==AXOM_NULLPTR, "supplied coordset group is null!" );
  SLIC_ERROR_IF( topology==AXOM_NULLPTR, "supplied topology group is null!" );
  SLIC_ERROR_IF( dim != extent->getDimension(),
                 "supplied extent does not match specified dimension!" );

  const char* dim_names[]     = { "dims/i", "dims/j", "dims/k" };
  const char* origin_names[]  = { "origin/x" , "origin/y", "origin/z" };
  const char* spacing_names[] = { "spacing/dx", "spacing/dy", "spacing/dz" };
  const char* topo_names[]    = { "elements/origin/i0",
                                  "elements/origin/j0",
                                  "elements/origin/k0"   };

  coordset->createView( "type" )->setString( "uniform" );
  for ( int i=0; i < dim; ++i )
  {
    coordset->createView( dim_names[ i ] )->setScalar( extent->size( i ) );
    coordset->createView( origin_names[ i ] )->setScalar( origin[ i ] );
    coordset->createView( spacing_names[ i ] )->setScalar( spacing[ i ] );

    // FIXME: cannot set int64 values in sidre
    int min = static_cast< int >( extent->min( i ) );
    topology->createView( topo_names[ i ] )->setScalar( min );
  } // END for

}

#else

static const char* msg = "mint::blueprint requires building with Sidre!";

bool validRootGroup( const sidre::Group* )
{
  SLIC_ERROR( msg );
  return ( false );
}

//------------------------------------------------------------------------------
bool validTopologyGroup( const sidre::Group* )
{
  SLIC_ERROR( msg );
  return ( false );
}

//------------------------------------------------------------------------------
bool validCoordsetGroup( const sidre::Group* )
{
  SLIC_ERROR( msg );
  return ( false );
}

//------------------------------------------------------------------------------
const sidre::Group* getTopologyGroup( const sidre::Group*, const std::string& )
{
  SLIC_ERROR( msg );
  return ( AXOM_NULLPTR );
}

//------------------------------------------------------------------------------
const sidre::Group* getCoordsetGroup( const sidre::Group*,
                                      const sidre::Group*  )
{
  SLIC_ERROR( msg );
  return ( AXOM_NULLPTR );
}

//------------------------------------------------------------------------------
void getMeshTypeAndDimension( int& , int& ,
                              sidre::Group* ,
                              const std::string&  )
{
  SLIC_ERROR( msg );
}

//------------------------------------------------------------------------------
void getMeshTypeAndDimension( int& , int& ,
                              const sidre::Group* )
{
  SLIC_ERROR( msg );
}

#endif

} /* namespace blueprint */
} /* namespace mint */
} /* namespace axom */
