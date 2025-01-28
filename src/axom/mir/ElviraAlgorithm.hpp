// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for internals.
//
// SPDX-License-Identifier: (BSD-3-Clause)
#ifndef AXOM_MIR_ELVIRA_ALGORITHM_HPP_
#define AXOM_MIR_ELVIRA_ALGORITHM_HPP_

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/mir.hpp"
#include "axom/slic.hpp"

// Include these directly for now.
#include "axom/mir/MIRAlgorithm.hpp"
#include "axom/mir/utilities/ExtractZones.hpp"
#include "axom/mir/utilities/MergeMeshes.hpp"
#include "axom/mir/utilities/NodeToZoneRelationBuilder.hpp"
#include "axom/mir/utilities/RecenterField.hpp"
#include "axom/mir/utilities/ZoneListBuilder.hpp"
#include "axom/mir/views/dispatch_coordset.hpp"
#include "axom/mir/views/MaterialView.hpp"

#include <conduit/conduit.hpp>

// RAJA
#if defined(AXOM_USE_RAJA)
  #include "RAJA/RAJA.hpp"
#endif

#include <algorithm>
#include <string>

// Uncomment to save inputs and outputs.
// #define AXOM_ELVIRA_DEBUG

#if defined(AXOM_ELVIRA_DEBUG)
  #include <conduit/conduit_relay_io_blueprint.hpp>
#endif

namespace axom
{
namespace mir
{
namespace elvira
{

enum class Direction
{
  VERTICAL = 0,
  HORIZONTAL = 1
};

struct Result2D
{
  int plane{2};
  int difference_used{0};
  double columns[3]{0., 0., 0.};
  double normal[3][2]{{0., 0.}, {0., 0.}, {0., 0.}};
};

/*!
 * \brief 2nd order 2-d mir function: returns the normal to the slope in 
 *        the direction away from the material.
 *
 * \param[out] result The result object.
 * \param vf   A view containing volume fractions for the current material
 *             arranged like this:
 *                       0 - lower left neighbor
 *                       1 - lower neighbor
 *                       2 - lower right neighbor
 *                       3 - left neighbor
 *                       4 - zone in question
 *                       5 - right neighbor
 *                       6 - upper left neighbor
 *                       7 - upper neighbor
 *                       8 - upper right neighbor
 * \param direction hv direction to use
 *
 * \return The result object will contain these values:
 *         normal  : normal (x,y) pointing away from the material
 *         err     : the L2 error of the answer returned
 *
 * \note Adapted from J. Grandy's Overlink code
 */
AXOM_HOST_DEVICE
void elvira2d(Result2D &result, const ArrayView<double> &vf, Direction direction)
{  
  constexpr int jtop[3] = {6,7,8} ; // i-direction
  constexpr int jmid[3] = {3,4,5} ;
  constexpr int jbot[3] = {0,1,2} ;

  constexpr int irit[3] = {2,5,8} ; // j-direction
  constexpr int imid[3] = {1,4,7} ;
  constexpr int ilft[3] = {0,3,6} ;

  const double jb = vf[0] + vf[1] + vf[2]; // bottom row 
  const double jm = vf[3] + vf[4] + vf[5]; // middle row 
  const double jt = vf[6] + vf[7] + vf[8]; // top row 

  const double il = vf[0] + vf[3] + vf[6]; // left column
  const double im = vf[1] + vf[4] + vf[7]; // middle column
  const double ir = vf[2] + vf[5] + vf[8]; // right column

  double slope[3]{0., 0., 0.}; // backward, cent, and forward differencing 
  double slopev = 0.;          // Sets overall sign of normal.    
  if (direction == Direction::VERTICAL)
  {
    slope[0] = (im - il) ;     // backward differencing 
    slope[1] = (ir - il)/2.0 ; // central differencing 
    slope[2] = (ir - im) ;     // forward differencing 
    slopev   = (jb - jt) ;

    result.columns[0] = il ;
    result.columns[1] = im ;
    result.columns[2] = ir ;
  }
  else
  {
    slope[0] = (jm - jb) ;     // backward differencing 
    slope[1] = (jt - jb)/2.0 ; // central differencing 
    slope[2] = (jt - jm) ;     // forward differencing 
    slopev   = (il - ir) ;

    result.columns[0] = jb ;
    result.columns[1] = jm ;
    result.columns[2] = jt ;
  }

  const double sign = ( slopev > 0.0 ) ? ( 1 ) : ( -1 ) ;

  // return the normal in the direction away from the material 
  // i = loop over three difference schemes.
  const int ihv = static_cast<int>(direction);
  for (int i=0; i<3; i++)
  {
    result.normal[i][1-ihv] = sign ;      // steep direction 
    result.normal[i][  ihv] =-slope[i] ;  // shallow  
  }
}

// THIS CODE LOOKS LIKE IT NEEDS THE STENCIL DATA TO BE 3D, GIVEN THE INDICES SUPPLIED.

// MAYBE I SHOULD JUST PUT ALL THE RELEVANT FUNCTIONS IN ONE FILE AND THEN ADAPT THEM.
/*!
 * \brief 2nd order 2-d mir function.
 *
 * \return normal (x,y,0) pointing away from the material
 */
AXOM_HOST_DEVICE
void elvira2xy(double n[3], const axom::ArrayView<double> &vf)
//    ( double *vf ) // (same as elvira3d)
{
  int stencil2d[3][9] = {{1,4,7,10,13,16,19,22,25},    /*yz*/
			 {3,12,21,4,13,22,5,14,23},    /*zx*/
			 {9,10,11,12,13,14,15,16,17}} ;  /*xy*/
  int ivf[9] =  {9,10,11,12,13,14,15,16,17} ;
  
  int diff ;
  double *n ;
  double chmin=0.0 ;
  double variance[2], vfs[27] ;
  double chfac[3] = {1.0, 0.25, 1.0} ;

  Result2D elv2d;

  MALLOC1 ( n, 3, double ) ;

  variance[0] = ( ( vf[13] - vf[10] ) * ( vf[13] - vf[10] ) +
                  ( vf[13] - vf[16] ) * ( vf[13] - vf[16] ) ) ;
  variance[1] = ( ( vf[13] - vf[12] ) * ( vf[13] - vf[12] ) +
                  ( vf[13] - vf[14] ) * ( vf[13] - vf[14] ) ) ;

  int small = ( variance[0] < variance[1] ) ? 0 : 1 ;
  int   ihv = 1 - small ;

  elv2d[0] = elvira2d( vf, (int*)stencil2d[2], ihv ) ;
  elv2d[0]->plane = 2 ;

		/* Choose difference schemes.   */
  for ( diff=0 ; diff < 3 ; diff++ ) {
    double d, chisq ;
    double* n3 = norm2d ( elv2d[0]->normal[diff] ) ;

    d = d_3cube ( n3, vf[13] ) ;
    vf_3cube    ( n3, d, vfs ) ;

		/* Compute diffs between computed, actual vf's.  */
    chisq = elvira_chisq ( vf, vfs, ivf, 9 ) ;
    chisq *= chfac[diff] ;  // Preference for central difference 

    if ( ( diff == 0 ) || ( chisq < chmin ) ) {
      chmin = chisq ;
      memcpy ( n, n3, 3*sizeof(double) ) ;
    }
    FREE1 ( n3 ) ;
  }
  FREE1(elv2d[0]) ;
  return( n ) ;
}
/*!
 * \brief Multiply the normal times the Jacobian and normalize.
 *
 * \param normal The normal to transform, nx,ny,nz
 * \param jac    The jacobian to multiply by
 *
 * \note Adapted from Overlink, J. Grandy
 */
AXOM_HOST_DEVICE
void transform(double *normal, const double jac[3][3]) 
{
  assert(normal != nullptr);

  double norm   = 0.0 ;
  double dfvsum = 0.0 ; // stores nx^2 + ny^2 + nz^2
  double dfv[3], delfv[3];

  dfv[0] = normal[2] ; // BJW: Why switch components? Rotation?
  dfv[1] = normal[1] ;
  dfv[2] = normal[0] ;

  for (int f=0; f < 3; f++) {

    delfv[f] = dfv[ 0 ] * jac[ 0 ][ f ] + 
               dfv[ 1 ] * jac[ 1 ][ f ] +
               dfv[ 2 ] * jac[ 2 ][ f ];

    norm   += delfv[f] * delfv[f]; 
    dfvsum += dfv[f] * dfv[f] ;
  } 

  // dfvsum tests small variations in vol frac ( see mira1_youngs )
  if ( ( dfvsum > 1.0e-24 ) && ( norm > 0.0 ) ) {
    
    norm = 1.0 / sqrt(norm) ;
    for (int f=0 ; f<3 ; f++) {   
      delfv[f] = delfv[f] * norm ; 
    }

  } else {
    delfv[0] = 1.0 ;   /* An arbitrary slope. */
    delfv[1] = 0.0 ;
    delfv[2] = 0.0 ;
  }

  /* Store the slope.  */
  normal[0] = delfv[0] ;
  normal[1] = delfv[1] ;
  normal[2] = delfv[2] ;
}

// modeled after MirOvl::mira1 in mira1.c:262
AXOM_HOST_DEVICE
void mira1(axom::IndexType zoneIndex,
           axom::IndexType numInterfaces,
           axom::IndexType offset,
           const axom::ArrayView<double> &zoneMatStencilView,
           axom::ArrayView<double> &fragmentVectorsView,
           int numVectorComponents,
           int iskip)
{
  // This is like mira1.c:262 MirOvl::mira1
  double jac[3][3];
  computeJacobian(jac);

  elvira(zoneIndex, numInterfaces, iskip, zoneMatStencilView, fragmentVectorsView);

  for(axom::IndexType i = 0; i < numInterfaces; i++)
  {
    //axom::ArrayView<double> stencilVFs(zoneMatStencilView.data() + ((offset + m) * StencilSize), StencilSize);
    double *normal = fragmentVectorsView.data() + ((offset + i) * numVectorComponents);
    transform(normal, jac);
  }

}

} // end namespace detail

/*!
 * \brief Implements Elvira algorithm for 2D structured meshes.
 */
template <typename ExecSpace, typename IndexPolicy, typename CoordsetView, typename MatsetView>
class ElviraAlgorithm2D : public axom::mir::MIRAlgorithm
{
  using reduce_policy = typename axom::execution_space<ExecSpace>::reduce_policy;

public:
  using TopologyView = axom::mir::views::StructuredTopologyView<IndexPolicy>;
  using ConnectivityType = typename TopologyView::ConnectivityType;

  // Make sure we only instantiate this algorithm for 2D meshes.
  static_assert(IndexPolicy::dimension() == 2, "IndexPolicy must be 2D");

  /*!
   * \brief Constructor
   *
   * \param topoView The topology view to use for the input data.
   * \param coordsetView The coordset view to use for the input data.
   * \param matsetView The matset view to use for the input data.
   */
  ElviraAlgorithm2D(const TopologyView &topoView,
                    const CoordsetView &coordsetView,
                    const MatsetView &matsetView)
    : axom::mir::MIRAlgorithm()
    , m_topologyView(topoView)
    , m_coordsetView(coordsetView)
    , m_matsetView(matsetView)
  { }

  /// Destructor
  virtual ~ElviraAlgorithm2D() = default;

  /*!
   * \brief Perform material interface reconstruction on a single domain.
   *
   * \param[in] n_topo The Conduit node containing the topology that will be used for MIR.
   * \param[in] n_coordset The Conduit node containing the coordset.
   * \param[in] n_fields The Conduit node containing the fields.
   * \param[in] n_matset The Conduit node containing the matset.
   * \param[in] n_options The Conduit node containing the options that help govern MIR execution.
   *
   * \param[out] n_newTopo A node that will contain the new clipped topology.
   * \param[out] n_newCoordset A node that will contain the new coordset for the clipped topology.
   * \param[out] n_newFields A node that will contain the new fields for the clipped topology.
   * \param[out] n_newMatset A Conduit node that will contain the new matset.
   * 
   */
  virtual void executeDomain(const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset) override
  {
#if 0
    namespace bputils = axom::mir::utilities::blueprint;
    using reduce_policy = axom::execution_space<ExecSpace>::reduce_policy;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    AXOM_ANNOTATE_SCOPE("ElviraAlgorithm2D");

    // Copy the options to make sure they are in the right memory space.
    conduit::Node n_options_copy;
    bputils::copy<ExecSpace>(n_options_copy, n_options);
    n_options_copy["topology"] = n_topo.name();

    // Get selected zones from the options.
    bputils::SelectedZones<ExecSpace> selectedZones(m_topologyView.numberOfZones(), n_options_copy);
    const auto selectedZonesView = selectedZones.view();

    // Partition the selected zones into clean, mixed lists.
    axom::Array<axom::IndexType> cleanZones, mixedZones;
    bputils::ZoneListBuilder<ExecSpace, TopologyView, MatsetView> zlb(m_topologyView, m_coordsetView);
    zlb.execute(selectedZonesView, cleanZones, mixedZones);

    if(cleanZones.size() > 0 && mixedZones.size() > 0)
    {
      // Some clean, some mixed.
    }
    else if(mixedZones.size() > 0)
    {
      // All mixed.
    }
    else if(cleanZones.size() > 0)
    {
      // All clean. We could just pass the inputs through.
      n_newTopo.set_external(n_topo);
      n_newCoordset.set_external(n_coordset);
      n_newFields.set_external(n_fields);
      n_newMatset.set_external(n_matset);
    }
#endif
  }
#if 0
  void processMixedZones(const axom::ArrayView<axom::IndexType> mixedZonesView,
                         const conduit::Node &n_topo,
                             const conduit::Node &n_coordset,
                             const conduit::Node &n_fields,
                             const conduit::Node &n_matset,
                             const conduit::Node &n_options,
                             conduit::Node &n_newTopo,
                             conduit::Node &n_newCoordset,
                             conduit::Node &n_newFields,
                             conduit::Node &n_newMatset)
  {
    AXOM_ANNOTATE_SCOPE("processMixedZones");
    using loop_policy = typename axom::execution_space<ExecSpace>::loop_policy;
    using reduce_policy =
      typename axom::execution_space<ExecSpace>::reduce_policy;
    const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

    // Count the number of fragments we'll make for the mixed zones.
    // We also 
    AXOM_ANNOTATE_BEGIN("counting");
    RAJA::ReduceSum<reduce_policy, axom::IndexType> num_reduce(0);

    const auto nzones = mixedZonesView.size();
    axom::Array<axom::IndexType> matCount(nzones, nzones, allocatorID);
    axom::Array<axom::IndexType> matZone(nzones, nzones, allocatorID);
    auto matCountView = matCount.view();
    auto matZoneView = matZone.view();

    const auto matsetView = m_matsetView;
    axom::for_all<ExecSpace>(
      mixedZonesView.size(),
      AXOM_LAMBDA(axom::IndexType szIndex) {
        const auto zoneIndex = mixedZonesView[szIndex];
        const auto nmats = matsetView.numberOfMaterials(zoneIndex);
        matCountView[szIndex] = nmats;
        matZoneView[szIndex] = zoneIndex;
        sum_reduce += nmats;
      });
    const auto numFragments = sum_reduce.get();
    AXOM_ANNOTATE_END("counting");

    // Sort the zones by the mat count. This should make adjacent zones in the
    // list more likely to have the same number of materials.
    AXOM_ANNOTATE_BEGIN("sorting");
    RAJA::sort_pairs<loop_policy>(RAJA::make_span(matCountView.data(), nzones),
                                  RAJA::make_span(matZoneView.data(), nzones));
    AXOM_ANNOTATE_END("sorting");

    AXOM_ANNOTATE_BEGIN("offsets");
    axom::Array<axom::IndexType> matOffset(nzones, nzones, allocatorID);
    auto matOffsetView = matOffset.view();
    axom::exclusive_scan<ExecSpace>(matCountView, matOffsetView);
    AXOM_ANNOTATE_END("offsets");

    // Retrieve stencil data for each zone material.
    AXOM_ANNOTATE_BEGIN("stencil");
    constexpr int StencilSize = 9; // for 2D.
    const auto stencilDataSize = numFragments * StencilSize;
    axom::Array<double> zoneMatStencil(stencilDataSize, stencilDataSize, allocatorID);
    auto zoneMatStencilView = zoneMatStencil.view();
    axom::Array<MatsetView::IndexType> sortedMaterialIds(numFragments, numFragments, allocatorID);
    auto sortedMaterialIdsView = sortedMaterialIds.view();

    // Traverse the selected zones based on how many materials there are in a zone.
    axom::for_all<ExecSpace>(matZoneView.size(), AXOM_LAMBDA(axom::IndexType szIndex)
    {
      // The selected zone index in the whole mesh.
      const auto zoneIndex = matZoneView[szIndex];
      const auto matCount = matCountView[szIndex];
      // Where to begin writing this zone's fragment data.
      const auto offset = matOffsetView[szIndex];

      // Get materials for this zone from the matset.
      MatsetView::IDList ids;
      MatsetView::VFList vfs;
      matsetView.zoneMaterials(zoneIndex, ids, vfs);

      // Sort the materials by the volume fraction. Save sorted ids in sortedMaterialIdsView.
      axom::utilities::sort_multiple(vfs.data(), ids.data());
      for(axom::IndexType m = 0; m < matCount; m++)
      {
        sortedMaterialIdsView[offset + m] = ids[m];
      }

      // Retrieve the stencil data from neighbor zones.
      auto logical = topologyView.indexing().IndexToLogicalIndex(zoneIndex);
      constexpr int stencil[9][2] = {{-1, -1}, {0, -1}, {1, -1},
                                     {-1, 0}, {0, 0}, {1, 0},
                                     {-1, 1}, {0, 1}, {1, 1}};
      for(int si = 0; si < 9; si++)
      {
        // Stencil neighbor
        TopologyView::LogicalIndex neighbor(logical);
        neighbor[0] += stencil[si][0];
        neighbor[1] += stencil[si][1];

        if(topologyView.indexing().contains(neighbor))
        {
          const auto neighborIndex = topologyView.indexing().LogicalIndexToIndex(neighbor);

          for(axom::IndexType m = 0; m < matCount; m++)
          {
            // Ask the neighbor zone for how much of the current material it contains.
            double vf = 0.;
            matsetView.zoneContainsMaterial(neighborIndex, ids[m], vf);

            // Store the vf into the stencil for the current material.
            const auto destIndex = (offset + m) * StencilSize + si;
            zoneMatStencilView[destIndex] = vf;
          }
        }
        else
        {
          // All of the material contributions for the neighbor are 0.
          for(axom::IndexType m = 0; m < matCount; m++)
          {
            // Store the vf into the stencil for the current material.
            const auto destIndex = (offset + m) * StencilSize + si;
            zoneMatStencilView[destIndex] = 0.;
          }
        }
      }
    });
    AXOM_ANNOTATE_END("stencil");

    AXOM_ANNOTATE_BEGIN("vectors");
    constexpr int numVectorComponents = 2;
    const auto vecSize = numFragments * numVectorComponents;
    axom::Array<double> fragmentVectors(vecSize, vecSize, allocatorID);
    auto fragmentVectorsView = fragmentVectors.view();

    axom::for_all<ExecSpace>(matZoneView.size(), AXOM_LAMBDA(axom::IndexType szIndex)
    {
      // The selected zone index in the whole mesh.
      const auto zoneIndex = matZoneView[szIndex];
      const auto matCount = matCountView[szIndex];
      // Where to begin writing this zone's fragment data.
      const auto offset = matOffsetView[szIndex];

      // This is like mira1.c:262 MirOvl::mira1
      double jac[3][3];
      computeJacobian(jac);

      for(axom::IndexType m = 0; m < matCount; m++)
      {
        axom::ArrayView<double> stencilVFs(zoneMatStencilView.data() + ((offset + m) * StencilSize), StencilSize);
        axom::ArrayView<double> matVector(fragmentVectorsView.data() + ((offset + m) * numVectorComponents), numVectorComponents);
        // TODO: Compute the normal for this material, store in matVectorView
        // elvira().

        // multiply matVector by jac.
        
      }
    });
    AXOM_ANNOTATE_END("vectors");

    // intermediate...

    // Now that we have all of the normals for each zone, clip the zones using those normals, making a fragment
    // shape that has the right VF. then save the shape, take the remaining piece and apply the next normal to it.

    // Add in clean zones - do we end up with a vector of polyhedral zones here that might be easier to pass into
    // the mapping stage? Or, do we make a BG topology for the polyhedra and make sure the mapper can eat that?

    // Make matset for combined mesh

  }
#endif
private:
  TopologyView m_topologyView;
  CoordsetView m_coordsetView;
  MatsetView m_matsetView;
};

}  // end namespace mir
}  // end namespace axom

#endif
