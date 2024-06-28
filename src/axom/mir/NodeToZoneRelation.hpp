
/**
 * \brief Build an o2m relation that lets us look up the zones for a node.
 */
template <typename ExecSpace>
class NodeToZoneRelationBuilder
{
public:

  /**
   * \brief This 
   */
  void execute(const conduit::Node &topo, conduit::Node &relation);

private:
  template <typename ViewType>
  void buildRelation(const ViewType &nodes_view, ViewType &zones_view, ViewType &offsets_view) const;

};

/**
 * \brief This function makes a unique array of values from an input list of keys.
 *
 * \tparam ExecSpace The execution space.
 * \tparam KeyType   The data type for the keys.
 *
 * \param keys_orig_view The input view that contains the input keys to be made unique.
 * \param[out] skeys     A sorted unique array of keys produced from keys_orig_view.
 * \param[out] sindices  An array of indices that indicate where in the original view the keys came from.
 *
 * \note This code is adapted from Ascent/DevilRay.
 *
 */
template <typename ExecSpace, typename KeyType>
void unique(const axom::ArrayView<KeyType> &keys_orig_view, axom::Array<KeyType> &skeys, axom::Array<axom::IndexType> &sindices)
{
  using DataType = typename KeyViewType::value_type;
  using for_policy = axom::execution_space<ExecSpace>::for_policy;
  using reduce_policy = axom::execution_space<ExecSpace>::reduce_policy;
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Make a copy of the keys and make original indices.
  const auto n = keys_orig_view.size();
  axom::Array<DataType> keys(n, n, allocatorID);
  axom::Array<axom::IndexType> indices(n, n, allocatorID);
  auto keys_view = keys.view();
  auto indices_view = indices.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    keys_view[i] = keys_orig_view[i];
    indices_view[i] = i;
  });

  // Sort the keys, indices in place.
  RAJA::sort_pairs<for_policy>(RAJA::make_span(keys_view, n),
                               RAJA::make_span(indices_view, n));

  // Make a mask array for where differences occur.
  axom::Array<axom::IndexType> mask(n, n, allocatorID);
  auto mask_view = mask.view();
  RAJA::ReduceSum<reduce_policy, axom::IndexType> mask_sum(0);
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    const axom::IndexType different = (keys_view[i] != keys_view[i - 1]) ? 1 : 0;
    const axom::IndexType m = (i >= 1) ? different : 1;
    mask_view[i] = m;
    mask_sum += m;
  });

  // Do a scan on the mask array to build an offset array.
  axom::Array<axom::IndexType> offsets(n, n, allocatorID);
  auto offsets_view = offsets.view();
  RAJA::exclusive_scan<for_policy>(RAJA::make_span(mask_view, n),
                                   RAJA::make_span(offsets_view, n),
                                   RAJA::operators::plus<axom::IndexType>{});

  // Allocate the output arrays.
  const axom::IndexType newsize = mask_sum.get();
  skeys = axom::Array<KeyType>(newsize, newsize, allocatorID);
  sindices = axom::Array<axom::IndexType>(newsize, newsize, allocatorID);

  // Iterate over the mask/offsets to store values at the right
  // offset in the new array.
  auto skeys_view = skeys.view();
  auto sindices_view = sindices.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    if(mask_view[i])
    {
      skeys_view[offsets_view[i]] = keys_view[i];
      sindices_view[offsets_view[i]] = indices_view[i];
    }
  });
}

/// SHIT! What does this relation mean when some nodes are not used in a zone as in the strided structured case?


/**
 * \brief Make an unstructured representation of a structured topology.
 */
template <typename ExecSpace>
void
to_unstructured(const conduit::Node &topo, const conduit::Node &coordset, const std::string &topoName, conduit::Node &mesh)
{
  const std::string type = topo.fetch_existing("type").as_string();
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  mesh["coordsets"][coordset.name()].set_external(coordset);
  conduit::Node &newtopo = mesh["topologies"][topoName];

  if(type == "unstructured")
  {
    newtopo.set_external(topo);
  }
  else
  {
    newtopo["type"] = "unstructured";
    conduit::Node &n_newconn = newtopo["elements/connectivity"];
    n_newconn.set_allocator(allocatorID);

    // Fill in the connectivity.
    views::dispatch_structured_topologies(topo, coordset, [&](const std::string &shape, auto &topoView)
    {
      const auto nzones = topoView.numberOfZones();
      int ptsPerZone = 2;
      if(shape == "quad")
        ptsPerZone = 4;
      else if(shape == "hex")
        ptsPerZone = 8;

      newtopo["elements/shape"] = shape;

      const auto connSize = nzones * ptsPerZone;
      n_newconn.set(conduit::DataType::index_t(connSize));
      axom::ArrayView<conduit::index_t> conn(reinterpret_cast<conduit::index_t *>(n_newconn.data_ptr()), connSize);
      auto conn_view = conn.view();
      topoView. template for_all_zones<ExecSpace>(AXOM_LAMBDA(auto zoneIndex, const auto &zone)
      {
        const auto start = zoneIndex * ptsPerZone;
        for(int i = 0; i < 4; i++)
          conn_view[start + i] = static_cast<conduit::index_t>(zone.getIds()[i]);
      });
    });
  }
}


/**
 * \brief Given views that contain the nodes and zones, sort the zones using the
 *        node numbers to produce a list of zones for each node and an offsets array
 *        that points to the start of each list of zones.
 * 
 * \param[in]    nodes_view A view that contains the set of all of the nodes in the topology (the connectivity)
 * \param[inout[ zones_view A view (same size as \a nodes_view) that contains the zone number of each node.
 * \param[out]   offsets_view A view that we fill with offsets so offsets_view[i] points to the start of the i'th list in \a zones_view.
 */
template <typename ViewType>
void
NodeToZoneRelationBuilder<ExecSpace>::buildRelation(const ViewType &nodes_view, ViewType &zones_view, ViewType &offsets_view) const
{
  assert(nodes_view.size() == zones_view.size());

  using for_policy = axom::execution_space<ExecSpace>::for_policy;
  using reduce_policy = axom::execution_space<ExecSpace>::reduce_policy;
  const int allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  // Make a copy of the nodes that we'll use as keys.
  const auto n = nodes_view.size();
  axom::Array<axom::IndexType> keys(n, n, allocatorID);
  auto keys_view = keys.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    keys_view[i] = nodes_view[i];
  });

  // Sort the keys, zones in place. This sorts the zones_view which we want for output.
  RAJA::sort_pairs<for_policy>(RAJA::make_span(keys_view, n),
                               RAJA::make_span(zones_view, n));

  // Make a mask array for where differences occur.
  axom::Array<axom::IndexType> mask(n, n, allocatorID);
  auto mask_view = mask.view();
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    const axom::IndexType different = (keys_view[i] != keys_view[i - 1]) ? 1 : 0;
    const axom::IndexType m = (i >= 1) ? different : 1;
    mask_view[i] = m;
  });

  // Do a scan on the mask array to build an offset array.
  axom::Array<axom::IndexType> dest_offsets(n, n, allocatorID);
  auto dest_offsets_view = dest_offsets.view();
  RAJA::exclusive_scan<for_policy>(RAJA::make_span(mask_view, n),
                                   RAJA::make_span(dest_offsets_view, n),
                                   RAJA::operators::plus<axom::IndexType>{});

  // Build the offsets to each node's zone ids.
  axom::for_all<ExecSpace>(offsets_view.size(), AXOM_LAMBDA(axom::IndexType i)
  {
    offsets_view[i] = 0;
  });
  axom::for_all<ExecSpace>(n, AXOM_LAMBDA(axom::IndexType i)
  {
    if(mask_view[i])
    {
      offsets_view[dest_offsets_view[i]] = i;
    }
  });
}

template <typename ExecSpace>
void
NodeToZoneRelationBuilder<ExecSpace>::execute(const conduit::Node &topo, conduit::Node &relation)
{
  using loop_policy = axom::execution_space<ExecSpace>::loop_policy;
  const std::string type = topo.fetch_existing("type").as_string();
  const auto allocatorID = axom::execution_space<ExecSpace>::allocatorID();

  conduit::Node &n_zones = relation["zones"];
  conduit::Node &n_sizes = relation["sizes"];
  conduit::Node &n_offsets = relation["offsets"];
  n_zones.set_allocator(allocatorID);
  n_sizes.set_allocator(allocatorID);
  n_offsets.set_allocator(allocatorID);

  if(type == "unstructured")
  {
    conduit::blueprint::mesh::utils::Shape shape(topo);
    const conduit::Node &n_connectivity = topo["elements/connectivity"];
    const auto intTypeId = n_connectivity.dtype().id()
    const auto connSize = n_connectivity.dtype().number_of_elements()

    // Use the coordset to get the number of nodes. Conduit should be able to do this using only metadata.
    const conduit::Node *coordset = conduit::blueprint::mesh::utils::find_reference_node(topo, "coordset");
    const auto nnodes = conduit::blueprint::mesh::utils::coordset::length(*coordset);

    if(shape.is_polyhedral())
    {
      const conduit::Node &n_topo_sizes = topo["elements/sizes"];
      const conduit::Node &n_topo_offsets = topo["elements/offsets"];
      // TODO: I could iterate over zones and store node ids into a big buffer but it would overcount zones in the end.
      //       Do something better.
    }
    else if(shape.is_polygonal())
    {
      const conduit::Node &n_topo_sizes = topo["elements/sizes"];
      const conduit::Node &n_topo_offsets = topo["elements/offsets"];

      const auto nzones = n_topo_sizes.dtype().number_of_elements();

      // Allocate Conduit arrays on the device in a data type that matches the connectivity.
      n_zones.set(conduit::DataType(intTypeId, connSize));
      n_sizes.set(conduit::DataType(intTypeId, nnodes));
      n_offsets.set(conduit::DataType(intTypeId, nnodes));

      // Make zones for each node
      views::IndexNode_to_ArrayView_same(n_zones, n_topo_sizes, n_topo_offsets, [&](auto zonesView, auto sizesView, auto offsetsView)
      {
        using DataType = typename decltype(zonesView)::value_type;
        axom::for_all<ExecSpace>(0, nzones, AXOM_LAMBDA(axom::IndexType zoneIndex)
        {
          for(DataType i = 0; i < sizesView[zoneIndex]; i++)
            zonesView[offsetsView[zoneIndex] + i] = zoneIndex;
        });
      });

      views::IndexNode_to_ArrayView_same(n_connectivity, n_zones, n_sizes, n_offsets, [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView)
      {
        // Make the relation, outputting into the zonesView and offsetsView.
        buildRelation(connectivityView, zonesView, offsetsView);
        
        // Compute sizes from offsets.
        axom::for_all<ExecSpace>(offsetsView.size(), AXOM_LAMBDA(axom::IndexType index)
        {
          sizes_view[i] = (i < offsetsView.size() - 1) ? (offsetsView[i + 1]  - offsetsView[i]) : (connSize - offsets_view[i]);
        });
      }
    }
    else
    {
      // Shapes are all the same size.
      const auto nodesPerShape = shape.indices;
      const auto nzones = connSize / nodesPerShape;

      // Allocate Conduit arrays on the device in a data type that matches the connectivity.
      n_zones.set(conduit::DataType(intTypeId, connSize));
      n_sizes.set(conduit::DataType(intTypeId, nnodes));
      n_offsets.set(conduit::DataType(intTypeId, nnodes));

      views::IndexNode_to_ArrayView_same(n_connectivity, n_zones, n_sizes, n_offsets, [&](auto connectivityView, auto zonesView, auto sizesView, auto offsetsView)
      {
        // Make zones for each node
        axom::for_all<ExecSpace>(0, connSize, AXOM_LAMBDA(axom::IndexType index)
        {
          zonesView[index] = index / nodesPerShape;
        });
        // Make the relation, outputting into the zonesView and offsetsView.
        buildRelation(connectivityView, zonesView, offsetsView);
        
        // Compute sizes from offsets.
        axom::for_all<ExecSpace>(offsetsView.size(), AXOM_LAMBDA(axom::IndexType index)
        {
          sizes_view[i] = (i < offsetsView.size() - 1) ? (offsetsView[i + 1]  - offsetsView[i]) : (connSize - offsets_view[i]);
        });
      }
    }
  }
  else
  {
    // These are all structured topos of some sort. Make an unstructured representation and recurse.

    conduit::Node mesh;
    to_unstructured<ExecSpace>(topo, *coordset, "newtopo", mesh);
        
    // Recurse using the unstructured mesh.
    execute(mesh["newtopo"], relation);
  }
}
