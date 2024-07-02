// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_MIR_VIEWS_MATERIAL_VIEW_HPP_
#define AXOM_MIR_VIEWS_MATERIAL_VIEW_HPP_

#include "axom/core.hpp"
#include "axom/mir/views/Shapes.hpp"

#include <conduit/conduit.hpp>

#include <string>
#include <map>
#include <vector>

namespace axom
{
namespace mir
{
namespace views
{

/**
 * \brief This object contains information about the materials as provided by a Conduit node.
 *
 * \note This would only be used on the host.
 */
class MaterialInformation
{
public:
  using StringIntMap = std::map<std::string, int>;
  using IntStringMap = std::map<int, std::string>;
  using IntVector = std::vector<int>;
  using StringVector = std::vector<std::string>;

  void set(const conduit::Node &matset)
  {
    if(matset.has_child("material_map"))
    {
      m_namesToIds.clear();
      m_idsToNames.clear();
      m_ids.clear();
      m_names.clear();

      const conduit::Node &mm = matset["material_map"];
      for(conduit::index_t i = 0; i < mm.number_of_children(); i++)
      {
        const auto id = mm[i].to_int32();
        m_namesToIds[mm[i].name()] = id;
        m_idsToNames[id] = mm[i].name();
        m_ids.push_back(id);
        m_names.push_back(mm[i].name());
      }
    }
  }

  const StringIntMap &getNamesToIds() const { return m_namesToIds; }
  const IntStringMap &getIdsToNames() const { return m_idsToNames; }
  const IntVector &getIds() const { return m_ids; }
  const StringVector &getNames() const { return m_names; }

private:
  StringIntMap m_namesToIds{};
  IntStringMap m_idsToNames{};
  IntVector m_ids{};
  StringVector m_names{};
};

//---------------------------------------------------------------------------
// Material views - These objects are meant to wrap Blueprint Matsets behind
//                  an interface that lets us query materials for a single
//                  zone. It is intended that these views will be used in
//                  device kernels.
//---------------------------------------------------------------------------

/**

matsets:
  matset:
    topology: topology
    material_map:
      a: 1
      b: 2
      c: 0
    material_ids: [0, 1, 2, 2, 2, 0, 1, 0]
    volume_fractions: [0, a0, b2, b1, b0, 0, a1, 0]
    sizes: [2, 2, 1]
    offsets: [0, 2, 4]
    indices: [1, 4, 6, 3, 2]

 */
template <typename IndexT, typename FloatType, size_t MAXMATERIALS>
class UnibufferMaterialView
{
public:
  using MaterialIndex = IndexT;
  using ZoneIndex = IndexT;
  using IndexType = IndexT;
  using IDList = StaticArray<MaterialIndex, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static size_t MaxMaterials = MAXMATERIALS;

  void set(const axom::ArrayView<IndexType> &material_ids,
           const axom::ArrayView<FloatType> &volume_fractions,
           const axom::ArrayView<IndexType> &sizes,
           const axom::ArrayView<IndexType> &offsets,
           const axom::ArrayView<IndexType> &indices)
  {
    m_material_ids = material_ids;
    m_volume_fractions = volume_fractions;
    m_sizes = sizes;
    m_offsets = offsets;
    m_indices = indices;
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfZones() const
  {
    return m_sizes.size();
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfMaterials(ZoneIndex zi) const
  {
    assert(zi < m_sizes.size());
    return m_sizes[zi];
  }

  AXOM_HOST_DEVICE
  void getZoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    assert(zi < m_sizes.size());

    ids.clear();
    vfs.clear();

    const auto sz = m_sizes[zi];
    const auto offset = m_offsets[zi];
    for(size_t i = 0; i < sz; i++)
    {
      const auto idx = m_indices[offset + i];

      ids.push_back(m_material_ids[idx]);
      vfs.push_back(m_volume_fractions[idx]);
    }
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat) const
  {
    assert(zi < m_sizes.size());
    const auto sz = m_sizes[zi];
    const auto offset = m_offsets[zi];
    for(size_t i = 0; i < sz; i++)
    {
      const auto idx = m_indices[offset + i];

      if(m_material_ids[idx] == mat)
        return true;
    }
    return false;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat, FloatType &vf) const
  {
    assert(zi < m_sizes.size());
    const auto sz = m_sizes[zi];
    const auto offset = m_offsets[zi];
    for(size_t i = 0; i < sz; i++)
    {
      const auto idx = m_indices[offset + i];

      if(m_material_ids[idx] == mat)
      {
        vf = m_volume_fractions[idx];
        return true;
      }
    }
    vf = 0;
    return false;
  }

private:
  axom::ArrayView<IndexType> m_material_ids;
  axom::ArrayView<FloatType> m_volume_fractions;
  axom::ArrayView<IndexType> m_sizes;
  axom::ArrayView<IndexType> m_offsets;
  axom::ArrayView<IndexType> m_indices;
};


/**

matsets:
  matset:
    topology: topology
    volume_fractions:
      a:
        values: [0, 0, 0, a1, 0, a0]
        indices: [5, 3]
      b:
        values: [0, b0, b2, b1, 0]
        indices: [1, 3, 2]
    material_map: # (optional)
      a: 0
      b: 1
 */

// NOTE: I'm not sure I 100% get this one.

template <typename IndexT, typename FloatType, size_t MAXMATERIALS>
class MultiBufferMaterialView
{
public:
  using MaterialIndex = IndexT;
  using ZoneIndex = IndexT;
  using IDList = StaticArray<IndexType, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static size_t MaxMaterials = MAXMATERIALS;

  void add(const axom::ArrayView<ZoneIndex> &ids, const axom::ArrayView<FloatType> &vfs)
  {
    assert(m_size + 1 < MaxMaterials);

    m_indices[m_size] = ids;
    m_values[m_size] = vfs;
    m_size++;
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfZones() const
  {
    size_t nzones = 0;
    for(int i = 0; i < m_size; i++)
       nzones = axom::utilities::max(nzones, m_indices[i].size());
    return nzones;
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfMaterials(ZoneIndex zi) const
  {
    size_t nmats = 0;
    for(size_t i = 0; i < m_size; i++)
    {
      if(zi < m_indices[i].size())
      {
        const auto idx = m_indices[zi];
        if(m_values[i][idx] > 0.)
          nmats++;
      }
    }

    return nmats;
  }

  AXOM_HOST_DEVICE
  void getZoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    ids.clear();
    vfs.clear();

    for(size_t i = 0; i < m_size; i++)
    {
      if(zi < m_indices[i].size())
      {
        const auto idx = m_indices[zi];
        if(m_values[i][idx] > 0.)
        {
          ids.push_back(i);
          vfs.push_back(m_values[i][idx]);
        }
      }
    }
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat) const
  {
    assert(mat < m_size);
    assert(zi < m_indices[mat].size());
    
    const auto idx = m_indices[mat][zi];
    return m_values[mat][zi] > 0.;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat, FloatType &vf) const
  {
    assert(mat < m_size);
    assert(zi < m_indices[mat].size());
    
    const auto idx = m_indices[mat][zi];
    vf = m_values[mat][zi];
    return vf > 0.;
  }
private:
  axom::StackArray<axom::ArrayView<FloatType>, MAXMATERIALS> m_values{};
  axom::StackArray<axom::ArrayView<ZoneIndex>, MAXMATERIALS> m_indices{};
  size_t m_size{0};
};

/**
matsets:
  matset:
    topology: topology
    volume_fractions:
      a: [a0, a1, 0]
      b: [b0, b1, b2]
      c: [0, 0, c2]
    material_map: # (optional)
      a: 0
      b: 1
      c: 2
 */
template <typename IndexT, typename FloatType, size_t MAXMATERIALS>
class ElementDominantMaterialView
{
public:
  using MaterialIndex = IndexT;
  using ZoneIndex = IndexT;
  using IDList = StaticArray<IndexType, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static size_t MaxMaterials = MAXMATERIALS;

  void add(const axom::ArrayView<FloatType> &vfs)
  {
    m_volume_fractions.push_back(vfs);
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfZones() const
  {
    return (m_volume_fractions.size() > 0) ? m_volume_fractions[0].size() : 0;
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfMaterials(ZoneIndex zi) const
  {
    size_t nmats = 0;
    if(m_volume_fractions.size() > 0)
    {
      assert(zi < m_volume_fractions[0].size());
      for(size_t i = 0; i < m_volume_fractions.size(); i++)
        nmats += m_volume_fractions[i][zi] > 0. ? 1 : 0;
    }
    return nmats;
  }

  AXOM_HOST_DEVICE
  void getZoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    ids.clear();
    vfs.clear();

    if(m_volume_fractions.size() > 0)
    {
      assert(zi < m_volume_fractions[0].size());
      for(size_t i = 0; i < m_volume_fractions.size(); i++)
      {
        if(m_volume_fractions[i][zi] > 0)
        {
          ids.push_back(i);
          vfs.push_back(m_volume_fractions[i][zi]);
        }
      }
    }
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat) const
  {
    bool contains = false;
    if(m_volume_fractions.size() > 0)
    {
      assert(zi < m_volume_fractions[0].size());
      contains = m_volume_fractions[mat][zi] > 0;
    }
    return contains;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat, FloatType &vf) const
  {
    bool contains = false;
    vf = 0;
    if(m_volume_fractions.size() > 0)
    {
      assert(zi < m_volume_fractions[0].size());
      vf = m_volume_fractions[mat][zi] > 0;
      contains = vf > 0;
    }
    return contains;
  }
private:
  axom::StaticArray<axom::ArrayView<FloatType>, MAXMATERIALS> m_volume_fractions{};
};

/**
matsets:
  matset:
    topology: topology
    volume_fractions:
      a: [a0, a1]
      b: [b0, b1, b2]
      c: [c2]
    element_ids:
      a: [0, 1]
      b: [0, 1, 2]
      c: [2]
    material_map: # (optional)
      a: 0
      b: 1
      c: 2
 */
/// NOTES: This matset type does not seem so GPU friendly since there is some work to do for some of the queries.

template <typename IndexT, typename FloatType, size_t MAXMATERIALS>
class MaterialDominantMaterialView
{
public:
  using MaterialIndex = IndexT;
  using ZoneIndex = IndexT;
  using IDList = StaticArray<IndexType, MAXMATERIALS>;
  using VFList = StaticArray<FloatType, MAXMATERIALS>;

  constexpr static size_t MaxMaterials = MAXMATERIALS;

  void add(const axom::ArrayView<ZoneIndex> &ids, const axom::ArrayView<FloatType> &vfs)
  {
    assert(m_size + 1 < MaxMaterials);

    m_element_ids[m_size] = ids;
    m_volume_fractions[m_size] = vfs;
    m_size++;
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfZones()
  {
    if(m_nzones == 0)
    {
      for(size_t mi = 0; mi < m_size; mi++)
      {
        const auto sz = m_element_ids[mi].size();
        for(size_t i = 0; i < sz; i++)
          m_nzones = axom::utilities::max(m_nzones, m_element_ids[mi][i]);
#if 0
        // host-only
        // Eh, do this.
        RAJA::ReduceMax rm(0);
        axom::forall<ExecSpace>(0, sz, AXOM_LAMBDA(int i)
        {
          rm.max(m_element_ids[mi][i]);
        }
        m_nzones = axom::utilties::max(m_nzones, rm.get());
#endif
      }
    }
    return m_nzones;
  }

  AXOM_HOST_DEVICE
  size_t getNumberOfMaterials(ZoneIndex zi) const
  {
    size_t nmats = 0;
    for(size_t mi = 0; mi < m_size; mi++)
    {
      const auto sz = m_element_ids[mi].size();
#if 1
      for(size_t i = 0; i < sz; i++)
      {
        if(m_element_ids[mi][i] == zi)
        {
          nmats++;
          break;
        }
      }
#else
      // host-only
      RAJA::ReduceMax rm(0);
      axom::forall<ExecSpace>(0, sz, AXOM_LAMBDA(int i)
      {
        rm.max((m_element_ids[mi][i] == zi) ? 1 : 0);
      }
      m_nzones += rm.get();
#endif
    }
    return nmats;
  }

  AXOM_HOST_DEVICE
  void getZoneMaterials(ZoneIndex zi, IDList &ids, VFList &vfs) const
  {
    ids.clear();
    vfs.clear();

    for(size_t mi = 0; mi < m_size; mi++)
    {
      const auto sz = m_element_ids[mi].size();
#if 1
      for(size_t i = 0; i < sz; i++)
      {
        if(m_element_ids[mi][i] == zi)
        {
          ids.push_back(mi);
          vfs.push_back(m_volume_fractions[mi][i]);
          break;
        }
      }
#else
      RAJA::ReduceMax rm(-1);
      axom::forall<ExecSpace>(0, sz, AXOM_LAMBDA(int i)
      {
        rm.max((m_element_ids[mi][i] == zi) ? i : -1);
      }
      const auto index = rm.get();
      if(index != -1)
      {
        ids.push_back(mi);
        vfs.push_back(m_volume_fractions[mi][index]);
      }
#endif
    }
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat) const
  {
    assert(mat < m_element_ids.size());

    bool found = false;
#if 1
    const auto element_ids = m_element_ids[mat];
    for(size_t i = 0; i < element_ids.size(); i++)
    {
      if(element_ids[i] == zi)
      {
        found = true;
        break;
      }
    }
#endif
    return found;
  }

  AXOM_HOST_DEVICE
  bool zoneContainsMaterial(ZoneIndex zi, MaterialIndex mat, FloatType &vf) const
  {
    assert(mat < m_element_ids.size());

    bool found = false;
#if 1
    const auto element_ids = m_element_ids[mat];
    for(size_t i = 0; i < element_ids.size(); i++)
    {
      if(element_ids[i] == zi)
      {
        found = true;
        vf = m_volume_fractions[mat][i];
        break;
      }
    }
#endif
    return found;
  }

private:
  StackArray<axom::ArrayView<IndexType>, MAXMATERIALS> m_element_ids{};
  StackArray<axom::ArrayView<FloatType>, MAXMATERIALS> m_volume_fractions{};
  size_t m_size{0};
  size_t m_nzones{0};
};

#if 0
//---------------------------------------------------------------------------
// Some host-algorithms on material views.
//---------------------------------------------------------------------------

/**
 */
template <typename ExecSpace>
axom::Array<int> makeMatsPerZone(const MaterialDominantMaterialView &view, size_t nzones)
{
  // Figure out the number of materials per zone.
  axom::Array<int> matsPerZone(nzones, nzones, axom::getAllocatorID<ExecSpace>());
  auto matsPerZone_view = matsPerZone.view();
  axom::forall<ExecSpace>(0, nzones, AXOM_LAMBDA(int i)
  {
    matsPerZone_view[i] = 0;
  });
  for(size_t mi = 0; mi < m_size; mi++)
  {
    auto element_ids_view = m_element_ids[mi].view();
    axom::forall<ExecSpace>(0, nzones, AXOM_LAMBDA(int i)
    {
      matsPerZone_view[element_ids_view[i]]++;
    });
  }
  return matsPerZone;
}

// NOTE: This needs to be a method of MaterialDominantMaterialView to access the view data.
template <typename ExecSpace, typename Predicate>
axom::Array<int> selectZones(const MaterialDominantMaterialView &view, Predicate &&pred)
{
  const auto nzones = view.getNumberOfZones();

  // Figure out the number of materials per zone.
  axom::Array<int> matsPerZone = makeMatsPerZone<ExecSpace>(view, nzones);
  auto matsPerZone_view = matsPerZone.view();

  // Count the clean zones.
  RAJA::ReduceSum num_selected(0);
  axom::forall<ExecSpace>(0, nzones, AXOM_LAMBDA(int i)
  {
    num_selected += pred(matsPerZone_view[i]) ? 1 : 0;
  });
  size_t outsize = num_selected.get();
  
  // Make an offset array that records where each thread can write its data.
  axom::Array<int, 1, ExecSpace> offsets(nzones);
  RAJA::inclusive_scan<ExecSpace>(RAJA::make_span(zones, zones.size()),
                                  RAJA::make_span(offsets, offsets.size()));

  // Make a list of the selected output zones.
  axom::Array<int> zonelist(outsize, outsize, axom::getAllocatorID<ExecSpace>());
  auto zonelist_view = zonelist.view();
  axom::forall<ExecSpace>(0, nzones, AXOM_LAMBDA(int zi)
  {
    if(pred(matsPerZone_view[zi]))
      zonelist_view[offset[zi]] = zi;
  });

  return zonelist;
}

template <typename ExecSpace>
axom::Array<int> selectZonesContainingMaterial(const MaterialDominantMaterialView &view, MaterialIndex mat)
{
  const auto zones_view = view.selectZonesContainingMaterial(mat);
  axom::Array<int> zones(zones_view.size(), zones_view.size(), axom::getAllocatorID<ExecSpace>());
  axom::copy(zones.data(), zones_view.data(), zones_view.size() * sizeof(int));
  return zones;
}

template <typename ExecSpace>
axom::Array<int> selectCleanZones(const MaterialDominantMaterialView &view)
{
  auto predicate = [](int nmats) -> bool { return nmats == 1; };
  return selectZones<ExecSpace, decltype(predicate)>(view, predicate);
}

template <typename ExecSpace>
axom::Array<int> selectMixedZones(const MaterialDominantMaterialView &view)
{
  auto predicate = [](int nmats) -> bool { return nmats > 1; };
  return selectZones<ExecSpace, decltype(predicate)>(view, predicate);
}

//---------------------------------------------------------------------------
template <typename ExecSpace, typename Predicate>
axom::Array<int> selectZones(const UnibufferMaterialView &view, MaterialIndex mat, Predicate &&pred) const
{
/**
 NOTE: I really do not like the code below because it forces the main Axom algorithm to use RAJA directly.
       In the case of the reducer, I'd prefer to do this:

       auto reducer = axom::execution_space<ExecSpace>::make_ReduceSum<int>(0);

       Then we could write the algorithm so we do RAJA things but we could make a reducer object for the serial non-RAJA case that does nothing.
 */
  const size_t nzones = view.getNumberOfZones();

  using REDUCE_POL = typename axom::execution_space<ExecSpace>::reduce_policy;
  RAJA::ReduceSum<REDUCE_POL, size_t> num_selected(0);
  axom::Array<int, 1, ExecutionPolicy> zones(nzones);
  auto zones_view = zones.view();
  axom::forall<ExecSpace>(0, zones.size(), AXOM_LAMBDA(int zi)
  {
    const int haveMat = pred(view, mat, zi) ? 1 : 0;
    zones_view[zi] = haveMat;
    num_selected += haveMat;
  });
  size_t outsize = num_selected.get();
  
  // Make an offset array that records where each thread can write its data.
  axom::Array<int, 1, ExecSpace> offsets(nzones);
  RAJA::inclusive_scan<ExecSpace>(RAJA::make_span(zones, zones.size()),
                                  RAJA::make_span(offsets, offsets.size()));

  // Make a list of the selected output zones.
  axom::Array<int, 1, ExecSpace> zonelist(outsize);
  auto zonelist_view = zonelist.view();
  axom::forall<ExecutionPolicy>(0, zones.size(), AXOM_LAMBDA(int zi)
  {
    if(zones[zi] > 0)
      zonelist_view[offset[zi]] = zi;
  });

  return zonelist;
}

template <typename ExecSpace>
axom::Array<int> selectZonesContainingMaterial(const UnibufferMaterialView &view, MaterialIndex mat)
{
  auto findMaterial = [](const UnibufferMaterialView &deviceView, MaterialIndex deviceMat, int zi)
  {
    return deviceView.zoneContainsMaterial(zi, deviceMat);
  };
  return selectZones<ExecSpace>(view, mat, findMaterial);
}

template <typename ExecSpace>
axom::Array<int> selectCleanZones(const UnibufferMaterialView &view)
{
  auto zoneIsClean = [](const UnibufferMaterialView &deviceView, MaterialIndex /*mat*/, int zi)
  {
    return view.getNumberOfMaterials(zi) == 1;
  };
  return selectZones<ExecSpace>(view, 0, zoneIsClean);
}

template <typename ExecSpace>
axom::Array<int> selectMixedZones(const UnibufferMaterialView &view)
{
  auto zoneIsMixed = [](const UnibufferMaterialView &deviceView, MaterialIndex /*mat*/, int zi)
  {
    return view.getNumberOfMaterials(zi) > 1;
  };
  return selectZones<ExecSpace>(view, 0, zoneIsClean);
}
#endif

} // end namespace views
} // end namespace mir
} // end namespace axom

#endif
