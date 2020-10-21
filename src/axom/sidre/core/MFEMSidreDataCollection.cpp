// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"

#ifdef AXOM_USE_MFEM

  #include <string>
  #include <iomanip>  // for setw, setfill
  #include <cstdio>   // for snprintf()
  #include <unordered_set>

  #include "conduit_blueprint.hpp"

  #include "MFEMSidreDataCollection.hpp"

using mfem::Array;
using mfem::GridFunction;
using mfem::Mesh;
using mfem::Ordering;

namespace axom
{
namespace sidre
{
// Constructor that will automatically create the sidre data store and necessary
// data groups for domain and global data.
MFEMSidreDataCollection::MFEMSidreDataCollection(const std::string& collection_name,
                                                 Mesh* the_mesh,
                                                 bool own_mesh_data)
  : mfem::DataCollection(collection_name, the_mesh)
  , m_owns_datastore(true)
  , m_owns_mesh_data(own_mesh_data)
  , m_meshNodesGFName("mesh_nodes")
{
  m_datastore_ptr = new sidre::DataStore();

  sidre::Group* global_grp =
    m_datastore_ptr->getRoot()->createGroup(collection_name + "_global");
  sidre::Group* domain_grp =
    m_datastore_ptr->getRoot()->createGroup(collection_name);

  m_bp_grp = domain_grp->createGroup("blueprint");
  // Currently only rank 0 adds anything to bp_index.
  m_bp_index_grp = global_grp->createGroup("blueprint_index/" + name);

  m_named_bufs_grp = domain_grp->createGroup("named_buffers");

  if(the_mesh)
  {
    MFEMSidreDataCollection::SetMesh(the_mesh);
  }
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  else
  {
    m_comm = MPI_COMM_NULL;
  }
  #endif
}

// Second constructor that allows external code to specify data groups to place
// domain and global data in.

// TODO - Conduit will have the capability to generate a blueprint index group
// in the future.  When this is available, all the blueprint index code can be
// removed from the data collection class.
MFEMSidreDataCollection::MFEMSidreDataCollection(const std::string& collection_name,
                                                 Group* bp_index_grp,
                                                 Group* domain_grp,
                                                 bool own_mesh_data)
  : mfem::DataCollection(collection_name)
  , m_owns_datastore(false)
  , m_owns_mesh_data(own_mesh_data)
  , m_meshNodesGFName("mesh_nodes")
  , m_datastore_ptr(nullptr)
  , m_bp_index_grp(bp_index_grp)
{
  m_bp_grp = domain_grp->createGroup("blueprint");

  m_named_bufs_grp = domain_grp->createGroup("named_buffers");

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  m_comm = MPI_COMM_NULL;
  #endif
}

MFEMSidreDataCollection::~MFEMSidreDataCollection()
{
  attr_map.DeleteData(true);
  if(m_owns_datastore)
  {
    delete m_datastore_ptr;
  }
  if(m_owns_mesh_obj)
  {
    delete mesh;
  }
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
void MFEMSidreDataCollection::SetComm(MPI_Comm comm)
{
  m_comm = comm;
  serial = false;
  appendRankToFileName = true;
  MPI_Comm_rank(m_comm, &myid);
  MPI_Comm_size(m_comm, &num_procs);
}
  #endif

// protected method
sidre::Group* MFEMSidreDataCollection::named_buffers_grp() const
{
  SLIC_ASSERT_MSG(
    m_named_bufs_grp != nullptr,
    "No group 'named_buffers' in data collection.  Verify that"
    " SetMesh was called to set the mesh in the data collection.");
  return m_named_bufs_grp;
}

// protected method
View* MFEMSidreDataCollection::alloc_view(Group* grp, const std::string& view_name)
{
  SLIC_ASSERT_MSG(grp, "Group pointer is nullptr");
  sidre::View* v = nullptr;

  if(!grp->hasView(view_name))
  {
    v = grp->createView(view_name);
    SLIC_ASSERT_MSG(v,
                    "error allocating View " << view_name << " in group "
                                             << grp->getPathName());
  }
  else
  {
    v = grp->getView(view_name);
  }

  return v;
}

// protected method
View* MFEMSidreDataCollection::alloc_view(Group* grp,
                                          const std::string& view_name,
                                          const DataType& dtype)
{
  SLIC_ASSERT_MSG(grp, "Group pointer is nullptr");
  sidre::View* v = nullptr;

  if(!grp->hasView(view_name))
  {
    v = grp->createView(view_name, dtype);
    SLIC_ASSERT_MSG(v,
                    "error allocating View " << view_name << " in group "
                                             << grp->getPathName());
  }
  else
  {
    v = grp->getView(view_name);
    SLIC_ASSERT_MSG(v->getSchema().dtype().equals(dtype), "");
  }
  return v;
}

// protected method
Group* MFEMSidreDataCollection::alloc_group(Group* grp,
                                            const std::string& group_name)
{
  SLIC_ASSERT_MSG(grp, "Group pointer is nullptr");
  sidre::Group* g = nullptr;

  if(!grp->hasGroup(group_name))
  {
    g = grp->createGroup(group_name);
    SLIC_ASSERT_MSG(g,
                    "error allocating Group " << group_name << " in group "
                                              << grp->getPathName());
  }
  else
  {
    g = grp->getGroup(group_name);
  }
  return g;
}

// protected method
std::string MFEMSidreDataCollection::get_file_path(const std::string& filename) const
{
  std::stringstream fNameSstr;

  // Note: If non-empty, prefix_path has a separator ('/') at the end
  fNameSstr << prefix_path << filename;

  if(GetCycle() >= 0)
  {
    fNameSstr << "_" << std::setfill('0') << std::setw(pad_digits_cycle)
              << GetCycle();
  }

  return fNameSstr.str();
}

View* MFEMSidreDataCollection::AllocNamedBuffer(const std::string& buffer_name,
                                                IndexType sz,
                                                TypeID type)
{
  sz = std::max(sz, sidre::IndexType(0));
  sidre::Group* f = named_buffers_grp();
  sidre::View* v = nullptr;

  if(!f->hasView(buffer_name))
  {
    // create a buffer view
    v = f->createViewAndAllocate(buffer_name, type, sz);
  }
  else
  {
    v = f->getView(buffer_name);
    SLIC_ASSERT_MSG(v->getTypeID() == type, "type does not match existing type");

    // Here v is the view holding the buffer in the named_buffers group, so
    // its size is the full size of the buffer.

    // check if we need to resize.
    if(!v->isApplied() || v->getNumElements() < sz)
    {
      // resize, even if the buffer has more than 1 View.
      // v->reallocate(sz); // this will not work for more than 1 view.
      sidre::DataType dtype(v->getSchema().dtype());
      dtype.set_number_of_elements(sz);
      f->destroyViewAndData(buffer_name);
      v = f->createViewAndAllocate(buffer_name, dtype);
    }
  }
  SLIC_ASSERT_MSG(v && v->isApplied(), "allocation failed");
  return v;
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintStubs(bool hasBP)
{
  if(!hasBP)
  {
    m_bp_grp->createGroup("state");
    m_bp_grp->createGroup("coordsets");
    m_bp_grp->createGroup("topologies");
    m_bp_grp->createGroup("fields");
  }

  // If rank is 0, set up blueprint index state group.
  if(myid == 0)
  {
    m_bp_index_grp->createGroup("state");
    m_bp_index_grp->createGroup("coordsets");
    m_bp_index_grp->createGroup("topologies");
    m_bp_index_grp->createGroup("fields");
  }
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintState(bool hasBP)
{
  if(!hasBP)
  {
    // Set up blueprint state group.
    m_bp_grp->createViewScalar("state/cycle", 0);
    m_bp_grp->createViewScalar("state/time", 0.);
    m_bp_grp->createViewScalar("state/domain", myid);
    m_bp_grp->createViewScalar("state/time_step", 0.);
  }

  // If rank is 0, set up blueprint index state group.
  if(myid == 0)
  {
    m_bp_index_grp->createViewScalar("state/cycle", 0);
    m_bp_index_grp->createViewScalar("state/time", 0.);
    m_bp_index_grp->createViewScalar("state/number_of_domains", num_procs);
  }
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintCoordset(bool hasBP)
{
  int dim = mesh->SpaceDimension();
  SLIC_ASSERT_MSG(dim >= 1 && dim <= 3, "invalid mesh dimension");

  // Assuming mfem::Vertex has the layout of a double array.
  const int NUM_COORDS = sizeof(mfem::Vertex) / sizeof(double);

  const int num_vertices = mesh->GetNV();
  const int coordset_len = NUM_COORDS * num_vertices;

  // Add blueprint if not present
  if(!hasBP)
  {
    m_bp_grp->createViewString("coordsets/coords/type", "explicit");

    sidre::DataType dtype = sidre::DataType::c_double(num_vertices);
    const size_t stride = dtype.stride();
    dtype.set_stride(stride * NUM_COORDS);

    // Set up views for x, y, z values
    sidre::View *vx, *vy = nullptr, *vz = nullptr;
    vx = m_bp_grp->createView("coordsets/coords/values/x", dtype);

    if(dim >= 2)
    {
      dtype.set_offset(dtype.offset() + stride);
      vy = m_bp_grp->createView("coordsets/coords/values/y", dtype);
    }
    if(dim >= 3)
    {
      dtype.set_offset(dtype.offset() + stride);
      vz = m_bp_grp->createView("coordsets/coords/values/z", dtype);
    }

    if(m_owns_mesh_data)
    {
      // Allocate buffer for coord values.
      sidre::Buffer* coordbuf =
        AllocNamedBuffer("vertex_coords", coordset_len)->getBuffer();

      vx->attachBuffer(coordbuf);
      if(dim >= 2)
      {
        vy->attachBuffer(coordbuf);
      }
      if(dim >= 3)
      {
        vz->attachBuffer(coordbuf);
      }
    }
    else
    {
      double* coordbuf = mesh->GetVertex(0);

      vx->setExternalDataPtr(coordbuf);
      if(dim >= 2)
      {
        vy->setExternalDataPtr(coordbuf);
      }
      if(dim >= 3)
      {
        vz->setExternalDataPtr(coordbuf);
      }
    }
  }

  // If rank 0, set up blueprint index for coordinate set.
  if(myid == 0)
  {
    m_bp_index_grp->createViewString(
      "coordsets/coords/path",
      m_bp_grp->getPathName() + "/coordsets/coords");

    m_bp_index_grp->getGroup("coordsets/coords")
      ->copyView(m_bp_grp->getView("coordsets/coords/type"));

    m_bp_index_grp->createViewString("coordsets/coords/coord_system/type",
                                     "cartesian");

    // These are empty views, their existence in the group tree is used to
    // define the number of dims
    m_bp_index_grp->createView("coordsets/coords/coord_system/axes/x");

    if(dim >= 2)
    {
      m_bp_index_grp->createView("coordsets/coords/coord_system/axes/y");
    }

    if(dim == 3)
    {
      m_bp_index_grp->createView("coordsets/coords/coord_system/axes/z");
    }
  }

  if(m_owns_mesh_data)
  {
    double* coord_values = GetNamedBuffer("vertex_coords")->getData();
    // Change ownership of the mesh vertex data to sidre
    mesh->ChangeVertexDataOwnership(coord_values, coordset_len, hasBP);
  }
}

// private method
void MFEMSidreDataCollection::createMeshBlueprintTopologies(
  bool hasBP,
  const std::string& mesh_name)
{
  const bool isBdry = (mesh_name == "boundary");

  const int num_elements = !isBdry ? mesh->GetNE() : mesh->GetNBE();

  const std::string mesh_topo_str = "topologies/" + mesh_name;
  const std::string mesh_attr_str = mesh_name + "_material_attribute";

  int num_indices = 0;
  int geom = 0;
  std::string eltTypeStr = "point";

  if(num_elements > 0)
  {
    mfem::Element* elem = !isBdry ? mesh->GetElement(0) : mesh->GetBdrElement(0);
    int element_size = elem->GetNVertices();

    num_indices = num_elements * element_size;

    // Find the element shape
    // Note: Assumes homogeneous elements, so only check the first element
    geom = elem->GetGeometryType();
    eltTypeStr = getElementName(elem->GetType());
  }

  // Create the blueprint "topology" group, if not present
  if(!hasBP)
  {
    sidre::Group* topology_grp = m_bp_grp->createGroup(mesh_topo_str);

    topology_grp->createViewString("type", "unstructured");
    topology_grp->createViewString("elements/shape", eltTypeStr);
    topology_grp->createViewAndAllocate("elements/connectivity",
                                        sidre::INT_ID,
                                        num_indices);
    topology_grp->createViewString("coordset", "coords");

    // If the mesh has nodes, set the name of the GridFunction holding the
    // mesh nodes in the blueprint group.
    if(!isBdry && mesh->GetNodes() != nullptr)
    {
      topology_grp->createViewString("grid_function", m_meshNodesGFName);
    }
  }

  // Add the mesh's attributes as an attribute field
  RegisterAttributeField(mesh_attr_str, isBdry);

  // Change ownership or copy the element arrays into Sidre
  if(num_elements > 0)
  {
    sidre::View* conn_view =
      m_bp_grp->getGroup(mesh_topo_str)->getView("elements/connectivity");

    // The MFEMSidreDataCollection always owns these arrays:
    mfem::Array<int> conn_array(conn_view->getData<int*>(), num_indices);
    mfem::Array<int>* attr_array = attr_map.Get(mesh_attr_str);
    if(!isBdry)
    {
      mesh->GetElementData(geom, conn_array, *attr_array);
    }
    else
    {
      mesh->GetBdrElementData(geom, conn_array, *attr_array);
    }
    SLIC_ASSERT_MSG(!conn_array.OwnsData(), "");
    SLIC_ASSERT_MSG(!attr_array->OwnsData(), "");
  }

  // If rank 0, set up blueprint index for topologies group
  if(myid == 0)
  {
    const std::string m_bp_grp_path = m_bp_grp->getPathName();

    if(isBdry)
    {
      // "Shallow" copy the m_bp_grp view into the m_bp_index_grp sub-group.
      // Note that the "topologies/mesh" sub-group has to exist, i.e. this
      // method should be called first with mesh_name = "mesh".
      m_bp_index_grp->getGroup("topologies/mesh")
        ->copyView(m_bp_grp->getView("topologies/mesh/boundary_topology"));
    }

    sidre::Group* bp_index_topo_grp = m_bp_index_grp->createGroup(mesh_topo_str);
    sidre::Group* topology_grp = m_bp_grp->getGroup(mesh_topo_str);

    bp_index_topo_grp->createViewString("path",
                                        m_bp_grp_path + "/" + mesh_topo_str);
    bp_index_topo_grp->copyView(topology_grp->getView("type"));
    bp_index_topo_grp->copyView(topology_grp->getView("coordset"));

    // If the mesh has nodes, set the name of the GridFunction holding the
    // mesh nodes in the blueprint_index group.
    if(!isBdry && mesh->GetNodes() != nullptr)
    {
      bp_index_topo_grp->copyView(topology_grp->getView("grid_function"));
    }
  }
}

  // private method
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
void MFEMSidreDataCollection::createMeshBlueprintAdjacencies(bool hasBP)
{
  mfem::ParMesh* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);

  const int GRP_SZ = 25;
  char group_str[GRP_SZ];

  // TODO(JRC): Separate this out into group hierarchy setup and data allocation
  // stages like all of the other "createMeshBlueprint*" functions.
  SLIC_WARNING_IF(hasBP, "The case hasBP == true is not supported yet!");

  sidre::Group* adjset_grp = nullptr;
  if(pmesh->GetNGroups() > 1)
  {
    adjset_grp = m_bp_grp->createGroup("adjsets/mesh");
    adjset_grp->createViewString("association", "vertex");
    adjset_grp->createViewString("topology", "mesh");

    if(myid == 0)
    {
      m_bp_index_grp->createGroup("adjsets");
    }
  }

  for(int gi = 1; gi < pmesh->GetNGroups(); ++gi)
  {
    int num_gneighbors = pmesh->gtopo.GetGroupSize(gi);
    int num_gvertices = pmesh->GroupNVertices(gi);

    // Skip creation of empty groups
    if(num_gneighbors > 1 && num_gvertices > 0)
    {
      std::snprintf(group_str,
                    GRP_SZ,
                    "groups/g%d_%d",
                    pmesh->gtopo.GetGroupMasterRank(gi),
                    pmesh->gtopo.GetGroupMasterGroup(gi));
      sidre::Group* group_grp = adjset_grp->createGroup(group_str);

      sidre::View* gneighbors_view =
        group_grp->createViewAndAllocate("neighbors",
                                         sidre::INT_ID,
                                         num_gneighbors - 1);
      int* gneighbors_data = gneighbors_view->getData<int*>();

      // skip local domain when adding Blueprint neighbors
      const int* gneighbors = pmesh->gtopo.GetGroup(gi);
      for(int ni = 0, noff = 0; ni < num_gneighbors; ++ni)
      {
        if(gneighbors[ni] == 0)
        {
          noff++;
        }
        else
        {
          gneighbors_data[ni - noff] =
            pmesh->gtopo.GetNeighborRank(gneighbors[ni]);
        }
      }

      sidre::View* gvertices_view =
        group_grp->createViewAndAllocate("values", sidre::INT_ID, num_gvertices);
      int* gvertices_data = gvertices_view->getData<int*>();

      for(int vi = 0; vi < num_gvertices; ++vi)
      {
        gvertices_data[vi] = pmesh->GroupVertex(gi, vi);
      }
    }
  }
}
  #endif

// private method
bool MFEMSidreDataCollection::verifyMeshBlueprint()
{
  conduit::Node mesh_node;
  m_bp_grp->createNativeLayout(mesh_node);

  conduit::Node verify_info;
  bool result = conduit::blueprint::mesh::verify(mesh_node, verify_info);
  // conduit::Node::to_string only available in latest version
  SLIC_WARNING_IF(!result,
                  "MFEMSidreDataCollection blueprint verification failed: " /*<<
                                                                               verify_info.to_string()
                                                                               */                          );
  return result;
}

bool MFEMSidreDataCollection::HasBoundaryMesh() const
{
  // check if this rank has any boundary elements
  int hasBndElts = mesh->GetNBE() > 0 ? 1 : 0;

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  // check if any rank has boundary elements
  mfem::ParMesh* pmesh = dynamic_cast<mfem::ParMesh*>(mesh);
  if(pmesh)
  {
    int hasBndElts_g;
    MPI_Allreduce(&hasBndElts, &hasBndElts_g, 1, MPI_INT, MPI_MAX, pmesh->GetComm());

    hasBndElts = hasBndElts_g;
  }
  #endif

  return hasBndElts > 0 ? true : false;
}

void MFEMSidreDataCollection::SetMesh(Mesh* new_mesh)
{
  DataCollection::SetMesh(new_mesh);

  // hasBP is used to indicate if the data currently in the blueprint should be
  // used to replace the data in the mesh.
  bool hasBP = m_bp_grp->getNumViews() > 0 || m_bp_grp->getNumGroups() > 0;
  bool has_bnd_elts = HasBoundaryMesh();

  createMeshBlueprintStubs(hasBP);
  createMeshBlueprintState(hasBP);
  createMeshBlueprintCoordset(hasBP);

  GridFunction* nodes = new_mesh->GetNodes();

  // register the "mesh" topology in the blueprint.
  createMeshBlueprintTopologies(hasBP, "mesh");

  if(has_bnd_elts)
  {
    // Set the "boundary_topology" of "mesh" to "boundary".
    m_bp_grp->createViewString("topologies/mesh/boundary_topology", "boundary");

    // register the "boundary" topology in the blueprint.
    createMeshBlueprintTopologies(hasBP, "boundary");
  }

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  mfem::ParMesh* new_pmesh = dynamic_cast<mfem::ParMesh*>(new_mesh);
  m_comm = new_pmesh ? new_pmesh->GetComm() : MPI_COMM_NULL;
  if(new_pmesh)
  {
    createMeshBlueprintAdjacencies(hasBP);
  }
  #endif

  if(nodes)
  {
    // See the comment at the definition of 'hasBP' above.
    if(hasBP)
    {
      // Get the bp mesh nodes name.
      sidre::View* v_bp_nodes_name =
        m_bp_grp->getView("topologies/mesh/grid_function");
      std::string bp_nodes_name(v_bp_nodes_name->getString());

      // Check that the names match, e.g. when loading the collection.
      SLIC_WARNING_IF(m_meshNodesGFName == bp_nodes_name,
                      "mismatch of requested and blueprint mesh nodes names");
      // Support renaming bp_nodes_name --> m_meshNodesGFName ?
    }

    if(m_owns_mesh_data)
    {
      // Make sure Sidre owns the data of the new_mesh's Nodes.
      if(!GetNamedBuffer(m_meshNodesGFName))
      {
        int sz = new_mesh->GetNodalFESpace()->GetVSize();
        double* gfData = AllocNamedBuffer(m_meshNodesGFName, sz)->getData();

        // See the comment at the definition of 'hasBP' above.
        if(!hasBP)
        {
          SLIC_ASSERT_MSG(nodes->Size() == sz, "");
          std::memcpy(gfData, nodes->GetData(), sizeof(double) * sz);
        }
      }
      // Since the named buffer for m_meshNodesGFName exists, the call to
      // RegisterField() below will replace the data of the nodes with the
      // data from the named buffer.
    }
    else
    {
      // Make sure Sidre does not have a named buffer for m_meshNodesGFName.
      SLIC_WARNING_IF(GetNamedBuffer(m_meshNodesGFName) != nullptr, "");
    }

    RegisterField(m_meshNodesGFName, nodes);

    if(own_data)
    {
      // Avoid double delete calls (for the nodes gf) when (own_data == true)
      // and the new_mesh owns its Nodes --> take ownership from new_mesh.
      // When new_mesh does not own its Nodes and (own_data == true), we can
      // not take ownership --> verify that does not happen.
      SLIC_WARNING_IF(!new_mesh->OwnsNodes(),
                      "mesh does not own its nodes, "
                      "can not take ownership");
      new_mesh->SetNodesOwner(false);
    }
  }
}

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
void MFEMSidreDataCollection::SetMesh(MPI_Comm comm, Mesh* new_mesh)
{
  // use MFEMSidreDataCollection's custom SetMesh, then set MPI info
  SetMesh(new_mesh);

  m_comm = comm;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &num_procs);
}
  #endif

void MFEMSidreDataCollection::SetGroupPointers(Group* bp_index_grp,
                                               Group* domain_grp)
{
  SLIC_WARNING_IF(!domain_grp->hasGroup("blueprint"),
                  "Domain group does not contain a blueprint group.");

  m_bp_grp = domain_grp->getGroup("blueprint");
  m_bp_index_grp = bp_index_grp;
  m_named_bufs_grp = domain_grp->getGroup("named_buffers");
}

void MFEMSidreDataCollection::Load(const std::string& path,
                                   const std::string& protocol)
{
  DataCollection::DeleteAll();
  // Reset DataStore?

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(m_comm != MPI_COMM_NULL)
  {
    IOManager reader(m_comm);
    // The conduit abstraction appears to automatically handle the ".root"
    // suffix, but the IOManager does not, so it gets added here
    reader.read(m_bp_grp->getDataStore()->getRoot(), path + ".root");
  }
  else
  #endif
  {
    m_bp_grp->load(path, protocol);
  }

  // If the data collection created the datastore, it knows the layout of where
  // the domain and global groups are, and can restore them after the Load().
  //
  // If the data collection did not create the datastore, the host code must
  // reset these pointers after the load operation and also reset the state
  // variables.
  if(m_owns_datastore)
  {
    // The commented-out line matches the logic used for the datacoll-created datastore
    // SetGroupPointers(m_datastore_ptr->getRoot()->getGroup(name + "_global/blueprint_index/" + name),
    SetGroupPointers(m_datastore_ptr->getRoot()->getGroup(name + "_global"),
                     m_datastore_ptr->getRoot()->getGroup(name));

    UpdateStateFromDS();
  }

  // Create a mesh from the datastore that was just read in
  reconstructMesh();

  // Create any fields from the datastore that was just read in
  reconstructFields();
}

void MFEMSidreDataCollection::LoadExternalData(const std::string& path)
{
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(m_comm != MPI_COMM_NULL)
  {
    IOManager reader(m_comm);
    reader.loadExternalData(m_bp_grp->getDataStore()->getRoot(), path);
  }
  else
  #endif
  {
    m_bp_grp->loadExternalData(path);
  }
}

void MFEMSidreDataCollection::UpdateStateFromDS()
{
  SetTime(m_bp_grp->getView("state/time")->getData<double>());
  SetCycle(m_bp_grp->getView("state/cycle")->getData<int>());
  SetTimeStep(m_bp_grp->getView("state/time_step")->getData<double>());
}

void MFEMSidreDataCollection::UpdateStateToDS()
{
  SLIC_ASSERT_MSG(
    mesh != NULL,
    "Need to set mesh before updating state in MFEMSidreDataCollection.");

  m_bp_grp->getView("state/cycle")->setScalar(GetCycle());
  m_bp_grp->getView("state/time")->setScalar(GetTime());
  m_bp_grp->getView("state/time_step")->setScalar(GetTimeStep());

  if(myid == 0)
  {
    m_bp_index_grp->getView("state/cycle")->setScalar(GetCycle());
    m_bp_index_grp->getView("state/time")->setScalar(time);
  }
}

void MFEMSidreDataCollection::PrepareToSave()
{
  verifyMeshBlueprint();
  UpdateStateToDS();
}

void MFEMSidreDataCollection::Save()
{
  std::string filename = name;
  std::string protocol = "sidre_hdf5";

  Save(filename, protocol);
}

void MFEMSidreDataCollection::Save(const std::string& filename,
                                   const std::string& protocol)
{
  PrepareToSave();

  create_directory(prefix_path, mesh, myid);

  std::string file_path = get_file_path(filename);

  sidre::Group* blueprint_indicies_grp = m_bp_index_grp->getParent();
  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  if(m_comm != MPI_COMM_NULL)
  {
    IOManager writer(m_comm);
    sidre::DataStore* datastore = m_bp_grp->getDataStore();
    writer.write(datastore->getRoot(), num_procs, file_path, protocol);
    if(myid == 0)
    {
      if(protocol == "sidre_hdf5")
      {
        writer.writeGroupToRootFile(blueprint_indicies_grp, file_path + ".root");
      }
      // Root file support only available in hdf5.
      else
      {
        writer.write(blueprint_indicies_grp, 1, file_path + ".root", protocol);
      }
    }
  }
  else
  #endif
  {
    // If serial, use sidre group writer.
    m_bp_grp->save(file_path, protocol);

    blueprint_indicies_grp->save(file_path + ".root", protocol);
  }
}

// private method
void MFEMSidreDataCollection::addScalarBasedGridFunction(
  const std::string& field_name,
  GridFunction* gf,
  const std::string& buffer_name,
  IndexType offset)
{
  sidre::Group* grp = m_bp_grp->getGroup("fields/" + field_name);
  SLIC_ASSERT_MSG(grp != nullptr, "field " << field_name << " does not exist");

  const int numDofs = gf->FESpace()->GetVSize();

  if(gf->GetData() == nullptr)
  {
    AllocNamedBuffer(buffer_name, offset + numDofs);
    // gf->data is set below.
  }

  /*
   *  Mesh blueprint for a scalar-based grid function is of the form
   *    /fields/field_name/basis
   *              -- string value is GridFunction's FEC::Name
   *    /fields/field_name/values
   *              -- array of size numDofs
   */

  // Make sure we have the View "values".
  sidre::View* vv = alloc_view(grp, "values");

  // Describe and apply the "values" View.
  // If the data store has buffer for field_name (e.g. AllocNamedBuffer was
  // called, or it was loaded from file), use that buffer.
  if(named_buffers_grp()->hasView(buffer_name))
  {
    sidre::View* bv = named_buffers_grp()->getView(buffer_name);
    SLIC_ASSERT_MSG(bv->hasBuffer() && bv->isDescribed(), "");

    // named buffers always have offset 0
    SLIC_ASSERT_MSG(bv->getSchema().dtype().offset() == 0, "");
    SLIC_ASSERT_MSG(bv->getNumElements() >= offset + numDofs, "");

    if(vv->isEmpty())
    {
      vv->attachBuffer(bv->getBuffer())->apply(sidre::DOUBLE_ID, numDofs, offset);
    }

    gf->NewDataAndSize(vv->getData(), numDofs);
  }
  else
  {
    // If we are not managing the grid function's data,
    // create a view with the external data
    vv->setExternalDataPtr(sidre::DOUBLE_ID, numDofs, gf->GetData());
  }
  SLIC_ASSERT_MSG((numDofs > 0 && vv->isApplied()) ||
                    (numDofs == 0 && vv->isEmpty() && vv->isDescribed()),
                  "invalid View state");
  SLIC_ASSERT_MSG(numDofs == 0 || vv->getData() == gf->GetData(),
                  "View data is different from GridFunction data");
  SLIC_ASSERT_MSG(vv->getNumElements() == numDofs,
                  "View size is different from GridFunction size");
}

// private method
void MFEMSidreDataCollection::addVectorBasedGridFunction(
  const std::string& field_name,
  GridFunction* gf,
  const std::string& buffer_name,
  IndexType offset)
{
  sidre::Group* grp = m_bp_grp->getGroup("fields/" + field_name);
  SLIC_ASSERT_MSG(grp != nullptr, "field " << field_name << " does not exist");

  const int FLD_SZ = 20;
  char fidxName[FLD_SZ];

  int vdim = gf->FESpace()->GetVDim();
  int ndof = gf->FESpace()->GetNDofs();
  Ordering::Type ordering = gf->FESpace()->GetOrdering();

  if(gf->GetData() == nullptr)
  {
    AllocNamedBuffer(buffer_name, offset + vdim * ndof);
    // gf->data is set below.
  }

  /*
   *  Mesh blueprint for a vector-based grid function is of the form
   *    /fields/field_name/basis
   *              -- string value is GridFunction's FEC::Name
   *    /fields/field_name/values/x0
   *    /fields/field_name/values/x1
   *    ...
   *    /fields/field_name/values/xn
   *              -- each coordinate is an array of size ndof
   */

  // Get/create the Group "values".
  sidre::Group* vg = alloc_group(grp, "values");

  // Create the Views "x0", "x1", etc inside the "values" Group, vg.
  // If we have a named buffer for field_name, attach it to the Views;
  // otherwise set the Views to use gf->GetData() as external data.
  sidre::DataType dtype = sidre::DataType::c_double(ndof);
  const int entry_stride = (ordering == Ordering::byNODES ? 1 : vdim);
  const int vdim_stride = (ordering == Ordering::byNODES ? ndof : 1);
  dtype.set_stride(dtype.stride() * entry_stride);

  if(named_buffers_grp()->hasView(buffer_name))
  {
    sidre::View* bv = named_buffers_grp()->getView(buffer_name);
    SLIC_ASSERT_MSG(bv->hasBuffer() && bv->isDescribed(), "");

    // named buffers always have offset 0
    SLIC_ASSERT_MSG(bv->getSchema().dtype().offset() == 0, "");
    dtype.set_offset(dtype.element_bytes() * offset);

    for(int d = 0; d < vdim; d++)
    {
      std::snprintf(fidxName, FLD_SZ, "x%d", d);
      sidre::View* xv = alloc_view(vg, fidxName, dtype);
      xv->attachBuffer(bv->getBuffer());
      dtype.set_offset(dtype.offset() + dtype.element_bytes() * vdim_stride);
    }

    gf->NewDataAndSize(bv->getData<double*>() + offset, vdim * ndof);
  }
  else
  {
    for(int d = 0; d < vdim; d++)
    {
      std::snprintf(fidxName, FLD_SZ, "x%d", d);
      sidre::View* xv = alloc_view(vg, fidxName, dtype);
      xv->setExternalDataPtr(gf->GetData());
      dtype.set_offset(dtype.offset() + dtype.element_bytes() * vdim_stride);
    }
  }

  #ifdef AXOM_DEBUG
  for(int d = 0; d < vdim; d++)
  {
    std::snprintf(fidxName, FLD_SZ, "x%d", d);
    sidre::View* xv = vg->getView(fidxName);
    SLIC_ASSERT_MSG((ndof > 0 && xv->isApplied()) ||
                      (ndof == 0 && xv->isEmpty() && xv->isDescribed()),
                    "invalid View state");
    SLIC_ASSERT_MSG(ndof == 0 || xv->getData() == gf->GetData() + d * vdim_stride,
                    "View data is different from GridFunction data");
    SLIC_ASSERT_MSG(xv->getNumElements() == ndof,
                    "View size is different from GridFunction size");
  }
  #endif
}

// private method
// Should only be called on mpi rank 0 ( or if serial problem ).
void MFEMSidreDataCollection::RegisterFieldInBPIndex(const std::string& field_name,
                                                     GridFunction* gf)
{
  sidre::Group* bp_field_grp = m_bp_grp->getGroup("fields/" + field_name);
  sidre::Group* bp_index_field_grp =
    m_bp_index_grp->createGroup("fields/" + field_name);

  bp_index_field_grp->createViewString("path", bp_field_grp->getPathName());
  bp_index_field_grp->copyView(bp_field_grp->getView("topology"));
  bp_index_field_grp->copyView(bp_field_grp->getView("basis"));

  // Note: The bp index requires GridFunction::VectorDim()
  //       since the GF might be scalar valued and have a vector basis
  //       (e.g. hdiv and hcurl spaces)
  const int number_of_components = gf->VectorDim();
  bp_index_field_grp->createViewScalar("number_of_components",
                                       number_of_components);
}

// private method
// Should only be called on mpi rank 0 ( or if serial problem ).
void MFEMSidreDataCollection::DeregisterFieldInBPIndex(const std::string& field_name)
{
  sidre::Group* fields_grp = m_bp_index_grp->getGroup("fields");
  SLIC_WARNING_IF(!fields_grp->hasGroup(field_name),
                  "No field exists in blueprint index with name " << name);

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  fields_grp->destroyGroup(field_name);
}

void MFEMSidreDataCollection::RegisterField(const std::string& field_name,
                                            GridFunction* gf,
                                            const std::string& buffer_name,
                                            IndexType offset)
{
  if(field_name.empty() || buffer_name.empty() || gf == nullptr ||
     gf->FESpace() == nullptr)
  {
    return;
  }

  SLIC_ASSERT_MSG(mesh != nullptr,
                  "Need to set mesh before registering attributes in "
                  "MFEMSidreDataCollection.");

  // Register field_name in the blueprint group.
  sidre::Group* f = m_bp_grp->getGroup("fields");

  if(f->hasGroup(field_name))
  {
    // There are two possibilities:
    // 1. If HasField(field_name) is true - we are overwriting a field that
    //    was previously registered.
    // 2. Otherwise, the field was loaded from a file, or defined outside of
    //    the data collection.
    if(HasField(field_name))
    {
  #ifdef AXOM_DEBUG
      // Warn about overwriting field.
      // Skip warning when re-registering the nodal grid function
      if(field_name != m_meshNodesGFName)
      {
        SLIC_WARNING("field with the name '"
                     << field_name
                     << "' is already "
                        "registered, overwriting the old field");
      }
  #endif
      DeregisterField(field_name);
    }
  }

  // This will return the existing field (if external), otherwise, a new group
  sidre::Group* grp = alloc_group(f, field_name);

  // Set the "basis" string using the gf's finite element space, overwrite if
  // necessary.
  sidre::View* v = alloc_view(grp, "basis");
  v->setString(gf->FESpace()->FEColl()->Name());

  // Set the topology of the GridFunction.
  // This is always 'mesh' except for a special case with the boundary material
  // attributes field.
  v = alloc_view(grp, "topology")->setString("mesh");

  SLIC_ASSERT_MSG(gf->Size() == gf->FESpace()->GetVSize(),
                  "GridFunction size does not match FiniteElementSpace size");

  // Set the data views of the grid function
  // e.g. the number of coefficients per DoF -- either scalar-valued or
  // vector-valued
  bool const isScalarValued = (gf->FESpace()->GetVDim() == 1);
  if(isScalarValued)
  {
    // Set the View "<m_bp_grp>/fields/<field_name>/values"
    addScalarBasedGridFunction(field_name, gf, buffer_name, offset);
  }
  else  // vector valued
  {
    // Set the Group "<m_bp_grp>/fields/<field_name>/values"
    addVectorBasedGridFunction(field_name, gf, buffer_name, offset);
  }

  // Register field_name in the blueprint_index group.
  if(myid == 0)
  {
    RegisterFieldInBPIndex(field_name, gf);
  }

  // Register field_name + gf in field_map.
  DataCollection::RegisterField(field_name, gf);
}

void MFEMSidreDataCollection::RegisterAttributeField(const std::string& attr_name,
                                                     bool is_bdry)
{
  SLIC_ASSERT_MSG(mesh != nullptr,
                  "Need to set mesh before registering attributes in "
                  "MFEMSidreDataCollection.");

  // Register attr_name in the blueprint group.
  sidre::Group* f = m_bp_grp->getGroup("fields");
  if(f->hasGroup(attr_name))
  {
    bool isAttr = attr_map.Has(attr_name);
    bool isFld = field_map.Has(attr_name);

    if(isAttr)
    {
      SLIC_WARNING("field with the name '"
                   << attr_name
                   << "' is already "
                      " registered as an attribute, overwriting old values.");
      DeregisterAttributeField(attr_name);
    }
    else if(isFld)
    {
      SLIC_WARNING("field with the name '"
                   << attr_name
                   << "' is already "
                      " registered as a field, skipping register attribute.");
      return;
    }
  }

  // Generate sidre views and groups for this mesh attribute and allocate space
  addIntegerAttributeField(attr_name, is_bdry);

  if(myid == 0)
  {
    RegisterAttributeFieldInBPIndex(attr_name);
  }

  // Register new attribute array with attr_map
  sidre::View* a =
    m_bp_grp->getGroup("fields")->getGroup(attr_name)->getView("values");
  mfem::Array<int>* attr =
    new mfem::Array<int>(a->getData<int*>(), a->getNumElements());

  attr_map.Register(attr_name, attr, true);
}

void MFEMSidreDataCollection::RegisterAttributeFieldInBPIndex(
  const std::string& attr_name)
{
  SLIC_ASSERT_MSG(m_bp_grp->getGroup("fields") != nullptr,
                  "Mesh blueprint does not have 'fields' group");
  SLIC_ASSERT_MSG(m_bp_index_grp->getGroup("fields") != nullptr,
                  "Mesh blueprint index does not have 'fields' group");

  // get the BP attr group
  sidre::Group* attr_grp = m_bp_grp->getGroup("fields")->getGroup(attr_name);

  // create blueprint index for this attribute
  sidre::Group* bp_index_attr_grp =
    m_bp_index_grp->getGroup("fields")->createGroup(attr_name);

  bp_index_attr_grp->createViewString("path", attr_grp->getPathName());
  bp_index_attr_grp->copyView(attr_grp->getView("association"));
  bp_index_attr_grp->copyView(attr_grp->getView("topology"));
  bp_index_attr_grp->createViewScalar("number_of_components", 1);
}

void MFEMSidreDataCollection::DeregisterAttributeField(const std::string& attr_name)
{
  attr_map.Deregister(name, true);

  sidre::Group* attr_grp = m_bp_grp->getGroup("fields");
  SLIC_WARNING_IF(!attr_grp->hasGroup(attr_name),
                  "No field exists in blueprint with name " << attr_name);

  // Delete attr_name from the blueprint group.

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  attr_grp->destroyGroup(attr_name);

  // Delete field_name from the blueprint_index group.
  if(myid == 0)
  {
    DeregisterAttributeFieldInBPIndex(attr_name);
  }

  // Delete field_name from the named_buffers group, if allocated.
  FreeNamedBuffer(attr_name);
}

void MFEMSidreDataCollection::DeregisterAttributeFieldInBPIndex(
  const std::string& attr_name)
{
  sidre::Group* fields_grp = m_bp_index_grp->getGroup("fields");
  SLIC_WARNING_IF(
    !fields_grp->hasGroup(attr_name),
    "No attribute exists in blueprint index with name " << attr_name);

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  fields_grp->destroyGroup(attr_name);
}

void MFEMSidreDataCollection::addIntegerAttributeField(const std::string& attr_name,
                                                       bool is_bdry)
{
  sidre::Group* fld_grp = m_bp_grp->getGroup("fields");
  SLIC_ASSERT_MSG(fld_grp != nullptr, "'fields' group does not exist");

  const int num_elem = is_bdry ? mesh->GetNBE() : mesh->GetNE();
  std::string topo_name = is_bdry ? "boundary" : "mesh";

  sidre::Group* attr_grp = fld_grp->createGroup(attr_name);
  attr_grp->createViewString("association", "element");
  attr_grp->createViewAndAllocate("values", sidre::INT_ID, num_elem);
  attr_grp->createViewString("topology", topo_name);
}

void MFEMSidreDataCollection::DeregisterField(const std::string& field_name)
{
  // Deregister field_name from field_map.
  DataCollection::DeregisterField(field_name);

  sidre::Group* fields_grp = m_bp_grp->getGroup("fields");
  SLIC_WARNING_IF(!fields_grp->hasGroup(field_name),
                  "No field exists in blueprint with name " << field_name);

  // Delete field_name from the blueprint group.

  // Note: This will destroy all orphaned views or buffer classes under this
  // group also.  If sidre owns this field data, the memory will be deleted
  // unless it's referenced somewhere else in sidre.
  fields_grp->destroyGroup(field_name);

  // Delete field_name from the blueprint_index group.
  if(myid == 0)
  {
    DeregisterFieldInBPIndex(field_name);
  }

  // Delete field_name from the named_buffers group, if allocated.
  FreeNamedBuffer(field_name);
}

// private method
std::string MFEMSidreDataCollection::getElementName(mfem::Element::Type elementEnum)
{
  // Note -- the mapping from Element::Type to string is based on
  //   enum Element::Type { POINT, SEGMENT, TRIANGLE, QUADRILATERAL,
  //                        TETRAHEDRON, HEXAHEDRON };
  // Note: -- the string names are from conduit's blueprint

  switch(elementEnum)
  {
  case mfem::Element::POINT:
    return "point";
  case mfem::Element::SEGMENT:
    return "line";
  case mfem::Element::TRIANGLE:
    return "tri";
  case mfem::Element::QUADRILATERAL:
    return "quad";
  case mfem::Element::TETRAHEDRON:
    return "tet";
  case mfem::Element::HEXAHEDRON:
    return "hex";
  case mfem::Element::WEDGE:
  default:;
  }

  return "unknown";
}

mfem::Geometry::Type MFEMSidreDataCollection::getElementTypeFromName(
  const std::string& name)
{
  if(name == "point")
    return mfem::Geometry::POINT;
  else if(name == "line")
    return mfem::Geometry::SEGMENT;
  else if(name == "tri")
    return mfem::Geometry::TRIANGLE;
  else if(name == "quad")
    // FIXME: Good enough to assume quad always means square?
    return mfem::Geometry::SQUARE;
  else if(name == "tet")
    return mfem::Geometry::TETRAHEDRON;
  else if(name == "hex")
    // FIXME: Is it good enough to assume that hex always means cube?
    return mfem::Geometry::CUBE;
  else
    return mfem::Geometry::INVALID;
}
  #include <unistd.h>
// Implemented because no current ParMesh constructor supports
// Virtual inheritance is no good here
class SidreParMeshWrapper : public mfem::ParMesh
{
  void CreateMesh(double* _vertices,
                  int num_vertices,
                  int* element_indices,
                  mfem::Geometry::Type element_type,
                  int* element_attributes,
                  int num_elements,
                  int* boundary_indices,
                  mfem::Geometry::Type boundary_type,
                  int* boundary_attributes,
                  int num_boundary_elements,
                  int dimension,
                  int space_dimension)
  {
    if(space_dimension == -1)
    {
      space_dimension = dimension;
    }

    InitMesh(dimension,
             space_dimension,
             /*num_vertices*/ 0,
             num_elements,
             num_boundary_elements);

    int element_index_stride = mfem::Geometry::NumVerts[element_type];
    int boundary_index_stride =
      num_boundary_elements > 0 ? mfem::Geometry::NumVerts[boundary_type] : 0;

    // assuming Vertex is POD
    vertices.MakeRef(reinterpret_cast<mfem::Vertex*>(_vertices), num_vertices);
    NumOfVertices = num_vertices;

    for(int i = 0; i < num_elements; i++)
    {
      elements[i] = NewElement(element_type);
      elements[i]->SetVertices(element_indices + i * element_index_stride);
      elements[i]->SetAttribute(element_attributes[i]);
    }
    NumOfElements = num_elements;

    for(int i = 0; i < num_boundary_elements; i++)
    {
      boundary[i] = NewElement(boundary_type);
      boundary[i]->SetVertices(boundary_indices + i * boundary_index_stride);
      boundary[i]->SetAttribute(boundary_attributes[i]);
    }
    NumOfBdrElements = num_boundary_elements;

    FinalizeTopology();
  }

  struct GroupInfo
  {
    // Avoid making a copy, could also write a basic array_view
    // std::vector<int> shared_verts;
    int* shared_verts;
    int num_shared_verts;
    std::vector<int> shared_edges;
    std::vector<int> shared_faces;
    std::string name;
    // We can make a copy of this, will be small
    std::vector<int> neighbors;
  };

  std::vector<GroupInfo> getGroupInfo(Group* groups_grp, const int dimension)
  {
    auto num_groups = groups_grp->getNumGroups();
    std::vector<GroupInfo> result(num_groups);

    // We keep a set of already-processed vertices,
    // if all vertices in an edge/face are in the set,
    // then we say that the edge/face is shared.
    // Iteration is done in backwards order, so the last
    // group will have only shared edges/faces if they are
    // in the set of shared vertices for that one last group.
    // But the penultimate group will check its own shared vertices
    // the shared vertices of the last group, et cetera.

    // We use a set structure for faster lookup
    std::unordered_set<int> shared_vertices;

    // Don't have duplicate entries of a shared edge or
    // face across multiple groups
    std::unordered_set<int> already_processed_edges;
    std::unordered_set<int> already_processed_faces;

    // We have to go backwards for some reason?
    int group_idx = num_groups - 1;

    std::vector<sidre::IndexType> groups_grp_idxs;
    for(auto idx = groups_grp->getFirstValidGroupIndex();
        sidre::indexIsValid(idx);
        idx = groups_grp->getNextValidGroupIndex(idx))
    {
      groups_grp_idxs.push_back(idx);
    }

    // Yes, this isn't needed, but there aren't that many groups anyways
    std::reverse(groups_grp_idxs.begin(), groups_grp_idxs.end());

    for(const auto& idx : groups_grp_idxs)
    {
      auto group_grp = groups_grp->getGroup(idx);
      auto& grp_info = result[group_idx];
      grp_info.name = group_grp->getName();

      // Copy the neighbors array
      auto group_neighbors = group_grp->getView("neighbors");
      int* neighbors_array = group_neighbors->getArray();
      std::size_t num_neighbors = group_neighbors->getNumElements();
      grp_info.neighbors.assign(neighbors_array, neighbors_array + num_neighbors);

      // Fill in the vertices
      auto group_values = group_grp->getView("values");
      grp_info.shared_verts = group_values->getArray();
      grp_info.num_shared_verts = group_values->getNumElements();

      // Fill in the edges, if applicable
      if(dimension >= 2)
      {
        // Convert to a hash table for faster lookup/search
        shared_vertices.insert(grp_info.shared_verts,
                               grp_info.shared_verts + grp_info.num_shared_verts);

        // Scratchpad array for storing vertices of an edge/face
        mfem::Array<int> verts;

        // Add all the shared edges
        for(int i = 0; i < GetNEdges(); i++)
        {
          GetEdgeVertices(i, verts);
          // Check if it hasn't been assigned yet and it's a shared edge
          if(!already_processed_edges.count(i) &&
             shared_vertices.count(verts[0]) && shared_vertices.count(verts[1]))
          {
            grp_info.shared_edges.push_back(i);
            already_processed_edges.insert(i);
          }
        }
        // Fill in the faces, if applicable
        if(dimension >= 3)
        {
          // Add all the shared faces
          for(int i = 0; i < GetNumFaces(); i++)
          {
            GetFaceVertices(i, verts);
            // Check if it hasn't been assigned yet and it's a shared face
            if(!already_processed_faces.count(i) &&
               std::all_of(verts.begin(),
                           verts.end(),
                           [&shared_vertices](const int vert) {
                             return shared_vertices.count(vert);
                           }))
            {
              grp_info.shared_faces.push_back(i);
              already_processed_faces.insert(i);
            }
          }
        }
      }
      group_idx--;  // Backwards?
    }
    return result;
  }

public:
  SidreParMeshWrapper(Group* bp_grp,
                      double* vertices,
                      int num_vertices,
                      int* element_indices,
                      mfem::Geometry::Type element_type,
                      int* element_attributes,
                      int num_elements,
                      int* boundary_indices,
                      mfem::Geometry::Type boundary_type,
                      int* boundary_attributes,
                      int num_boundary_elements,
                      int dimension,
                      int space_dimension = -1)
  {
    MyComm = MPI_COMM_WORLD;  // FIXME parameter
    gtopo.SetComm(MPI_COMM_WORLD);
    CreateMesh(vertices,
               num_vertices,
               element_indices,
               element_type,
               element_attributes,
               num_elements,
               boundary_indices,
               boundary_type,
               boundary_attributes,
               num_boundary_elements,
               dimension,
               space_dimension);

    MPI_Comm_size(MyComm, &NRanks);
    MPI_Comm_rank(MyComm, &MyRank);

    ReduceMeshGen();  // determine the global 'meshgen'

    // FIXME: The word "group" is ridiculously overloaded in this context
    // Can we do any better with variable naming?
    // A group can refer to a communication group or a Sidre group

    if(!bp_grp->hasGroup("adjsets"))
    {
      SLIC_ERROR("Does not contain any adjacency sets");
    }
    // {
    // volatile int i = 0;
    // char hostname[256];
    // printf("PID %d ready for attach\n", getpid());
    // fflush(stdout);
    // while (0 == i)
    //     sleep(5);
    // }

    auto groups_grp = bp_grp->getGroup("adjsets/mesh/groups");
    auto num_groups = groups_grp->getNumGroups();

    auto groups_info = getGroupInfo(groups_grp, dimension);

    // TODO: Initialize the GroupTopology object
    int group_master_rank;   // The rank of the group master for a given group
    int group_master_group;  // The group number in the master for a given group
    mfem::ListOfIntegerSets group_integer_sets;

    // FIXME: This logic is probably completely wrong

    // The first group always contains only the current rank
    mfem::IntegerSet integer_set_0;
    mfem::Array<int>& array_0 = integer_set_0;  // Bind a ref so we can modify it
    array_0.Append(MyRank);
    group_integer_sets.Insert(integer_set_0);

    for(const auto& group_info : groups_info)
    {
      const int rc = std::sscanf(group_info.name.c_str(),
                                 "g%d_%d",
                                 &group_master_rank,
                                 &group_master_group);
      SLIC_ERROR_IF(rc != 2, "SidreDC: Failed to parse adjacency group name");

      mfem::IntegerSet integer_set;
      mfem::Array<int>& array = integer_set;  // Bind a ref so we can modify it

      // On a given rank, the group ID corresponding to the current rank
      // is a member of each group
      array.Append(MyRank);

      // TODO: Determine if the number of ranks used to build the blueprint differs
      // from the current number of ranks
      for(const auto neighbor : group_info.neighbors)
      {
        array.Append(neighbor);
      }
      group_integer_sets.Insert(integer_set);
    }

    gtopo.Create(group_integer_sets, 823);  // 822 or 823?

    std::size_t total_shared_vertices = 0;
    std::size_t total_shared_edges = 0;
    std::size_t total_shared_faces = 0;

    for(const auto& group_info : groups_info)
    {
      total_shared_vertices += group_info.num_shared_verts;
      total_shared_edges += group_info.shared_edges.size();
      total_shared_faces += group_info.shared_faces.size();
    }

    // Resizing for shared vertex data
    svert_lvert.SetSize(total_shared_vertices);
    group_svert.SetDims(num_groups, total_shared_vertices);

    if(dimension >= 2)
    {
      // Resizing for shared edge data
      sedge_ledge.SetSize(total_shared_edges);
      shared_edges.SetSize(total_shared_edges);
      group_sedge.SetDims(num_groups, total_shared_edges);
    }
    if(dimension >= 3)
    {
      // Resizing for shared face data
      sface_lface.SetSize(total_shared_faces);
      group_stria.MakeI(num_groups);
      group_squad.MakeI(num_groups);
    }
    else
    {
      group_stria.SetSize(num_groups, 0);  // create empty group_stria
      group_squad.SetSize(num_groups, 0);  // create empty group_squad
    }
    // Do we need to subtract 1 here?

    std::size_t group_idx = 1;  // The ID of the current group
    std::size_t global_shared_vert_idx = 0;
    std::size_t global_shared_edge_idx = 0;

    for(const auto& group_info : groups_info)
    {
      // Add shared vertices
      auto n_shared_verts = group_info.num_shared_verts;
      n_shared_verts +=
        global_shared_vert_idx;  // Get the index into the global table
      SLIC_ERROR_IF(n_shared_verts > group_svert.Size_of_connections(),
                    "incorrect number of total_shared_vertices");
      group_svert.GetI()[group_idx] = n_shared_verts;

      for(int i = 0; i < group_info.num_shared_verts; i++)
      {
        group_svert.GetJ()[global_shared_vert_idx] = global_shared_vert_idx;
        svert_lvert[global_shared_vert_idx] = group_info.shared_verts[i];
        global_shared_vert_idx++;
      }

      if(dimension >= 2)
      {
        // Scratchpad array for storing vertices of an edge/face
        mfem::Array<int> verts;

        if(MyRank == 0)
        {
          SLIC_INFO("Number of shared edges: " << group_info.shared_edges.size());
        }

        // Add all the shared edges
        auto n_shared_edges = group_info.shared_edges.size();
        n_shared_edges += global_shared_edge_idx;
        SLIC_ERROR_IF(n_shared_edges > static_cast<std::size_t>(
                                         group_sedge.Size_of_connections()),
                      "incorrect number of total_shared_edges");
        group_svert.GetI()[group_idx] = n_shared_edges;

        for(const auto shared_edge_idx : group_info.shared_edges)
        {
          GetEdgeVertices(shared_edge_idx, verts);
          group_sedge.GetJ()[global_shared_edge_idx] = global_shared_edge_idx;
          shared_edges[global_shared_edge_idx] =
            new mfem::Segment(verts[0], verts[1], 1);
          global_shared_edge_idx++;
          if(MyRank == 0)
          {
            SLIC_INFO("Shared edge with verts: " << verts[0] << " " << verts[1]);
          }
        }

        if(dimension >= 3)
        {
          // Add all the shared faces
          std::size_t triangles_added = 0;
          std::size_t quadrilaterals_added = 0;

          if(MyRank == 0)
          {
            SLIC_INFO(
              "Number of shared faces: " << group_info.shared_faces.size());
          }

          for(const auto shared_face_idx : group_info.shared_faces)
          {
            GetFaceVertices(shared_face_idx, verts);
            switch(GetFaceBaseGeometry(shared_face_idx))
            {
            case mfem::Geometry::TRIANGLE:
              shared_trias.Append({verts[0], verts[1], verts[2]});
              if(MyRank == 0)
              {
                SLIC_INFO("Shared triangle with verts: "
                          << verts[0] << " " << verts[1] << " " << verts[2]);
              }
              triangles_added++;
              break;
            case mfem::Geometry::SQUARE:
              shared_quads.Append({verts[0], verts[1], verts[2], verts[3]});
              quadrilaterals_added++;
              break;
            default:
              SLIC_ERROR("Face " << shared_face_idx << " had invalid geometry");
            }
          }
          group_stria.AddColumnsInRow(group_idx - 1, triangles_added);
          group_squad.AddColumnsInRow(group_idx - 1, quadrilaterals_added);
        }
      }

      group_idx++;
    }

    if(dimension >= 3)
    {
      SLIC_ERROR_IF(
        shared_trias.Size() + shared_quads.Size() != sface_lface.Size(),
        "incorrect number of total_shared_faces");
      // Define the J arrays of group_stria and group_squad -- they just contain
      // consecutive numbers starting from 0 up to shared_trias.Size()-1 and
      // shared_quads.Size()-1, respectively.
      group_stria.MakeJ();
      for(int i = 0; i < shared_trias.Size(); i++)
      {
        group_stria.GetJ()[i] = i;
      }
      group_squad.MakeJ();
      for(int i = 0; i < shared_quads.Size(); i++)
      {
        group_squad.GetJ()[i] = i;
      }
    }

    const bool fix_orientation = false;
    const bool refine = false;
    Finalize(refine, fix_orientation);
  }
};

// private method
void MFEMSidreDataCollection::reconstructMesh()
{
  // assumes m_own_mesh_data

  // Sufficient to always grab the underlying data to the x-coords
  // as the layout of the internal buffer is identical to what MFEM
  // expects, i.e., x0 y0 z0 x1 y1 z1 regardless of dimension
  double* vertices = m_bp_grp->getView("coordsets/coords/values/x")->getData();

  // Use the x to get the number of vertices
  int num_vertices =
    m_bp_grp->getView("coordsets/coords/values/x")->getNumElements();

  // Element indices
  int* element_indices =
    m_bp_grp->getView("topologies/mesh/elements/connectivity")->getData<int*>();

  // Name of the element type - convert to mfem::Geometry::Type
  std::string element_name =
    m_bp_grp->getView("topologies/mesh/elements/shape")->getString();

  // Element attributes
  int* element_attributes =
    m_bp_grp->getView("fields/mesh_material_attribute/values")->getData<int*>();

  // Number of elements
  int num_elements =
    m_bp_grp->getView("fields/mesh_material_attribute/values")->getNumElements();

  // Indices of boundary elements
  int* boundary_indices =
    m_bp_grp->getView("topologies/boundary/elements/connectivity")->getData<int*>();

  // Name of the element type - convert to mfem::Geometry::Type
  std::string bdr_element_name =
    m_bp_grp->getView("topologies/boundary/elements/shape")->getString();

  // Boundary attributes
  int* boundary_attributes =
    m_bp_grp->getView("fields/boundary_material_attribute/values")->getData<int*>();

  // Number of boundariy elements
  int num_boundary_elements =
    m_bp_grp->getView("fields/boundary_material_attribute/values")->getNumElements();

  int dimension = 1;

  if(m_bp_grp->hasView("coordsets/coords/values/z"))
  {
    dimension = 3;
  }

  else if(m_bp_grp->hasView("coordsets/coords/values/y"))
  {
    dimension = 2;
  }

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  // If it has an adjacencies group, the reloaded state was a ParMesh

  // FIXME: This does not work, the mesh read in on a given node is already
  // only the local part of the mesh, attempting to create a ParMesh with it
  // will result in the local partition being repartitioned
  if(m_bp_grp->hasGroup("adjsets/mesh"))
  {
    // I don't think it's possible to save an MPI communicator in the datastore
    // Though typically they are just integers
    SLIC_ERROR_IF(m_comm == MPI_COMM_NULL,
                  "Must set the communicator with SetComm before a ParMesh can "
                  "be reconstructed");

    mesh = new SidreParMeshWrapper(m_bp_grp,
                                   vertices,
                                   num_vertices,
                                   element_indices,
                                   getElementTypeFromName(element_name),
                                   element_attributes,
                                   num_elements,
                                   boundary_indices,
                                   getElementTypeFromName(bdr_element_name),
                                   boundary_attributes,
                                   num_boundary_elements,
                                   dimension);
  }
  else
  #endif
  {
    mesh = new mfem::Mesh(vertices,
                          num_vertices,
                          element_indices,
                          getElementTypeFromName(element_name),
                          element_attributes,
                          num_elements,
                          boundary_indices,
                          getElementTypeFromName(bdr_element_name),
                          boundary_attributes,
                          num_boundary_elements,
                          dimension);
  }

  // The DataCollection dtor is now responsible for deleting the mesh
  // We can't use owns_data in the base class here because that would
  // result in the deletion of the GridFunction objects which as of
  // now can only be externally managed (non-owning pointers are
  // obtained through RegisterField).
  m_owns_mesh_obj = true;
}

void MFEMSidreDataCollection::reconstructFields()
{
  sidre::Group* f = m_bp_grp->getGroup("fields");
  for(auto idx = f->getFirstValidGroupIndex(); sidre::indexIsValid(idx);
      idx = f->getNextValidGroupIndex(idx))
  {
    Group* field_grp = f->getGroup(idx);
    // Filter out the non-user-registered attribute fields
    if(!field_grp->hasView("association"))
    {
      // Build the new gridfunction

      // FiniteElementCollection
      auto basis_name = field_grp->getView("basis")->getString();
      m_fecolls.emplace_back(mfem::FiniteElementCollection::New(basis_name));

      // FiniteElementSpace - mesh ptr and FEColl ptr
      m_fespaces.emplace_back(
        new mfem::FiniteElementSpace(mesh, m_fecolls.back().get()));

      double* values = nullptr;
      // Scalar grid function
      if(field_grp->hasView("values"))
      {
        values = field_grp->getView("values")->getData();
      }

      // Vector grid function
      if(field_grp->hasGroup("values"))
      {
        // Sufficient to use address of first component as data is interleaved
        values = field_grp->getGroup("values")->getView("x0")->getData();
      }

      m_sidre_owned_gfs.emplace_back(
        new mfem::GridFunction(m_fespaces.back().get(), values));

      // Register a non-owning pointer with the base subobject
      DataCollection::RegisterField(field_grp->getName(),
                                    m_sidre_owned_gfs.back().get());
    }
  }
}

} /* namespace sidre */
} /* namespace axom */

#endif  // AXOM_USE_MFEM
