// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef MFEM_SIDRE_DATACOLLECTION_HPP_
#define MFEM_SIDRE_DATACOLLECTION_HPP_

#include "axom/config.hpp"

#ifdef AXOM_USE_MFEM

  #include "mfem.hpp"

  #include "axom/sidre/core/sidre.hpp"

namespace axom
{
namespace sidre
{
/** @brief Data collection with Sidre routines following the Conduit mesh
    blueprint specification. */
/** MFEMSidreDataCollection provides an HDF5-based file format for visualization
   or
    restart capability.  This functionality is aimed primarily at customers of
    LLNL's axom project that run problems at extreme scales.

    For more information, see:
    - LLNL conduit/blueprint library, https://github.com/LLNL/conduit
    - HDF5 library, https://support.hdfgroup.org/HDF5

    The layout created in the Sidre DataStore is: (`"──"` denote groups,
    `"─•"` denote views, `"─>"` denote links, i.e. shallow-copy view)

        <root>
         ├── <collection-name>_global          (global group)
         │    └── blueprint_index
         │         └── <collection-name>       (bp_index group)
         │              ├── state
         │              │    ├─• cycle
         │              │    ├─• time
         │              │    └─• number_of_domains = <mesh-mpi-comm-size>
         │              ├── coordsets
         │              │    └── coords
         │              │         ├─• path = "<bp-path>/coordsets/coords"
         │              │         ├─• type ─> <bp-grp>/coordsets/coords/type =
            "explicit"
         │              │         └─• coord_system = "x"|"xy"|"xyz"
         │              ├── topologies
         │              │    ├── mesh
         │              │    │    ├─• path = "<bp-path>/topologies/mesh"
         │              │    │    ├─• type              ─>
            <bp-grp>/topologies/mesh/type = "unstructured"
         │              │    │    ├─• coordset          ─>
            <bp-grp>/topologies/mesh/coordset = "coords"
         │              │    │    ├─• grid_function     ─>
            <bp-grp>/topologies/mesh/grid_function = "<nodes-field-name>"
         │              │    │    └─• boundary_topology ─>
            <bp-grp>/topologies/mesh/boundary_topology = "boundary"
         │              │    └── boundary
         │              │         ├─• path = "<bp-path>/topologies/mesh"
         │              │         ├─• type     ─>
            <bp-grp>/topologies/boundary/type = "unstructured"
         │              │         └─• coordset ─>
            <bp-grp>/topologies/boundary/coordset = "coords"
         │              └── fields
         │                   ├── mesh_material_attribute
         │                   │    ├─• path =
            "<bp-path>/fields/mesh_material_attribute"
         │                   │    ├─• association ─>
            <bp-grp>/fields/mesh_material_attribute/association = "element"
         │                   │    ├─• topology    ─>
            <bp-grp>/fields/mesh_material_attribute/topology = "mesh"
         │                   │    └─• number_of_components = 1
         │                   ├── boundary_material_attribute
         │                   │    ├─• path =
            "<bp-path>/fields/boundary_material_attribute"
         │                   │    ├─• association ─>
            <bp-grp>/fields/boundary_material_attribute/association = "element"
         │                   │    ├─• topology    ─>
            <bp-grp>/fields/boundary_material_attribute/topology = "boundary"
         │                   │    └─• number_of_components = 1
         │                   ├── grid-function-1
         │                   │    ├─• path = "<bp-path>/fields/grid-function-1"
         │                   │    ├─• basis    ─>
            <bp-grp>/fields/grid-function-1/basis = "<fe-coll-name>"
         │                   │    ├─• topology ─>
            <bp-grp>/fields/grid-function-1/topology = "mesh"
         │                   │    └─• number_of_components = gf1->VectorDim()
         │                   ├── grid-function-2
         │                   │    ├─• path = "<bp-path>/fields/grid-function-2"
         │                   │    ├─• basis    ─>
            <bp-grp>/fields/grid-function-2/basis = "<fe-coll-name>"
         │                   │    ├─• topology ─>
            <bp-grp>/fields/grid-function-2/topology = "mesh"
         │                   │    └─• number_of_components = gf2->VectorDim()
         │                   ├── ...
         │                  ...
         └── <collection-name>                 (domain group)
              ├── blueprint                    (blueprint group)
              │    ├── state
              │    │    ├─• cycle
              │    │    ├─• time
              │    │    ├─• domain = <mesh-mpi-rank>
              │    │    └─• time_step
              │    ├── coordsets
              │    │    └── coords
              │    │         ├─• type = "explicit"
              │    │         └── values
              │    │              ├─• x = view in
                 <vertex-coords-buffer>/<ext-double-data>
              │    │              ├─• y = view in
                 <vertex-coords-buffer>/<ext-double-data>
              │    │              └─• z = view in
                 <vertex-coords-buffer>/<ext-double-data>
              │    ├── topologies
              │    │    ├── mesh
              │    │    │    ├─• type = "unstructured"
              │    │    │    ├── elements
              │    │    │    │    ├─• shape = "points"|"lines"|...
              │    │    │    │    └─• connectivity = <vert-idx-array>
              │    │    │    ├─• coordset = "coords"
              │    │    │    ├─• grid_function = "<nodes-field-name>"
              │    │    │    └─• boundary_topology = "boundary"
              │    │    └── boundary
              │    │         ├─• type = "unstructured"
              │    │         ├── elements
              │    │         │    ├─• shape = "points"|"lines"|...
              │    │         │    └─• connectivity = <vert-idx-array>
              │    │         └─• coordset = "coords"
              │    └── fields
              │         ├── mesh_material_attribute
              │         │    ├─• association = "element"
              │         │    ├─• topology = "mesh"
              │         │    └─• values = <attr-array>
              │         ├── boundary_material_attribute
              │         │    ├─• association = "element"
              │         │    ├─• topology = "boundary"
              │         │    └─• values = <attr-array>
              │         ├── grid-function-1   (name can include path)
              │         │    ├─• basis = "<fe-coll-name>"
              │         │    ├─• topology = "mesh"
              │         │    └─• values = <ext-double-array>/<named-buffer>
                 (vdim == 1)
              │         ├── grid-function-2   (name can include path)
              │         │    ├─• basis = "<fe-coll-name>"
              │         │    ├─• topology = "mesh"
              │         │    └── values   (vdim > 1)
              │         │         ├─• x0 = view into
                 <ext-double-array>/<named-buffer>
              │         │         ├─• x1 = view into
                 <ext-double-array>/<named-buffer>
              │         │         └─• x2 = view into
                 <ext-double-array>/<named-buffer>
              │         ├── ...
              │        ...
              └── named_buffers                (named_buffers group)
                   ├─• vertex_coords = <double-array>
                   ├─• grid-function-1 = <double-array>
                   ├─• grid-function-2 = <double-array>
                  ...

    @note blueprint_index is used both in serial and in parallel. In parallel,
    only rank 0 will add entries to the blueprint index.

    @note QuadratureFunction%s (q-fields) are not supported.

    @note MFEMSidreDataCollection will attempt to reconstruct meshes and fields
    on a restart, and will automatically register them.  This functionality is
    experimental and only applies to meshes and fields written to a file by an
    instance of this class.

    @warning This class is still _experimental_, meaning that in future
    releases, it may not be backward compatible, and the output files generated
    by the current version may become unreadable.
 */
class MFEMSidreDataCollection : public mfem::DataCollection
{
public:
  using AttributeFieldMap = mfem::NamedFieldsMap<mfem::Array<int>>;
  AttributeFieldMap attr_map;

public:
  /// Constructor that allocates and initializes a Sidre DataStore.
  /**
      @param[in] collection_name  Name of the collection used as a file name
                                  when saving
      @param[in] the_mesh         Mesh shared by all grid functions in the
                                  collection (can be nullptr)
      @param[in] owns_mesh_data   Does the SidreDC own the mesh vertices?

      With this constructor, the MFEMSidreDataCollection owns the allocated
         Sidre
      DataStore.
   */
  explicit MFEMSidreDataCollection(const std::string& collection_name,
                                   mfem::Mesh* the_mesh = nullptr,
                                   bool owns_mesh_data = false);

  /// Constructor that links to an external Sidre DataStore.
  /** Specifically, the global and domain groups can be at arbitrary paths.

      @param[in] collection_name  Name of the collection used as a file name
                                  when saving
      @param[in] bp_index_grp     Pointer to the blueprint index group in the
                                  datastore, see the above schematic
      @param[in] domain_grp       Pointer to the domain group in the datastore,
                                  see the above schematic
      @param[in] owns_mesh_data   Does the SidreDC own the mesh vertices?

      With this constructor, the MFEMSidreDataCollection does not own the Sidre
      DataStore.
      @note No mesh or fields are read from the given Groups. The mesh has
      to be set with SetMesh() and fields registered with RegisterField().
   */
  MFEMSidreDataCollection(const std::string& collection_name,
                          Group* bp_index_grp,
                          Group* domain_grp,
                          bool owns_mesh_data = false);

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  /// Associate an MPI communicator with the collection.
  /** If no mesh was associated with the collection, this method should be
      called before using any of the Load() methods to read parallel data. */
  void SetComm(MPI_Comm comm);
  #endif

  /// Register a GridFunction in the Sidre DataStore.
  /** This method is a shortcut for the call
      `RegisterField(field_name, gf, field_name, 0)`.
   */
  virtual void RegisterField(const std::string& field_name, mfem::GridFunction* gf)
  {
    RegisterField(field_name, gf, field_name, 0);
  }

  /// Register a GridFunction in the Sidre DataStore.
  /** The registration procedure is as follows:
      - if (@a gf's data is nullptr), allocate named buffer with the name
        @a buffer_name with size _offset + gf->FESpace()->GetVSize()_ and use
        its data (plus the given @a offset) to set @a gf's data;
      - else, if (DataStore has a named buffer @a buffer_name), replace @a gf's
        data array with that named buffer plus the given @a offset;
      - else, use @a gf's data as external data associated with @a field_name
        in the DataStore;
      - register @a field_name in #field_map.

      Both the @a field_name and @a buffer_name can contain a path prefix.
      @note If @a field_name or @a buffer_name is empty, the method does
      nothing.
      @note If the GridFunction pointer @a gf or it's FiniteElementSpace
      pointer are nullptr, the method does nothing.
   */
  void RegisterField(const std::string& field_name,
                     mfem::GridFunction* gf,
                     const std::string& buffer_name,
                     IndexType offset);

  /// Registers an attribute field in the Sidre DataStore
  /** The registration process is similar to that of RegisterField()
      The attribute field is associated with the elements of the mesh
      when @a is_bdry is false, and with the boundary elements, when
      @a is_bdry is true.
      @sa RegisterField()  */
  void RegisterAttributeField(const std::string& name, bool is_bdry);
  void DeregisterAttributeField(const std::string& name);

  /** Returns a pointer to the attribute field associated with
      @a field_name, or nullptr when there is no associated field */
  mfem::Array<int>* GetAttributeField(const std::string& field_name) const
  {
    return attr_map.Get(field_name);
  }

  /** Checks if there is an attribute field associated with @a field_name */
  bool HasAttributeField(const std::string& field_name) const
  {
    return attr_map.Has(field_name);
  }

  /** Checks if any rank in the mesh has boundary elements */
  bool HasBoundaryMesh() const;

  /// Set the name of the mesh nodes field.
  /** This name will be used by SetMesh() to register the mesh nodes, if not
      already registered. Also, this method should be called if the mesh nodes
      GridFunction was or will be registered directly by the user. The default
      value for the name is "mesh_nodes". */
  void SetMeshNodesName(const std::string& nodes_name)
  {
    if(!nodes_name.empty())
    {
      m_meshNodesGFName = nodes_name;
    }
  }

  /// De-register @a field_name from the MFEMSidreDataCollection.
  /** The field is removed from the #field_map and the DataStore, including
      deleting it from the named_buffers group, if allocated. */
  virtual void DeregisterField(const std::string& field_name);

  /// Delete all owned data.
  virtual ~MFEMSidreDataCollection();

  /// Set/change the mesh associated with the collection
  /** Uses the field name "mesh_nodes" or the value set by SetMeshNodesName()
      to register the mesh nodes GridFunction, if the mesh uses nodes. */
  virtual void SetMesh(mfem::Mesh* new_mesh);

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  /// Set/change the mesh associated with the collection
  /** Uses the field name "mesh_nodes" or the value set by SetMeshNodesName()
      to register the mesh nodes GridFunction, if the mesh uses nodes. */
  virtual void SetMesh(MPI_Comm comm, mfem::Mesh* new_mesh);
  #endif

  /// Reset the domain and global datastore group pointers.
  /** These are set in the constructor, but if a host code changes the
      datastore contents ( such as wiping out the datastore and loading in new
      contents from a file, i.e. a restart ) these pointers will need to be
      reset to valid groups in the datastore.
      @sa Load(const std::string &path, const std::string &protocol).
   */
  void SetGroupPointers(Group* global_grp, Group* domain_grp);

  Group* GetBPGroup() { return m_bp_grp; }
  Group* GetBPIndexGroup() { return m_bp_index_grp; }

  /// Prepare the DataStore for writing
  virtual void PrepareToSave();

  /// Save the collection to file.
  /** This method calls `Save(collection_name, "sidre_hdf5")`. */
  virtual void Save();

  /// Save the collection to @a filename.
  /** The collection path prefix is prepended to the @a filename and the
      current cycle is appended, if cycle >= 0. */
  void Save(const std::string& filename, const std::string& protocol);

  /// Load the Sidre DataStore from file.
  /** The mesh and fields will be read from the loaded datastore

      If the data collection created the datastore, it knows the layout of
      where the domain and global groups are, and can restore them after the
      Load().

      If, however, the data collection does not own the datastore (e.g. it did
      not create the datastore), the host code must reset these pointers after
      the load operation, using SetGroupPointers(), and also reset the state
      variables, using UpdateStateFromDS().
   */
  void Load(const std::string& path, const std::string& protocol);

  /// Load MFEMSidreDataCollection from file.
  /** The used file path is based on the current prefix path, collection name,
      and the given @a cycle_. The protocol is "sidre_hdf5".
      @sa Load(const std::string &path, const std::string &protocol).
   */
  virtual void Load(int cycle_ = 0)
  {
    SetCycle(cycle_);
    Load(get_file_path(name), "sidre_hdf5");
  }

  /// Load external data after registering externally owned fields.
  void LoadExternalData(const std::string& path);

  /** @brief Updates the DataCollection's cycle, time, and time-step variables
      with the values from the data store. */
  void UpdateStateFromDS();

  /** @brief Updates the data store's cycle, time, and time-step variables with
      the values from the MFEMSidreDataCollection. */
  void UpdateStateToDS();

  /** @name Methods for named buffer access and manipulation. */
  ///@{

  /** @brief Get a pointer to the sidre::View holding the named buffer for
      @a buffer_name. */
  /** If such named buffer is not allocated, the method returns nullptr.
      @note To access the underlying pointer, use View::getData().
      @note To query the size of the buffer, use View::getNumElements().
   */
  View* GetNamedBuffer(const std::string& buffer_name) const
  {
    return named_buffers_grp()->hasView(buffer_name)
      ? named_buffers_grp()->getView(buffer_name)
      : nullptr;
  }

  /// Return newly allocated or existing named buffer for @a buffer_name.
  /** The buffer is stored in the named_buffers group. If the currently
      allocated buffer size is smaller than @a sz, then the buffer is
      reallocated with size @a sz, destroying its contents.
      @note To access the underlying pointer, use View::getData().
   */
  View* AllocNamedBuffer(const std::string& buffer_name,
                         IndexType sz,
                         TypeID type = DOUBLE_ID);

  /// Deallocate the named buffer @a buffer_name, if allocated.
  void FreeNamedBuffer(const std::string& buffer_name)
  {
    if(named_buffers_grp()->hasView(buffer_name))
    {
      named_buffers_grp()->destroyViewAndData(buffer_name);
    }
  }

  /** @brief Updates the DataCollection's mesh and registered fields
      with the values from the data store. */
  void UpdateMeshAndFieldsFromDS()
  {
    reconstructMesh();
    reconstructFields();
  }

  /// Verifies that the contents of the mesh blueprint data is valid.
  bool verifyMeshBlueprint();

  ///@}

private:
  // Used if the Sidre data collection is providing the datastore itself.
  const bool m_owns_datastore;

  // TODO - Need to evaluate if this bool member can be combined with own_data
  // in parent data collection class. m_owns_mesh_data indicates whether the
  // Sidre dc owns the mesh element data and node positions gf. The DC base
  // class own_data indicates if the dc owns the mesh object pointer itself and
  // GF objects. Can we use one flag and just have DC own all objects vs none?
  const bool m_owns_mesh_data;

  // Name to be used for registering the mesh nodes in the
  // MFEMSidreDataCollection.
  // This name is used by SetMesh() and can be overwritten by the method
  // SetMeshNodesName().
  // Default value: "mesh_nodes".
  std::string m_meshNodesGFName;

  // If the data collection owns the datastore, it will store a pointer to it.
  // Otherwise, this pointer is nullptr.
  DataStore* m_datastore_ptr;

protected:
  Group* named_buffers_grp() const;

  View* alloc_view(Group* grp, const std::string& view_name);

  View* alloc_view(Group* grp, const std::string& view_name, const DataType& dtype);

  Group* alloc_group(Group* grp, const std::string& group_name);

  // return the filename based on prefix_path, collection name and cycle.
  std::string get_file_path(const std::string& filename) const;

private:
  // If the data collection does not own the datastore, it will need pointers
  // to the blueprint and blueprint index group to use.
  Group* m_bp_grp;
  Group* m_bp_index_grp;

  // This is stored for convenience.
  Group* m_named_bufs_grp;

  // Used to retain ownership of components of reconstructed Meshes and GridFuncs
  // Instead of using flags to keep track of ownership between this class
  // and the mfem::DataCollection subobject, always have the subobject
  // retain a non-owning pointer and manage memory via unique_ptr in
  // this class for consistency
  std::unique_ptr<mfem::Mesh> m_owned_mesh;
  std::vector<std::unique_ptr<mfem::FiniteElementCollection>> m_fecolls;
  std::vector<std::unique_ptr<mfem::FiniteElementSpace>> m_fespaces;
  std::vector<std::unique_ptr<mfem::GridFunction>> m_owned_gridfuncs;

  // Private helper functions

  void RegisterFieldInBPIndex(const std::string& field_name,
                              mfem::GridFunction* gf);
  void DeregisterFieldInBPIndex(const std::string& field_name);

  void RegisterAttributeFieldInBPIndex(const std::string& attr_name);
  void DeregisterAttributeFieldInBPIndex(const std::string& attr_name);

  /** @brief Return a string with the conduit blueprint name for the given
      Element::Type. */
  std::string getElementName(mfem::Element::Type elementEnum);

  /** @brief Return an mfem::Geometry::Type for the given
      string with the conduit blueprint name. */
  // Why is there mfem::Element::Type and mfem::Geometry::Type? They look the same
  mfem::Geometry::Type getElementTypeFromName(const std::string& name);

  // Reconstructs a mesh using the current contents of the datastore
  // Used as part of Load()
  void reconstructMesh();

  // Reconstructs all non-mesh-related fields using the current contents
  // of the datastore, used as part of Load()
  void reconstructFields();

  /**
   * \brief A private helper function to set up the views associated with the
      data of a scalar valued grid function in the blueprint style.
   * \pre gf is not null
   * \note This function is expected to be called by RegisterField()
   * \note Handles cases where hierarchy is already set up,
   *      where the data was allocated by this data collection
   *      and where the gridfunction data is external to Sidre
   */
  void addScalarBasedGridFunction(const std::string& field_name,
                                  mfem::GridFunction* gf,
                                  const std::string& buffer_name,
                                  IndexType offset);

  /**
   * \brief A private helper function to set up the views associated with the
      data of a vector valued grid function in the blueprint style.
   * \pre gf is not null
   * \note This function is expected to be called by RegisterField()
   * \note Handles cases where hierarchy is already set up,
   *      where the data was allocated by this data collection
   *      and where the gridfunction data is external to Sidre
   */
  void addVectorBasedGridFunction(const std::string& field_name,
                                  mfem::GridFunction* gf,
                                  const std::string& buffer_name,
                                  IndexType offset);

  /** @brief A private helper function to set up the Views associated with
      attribute field named @a field_name */
  void addIntegerAttributeField(const std::string& field_name, bool is_bdry);

  /// Sets up the four main mesh blueprint groups.
  /**
   * \param hasBP Indicates whether the blueprint has already been set up.
   */
  void createMeshBlueprintStubs(bool hasBP);

  /// Sets up the mesh blueprint 'state' group.
  /**
   * \param hasBP Indicates whether the blueprint has already been set up.
   */
  void createMeshBlueprintState(bool hasBP);

  /// Sets up the mesh blueprint 'coordsets' group.
  /**
   * \param hasBP Indicates whether the blueprint has already been set up.
   */
  void createMeshBlueprintCoordset(bool hasBP);

  /// Sets up the mesh blueprint 'topologies' group.
  /**
   * This method is called from SetMesh().
   * \param hasBP Indicates whether the blueprint has already been set up.
   * \param mesh_name The name of the topology.
   * \note Valid values for @a mesh_name are "mesh" and "boundary" and the
           former has to be created with this method before the latter.
   */
  void createMeshBlueprintTopologies(bool hasBP, const std::string& mesh_name);

  #if defined(AXOM_USE_MPI) && defined(MFEM_USE_MPI)
  /// Sets up the mesh blueprint 'adjacencies' group.
  /**
   * \param hasBP Indicates whether the blueprint has already been set up.
   * \note Only valid when using parallel meshes
   */
  void createMeshBlueprintAdjacencies(bool hasBP);
  #endif

  // /// Verifies that the contents of the mesh blueprint data is valid.
  // void verifyMeshBlueprint();

  // The names for the mesh and boundary topologies in the blueprint group,
  // and the suffix used to store their attributes (as fields)
  static const std::string s_mesh_topology_name;
  static const std::string s_boundary_topology_name;
  static const std::string s_attribute_suffix;
  static const std::string s_coordset_name;
};

} /* namespace sidre */
} /* namespace axom */

#endif  // AXOM_USE_MFEM

#endif
