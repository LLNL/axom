// Copyright (c) 2017-2021, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_SPIN_BVHTREE_HPP_
#define AXOM_SPIN_BVHTREE_HPP_

// axom includes
#include "axom/core/Macros.hpp"  // for DISABLE_COPY_AND_ASSIGNMENT macro
#include "axom/core/Types.hpp"   // for nullptr

#include "axom/primal/geometry/BoundingBox.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/operators/squared_distance.hpp"

#include "axom/slic/interface/slic.hpp"  // for logging

// C/C++ includes
#include <fstream>  // for file stream
#include <queue>    // for std::queue
#include <sstream>  // std::ostringstream
#include <vector>   // for std::vector

namespace axom
{
namespace spin
{
/*!
 * \class BVHTree
 *
 * \brief The BVHTree class provides functionality for generating a
 *  Bounding Volume Hierarchy (BVH) of axis-aligned bounding boxes, i.e.,
 *  buckets that spatially partitions a set of objects. Each object is defined
 *  by its bounding box and user-supplied data-structure specified as a template
 *  argument. Once generated, the BVHTree can then be used as an
 *  acceleration structure to speed up the performance of point queries on the
 *  given objects, e.g., minimum distance, closest point, orientation, etc.
 *
 * \tparam T the type of data inserted in a BVH bin
 * \tparam NDIMS the number of dimensions
 *
 * \see BoundingBox, Point
 */
template <typename T, int NDIMS>
class BVHTree
{
public:
  using BoxType = primal::BoundingBox<double, NDIMS>;
  using PointType = primal::Point<double, NDIMS>;

public:
  /*!
   * \brief Default constructor, creates an empty BVHTree instance.
   * \param [in] maxNumLevel maximum number of subdivision levels. Default is 5.
   */
  BVHTree(int maxNumLevel = 5);

  /*!
   * \brief Creates a BVHTree instance that can grow up to the specified
   *  number of levels and can hold approximately the given number of objects.
   * \param [in] estNumObjects user-supplied estimated number of objects.
   * \param [in] maxNumLevels the maximum number of levels.
   * \note It is recommended to use this constructor for best-performance.
   */
  BVHTree(int estNumObjects, int maxNumLevels);

  /*!
   * \brief Destructor.
   */
  ~BVHTree();

  /*!
   * \brief Returns the maximum number of levels in the tree possible.
   * \return N the maximum number of levels in the tree.
   */
  int getMaxNumLevels() const { return m_maxNumLevels; };

  /*!
   * \brief Returns the number of levels in the given BVHTree instance.
   * \return N the number of levels in the BVHTree.
   */
  int getNumLevels() const { return m_numLevels; };

  /*!
   * \brief Inserts an object, identified by its bounding box and user-supplied
   *  data, to the bucket tree.
   * \param [in] box Axis-Aligned bounding box enclosing the object.
   * \param [in] data user-supplied data associated with the object.
   */
  void insert(const BoxType& box, const T& data);

  /*!
   * \brief Returns the total number of inserted objects in the tree.
   * \return N the total number of objects.
   * \post N >= 0.
   */
  int getNumberOfObjects() const { return static_cast<int>(m_objects.size()); };

  /*!
   * \brief Builds a hierarchical decomposition of axis-aligned bounding boxes
   *  that partitions the set of objects in different buckets.
   * \param [in] threshold subdivision threshold, buckets consisting of a number
   *  of objects greater than the given threshold will be subdivided.
   */
  void build(int threshold = 25);

  /*!
   * \brief Checks if the tree is empty.
   * \return status true if the tree is empty, otherwise, false.
   */
  bool empty() const { return m_tree.empty(); };

  /*!
   * \brief Deletes the tree and all associated objects.
   * \post this->empty() == true.
   * \post this->getNumberOfObjects() == 0.
   */
  void clear();

  /// \name Point Query API
  /// @{

  /*!
   * \brief Checks if the point is inside the BVHTree.
   * \param [in] pt the point in query.
   * \return status true if the point is inside, else false.
   */
  bool contains(const PointType& pt) const;

  /*!
   * \brief Finds a subset of candidate buckets for the given point query.
   * \param [in] pt the point in query.
   * \param [out] candidate_buckets list of bucket IDs for the given point
   * query.
   * \pre this->empty()==false.
   * \pre candidate_buckets.size()==0
   */
  void find(const PointType& pt, std::vector<int>& candidate_buckets) const;

  /// @}

  /// \name BVHTree Access methods
  /// @{

  /*!
   * \brief Returns a const reference to the bounding box of the given bucket.
   * \param [in] bucketIdx the index of the bucket in query.
   * \return Box the bounding box of the bucket.
   * \pre bucketIdx >= 0 && bucketIdx < m_tree.size()
   * \post Box.isValid() == true.
   */
  const BoxType& getBucketBox(int bucketIdx) const;

  /*!
   * \brief Returns the number of object within the given bucket.
   * \param [in] bucketIdx the index of the bucket in query.
   * \return N the number of object within the given bucket.
   * \pre bucketIdx >= 0 && bucketIdx < m_tree.size()
   * \post N >= 0
   */
  int getBucketNumObjects(int bucketIdx) const;

  /*!
   * \brief Return const pointer to the object array of the given bucket.
   * \param [in] bucketIdx the index of the bucket in query.
   * \return arrayPtr pointer to the buckets object array
   * \pre bucketIdx >= 0 && bucketIdx < m_tree.size()
   * \post arryPtr != nullptr
   */
  const int* getBucketObjectArray(int bucketIdx) const;

  /*!
   * \brief Returns a const reference to the bounding box of the given object.
   * \param [in] objIdx the ID of the object in query.
   * \return Box the bounding box of the given object.
   * \post Box.isValid() == true.
   */
  const BoxType& getObjectBox(int objIdx) const;

  /*!
   * \brief Returns const reference to the data associated with the given
   * object.
   * \param objIdx the ID of the object in query.
   * \return data the data associated with the object.
   * \note The return data is the data supplied when the object was inserted.
   * \pre objIdx >= 0 && objIdx < this->getNumberOfObjects()
   */
  const T& getObjectData(int objIdx) const;

  /*!
   * \brief Returns the bucket index of the given object.
   * \param [in] objIdx the ID of the object in query.
   * \return bucketIdx the bucket index.
   * \pre objIdx >= 0 && objIdx < this->getNumberOfObjects()
   * \post The return value, butcketIdx will be:
   *  <ul>
   *   <li> bucketIdx < 0, iff the tree is not built </li>
   *   <li> bucketIdx >= 0 && bucektIdx < this->getNumberO</li>
   *  </ul>
   */
  int getObjectBucketIndex(int objIdx) const;

  /// @}

  /*!
   * \brief Writes the hierarchical decomposition in a VTK file.
   * \param [in] fileName the name of the VTK formatted file to generate.
   * \param [in] include_objects optional parameter, when true it also writes
   *  out the bounding boxes of the supplied objects.
   */
  void writeVtkFile(const std::string& fileName,
                    bool include_objects = false) const;

  /*!
   * \brief Writes the user-supplied set of bins in a VTK file.
   * \param [in] fileName the name of the VTK formatted file to generate.
   * \param [in] bins array of bins to dump in the VTK file.
   * \param [in] nbins the number of bins.
   * \note Primarily used for debugging.
   */
  void writeVtkFile(const std::string& fileName, const int* bins, int nbins) const;

private:
  /// \name Internal Helper Routines for Binary Tree traversal
  /// @{

  /*!
   * \brief Given a flat index of a bucket, find its parent index in the tree.
   * \param [in] i the index of the bucket in query.
   * \return p the flat index of the parent bucket.
   * \post p==-1 if i corresponds to the root bucket, otherwise, p >= 0.
   * \pre i >= 0
   */
  static int parent(int i)
  {
    SLIC_ASSERT(i >= 0);
    return (((i - 1) >> 1));
  };

  /*!
   * \brief Given a flat index of a bucket, find its right child tree index.
   * \param [in] i the index of the bucket in query.
   * \return r the flat index of the right child bucket.
   * \pre i >= 0
   * \post r > 0
   */
  static int right_child(int i)
  {
    SLIC_ASSERT(i >= 0);
    return ((i << 1) + 2);
  };

  /*!
   * \brief Given the flat index of a bucket, find its left child tree index.
   * \param [in] i the index of the bucket in query.
   * \return l the flat index of the left child bucket.
   * \pre i >= 0
   * \post l > 0
   */
  static int left_child(int i)
  {
    SLIC_ASSERT(i >= 0);
    return ((i << 1) + 1);
  };

  /*!
   * \brief Computes the value of 2 raised to the nth power.
   * \param n user-supplied power.
   * \note Uses bit shifting to avoid potentially expensive calls to std::pow().
   * \return result the value of \f$ 2^n \f$
   * \pre n >= 0
   */
  static int pow2(int n)
  {
    SLIC_ASSERT(n >= 0);
    return (1 << n);
  };

  /*!
   * \brief Computes the linear index of the 1st bucket within the level.
   * \param [in] l the level index in query (starting from 0).
   * \return idx the index of the 1st bucket within the level.
   * \pre l >= 0
   */
  static int level_begin(int l)
  {
    SLIC_ASSERT(l >= 0);
    return (BVHTree::pow2(l) - 1);
  }

  /*!
   * \brief Computes the linear index of the last bucket within the level.
   * \param [in] l the level index in query (starting from 0).
   * \return idx the index of the last bucket within the level.
   * \pre l >= 0
   */
  static int level_end(int l)
  {
    SLIC_ASSERT(l >= 0);
    return (BVHTree::pow2(l + 1) - 2);
  }

  /*!
   * \brief Writes the box to the given C++ stream.
   * \param [in,out] oss the output string stream to write the box coordinates.
   * \param [in] box a 2-D bounding box instance.
   * \note This is just a helper routine used when writing out the tree.
   */
  static int write_box(std::ostringstream& oss,
                       const primal::BoundingBox<double, 2>& box)
  {
    const PointType& min = box.getMin();
    const PointType& max = box.getMax();

    oss << min[0] << " " << min[1] << " 0.0\n";
    oss << max[0] << " " << min[1] << " 0.0\n";
    oss << max[0] << " " << max[1] << " 0.0\n";
    oss << min[0] << " " << max[1] << " 0.0\n";
    oss << std::endl;

    return 4;
  }

  /*!
   * \brief Writes the box to the given C++ stream.
   * \param [in,out] oss the output string stream to write the box coordinates.
   * \param [in] box a 3-D bounding box instance.
   * \note This is just a helper routine used when writing out the tree.
   */
  static int write_box(std::ostringstream& oss,
                       const primal::BoundingBox<double, 3>& box)
  {
    const PointType& min = box.getMin();
    const PointType& max = box.getMax();

    oss << min[0] << " " << min[1] << " " << min[2] << std::endl;
    oss << max[0] << " " << min[1] << " " << min[2] << std::endl;
    oss << max[0] << " " << max[1] << " " << min[2] << std::endl;
    oss << min[0] << " " << max[1] << " " << min[2] << std::endl;

    oss << min[0] << " " << min[1] << " " << max[2] << std::endl;
    oss << max[0] << " " << min[1] << " " << max[2] << std::endl;
    oss << max[0] << " " << max[1] << " " << max[2] << std::endl;
    oss << min[0] << " " << max[1] << " " << max[2] << std::endl;

    return 8;
  }

  /// @}

  /*!
   * \brief Helper routine used to determine if bucket should be refined.
   * \param [in] idx index of the bucket in the tree.
   * \param [in] threshold allowable number of objects in a bucket.
   * \return status true if the bucket should be refined, otherwise false.
   * \note A bucket is refined if the following 2 criteria are satisfied:
   *  <ol>
   *   <li> numObjectsInBucket > threshold </li>
   *   <li> bucketLevel < maxNumLevels-1 </li>
   *  </ol>
   * \note Used in BVHTree::build()
   * \pre idx >= 0 && idx < m_tree.size()
   */
  bool shouldRefine(int idx, int threshold);

  /*!
   * \brief Helper routine used to refine a bucket in the tree.
   * \param [in] idx index of the bucket in the tree.
   * \note Used in BVHTree::build()
   * \pre idx >= 0 && idx  < m_tree.size()
   */
  void refine(int idx);

  /*!
   * \brief Resizes this BVHTree instance (as needed) to store the level.
   * \param [in] level the level to resize the tree to.
   * \note Helper method used in refine().
   * \note Method resizes as needed, otherwise, this method returns immediately.
   */
  void resizeToLevel(int level);

  /*!
   * \brief Propagates the data from the parent bucket to the children.
   * \param [in] parent the ID of the parent bucket.
   * \param [in] rChild the ID of the right child.
   * \param [in] lChild the ID of the left child.
   * \note Helper method used in refine()
   * \pre parent >= 0 && parent < m_tree.size()
   * \pre rChild >= 0 && rChild < m_tree.size()
   * \pre lChild >= 0 && lChild < m_tree.size()
   * \pre m_tree[ rChild ].Level > m_tree[ parent ].Level
   * \pre m_tree[ lChild ].Level > m_tree[ parent ].Level
   */
  void percolateDown(int parent, int rChild, int lChild);

  /*!
   * \brief Computes the min squared distance of a point to the given bucket.
   * \param [in] bucketIdx index of the bucket in query.
   * \param [in] pt the user-supplied point in query.
   * \return min the minimum squared distance to the bucket.
   * \pre bucketIdx >= 0 && bucketIdx < m_tree.size()
   */
  double getMinSqDistanceToBucket(int bucketIdx, const PointType& pt) const;

  /*!
   * \brief Computes the max squared distance of a point to the given bucket.
   * \param [in] bucketIdx index of the bucket in query.
   * \param [in] pt the user-supplied point in query.
   * \return max the maximum squared distance to the bucket.
   * \pre bucketIdx >= 0 && bucketIdx < m_tree.size()
   */
  double getMaxSqDistanceToBucket(int bucketIdx, const PointType& pt) const;

  /// \name Internal Data-Structure Definitions
  /// @{

  struct Object
  {
    BoxType Box;
    T Data;
    int BucketIdx;
  };

  struct Bucket
  {
    int Level;
    bool Void;
    bool Refined;
    BoxType Box;
    std::vector<int> ObjectArray;
  };

  /// @}

  std::vector<Object> m_objects; /*!< List of user-supplied objects   */
  std::vector<Bucket> m_tree;    /*!< Flat binary tree decomposition  */

  int m_numLevels;    /*!< current number of levels        */
  int m_maxNumLevels; /*!< maximum level of the tree       */

  DISABLE_COPY_AND_ASSIGNMENT(BVHTree);
  DISABLE_MOVE_AND_ASSIGNMENT(BVHTree);
};

} /* namespace spin */
}  // namespace axom

//------------------------------------------------------------------------------
//           BVHTree Public API Implementation
//------------------------------------------------------------------------------
namespace axom
{
namespace spin
{
template <typename T, int NDIMS>
BVHTree<T, NDIMS>::BVHTree(int maxNumLevels)
  : m_numLevels(0)
  , m_maxNumLevels(maxNumLevels)

{
  // pre-allocate internal data-structures.
  m_objects.reserve(100);

  int numBuckets = BVHTree::pow2(m_maxNumLevels) - 1;
  m_tree.reserve(numBuckets);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BVHTree<T, NDIMS>::BVHTree(int estNumObjects, int maxNumLevels)
  : m_numLevels(0)
  , m_maxNumLevels(maxNumLevels)

{
  // pre-allocate internal data-structures.
  m_objects.reserve(estNumObjects);

  int numBuckets = BVHTree::pow2(m_maxNumLevels) - 1;
  m_tree.reserve(numBuckets);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
BVHTree<T, NDIMS>::~BVHTree()
{
  this->clear();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::insert(const BoxType& box, const T& data)
{
  // Sanity checks
  SLIC_ASSERT(box.isValid());

  // STEP 0: create a BVHTree Object instance.
  Object obj;
  obj.Box = box;
  obj.Data = data;
  obj.BucketIdx = 0;
  SLIC_ASSERT(box == obj.Box);

  // STEP 1: insert the object to the global object list
  m_objects.push_back(obj);
  const int idx = static_cast<int>(m_objects.size()) - 1;

  // STEP 2: create the tree if empty
  if(this->empty())
  {
    ++m_numLevels;

    m_tree.resize(1);

    m_tree[0].Level = 0;
    m_tree[0].Refined = false;
    m_tree[0].Void = false;
    m_tree[0].ObjectArray.reserve(m_objects.size());
  }  // END if

  // STEP 3: Update bounds of the root bounding box to encapsulate the
  // bounding box of the object.
  m_tree[0].Box.addBox(obj.Box);

  // STEP 4: push back the object to the object array
  m_tree[0].ObjectArray.push_back(idx);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::resizeToLevel(int level)
{
  if(level < m_numLevels)
  {
    /* short-circuit */
    return;
  }

  const int begin = BVHTree::level_begin(level);
  const int end = BVHTree::level_end(level);
  if(end >= static_cast<int>(m_tree.size()))
  {
    const int new_size = end + 1;
    m_tree.resize(new_size);

    // initialize the buckets at the new level
    for(int i = begin; i < new_size; ++i)
    {
      m_tree[i].Level = level;
      m_tree[i].Void = true;
      m_tree[i].Refined = false;
    }  // END for

    ++m_numLevels;

  }  // END if
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::percolateDown(int parent, int rChild, int lChild)
{
  // Sanity Checks
  SLIC_ASSERT(parent >= 0 && parent < static_cast<int>(m_tree.size()));
  SLIC_ASSERT(rChild >= 0 && rChild < static_cast<int>(m_tree.size()));
  SLIC_ASSERT(lChild >= 0 && lChild < static_cast<int>(m_tree.size()));
  SLIC_ASSERT(m_tree[parent].Level < m_tree[rChild].Level);
  SLIC_ASSERT(m_tree[parent].Level < m_tree[lChild].Level);
  SLIC_ASSERT(!m_tree[parent].Void);

  // STEP 0: Get the bounding boxes of the two children
  const BoxType rightBox = m_tree[rChild].Box;
  const BoxType leftBox = m_tree[lChild].Box;

  // STEP 1: pre-allocated space for object at the children
  const int numObjects = static_cast<int>(m_tree[parent].ObjectArray.size());
  const int estChildSize = static_cast<int>(0.5 * numObjects);
  m_tree[rChild].ObjectArray.reserve(estChildSize);
  m_tree[lChild].ObjectArray.reserve(estChildSize);

  // STEP 2: loop over objects and place on right or left bucket.
  for(int i = 0; i < numObjects; ++i)
  {
    // STEP A: get index in to the object list
    const int objIdx = m_tree[parent].ObjectArray[i];
    SLIC_ASSERT((objIdx >= 0) && (objIdx < static_cast<int>(m_objects.size())));

    // STEP B: get bounding box of the object.
    const BVHTree::BoxType objBox = m_objects[objIdx].Box;
    SLIC_ASSERT(objBox.isValid());

    // STEP C: get centroid of the object's bounding box
    const BVHTree::PointType centroid = objBox.getCentroid();

    int bucketIdx = -1;
    if(rightBox.contains(centroid))
    {
      bucketIdx = rChild;
    }
    else if(leftBox.contains(centroid))
    {
      bucketIdx = lChild;
    }
    else
    {
      SLIC_ERROR("Failed to place object in bucket!");
    }

    SLIC_ASSERT((bucketIdx >= 0) && (bucketIdx < static_cast<int>(m_tree.size())));

    m_tree[bucketIdx].Box.addBox(objBox);
    m_tree[bucketIdx].ObjectArray.push_back(objIdx);
    m_objects[objIdx].BucketIdx = bucketIdx;
  }  // END for all objects

  // STEP 3: unmark children as void
  if(m_tree[rChild].ObjectArray.size() > 0)
  {
    m_tree[rChild].Void = false;
  }

  if(m_tree[lChild].ObjectArray.size() > 0)
  {
    m_tree[lChild].Void = false;
  }
  m_tree[parent].ObjectArray.clear();
  m_tree[parent].Refined = true;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::refine(int idx)
{
  // Sanity checks
  SLIC_ASSERT(idx >= 0 && idx < static_cast<int>(m_tree.size()));
  SLIC_ASSERT(m_tree[idx].Box.isValid());
  SLIC_ASSERT(!m_tree[idx].Refined);
  SLIC_ASSERT(!m_tree[idx].Void);

  // STEP 0: get current & next level
  const int currentLevel = m_tree[idx].Level;
  const int nextLevel = currentLevel + 1;
  SLIC_ASSERT(nextLevel < m_maxNumLevels);

  // STEP 1: if needed, resize the m_tree for the next level.
  this->resizeToLevel(nextLevel);

  // STEP 2: get bucket index of right & left children and mark them valid
  const int rightChild = BVHTree::right_child(idx);
  const int leftChild = BVHTree::left_child(idx);
  SLIC_ASSERT(rightChild < static_cast<int>(m_tree.size()));
  SLIC_ASSERT(leftChild < static_cast<int>(m_tree.size()));

  // STEP 3: split bounding box
  // TODO: think about better strategies for figuring out the location of the
  // cut-plane to bisect, e.g., center-of-mass density, etc.
  m_tree[idx].Box.bisect(m_tree[rightChild].Box, m_tree[leftChild].Box);
  SLIC_ASSERT(m_tree[rightChild].Box.isValid());
  SLIC_ASSERT(m_tree[leftChild].Box.isValid());

  // STEP 5: percolate the data from the parent to the children
  this->percolateDown(idx, rightChild, leftChild);
  SLIC_ASSERT(this->getBucketNumObjects(idx) == 0);
  SLIC_ASSERT(m_tree[idx].Refined);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
bool BVHTree<T, NDIMS>::shouldRefine(int idx, int threshold)
{
  SLIC_ASSERT(idx >= 0 && idx < static_cast<int>(m_tree.size()));
  const int numObjects = static_cast<int>(m_tree[idx].ObjectArray.size());
  const int level = m_tree[idx].Level;

  if((numObjects > threshold) && (level < m_maxNumLevels - 1))
  {
    return true;
  }

  return false;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::build(int threshold)
{
  SLIC_ASSERT(threshold >= 1);

  // STEP 0: check if the tree is empty
  if(this->empty())
  {
    /* short-circuit, tree is empty */
    return;
  }

  // STEP 1: put root bucket in the queue
  std::queue<int> workqueue;
  workqueue.push(0);

  //  int count = 0;
  //  std::ostringstream oss;
  //  oss << "bucket_tree_" << count << ".vtk";
  //  this->writeLegacyVtkFile( oss.str() );
  //  ++count;

  // STEP 2: Iteratively subdivide until threshold criteria is satisfied
  while(!workqueue.empty())
  {
    // STEP A: pop next bucket from the queue
    int bucketIdx = workqueue.front();
    workqueue.pop();

    // STEP B: check if we should refine this bucket
    if(this->shouldRefine(bucketIdx, threshold))
    {
      // STEP C: refine bucket
      this->refine(bucketIdx);
      SLIC_ASSERT(m_tree[bucketIdx].Refined);

      // STEP D: update workqueue accordingly
      int leftIdx = BVHTree::left_child(bucketIdx);
      int rightIdx = BVHTree::right_child(bucketIdx);

      if(this->shouldRefine(leftIdx, threshold))
      {
        workqueue.push(leftIdx);
      }

      if(this->shouldRefine(rightIdx, threshold))
      {
        workqueue.push(rightIdx);
      }
    }  // END if > threshold

    //      oss.clear();
    //      oss.str("");
    //      oss << "bucket_tree_" << count << ".vtk";
    //      this->writeLegacyVtkFile( oss.str() );
    //      ++count;

  }  // END while
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::clear()
{
  m_objects.clear();
  m_tree.clear();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
bool BVHTree<T, NDIMS>::contains(const PointType& pt) const
{
  if(this->empty())
  {
    return false;
  }

  return (m_tree[0].Box.contains(pt));
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline double BVHTree<T, NDIMS>::getMinSqDistanceToBucket(int bucketIdx,
                                                          const PointType& pt) const
{
  SLIC_ASSERT(bucketIdx >= 0 && bucketIdx < static_cast<int>(m_tree.size()));

  double dist = std::numeric_limits<double>::max();

  if(!m_tree[bucketIdx].Void)
  {
    dist = squared_distance(pt, m_tree[bucketIdx].Box);
  }

  return (dist);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
inline double BVHTree<T, NDIMS>::getMaxSqDistanceToBucket(int bucketIdx,
                                                          const PointType& pt) const
{
  SLIC_ASSERT(bucketIdx >= 0 && bucketIdx < static_cast<int>(m_tree.size()));

  double dist = std::numeric_limits<double>::min();

  if(!m_tree[bucketIdx].Void)
  {
    std::vector<PointType> pnts;
    BoxType::getPoints(m_tree[bucketIdx].Box, pnts);
    SLIC_ASSERT(pnts.size() == 4 || pnts.size() == 8);

    const int npoints = static_cast<int>(pnts.size());
    for(int i = 0; i < npoints; ++i)
    {
      dist = std::max(dist, squared_distance(pt, pnts[i]));
    }  // END for all points
  }

  return dist;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::find(const PointType& pt,
                             std::vector<int>& candidate_buckets) const
{
  /* Sanity checks */
  SLIC_ASSERT(candidate_buckets.size() == 0);

  // STEP 0: if empty, or just a single level return immediately
  if(this->empty())
  {
    /* short-circuit */
    return;
  }
  else if(m_numLevels == 1)
  {
    /* short-circuit */
    candidate_buckets.push_back(0);
    return;
  }

  const double TOL = 1.0e-9;

  int globalClosestBucketIdx = -1;
  double estLowerBound = std::numeric_limits<double>::max();
  double estUpperBound = std::numeric_limits<double>::min();
  std::vector<int> first_set_candidates;
  first_set_candidates.reserve(BVHTree::pow2(m_numLevels - 1));

  // STEP 1: pre-allocate buffer storage for buckets to check at each level.
  std::vector<int> buckets_to_check;
  buckets_to_check.reserve(BVHTree::pow2(m_numLevels - 1));

  // STEP 2: pre-populate buffer with the two buckets in level one
  buckets_to_check.push_back(1);
  buckets_to_check.push_back(2);

  // STEP 3: iteratively descend down in the hierarchy, pruning away
  // buckets that are too far away from the query point, pt, until a
  // first set of candidate leaf buckets is attained. This also, computes
  // the estLowerBound and estUpperBound
  for(int level = 1; level < m_numLevels; ++level)
  {
    const int nbuckets = static_cast<int>(buckets_to_check.size());
    if(nbuckets == 0)
    {
      /* short-circuit */
      break;
    }

    // STEP 3.A: find closest bucket, and lowerDistBound
    int closest_bucket = -1;
    double lowerDistBound = std::numeric_limits<double>::max();
    for(int j = 0; j < nbuckets; ++j)
    {
      const int bucketIdx = buckets_to_check[j];
      SLIC_ASSERT(m_tree[bucketIdx].Level == level);

      const double d2b = this->getMinSqDistanceToBucket(bucketIdx, pt);
      if(d2b < lowerDistBound)
      {
        lowerDistBound = d2b;
        closest_bucket = bucketIdx;
      }
    }  // END for all buckets

    SLIC_ASSERT(closest_bucket >= 0 &&
                closest_bucket < static_cast<int>(m_tree.size()));

    // STEP 3.B: get upper distance bound -- max distance to closest bucket
    double upperDistBound = this->getMaxSqDistanceToBucket(closest_bucket, pt);

    // STEP 3.C: find candidate buckets & prune buckets that are too far
    for(int j = 0; j < nbuckets; ++j)
    {
      const int bucketIdx = buckets_to_check[j];
      bool keep = false;

      if(bucketIdx == closest_bucket)
      {
        keep = true;
      }
      else
      {
        const double d2b = this->getMinSqDistanceToBucket(bucketIdx, pt);
        if(d2b <= (upperDistBound + TOL))
        {
          keep = true;
        }

      }  // END else

      if(keep && m_tree[bucketIdx].Refined)
      {
        buckets_to_check.push_back(BVHTree::left_child(bucketIdx));
        buckets_to_check.push_back(BVHTree::right_child(bucketIdx));

      }  // END if
      else if(keep && !m_tree[bucketIdx].Void)
      {
        first_set_candidates.push_back(bucketIdx);

        double minDist = this->getMinSqDistanceToBucket(bucketIdx, pt);
        if(minDist < estLowerBound)
        {
          globalClosestBucketIdx = bucketIdx;
          estLowerBound = minDist;
          estUpperBound = this->getMaxSqDistanceToBucket(bucketIdx, pt);
        }

      }  // END else if

    }  // END for all buckets

    // STEP 3.D: remove the first nbuckets from the list of buckets to check
    buckets_to_check.erase(buckets_to_check.begin(),
                           buckets_to_check.begin() + nbuckets);

  }  // END for all levels

  SLIC_ASSERT(globalClosestBucketIdx >= 0 &&
              globalClosestBucketIdx < static_cast<int>(m_tree.size()));

  // STEP 4: further filter first_set_candidate buckets into the final set
  // of candidate buckets to pass back to the caller.
  const int numCandidates = static_cast<int>(first_set_candidates.size());
  candidate_buckets.clear();
  candidate_buckets.reserve(numCandidates);
  candidate_buckets.push_back(globalClosestBucketIdx);

  for(int i = 0; i < numCandidates; ++i)
  {
    const int bucketIdx = first_set_candidates[i];
    if(bucketIdx != globalClosestBucketIdx)
    {
      double minDist = this->getMinSqDistanceToBucket(bucketIdx, pt);
      if(minDist <= (estUpperBound + TOL))
      {
        candidate_buckets.push_back(bucketIdx);
      }
    }

  }  // END for all first set candidates
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
const primal::BoundingBox<double, NDIMS>& BVHTree<T, NDIMS>::getBucketBox(
  int bucketIdx) const
{
  SLIC_ASSERT(bucketIdx >= 0 && bucketIdx < static_cast<int>(m_tree.size()));
  return m_tree[bucketIdx].Box;
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
int BVHTree<T, NDIMS>::getBucketNumObjects(int bucketIdx) const
{
  SLIC_ASSERT(bucketIdx >= 0 && bucketIdx < static_cast<int>(m_tree.size()));
  return (static_cast<int>(m_tree[bucketIdx].ObjectArray.size()));
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
const int* BVHTree<T, NDIMS>::getBucketObjectArray(int bucketIdx) const
{
  SLIC_ASSERT(bucketIdx >= 0 && bucketIdx < static_cast<int>(m_tree.size()));
  return (&m_tree[bucketIdx].ObjectArray[0]);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
const primal::BoundingBox<double, NDIMS>& BVHTree<T, NDIMS>::getObjectBox(
  int objIdx) const
{
  SLIC_ASSERT(objIdx >= 0 && objIdx < static_cast<int>(m_objects.size()));
  return (m_objects[objIdx].Box);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
const T& BVHTree<T, NDIMS>::getObjectData(int objIdx) const
{
  SLIC_ASSERT(objIdx >= 0 && objIdx < static_cast<int>(m_objects.size()));
  return (m_objects[objIdx].Data);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
int BVHTree<T, NDIMS>::getObjectBucketIndex(int objIdx) const
{
  SLIC_ASSERT(objIdx >= 0 && objIdx < static_cast<int>(m_objects.size()));
  return (m_objects[objIdx].BucketIdx);
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::writeVtkFile(const std::string& fileName,
                                     const int* bins,
                                     int nbins) const
{
  // STEP 0: Write VTK header
  std::ofstream ofs;
  ofs.open(fileName.c_str());
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " BVHTree \n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  std::ostringstream coordinates;
  std::ostringstream cells;

  int id = 0;
  for(int ibin = 0; ibin < nbins; ++ibin)
  {
    const int bucketIdx = bins[ibin];
    const int nnodes = BVHTree::write_box(coordinates, m_tree[bucketIdx].Box);
    coordinates << std::endl;

    cells << nnodes << " ";
    for(int i = 0; i < nnodes; ++i, ++id)
    {
      cells << id << " ";
    }
    cells << std::endl;
  }

  ofs << "POINTS " << id << " double\n";
  ofs << coordinates.str() << std::endl;

  int nnodes = (NDIMS == 2) ? 4 : 8;
  ofs << "CELLS " << nbins << " " << nbins * (nnodes + 1) << std::endl;
  ofs << cells.str() << std::endl;

  ofs << "CELL_TYPES " << nbins << std::endl;
  int cellType = (NDIMS == 2) ? 9 : 12;
  for(int i = 0; i < nbins; ++i)
  {
    ofs << cellType << std::endl;
  }  // END for all cells
  ofs << std::endl;

  ofs.close();
}

//------------------------------------------------------------------------------
template <typename T, int NDIMS>
void BVHTree<T, NDIMS>::writeVtkFile(const std::string& fileName,
                                     bool AXOM_NOT_USED(include_objects)) const
{
  // STEP 0: Write VTK header
  std::ofstream ofs;
  ofs.open(fileName.c_str());
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << " BVHTree \n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  std::ostringstream coordinates;
  std::ostringstream cells;
  std::ostringstream cell_level_flag;
  std::ostringstream cell_index_flag;

  // Loop through all levels
  int id = 0;
  int numCells = 0;

  for(int i = 0; i < m_numLevels; ++i)
  {
    int start = BVHTree::level_begin(i);
    int end = BVHTree::level_end(i);
    SLIC_ASSERT(start >= 0 && start < static_cast<int>(m_tree.size()));
    SLIC_ASSERT(end >= 0 && end < static_cast<int>(m_tree.size()));

    for(int j = start; j <= end; ++j)
    {
      if(m_tree[j].Void)
      {
        continue;
      }

      const int nnodes = BVHTree::write_box(coordinates, m_tree[j].Box);
      coordinates << std::endl;

      cells << nnodes << " ";
      for(int k = 0; k < nnodes; ++k)
      {
        cells << id << " ";
        ++id;
      }

      cells << std::endl << std::endl;

      cell_level_flag << i << std::endl;
      cell_index_flag << numCells << std::endl;
      ++numCells;

    }  // END for all buckets within this level

  }  // END for all levels

  ofs << "POINTS " << id << " double\n";
  ofs << coordinates.str() << std::endl;

  int nnodes = (NDIMS == 2) ? 4 : 8;
  ofs << "CELLS " << numCells << " " << numCells * (nnodes + 1) << std::endl;
  ofs << cells.str() << std::endl;

  ofs << "CELL_TYPES " << numCells << std::endl;
  int cellType = (NDIMS == 2) ? 9 : 12;
  for(int i = 0; i < numCells; ++i)
  {
    ofs << cellType << std::endl;

  }  // END for all cells
  ofs << std::endl;

  ofs << "CELL_DATA " << numCells << std::endl;

  ofs << "SCALARS level int\n";
  ofs << "LOOKUP_TABLE default\n";
  ofs << cell_level_flag.str() << std::endl;
  ofs << std::endl;

  ofs << "SCALARS ID int\n";
  ofs << "LOOKUP_TABLE default\n";
  ofs << cell_index_flag.str() << std::endl;
  ofs << std::endl;

  ofs.close();
}

}  // namespace spin
}  // namespace axom

#endif  // AXOM_SPIN_BVHTREE_HPP_
