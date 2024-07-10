// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef SINA_DATAHOLDER_HPP
#define SINA_DATAHOLDER_HPP

/**
 * @file
 *
 * Contains the definition of the DataHolder class. 
 */

#include <string>
#include <memory>
#include <unordered_map>

#include "conduit.hpp"

#include "axom/sina/include/Datum.hpp"
#include "axom/sina/include/CurveSet.hpp"

namespace axom
{
namespace sina
{

/**
 * A DataHolder is a basic container for certain types of information.
 *
 * DataHolders contain curves, libraries, and data (Datum), and represent
 * all the information a library can have associated with it. Records expand
 * on DataHolders to contain additional info.
 *
 * \see Record
 * \see LibraryDataMap
 */
class DataHolder {
public:
    /**
     * An unordered map of Datum objects.
     */
    using DatumMap = std::unordered_map<std::string, Datum>;

    /**
     * An unordered map of CurveSet objects.
     */
    using CurveSetMap = std::unordered_map<std::string, CurveSet>;

    /**
     * An unordered map of shared pointers to DataHolder objects.
     */
    using LibraryDataMap = std::unordered_map<std::string, std::shared_ptr<DataHolder>>;

    /**
     * Construct an empty DataHolder.
     */
    DataHolder() = default;

    /**
     * Virtual destructor to automatically clean up resources held by an instance of the DataHolder class.
     */
    virtual ~DataHolder() = default;

    /**
     * Copy constructor that disallows this constructor type.
     */
    DataHolder(DataHolder const &) = delete;

    /**
     * Disable copy assignment.
     */
    DataHolder &operator=(DataHolder const &) = delete;

    /**
     * Construct a DataHolder from its conduit Node representation.
     *
     * @param asNode the DataHolder as a Node
     */
    explicit DataHolder(conduit::Node const &asNode);

    /**
     * Get the DataHolder's data.
     *
     * @return the DataHolder's data
     */
    DatumMap const &getData() const noexcept {
        return data;
    }

    /**
     * Add a Datum to this DataHolder.
     *
     * @param name the key for the Datum to add
     * @param datum the Datum to add
     */
    void add(std::string name, Datum datum);

    /**
     * Add a CurveSet to this DataHolder.
     *
     * @param curveSet the CurveSet to add
     */
    void add(CurveSet curveSet);

    /**
     * Get the curve sets associated with this DataHolder.
     *
     * @return the dataholder's curve sets
     */
    CurveSetMap const &getCurveSets() const noexcept {
        return curveSets;
    }

    /**
     * Add a new library to this DataHolder.
     *
     * If you try to add a library with a name that already exists, the old
     * library will be replaced.
     *
     * @return a pointer to a new DataHolder for a library
     * of the given name.
     */
    std::shared_ptr<DataHolder> addLibraryData(std::string const &name);

    /**
     * Add a new library to this DataHolder with existing library data.
     *
     * @return a pointer to a new DataHolder for a library of the given name.
     */
    std::shared_ptr<DataHolder> addLibraryData(std::string const &name, conduit::Node existingLibraryData);

    /**
     * Get all library data associated with this DataHolder.
     *
     * @return the dataholder's library data
     */
    LibraryDataMap const &getLibraryData() const noexcept {
        return libraryData;
    }

    /**
     * Get a specific library associated with this DataHolder.
     *
     * @return the dataholder's library data
     */
    std::shared_ptr<DataHolder> getLibraryData(std::string const &libraryName) {
      return libraryData.at(libraryName);
    }

    /**
     * Get a specific library associated with this DataHolder.
     * Used when the object on which this function is called is a const.
     *
     * @return the dataholder's library data
     */
    std::shared_ptr<DataHolder> const getLibraryData(std::string const &libraryName) const {
      return libraryData.at(libraryName);
    }

    /**
     * Get the user-defined content of the object.
     *
     * @return the user-defined content
     */
    conduit::Node const &getUserDefinedContent() const noexcept {
        return userDefined;
    }

    /**
     * Get the user-defined content of the object.
     *
     * @return the user-defined content
     */
    conduit::Node &getUserDefinedContent() noexcept {
        return userDefined;
    }

    /**
     * Set the user-defined content of the object.
     *
     * @param userDefined the user-defined content. Must be an object (key/value pairs)
     */
    void setUserDefinedContent(conduit::Node userDefined);

    /**
     * Convert this DataHolder to its conduit Node representation.
     *
     * @return the Node representation of this DataHolder.
     */
    virtual conduit::Node toNode() const;


private:
    CurveSetMap curveSets;
    DatumMap data;
    LibraryDataMap libraryData;
    conduit::Node userDefined;
};

}  // end sina namespace
}  // end axom namespace

#endif //SINA_DATAHOLDER_HPP
