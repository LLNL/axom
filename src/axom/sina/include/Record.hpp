// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)


#ifndef SINA_RECORD_HPP
#define SINA_RECORD_HPP

/// @file

#include <functional>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "conduit.hpp"

#include "axom/sina/include/ID.hpp"
#include "axom/sina/include/DataHolder.hpp"
#include "axom/sina/include/CurveSet.hpp"
#include "axom/sina/include/Datum.hpp"
#include "axom/sina/include/File.hpp"

namespace axom
{
namespace sina
{

/**
 * FileEqualByURI is used to store files in a Record.
 */
struct FileEqualByURI {
    bool operator()(const File &file1, const File &file2) const {
        return file1.getUri() == file2.getUri();
    }
};

/**
 * FileHashByURI is used to store files in a Record. Files are stored according
 * to the hash of their URI.
 */
struct FileHashByURI {
    size_t operator()(const File &file) const {
        return std::hash<std::string>()(file.getUri());
    }
};

/**
 * The Record class represents an entry in a Document's Record list. Records represent the data to be stored
 * (as opposed to the relationships between data)--natural scopes for Records include things like a single run
 * of an application, an msub job, a cluster of runs that has some metadata attached to the cluster (this Record
 * might have a "contains" Relationship for all the runs within it), etc.
 *
 * Each Record must have a type and an id. Each Record can also have a list of
 * File objects and a map of Datum objects.
 *
 * \code
 * axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
 * std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};
 * std::vector<std::string> myTags{"input"};
 * axom::sina::Datum myDatum{12, myTags};
 * myRecord->add("my_scalar",std::move(myDatum));
 * std::cout << myRecord->toNode().to_json() << std::endl;
 * \endcode
 *
 * The output would be:
 * \code{.json}
 * {"local_id":"my_record","type":"my_type","data":{"my_scalar":{"tags":["input"],"value":12.0}}}
 * \endcode
 */
class Record : public DataHolder {
public:
    /**
     * An unordered set of File objects.
     */
    using FileSet = std::unordered_set<File, FileHashByURI, FileEqualByURI>;

    /**
     * Construct a new Record.
     *
     * @param id the ID of the record
     * @param type the type of the record
     */
    Record(ID id, std::string type);

    /**
     * Construct a Record from its conduit Node representation.
     *
     * @param asNode the Record as a Node
     */
    explicit Record(conduit::Node const &asNode);

    /**
     * Disable the copy constructor.
     */
    Record(Record const &) = delete;

    /**
     * Disable copy assignment.
     */
    Record &operator=(Record const &) = delete;

    /**
     * Get the Record's ID.
     *
     * @return the ID
     */
    ID const &getId() const noexcept {
        return id.getID();
    }

    /**
     * Get the Record's type.
     *
     * @return the Record's type
     */
    std::string const &getType() const noexcept {
        return type;
    }

    /**
     * Remove a File from this record.
     *
     * @param file the File to remove
     */
    void remove(File const& file);



    using DataHolder::add;
    /**
     * Add a File to this record.
     *
     * @param file the File to add
     */
    void add(File file);


    /**
     * Get the files associated with this record.
     *
     * @return the record's files
     */
    FileSet const &getFiles() const noexcept {
        return files;
    }

    /**
     * Convert this record to its conduit Node representation.
     *
     * @return the Node representation of this record.
     */
    conduit::Node toNode() const override;

   /**
    * Add another record to this one as library data.
    *
    * Useful for libraries that can run in standalone mode; the host
    * simply calls this method on the record the library produces.
    * Merges file lists.
    */
    void addRecordAsLibraryData(Record const &childRecord, std::string const &name);

private:
    internal::IDField id;
    std::string type;
    FileSet files;
};


/**
 * A RecordLoader is used to convert conduit::Node instances which represent
 * Sina Records into instances of their corresponding axom::sina::Record
 * subclasses. For convenience, a RecordLoader capable of handling Records of all known
 * types can be created using createRecordLoaderWithAllKnownTypes:
 *
 * \code
 * axom::sina::Document myDocument = axom::sina::Document(jObj, axom::sina::createRecordLoaderWithAllKnownTypes());
 * \endcode
 */
class RecordLoader {
public:
    /**
     * A TypeLoader is a function which converts records of a specific type
     * to their corresponding sub classes.
     */
    using TypeLoader = std::function<std::unique_ptr<Record>(
            conduit::Node const &)>;

    /**
     * Add a function for loading records of the specified type.
     *
     * @param type the type of records this function can load
     * @param loader the function which can load the records
     */
    void addTypeLoader(std::string const &type, TypeLoader loader);

    /**
     * Load a axom::sina::Record from its conduit Node representation.
     *
     * @param recordAsNode the Record as a Node
     * @return the Record
     */
    std::unique_ptr<Record> load(conduit::Node const &recordAsNode) const;

    /**
     * Check whether this loader can load records of the given type.
     *
     * @param type the type of the records to check
     * @return whether records of the given type can be loaded
     */
    bool canLoad(std::string const &type) const;

private:
    std::unordered_map<std::string, TypeLoader> typeLoaders;
};

/**
 * Create a RecordLoader which can load records of all known types.
 *
 * @return the newly-created loader
 */
RecordLoader createRecordLoaderWithAllKnownTypes();

}  // end sina namespace
}  // end axom namespace

#endif //SINA_RECORD_HPP
