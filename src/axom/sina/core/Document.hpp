// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef SINA_DOCUMENT_HPP
#define SINA_DOCUMENT_HPP

/*!
 ******************************************************************************
 *
 * \file Document.hpp
 *
 * \brief   Header file for Sina Document class
 *
 ******************************************************************************
 */

#include <memory>
#include <vector>

#include "conduit.hpp"

#include "axom/sina/core/Record.hpp"
#include "axom/sina/core/Relationship.hpp"

#define SINA_FILE_FORMAT_VERSION_MAJOR 1
#define SINA_FILE_FORMAT_VERSION_MINOR 0

namespace axom
{
namespace sina
{

/**
 * \brief An object representing the top-level object of a Sina JSON file
 *
 * A Document represents the top-level object of a JSON file conforming to the
 * Sina schema. When serialized, these documents can be ingested into a
 * Sina database and used with the Sina tool.
 *
 * Documents contain at most two objects: a list of Records and a list of Relationships. A simple, empty document:
 * \code{.json}
 *   {
 *      "records": [],
 *      "relationships": []
 *   }
 * \endcode
 *
 * The "records" list can contain Record objects and their inheriting types, such as Run (for a full list, please see
 * the inheritance diagram in the Record documentation). The "relationships" list can contain Relationship objects.
 *
 * Documents can be assembled programatically and/or generated from existing JSON. An example of an assembled
 * Document is provided on the main page. To load a Document from an existing JSON file:
 * \code
 *   axom::sina::Document myDocument = axom::sina::loadDocument("path/to/infile.json");
 * \endcode
 *
 * To generate a Document from a JSON string and vice versa:
 * \code
 *   std::string my_json = "{\"records\":[{\"type\":\"run\",\"id\":\"test\"}],\"relationships\":[]}";
 *   axom::sina::Document myDocument = axom::sina::Document(my_json, axom::sina::createRecordLoaderWithAllKnownTypes());
 *   std::cout << myDocument.toJson() << std::endl;);
 * \endcode
 *
 * You can add further entries to the Document using add():
 * \code
 *   std::unique_ptr<sina::Record> myRun{new axom::sina::Run{someID, "My Sim Code", "1.2.3", "jdoe"}};
 *   axom::sina::Relationship myRelationship{someID, "comes before", someOtherID};
 *   myDocument.add(myRun);
 *   myDocument.add(myRelationship);
 * \endcode
 *
 * You can also export your Document to file:
 * \code
 *   axom::sina::saveDocument(myDocument, "path/to/outfile.json")
 * \endcode
 *
 * Check the Sina file format version with:
 * \code
 *   axom::sina::getSinaFileFormatVersion();
 * \endcode
 *
 */
class Document
{
public:
  /**
     * A vector of pointers to Record objects.
     */
  using RecordList = std::vector<std::unique_ptr<Record>>;

  /**
     * A vector of Relationship objects.
     */
  using RelationshipList = std::vector<Relationship>;

  /**
     * Construct an empty Document.
     */
  Document() = default;

  /**
     * Disable copying Document objects. We must do this since we hold
     * pointers to polymorphic objects.
     */
  Document(Document const &) = delete;

  /**
     * Disabling copy assignment.
     */
  Document &operator=(Document const &) = delete;

  /**
     * Move constructor which should be handled by the compiler.
     */
  Document(Document &&) = default;

  /**
     * Move assignment which should be handled by the compiler.
     */
  Document &operator=(Document &&) = default;

  /**
     * \brief Create a Document from its Conduit Node representation
     *
     * \param asNode the Document as a Node
     * \param recordLoader an RecordLoader to use to load the different
     *                     types of records which may be in the document
     */
  Document(conduit::Node const &asNode, RecordLoader const &recordLoader);

  /**
     * \brief Create a Document from a JSON string representation
     *
     * \param asJson the Document as a JSON string
     * \param recordLoader an RecordLoader to use to load the different
     *                     types of records which may be in the document
     */
  Document(std::string const &asJson, RecordLoader const &recordLoader);

  /**
     * \brief Add the given record to this document.
     *
     * \param record the record to add
     */
  void add(std::unique_ptr<Record> record);

  /**
     * \brief Get the list of records currently in this document.
     *
     * \return the list of records
     */
  RecordList const &getRecords() const noexcept { return records; }

  /**
     * \brief Add a relationship to this document
     *
     * \param relationship the relationship to add
     */
  void add(Relationship relationship);

  /**
     * \brief Get the list of relationships in this document.
     *
     * \return the list of relationships
     */
  RelationshipList const &getRelationships() const noexcept
  {
    return relationships;
  }

  /**
     * \brief Convert this document to a conduit Node.
     *
     * \return the contents of the document as a Node
     */
  conduit::Node toNode() const;

  /**
     * \brief Convert this document to a JSON string.
     *
     * \return the contents of the document as a JSON string
     */
  std::string toJson(conduit::index_t indent = 0,
                     conduit::index_t depth = 0,
                     const std::string &pad = "",
                     const std::string &eoe = "") const;

private:
  /**
     * Constructor helper method, extracts info from a conduit Node.
     */
  void createFromNode(conduit::Node const &asNode,
                      RecordLoader const &recordLoader);
  RecordList records;
  RelationshipList relationships;
};

/**
 * \brief Save the given Document to the specified location. If the given file exists,
 *        it will be overwritten.
 *
 * \param document the Document to save
 * \param fileName the location to which to save the file
 * \throws std::ios::failure if there are any IO errors
 */
void saveDocument(Document const &document, std::string const &fileName);

/**
 * \brief Get the current file format version.
 *
 * \return A string representing the file format version.
 */
inline std::string getSinaFileFormatVersion()
{
  return std::to_string(SINA_FILE_FORMAT_VERSION_MAJOR) + "." +
    std::to_string(SINA_FILE_FORMAT_VERSION_MINOR);
}

/**
 * \brief Load a document from the given path. Only records which this library
 *        knows about will be able to be loaded.
 *
 * \param path the file system path from which to load the document
 * \return the loaded Document
 */
Document loadDocument(std::string const &path);

/**
 * \brief Load a document from the given path.
 *
 * \param path the file system path from which to load the document
 * \param recordLoader the RecordLoader to use to load the different types
 *                     of records
 * \return the loaded Document
 */
Document loadDocument(std::string const &path, RecordLoader const &recordLoader);

}  // namespace sina
}  // namespace axom

#endif  //SINA_DOCUMENT_HPP
