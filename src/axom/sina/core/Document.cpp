// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

/*!
 ******************************************************************************
 *
 * \file Document.cpp
 *
 * \brief   Implementation file for Sina Document class
 *
 ******************************************************************************
 */

#include "axom/sina/core/Document.hpp"
#include "axom/sina/core/CurveSet.hpp"
#include "axom/sina/core/Curve.hpp"
#include "axom/sina/core/Record.hpp"

#include <cstdio>
#include <fstream>
#include <ios>
#include <iostream>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <nlohmann/json.hpp>
#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp"
#include <algorithm>
#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp"


using json = nlohmann::json;

namespace axom
{
namespace sina
{

namespace
{
char const RECORDS_KEY[] = "records";
char const RELATIONSHIPS_KEY[] = "relationships";
char const SAVE_TMP_FILE_EXTENSION[] = ".sina.tmp";
}

std::string protocol_set(std::string protocol, std::string const &name) {
  size_t pos = name.rfind('.');
  std::string toReturn = name;

  if (pos != std::string::npos) {
        std::string found = name.substr(pos+1);

        if (("." + found) != protocol) {
          for (const std::string& file_type : supported_types) {
            std::string lower_case_type = file_type; 
            std::transform(lower_case_type.begin(), lower_case_type.end(), lower_case_type.begin(), [](unsigned char c) { return std::tolower(c); });

            if ((lower_case_type != protocol) && (lower_case_type == found)) {
              std::cout << "|| WARNING: INCORRECT FILE EXTENSION FOUND (FOUND: ." << found << "; EXPECTED: " << protocol << "),  SINA WILL BE REPLACING IT TO THE " << protocol << " FILE EXTENSION ||";
              toReturn = name.substr(0, pos);
              toReturn = toReturn.append(protocol);
              return toReturn;
            }
          }
          std::cout << "|| WARNING: BROKEN FILE EXTENSION FOUND, SINA WILL BE APPENDING THE " << protocol << " FILE EXTENSION ||";
          toReturn = toReturn.append(protocol);
          return toReturn;
          
        } else {
          return name;
        }
  } else {
    std::cout << "|| WARNING: NO FILE EXTENSION FOUND, SINA WILL BE ADDING THE " << protocol << " FILE EXTENSION ||";
    toReturn = toReturn.append(protocol);
    return toReturn;
  }
}

void removeSlashes(const conduit::Node& originalNode, conduit::Node& modifiedNode)
{
    for (auto it = originalNode.children(); it.has_next();)
    {
        it.next();
        std::string key = it.name();
        std::string modifiedKey = key;
        //modifiedKey.erase(std::remove(modifiedKey.begin(), modifiedKey.end(), '/'), modifiedKey.end());

        std::string toReplace = "/";
        std::string replacement = "__SINA_SLASHREPLACE__";

        size_t pos = 0;
        // Find and replace all occurrences of "/"
        while ((pos = modifiedKey.find(toReplace, pos)) != std::string::npos) {
          modifiedKey.replace(pos, toReplace.length(), replacement);
          pos += replacement.length(); // Move past the replaced substring
        }

        modifiedNode[modifiedKey] = it.node();
        
        if (it.node().dtype().is_object())
        {
            conduit::Node nestedNode;
            removeSlashes(it.node(), nestedNode);
            modifiedNode[modifiedKey].set(nestedNode); 
        }
    }
}

void restoreSlashes(const conduit::Node& modifiedNode, conduit::Node& restoredNode)
{
    // Check if List or Object, if its a list the else statement would turn it into an object
    // which breaks the Document

    if (modifiedNode.dtype().is_list())
    {
        // If its empty with no children it's the end of a tree

        for (auto it = modifiedNode.children(); it.has_next();)
        {
            it.next();
            conduit::Node& newChild = restoredNode.append();

            // Leaves empty nodes empty, if null data is set the
            // Document breaks

            if (it.node().dtype().is_string() || it.node().dtype().is_number())
            {
                newChild.set(it.node()); // Lists need .set
            }
            
            // Recursive Call
            if (it.node().number_of_children() > 0)
            {
                restoreSlashes(it.node(), newChild);
            }
        }
    }
    else
    {
        for (auto it = modifiedNode.children(); it.has_next();)
        {
            it.next();
            std::string key = it.name();
            std::string restoredKey = key;
            std::string toReplace = "__SINA_SLASHREPLACE__";
            std::string replacement = "/";

            size_t pos = 0;
            // Find and replace all occurrences of "__SLASH__"

            while ((pos = restoredKey.find(toReplace, pos)) != std::string::npos) {
                restoredKey.replace(pos, toReplace.length(), replacement);
                pos += replacement.length();
            }


            // Initialize a new node for the restored key
            conduit::Node& newChild = restoredNode.add_child(restoredKey);

            // Leaves empty keys empty but continues recursive call if its a list
            if (it.node().dtype().is_string() || it.node().dtype().is_number() || it.node().dtype().is_object())
            {
                newChild.set(it.node());
            }
            else if (it.node().dtype().is_list())
            {
                restoreSlashes(it.node(), newChild);  // Handle nested lists
            }

            // If the node has children, recursively restore them
            if (it.node().number_of_children() > 0)
            {
                conduit::Node nestedNode;
                restoreSlashes(it.node(), nestedNode);
                newChild.set(nestedNode);
            }
        }
    }
}

void Document::add(std::unique_ptr<Record> record)
{
  records.emplace_back(std::move(record));
}

void Document::add(Relationship relationship)
{
  relationships.emplace_back(std::move(relationship));
}

conduit::Node Document::toNode() const
{
  conduit::Node document(conduit::DataType::object());
  document[RECORDS_KEY] = conduit::Node(conduit::DataType::list());
  document[RELATIONSHIPS_KEY] = conduit::Node(conduit::DataType::list());
  for(auto &record : records)
  {
    auto &list_entry = document[RECORDS_KEY].append();
    list_entry.set_node(record->toNode());
  }
  for(auto &relationship : relationships)
  {
    auto &list_entry = document[RELATIONSHIPS_KEY].append();
    list_entry = relationship.toNode();
  }
  return document;
}

void Document::createFromNode(conduit::Node const &asNode,
                              RecordLoader const &recordLoader)
{
  if(asNode.has_child(RECORDS_KEY))
  {
    conduit::Node record_nodes = asNode[RECORDS_KEY];
    if(record_nodes.dtype().is_list())
    {
      auto recordIter = record_nodes.children();
      while(recordIter.has_next())
      {
        auto record = recordIter.next();
        add(recordLoader.load(record));
      }
    }
    else
    {
      std::ostringstream message;
      message << "The '" << RECORDS_KEY
              << "' element of a document must be an array";
      throw std::invalid_argument(message.str());
    }
  }

  if (asNode.has_child(RELATIONSHIPS_KEY))
    {
        conduit::Node relationship_nodes = asNode[RELATIONSHIPS_KEY];
        if (relationship_nodes.number_of_children() == 0)
        {
            relationship_nodes.set(conduit::DataType::list());
        }
        else if (!relationship_nodes.dtype().is_list())
        {
            std::ostringstream message;
            message << "The '" << RELATIONSHIPS_KEY << "' element of a document must be an array";
            throw std::invalid_argument(message.str());
        }

        auto relationshipsIter = relationship_nodes.children();
        while (relationshipsIter.has_next())
        {
            auto &relationship = relationshipsIter.next();
            add(Relationship{relationship});
        }
    }
}

Document::Document(conduit::Node const &asNode, RecordLoader const &recordLoader)
{
  this->createFromNode(asNode, recordLoader);
}

Document::Document(std::string const &asJson, RecordLoader const &recordLoader)
{
  conduit::Node asNode;
  asNode.parse(asJson, "json");
  this->createFromNode(asNode, recordLoader);
}

void Document::toHDF5(const std::string &filename) const 
{
    conduit::Node node;
    conduit::Node &recordsNode = node["records"];
    conduit::Node &relationshipsNode = node["relationships"];

    for (const auto& record : getRecords())
    {
        conduit::Node recordNode = record->toNode();
        conduit::Node modifiedRecordNode;

        removeSlashes(recordNode, modifiedRecordNode);

        recordsNode.append() = modifiedRecordNode; 
    }

    // Process relationships
    for (const auto& relationship : getRelationships())
    {
        conduit::Node relationshipNode = relationship.toNode();
        //TODO:: Remove Slashes
        relationshipsNode.append() = relationshipNode;
    }

    conduit::relay::io::save(node, filename, "hdf5");
}

//

std::string Document::toJson(conduit::index_t indent,
                             conduit::index_t depth,
                             const std::string &pad,
                             const std::string &eoe) const
{
  return this->toNode().to_json("json", indent, depth, pad, eoe);
}

void saveDocument(Document const &document, std::string const &fileName, Protocol protocol)
{
  std::string tmpFileName = fileName + SAVE_TMP_FILE_EXTENSION;
  std::string name;

  try
  {
      if (protocol == Protocol::JSON)
      {
          name = protocol_set(".json", fileName);
          auto asJson = document.toJson();
          std::ofstream fout {tmpFileName};
          fout.exceptions(std::ostream::failbit | std::ostream::badbit);
          fout << asJson;
          fout.close();
      }
      else if (protocol == Protocol::HDF5)
      {
          name = protocol_set(".hdf5", fileName);
          document.toHDF5(tmpFileName);
      }
      else
      {
          throw std::invalid_argument("Invalid format choice. Please enter 'json' or 'hdf5'.");
      }

      if (rename(tmpFileName.c_str(), name.c_str()) != 0)
      {
          std::string message {"Could not save to '"};
          message += fileName;
          message += "'";
          throw std::ios::failure {message};
      }
    }
    catch (const std::exception &e)
    {
        std::cerr << "An error occurred: " << e.what() << "\n";
        std::remove(tmpFileName.c_str());
        throw;
    }

  // It is a common use case for users to want to overwrite their files as
  // the simulation progresses. However, this operation should be atomic so
  // that if a write fails, the old file is left intact. For this reason,
  // we write to a temporary file first and then move the file. The temporary
  // file is in the same directory to ensure that it is part of the same
  // file system as the destination file so that the move operation is
  // atomic.
}

Document loadDocument(std::string const &path, Protocol protocol)
{
  return loadDocument(path, createRecordLoaderWithAllKnownTypes(), protocol);
}

Document loadDocument(std::string const &path, RecordLoader const &recordLoader, Protocol protocol)
{
    conduit::Node node;
    
    // Load the file depending on the protocol
    if (protocol == Protocol::JSON) {
        std::ifstream file_in {path};
        std::ostringstream file_contents;
        file_contents << file_in.rdbuf();
        file_in.close();
        node.parse(file_contents.str(), "json");
        return Document {node, recordLoader};
    } else if (protocol == Protocol::HDF5) {
        conduit::Node modifiedNode;
        conduit::relay::io::load(path, "hdf5", node);
        restoreSlashes(node, modifiedNode);
        return Document {modifiedNode, recordLoader};
    }
}

}  // namespace sina
}  // namespace axom
