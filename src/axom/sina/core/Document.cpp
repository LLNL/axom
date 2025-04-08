// Copyright (c) 2017-2025, Lawrence Livermore National Security, LLC and
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
#include "axom/json.hpp"
#include "axom/config.hpp"
#include "axom/core/Path.hpp"

#include "axom/config.hpp"
#include "axom/core/Path.hpp"
#include "axom/core/utilities/StringUtilities.hpp"

#include "conduit.hpp"
#ifdef AXOM_USE_HDF5
  #include "conduit_relay.hpp"
  #include "conduit_relay_io/utilities/StringUtilities.hpp"
#endif

#include "conduit.hpp"
#ifdef AXOM_USE_HDF5
  #include "conduit_relay.hpp"
  #include "conduit_relay_io.hpp"
#endif

#include <functional>
#include <set>
#include <string>
#include <cstdio>
#include <fstream>
#include <ios>
#include <iostream>
#include <utility>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp"
#include "conduit_relay_io_hdf5.hpp"


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
}  // namespace

void protocolWarn(std::string const protocol, std::string const &name)
{
  std::unordered_map<std::string, std::string> protocolMessages = {
    {".json", ".json extension not found, did you mean to save to this format?"},
    {".hdf5",
     ".hdf5 extension not found, did you use one of its other supported types? "
     "(h5, hdf, ...)"}};

  Path path(name, '.');

  if(protocol != '.' + path.baseName())
  {
    auto messageIt = protocolMessages.find(protocol);
    if(messageIt != protocolMessages.end())
    {
      std::cerr << messageIt->second;
    }
  }
}

std::string get_supported_file_types()
{
  std::string types = "[";
  for(size_t i = 0; i < supported_types.size(); ++i)
  {
    types += supported_types[i];
    if(i < supported_types.size() - 1)
    {
      types += ", ";
    }
  }
  types += "]";
  return types;
}

void Document::add(std::unique_ptr<Record> record) { records.emplace_back(std::move(record)); }

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

void Document::createFromNode(const conduit::Node &asNode, const RecordLoader &recordLoader)
{
  conduit::Node nodeCopy = asNode;

  auto processChildNodes = [&](const char *key, std::function<void(conduit::Node &)> addFunc) {
    if(nodeCopy.has_child(key))
    {
      conduit::Node &childNodes = nodeCopy[key];

      // -- 1. Check if this node is a primitive leaf (throw immediately if so)
      // Customize these checks to match exactly what you consider "primitive."
      if(childNodes.dtype().is_number() || childNodes.dtype().is_char8_str() ||
         childNodes.dtype().is_string())
      {
        std::ostringstream message;
        message << "The '" << key << "' element of a document cannot be a primitive value.";
        throw std::invalid_argument(message.str());
      }

      // -- 2. Not a primitive. Check if it has no children.
      if(childNodes.number_of_children() == 0)
      {
        // Turn it into an empty list
        childNodes.set(conduit::DataType::list());
      }

      // -- 3. If it's still not a list, throw
      if(!childNodes.dtype().is_list())
      {
        std::ostringstream message;
        message << "The '" << key << "' element of a document must be an array/list.";
        throw std::invalid_argument(message.str());
      }

      // -- 4. Now it's guaranteed to be a list, so iterate
      auto childIter = childNodes.children();
      while(childIter.has_next())
      {
        conduit::Node child = childIter.next();
        addFunc(child);
      }
    }
  };
  processChildNodes(RECORDS_KEY, [&](conduit::Node &record) { add(recordLoader.load(record)); });

  processChildNodes(RELATIONSHIPS_KEY,
                    [&](conduit::Node &relationship) { add(Relationship {relationship}); });
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

#ifdef AXOM_USE_HDF5
void removeSlashes(const conduit::Node &originalNode, conduit::Node &modifiedNode)
{
  for(auto it = originalNode.children(); it.has_next();)
  {
    it.next();
    std::string key = it.name();
    std::string modifiedKey = axom::utilities::string::replaceAllInstances(key, "/", slashSubstitute);

    modifiedNode[modifiedKey] = it.node();

    if(it.node().dtype().is_object())
    {
      conduit::Node nestedNode;
      removeSlashes(it.node(), nestedNode);
      modifiedNode[modifiedKey].set(nestedNode);
    }
  }
}

void restoreSlashes(const conduit::Node &modifiedNode, conduit::Node &restoredNode)
{
  // Check if List or Object, if its a list the else statement would turn it into an object
  // which breaks the Document

  if(modifiedNode.dtype().is_list())
  {
    // If its empty with no children it's the end of a tree

    for(auto it = modifiedNode.children(); it.has_next();)
    {
      it.next();
      conduit::Node &newChild = restoredNode.append();
      auto data_type = it.node().dtype();

      // Leaves empty nodes empty, if null data is set the
      // Document breaks

      if(data_type.is_string() || data_type.is_number())
      {
        newChild.set(it.node());  // Lists need .set
      }

      // Recursive Call
      if(it.node().number_of_children() > 0)
      {
        restoreSlashes(it.node(), newChild);
      }
    }
  }
  else
  {
    for(auto it = modifiedNode.children(); it.has_next();)
    {
      it.next();
      std::string key = it.name();
      std::string restoredKey =
        axom::utilities::string::replaceAllInstances(key, slashSubstitute, "/");

      // Initialize a new node for the restored key
      conduit::Node &newChild = restoredNode.add_child(restoredKey);
      auto data_type = it.node().dtype();

      // Leaves empty keys empty but continues recursive call if its a list
      if(data_type.is_string() || data_type.is_number() || data_type.is_object())
      {
        newChild.set(it.node());
      }
      else if(data_type.is_list())
      {
        restoreSlashes(it.node(), newChild);  // Handle nested lists
      }

      // If the node has children, recursively restore them
      if(it.node().number_of_children() > 0)
      {
        conduit::Node nestedNode;
        restoreSlashes(it.node(), nestedNode);
        newChild.set(nestedNode);
      }
    }
  }
}

void Document::toHDF5(const std::string &filename) const
{
  conduit::Node node;
  conduit::Node &recordsNode = node["records"];
  conduit::Node &relationshipsNode = node["relationships"];

  for(const auto &record : getRecords())
  {
    conduit::Node recordNode = record->toNode();

    removeSlashes(recordNode, recordsNode.append());
  }

  // Process relationships
  for(const auto &relationship : getRelationships())
  {
    conduit::Node relationshipNode = relationship.toNode();

    removeSlashes(relationshipNode, relationshipsNode.append());
  }

  conduit::relay::io::save(node, filename, "hdf5");
}
#endif

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
  // It is a common use case for users to want to overwrite their files as
  // the simulation progresses. However, this operation should be atomic so
  // that if a write fails, the old file is left intact. For this reason,
  // we write to a temporary file first and then move the file. The temporary
  // file is in the same directory to ensure that it is part of the same
  // file system as the destination file so that the move operation is
  // atomic.

  std::string tmpFileName = fileName + SAVE_TMP_FILE_EXTENSION;

  switch(protocol)
  {
  case Protocol::JSON:
  {
    protocolWarn(".json", fileName);
    auto asJson = document.toJson();
    std::ofstream fout {tmpFileName};
    fout.exceptions(std::ostream::failbit | std::ostream::badbit);
    fout << asJson;
    fout.close();
  }
  break;
#ifdef AXOM_USE_HDF5
  case Protocol::HDF5:
    protocolWarn(".hdf5", fileName);
    document.toHDF5(tmpFileName);
    break;
#endif
  default:
  {
    std::ostringstream message;
    message << "Invalid format choice. Please choose from one of the supported "
               "protocols: "
            << get_supported_file_types();
    throw std::invalid_argument(message.str());
  }
  }

  if(rename(tmpFileName.c_str(), fileName.c_str()) != 0)
  {
    std::string message {"Could not save to '"};
    message += fileName;
    message += "'";
    throw std::ios::failure {message};
  }
}

Document loadDocument(std::string const &path, Protocol protocol)
{
  return loadDocument(path, createRecordLoaderWithAllKnownTypes(), protocol);
}

Document loadDocument(std::string const &path, RecordLoader const &recordLoader, Protocol protocol)
{
  conduit::Node node, modifiedNode;
  std::ostringstream file_contents;
  std::ifstream file_in {path};

  // Load the file depending on the protocol
  switch(protocol)
  {
  case Protocol::JSON:
    file_contents << file_in.rdbuf();
    file_in.close();
    node.parse(file_contents.str(), "json");
    return Document {node, recordLoader};
#ifdef AXOM_USE_HDF5
  case Protocol::HDF5:
    file_in.close();
    conduit::relay::io::load(path, "hdf5", node);
    restoreSlashes(node, modifiedNode);
    return Document {modifiedNode, recordLoader};
#endif
  default:
    std::ostringstream message;
    message << "Invalid format choice. Please choose from one of the supported "
               "protocols: "
            << get_supported_file_types();
    throw std::invalid_argument(message.str());
    break;
  }
}

// Unified helper function that validates a set of curves (either dependent or independent).
// Parameters:
//   new_curves     : a map of new curves (e.g. std::map<std::string, Curve>).
//   existing_keys  : a set of keys that are present in the existing data.
//   getExistingSize: a callable that takes a curve key and returns the size (int)
//                    from the existing data or -1 if not found.
//   curveType      : "dependent" or "independent" (for erro messages).
//   recordId       : identifier for the record (for error messages).
//   curveSetId     : identifier for the curve set (for error messages).
//   baseline       : reference to an int that will hold the computed baseline value
//                    (initialize to -1 to have the function set baseline).
//
// Templated on the container type to avoid conversion issues (e.g., unordered_map vs. map).
template<typename CurveMap>
bool validate_curves_unified(
    const CurveMap &new_curves,
    const std::set<std::string>& existing_keys,
    std::function<int(const std::string&)> getExistingSize,
    const std::string &curveType,
    const std::string &recordId,
    const std::string &curveSetId,
    int baseline)
{
    std::set<std::string> unionKeys = existing_keys;
    for (const auto &pair : new_curves) {
        unionKeys.insert(pair.first);
    }

    for (const auto &key : unionKeys) {
        int newSize = 0;
        auto newItr = new_curves.find(key);
        if (newItr != new_curves.end()) {
            newSize = static_cast<int>(newItr->second.getValues().size());
        }
        int existingSize = getExistingSize(key);

        // Get total size but ignore -1 returns
        int total = (existingSize >= 0 ? newSize + existingSize : newSize);

        if (baseline < 0) {
            baseline = total;
        } else if (baseline != total) {
            std::cerr << "Error validating " << curveType << ": Record " << recordId
                      << ", Curve Set " << curveSetId << ", Curve " << key
                      << " size mismatch (expected " << baseline << ", got " << total << ")."
                      << std::endl;
            return false;
        }
    }
    return true;
}

bool validate_curve_sets_json(const DataHolder::CurveSetMap new_curve_sets,
                              const nlohmann::json &existing_curve_sets,
                              const std::string recordId)
{
    for (auto & [curveSetId, existing_curve_set] : existing_curve_sets.items()) {
        // Create an alias to avoid capturing the structured binding directly.
        auto &ecs = existing_curve_set;

        if (new_curve_sets.find(curveSetId) != new_curve_sets.end()) {
            const auto &new_curve_set = new_curve_sets.at(curveSetId);

            // Validate dependent curves.
            std::set<std::string> existingDepKeys;
            if (ecs.contains("dependent")) {
                for (auto & [key, val] : ecs["dependent"].items()) {
                    existingDepKeys.insert(key);
                }
            }
            
            auto getExistingDepSize = [&ecs](const std::string &key) -> int {
                if (ecs.contains("dependent") &&
                    ecs["dependent"].contains(key)) {
                    return static_cast<int>(ecs["dependent"][key]["value"].size());
                }
                return -1;
            };

            if (!new_curve_set.getDependentCurves().empty() &&
                !validate_curves_unified(new_curve_set.getDependentCurves(),
                                         existingDepKeys,
                                         getExistingDepSize,
                                         "dependent",
                                         recordId,
                                         curveSetId,
                                         -1))
            {
                return false;
            }

            // Validate independent curves.
            std::set<std::string> existingIndepKeys;
            if (ecs.contains("independent")) {
                for (auto & [key, val] : ecs["independent"].items()) {
                    existingIndepKeys.insert(key);
                }
            }

            auto getExistingIndepSize = [&ecs](const std::string &key) -> int {
                if (ecs.contains("independent") &&
                    ecs["independent"].contains(key)) {
                    return static_cast<int>(ecs["independent"][key]["value"].size());
                }
                return -1;
            };

            if (!new_curve_set.getIndependentCurves().empty() &&
                !validate_curves_unified(new_curve_set.getIndependentCurves(),
                                         existingIndepKeys,
                                         getExistingIndepSize,
                                         "independent",
                                         recordId,
                                         curveSetId,
                                         -1))
            {
                return false;
            }
        } else {
            std::cerr << "Curve set " << curveSetId
                      << " not found in new data for record " << recordId << std::endl;
            return false;
        }
    }
    return true;
}

bool append_to_json(const std::string& jsonFilePath,
                 Document const &newData,
                 const int data_protocol,
                 const int udc_protocol)
{
 std::ifstream file(jsonFilePath);
 if (!file.is_open()) {
     std::cerr << "Error opening file: " << jsonFilePath << std::endl;
     return false;
 }

 nlohmann::json j;
 try {
     file >> j;
 } catch (const nlohmann::json::parse_error& e) {
     std::cerr << "JSON parsing error: " << e.what() << std::endl;
     return false;
 }
 file.close();

 if (j["records"].size() != newData.getRecords().size()) {
     std::cerr << "Mismatch in the number of records." << std::endl;
     return false;
 }

 std::vector<std::function<void()>> write_queue;

 for (auto& new_record : newData.getRecords()) {
     for (auto& existing_record : j["records"]) {
         if (!validate_curve_sets_json(new_record->getCurveSets(), existing_record["curve_sets"], existing_record["id"])) {
             return false;
         }
     }
 }

 for (auto& new_record : newData.getRecords()) {
     bool found = false;
     for (auto& existing_record : j["records"]) {
         if (new_record->getId().getId() == existing_record["id"]) {
             found = true;
             // ----------------------
             // Queue update of CURVE SETS.
             // ----------------------
             if ((new_record->getCurveSets().size() > 0) &&
                 existing_record.contains("curve_sets"))
             {
                 nlohmann::json& existing_curve_sets = existing_record["curve_sets"];
                 auto& new_curve_sets = new_record->getCurveSets();
                 for (auto& new_curve_set : new_curve_sets) {
                     auto& [new_cs_key, new_cs_values] = new_curve_set;
                     auto cs_key = new_cs_key;
                     
                     for (auto& new_dependent : new_cs_values.getDependentCurves()) {
                         auto& [new_dep_key, new_dep_values] = new_dependent;
                         auto dep_key = new_dep_key; // local copy for lambda capture
                         for (const auto& value : new_dep_values.getValues()) {
                             write_queue.push_back(
                                 [ &existing_curve_sets, cs_key, dep_key, value ]() {
                                     existing_curve_sets[cs_key]["dependent"][dep_key]["value"]
                                         .push_back(static_cast<double>(value));
                                 }
                             );
                         }
                     }
                     
                     // For independent curves.
                     for (auto& new_independent : new_cs_values.getIndependentCurves()) {
                         auto& [new_indep_key, new_indep_values] = new_independent;
                         auto indep_key = new_indep_key; // local copy for lambda capture
                         for (const auto& value : new_indep_values.getValues()) {
                             write_queue.push_back(
                                 [ &existing_curve_sets, cs_key, indep_key, value ]() {
                                     existing_curve_sets[cs_key]["independent"][indep_key]["value"]
                                         .push_back(static_cast<double>(value));
                                 }
                             );
                         }
                     }
                 }
             }
             
             // ----------------------
             // Queue update of DATA VALUES.
             // ----------------------
             if ((new_record->getData().size() > 0) &&
                 existing_record.contains("data"))
             {
                 nlohmann::json& existing_data_sets = existing_record["data"];
                 auto& new_data_sets = new_record->getData();
                 for (auto& new_data : new_data_sets) {
                     auto& [new_data_key, new_data_pair] = new_data;
                     auto data_key = new_data_key; // local copy for capture
                     nlohmann::json obj = nlohmann::json::parse(new_data_pair.toNode().to_json());
                     if (existing_data_sets.contains(data_key)) {
                         // Duplicate Handling
                         switch(data_protocol) {
                             case 1:
                                 write_queue.push_back(
                                     [ &existing_data_sets, data_key, obj ]() {
                                         existing_data_sets[data_key] = obj;
                                     }
                                 );
                                 break;
                             case 2:
                                 break;
                             case 3:
                                 std::cerr << "Found a duplicate data entry, protocol 3 dictates append cancellation." << std::endl;
                                 return false;
                             default:
                                 std::cerr << "Invalid Data Protocol Entry. Append cancelled." << std::endl;
                                 return false;
                         }
                     } else {
                         write_queue.push_back(
                             [ &existing_data_sets, data_key, obj ]() {
                                 existing_data_sets[data_key] = obj;
                             }
                         );
                     }
                 }
             }
             
             // ----------------------
             // Queue update of USER DEFINED CONTENT.
             // ----------------------
             if ((!new_record->getUserDefinedContent().dtype().is_empty()) &&
                 existing_record.contains("user_defined"))
             {
                 nlohmann::json& existing_udc = existing_record["user_defined"];
                 auto& new_udc = new_record->getUserDefinedContent();
                 for (auto& udc : new_udc.children()) {
                     std::string udc_name = udc.name();
                     if (existing_record["user_defined"].contains(udc_name)) {
                         switch(udc_protocol) {
                             case 1:
                                 write_queue.push_back(
                                     [ &existing_udc, udc_name, udc ]() {
                                         existing_udc[udc_name] = udc.as_string();
                                     }
                                 );
                                 break;
                             case 2:
                                 break;
                             case 3:
                                 std::cerr << "Found duplicate UDC, protocol 3 dictates append cancellation." << std::endl;
                                 return false;
                             default:
                                 std::cerr << "Invalid UDC Protocol Entry. Append cancelled." << std::endl;
                                 return false;
                         }
                     } else {
                         write_queue.push_back(
                             [ &existing_udc, udc_name, udc ]() {
                                 existing_udc[udc_name] = udc.as_string();
                             }
                         );
                     }
                 }
             }
         }
     }
     if (!found) {
         nlohmann::json obj = nlohmann::json::parse(new_record->toNode().to_json());
         write_queue.push_back(
             [&j, obj]() {
                 j["records"][j["records"].size()] = obj;
             }
         );
     }
 }
 
 // Execute all queued operations.
 for (auto &op : write_queue) {
     try {
         op();
     } catch (const std::exception &e) {
         std::cerr << "Error executing queued operation: " << e.what() << std::endl;
         return false;
     }
 }
 
 std::ofstream outFile(jsonFilePath);
 if (!outFile.is_open()) {
     std::cerr << "Error saving updated JSON to file!" << std::endl;
     return false;
 }
 outFile << j.dump(4);
 outFile.close();

 return true;
}

#ifdef AXOM_USE_HDF5
// Top-level HDF5 validation function.
bool validate_curve_sets_hdf5(const Document &newData,
                              conduit::relay::io::IOHandle &existing_file,
                              const conduit::Node &info)
{
    std::vector<std::string> record_list;
    existing_file.list_child_names("records", record_list);

    for (auto &record : newData.getRecords()) {
        for (const auto &rec_name : record_list) {
            std::string id_path = "records/" + rec_name + "/id/";
            conduit::Node id;
            existing_file.read(id_path, id);
            if (id.to_string() == "\"" + record->getId().getId() + "\"") {
                try {
                    std::string cs_path = "records/" + rec_name + "/curve_sets/";
                    std::vector<std::string> curve_sets_list;
                    existing_file.list_child_names(cs_path, curve_sets_list);
                    auto &new_curve_sets = record->getCurveSets();

                    for (const auto &cs_name : curve_sets_list) {
                        auto cs_itr = new_curve_sets.find(cs_name);
                        if (cs_itr == new_curve_sets.end())
                            continue;

                        const auto &new_curve_set = cs_itr->second;
                        std::set<std::string> existingDepKeys;
                        if (info.has_path(cs_path + cs_name + "/dependent/")) {
                          for (const auto &key : info[cs_path + cs_name + "/dependent/"].child_names()) {
                            existingDepKeys.insert(key);
                          }
                        }

                        // Validate dependent curves.
                        auto getExistingDepSize = [&info, cs_path, cs_name](const std::string &key) -> int {
                            std::string fullPath = cs_path + cs_name + "/dependent/" + key + "/value";
                            if (info.has_path(fullPath)) {
                                return info[fullPath]["num_elements"].to_int();
                            }
                            return -1;
                        };

                        if (!new_curve_set.getDependentCurves().empty() &&
                            !validate_curves_unified(new_curve_set.getDependentCurves(),
                                                     existingDepKeys,
                                                     getExistingDepSize,
                                                     "dependent",
                                                     record->getId().getId(),
                                                     cs_name,
                                                     -1))
                        {
                            return false;
                        }

                        // Validate independent curves.
                        std::set<std::string> existingIndepKeys;
                        if (info.has_path(cs_path + cs_name + "/independent/")) {
                          for (const auto &key : info[cs_path + cs_name + "/independent/"].child_names()) {
                            existingIndepKeys.insert(key);
                          }
                        }

                        auto getExistingIndepSize = [&info, cs_path, cs_name](const std::string &key) -> int {
                            std::string fullPath = cs_path + cs_name + "/independent/" + key + "/value";
                            if (info.has_path(fullPath)) {
                                return info[fullPath]["num_elements"].to_int();
                            }
                            return -1;
                        };

                        if (!new_curve_set.getIndependentCurves().empty() &&
                            !validate_curves_unified(new_curve_set.getIndependentCurves(),
                                                     existingIndepKeys,
                                                     getExistingIndepSize,
                                                     "independent",
                                                     record->getId().getId(),
                                                     cs_name,
                                                     -1))
                        {
                            return false;
                        }
                    }
                } catch (const std::exception &e) {
                    std::cerr << "Error validating curve sets: " << e.what() << std::endl;
                    return false;
                }
                break; // Found matching record; proceed to next.
            }
        }
    }
    return true;
}

bool append_to_hdf5(const std::string &hdf5FilePath,
                 const Document &newData,
                 const int data_protocol,
                 const int udc_protocol)
{
 conduit::relay::io::IOHandle existing_file;
 conduit::Node to_load;
 conduit::Node to_set;
 conduit::Node opts;
 conduit::Node info;
 opts["stride"] = 1;
 
 // Open the existing HDF5 file.
 try {
     existing_file.open(hdf5FilePath);
     conduit::relay::io::hdf5_read_info(hdf5FilePath, info);
 } catch (const std::exception &e) {
     std::cerr << "Error loading HDF5 file: " << e.what() << std::endl;
     return false;
 }
 
 std::string cs_path, data_path, udc_path;
 int extension = info["records"].number_of_children();

 if (!validate_curve_sets_hdf5(newData, existing_file, info)) {
     return false;
 }
 
 std::vector<std::string> record_list;
 existing_file.list_child_names("records", record_list);
 
 // Queue for write and remove operations.
 std::vector<std::function<void()>> write_queue;
 
 for (auto &record : newData.getRecords()) {
     bool found = false;
     for (const auto &rec_name : record_list) {
         std::string id_path = "records/" + rec_name + "/id/";
         conduit::Node id;
         existing_file.read(id_path, id);
         if (id.to_string() == "\"" + record->getId().getId() + "\"") {
             found = true;
             std::string record_path = "records/" + rec_name;
             
             // ----------------------
             // Queue update of DATA VALUES.
             // ----------------------
             try {
                 data_path = "records/" + rec_name + "/data/";
                 std::vector<std::string> data_keys;
                 existing_file.list_child_names(data_path, data_keys);
                 auto &new_data_sets = record->getData();
                 for (auto &new_data : new_data_sets) {
                     auto new_data_key = new_data.first;
                     auto new_data_pair = new_data.second;
                     if (std::find(data_keys.begin(), data_keys.end(), new_data_key) != data_keys.end()) {
                         // Duplicate Handling
                         switch(data_protocol) {
                             case 1:
                                 write_queue.push_back([&existing_file, data_path, new_data_key, new_data_pair]() {
                                     existing_file.remove(data_path + new_data_key);
                                     existing_file.write(new_data_pair.toNode(), data_path + new_data_key);
                                 });
                                 break;
                             case 2:
                                 break;
                             case 3:
                                 std::cerr << "Found a duplicate data entry, protocol 3 dictates append cancellation." << std::endl;
                                 return false;
                             default:
                                 std::cerr << "Invalid Data Protocol Entry. Append cancelled." << std::endl;
                                 return false;
                         }
                     } else {
                         write_queue.push_back([&existing_file, data_path, new_data_key, new_data_pair]() {
                             existing_file.write(new_data_pair.toNode(), data_path + new_data_key);
                         });
                     }
                 }
             } catch (const std::exception &e) {
                 std::cerr << "Error updating data values: " << e.what() << std::endl;
             }
                             
             // ----------------------
             // Queue update of USER DEFINED CONTENT.
             // ----------------------
             try {
                 udc_path = "records/" + rec_name + "/user_defined/";
                 std::vector<std::string> udc_list;
                 existing_file.list_child_names(udc_path, udc_list);
                 auto &new_udc_sets = record->getUserDefinedContent();
                 for (auto &new_udc : new_udc_sets.children()) {
                     std::string udc_name = new_udc.name();
                     if (std::find(udc_list.begin(), udc_list.end(), udc_name) != udc_list.end()) {
                         // Duplicate Handling
                         switch(udc_protocol) {
                             case 1:
                                 write_queue.push_back([&existing_file, udc_path, udc_name, new_udc]() {
                                     existing_file.remove(udc_path + udc_name);
                                     existing_file.write(new_udc, udc_path + udc_name);
                                 });
                                 break;
                             case 2:
                                 break;
                             case 3:
                                 std::cerr << "Found duplicate UDC, protocol 3 dictates append cancellation." << std::endl;
                                 return false;
                             default:
                                 std::cerr << "Invalid UDC Protocol Entry. Append cancelled." << std::endl;
                                 return false;
                         }
                     } else {
                         write_queue.push_back([&existing_file, udc_path, udc_name, new_udc]() {
                             existing_file.write(new_udc, udc_path + udc_name);
                         });
                     }
                 }
             } catch (const std::exception &e) {
                 std::cerr << "Error updating user defined content: " << e.what() << std::endl;
             }
             
             // ----------------------
             // Queue update of CURVE SETS.
             // ----------------------
             try {
                 cs_path = "records/" + rec_name + "/curve_sets/";
                 std::vector<std::string> curve_sets_list;
                 existing_file.list_child_names("records/" + rec_name + "/curve_sets", curve_sets_list);
                 auto &new_curve_sets = record->getCurveSets();
                 
                 for (const auto &curve_set_pair : new_curve_sets)
                 {
                     const std::string &cs_name = curve_set_pair.first;
                     const auto &cs_data = curve_set_pair.second;
                     
                     bool cs_exists = (std::find(curve_sets_list.begin(),
                                                 curve_sets_list.end(), cs_name)
                                                 != curve_sets_list.end());

                     if (cs_exists) {
                         const auto &new_dependents = cs_data.getDependentCurves();
                         for (const auto &dep_pair : new_dependents)
                         {
                             const std::string &dep_name = dep_pair.first;
                             std::string curve_set_path = record_path + "/curve_sets/" + cs_name +
                                                         "/dependent/" + dep_name + "/value";
                             
                             if (info.has_path(curve_set_path) &&
                                 info[curve_set_path].has_child("num_elements"))
                             {
                                 int current_size = info[curve_set_path]["num_elements"].to_int();
                                 opts["offset"] = current_size;
                             }
                             else
                             {
                                 opts["offset"] = 0;
                             }
                             
                             const auto &new_dep_values = dep_pair.second.getValues();
                             std::vector<double> current_array(new_dep_values.begin(), new_dep_values.end());
                             to_set.set(current_array);
                             
                             write_queue.push_back([curve_set_path, hdf5FilePath, opts, to_set]() mutable {
                                 conduit::relay::io::hdf5_write(to_set, hdf5FilePath,
                                                                 curve_set_path, opts, true);
                             });
                         }
                         
                         const auto &new_independents = cs_data.getIndependentCurves();
                         for (const auto &indep_pair : new_independents)
                         {
                             const std::string &indep_name = indep_pair.first;
                             std::string curve_set_path = record_path + "/curve_sets/" + cs_name +
                                                         "/independent/" + indep_name + "/value";
                             
                             if (info.has_path(curve_set_path) &&
                                 info[curve_set_path].has_child("num_elements"))
                             {
                                 int current_size = info[curve_set_path]["num_elements"].to_int();
                                 opts["offset"] = current_size;
                             }
                             else
                             {
                                 opts["offset"] = 0;
                             }
                             
                             const auto &new_indep_values = indep_pair.second.getValues();
                             std::vector<double> current_array(new_indep_values.begin(), new_indep_values.end());
                             to_set.set(current_array);
                             
                             write_queue.push_back([curve_set_path, hdf5FilePath, opts, to_set]() mutable {
                                 conduit::relay::io::hdf5_write(to_set, hdf5FilePath,
                                                                 curve_set_path, opts, true);
                             });
                         }
                     } else {
                         conduit::Node cs_node;
                         const auto &new_dependents = cs_data.getDependentCurves();
                         for (const auto &dep_pair : new_dependents)
                         {
                             const std::string &dep_name = dep_pair.first;
                             const auto &new_dep_values = dep_pair.second.getValues();
                             std::vector<double> dep_array(new_dep_values.begin(), new_dep_values.end());
                             cs_node["dependent"][dep_name]["value"].set(dep_array);
                         }
                         
                         // Process independent curves.
                         const auto &new_independents = cs_data.getIndependentCurves();
                         for (const auto &indep_pair : new_independents)
                         {
                             const std::string &indep_name = indep_pair.first;
                             const auto &new_indep_values = indep_pair.second.getValues();
                             std::vector<double> indep_array(new_indep_values.begin(), new_indep_values.end());
                             cs_node["independent"][indep_name]["value"].set(indep_array);
                         }
                         
                         write_queue.push_back([record_path, cs_name, hdf5FilePath, cs_node]() mutable {
                             conduit::relay::io::hdf5_write(cs_node, hdf5FilePath,
                                                             record_path + "/curve_sets/" + cs_name,
                                                             conduit::Node(), true);
                         });
                     }
                 }
             } catch (const std::exception &e) {
                 std::cerr << "Error updating curve sets: " << e.what() << std::endl;
             }
         }
     }
     if (!found) {
         std::string record_path = "records/" + std::to_string(extension++) + "/";
         write_queue.push_back([&existing_file, record_path, rec_ptr = record.get()]() {
             existing_file.write(rec_ptr->toNode(), record_path);
         });
     }
 }
 
 // Execute all queued operations.
 for (auto &write_op : write_queue) {
     try {
         write_op();
     } catch (const std::exception &e) {
         std::cerr << "Error executing queued write operation: " << e.what() << std::endl;
         return false;
     }
 }
 
 return true;
}
#endif

}  // namespace sina
}  // namespace axom
