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
#include "axom/core/utilities/StringUtilities.hpp"

#include "conduit.hpp"
#ifdef AXOM_USE_HDF5
  #include "conduit_relay.hpp"
  #include "conduit_relay_io.hpp"
#endif

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
#include <algorithm>
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

void Document::createFromNode(const conduit::Node &asNode,
                              const RecordLoader &recordLoader)
{
  conduit::Node nodeCopy = asNode;

  auto processChildNodes = [&](const char *key,
                               std::function<void(conduit::Node &)> addFunc) {
    if(nodeCopy.has_child(key))
    {
      conduit::Node &childNodes = nodeCopy[key];

      // -- 1. Check if this node is a primitive leaf (throw immediately if so)
      // Customize these checks to match exactly what you consider "primitive."
      if(childNodes.dtype().is_number() || childNodes.dtype().is_char8_str() ||
         childNodes.dtype().is_string())
      {
        std::ostringstream message;
        message << "The '" << key
                << "' element of a document cannot be a primitive value.";
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
        message << "The '" << key
                << "' element of a document must be an array/list.";
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
  processChildNodes(RECORDS_KEY, [&](conduit::Node &record) {
    add(recordLoader.load(record));
  });

  processChildNodes(RELATIONSHIPS_KEY, [&](conduit::Node &relationship) {
    add(Relationship {relationship});
  });
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
    std::string modifiedKey =
      axom::utilities::string::replaceAllInstances(key, "/", slashSubstitute);

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

      // Leaves empty nodes empty, if null data is set the
      // Document breaks

      if(it.node().dtype().is_string() || it.node().dtype().is_number())
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

      // Leaves empty keys empty but continues recursive call if its a list
      if(it.node().dtype().is_string() || it.node().dtype().is_number() ||
         it.node().dtype().is_object())
      {
        newChild.set(it.node());
      }
      else if(it.node().dtype().is_list())
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

std::string Document::toJson(conduit::index_t indent,
                             conduit::index_t depth,
                             const std::string &pad,
                             const std::string &eoe) const
{
  return this->toNode().to_json("json", indent, depth, pad, eoe);
}

void saveDocument(Document const &document,
                  std::string const &fileName,
                  Protocol protocol)
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

Document loadDocument(std::string const &path,
                      RecordLoader const &recordLoader,
                      Protocol protocol)
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

bool validate_curve_sets_json(const DataHolder::CurveSetMap new_curve_sets, const nlohmann::json& existing_curve_sets, const std::string id) {
  int largest = 0;
  for (auto& [existing_key, existing_curve_set] : existing_curve_sets.items()) {
      bool match_found = true;

      if (new_curve_sets.find(existing_key) != new_curve_sets.end()) {
          bool dep_match = true;
          bool indep_match = true;
          const auto& new_curve_set = new_curve_sets.at(existing_key);
          if (((new_curve_set.getDependentCurves().size() > 0) && existing_curve_set.contains("dependent"))) {
              for (auto& [dep_key, dep_item] : existing_curve_set["dependent"].items()) {
                  auto& dependents = new_curve_set.getDependentCurves();
                  if (dependents.find(dep_key) != dependents.end()) {
                      size_t new_size = dependents.at(dep_key).getValues().size();
                      if (largest == 0) {
                          largest = dep_item["value"].size() + static_cast<int>(new_size);
                      } else if (largest != 0 && largest != ((int)dep_item["value"].size() + static_cast<int>(new_size))) {
                          std::cerr << "Error validating dependents: Record " << id << ", Curve Set " << existing_key << ", Dependent " << dep_key
                                          << "'s size after append will mismatch with an earlier curve post-append." << std::endl;
                          return false;
                      }
                  } else {
                      size_t new_size = dep_item["value"].size();
                      if (largest == 0) {
                          largest = static_cast<int>(new_size);
                      } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                          std::cerr << "Error validating dependents: Record " << id << ", Curve Set " << existing_key << ", Dependent " << dep_key
                                         << " is not being appended to and will mismatch with curves that were." << std::endl;
                          return false;
                      }
                  }
              }

              for (auto& [new_dep_key, new_dep_item] : new_curve_set.getDependentCurves()) {
                  auto& dependents = existing_curve_set["dependent"];
                  if (!dependents.contains(new_dep_key)) {
                      size_t new_size = new_dep_item.getValues().size();
                      if (largest == 0) {
                          largest = static_cast<int>(new_size);
                      } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                          std::cerr << "Error validating dependents: Record " << id << ", Curve Set " << existing_key << ", Dependent " << new_dep_key
                                         << " will mistmatch with earlier curves in its curve_set if appended to." << std::endl;
                          return false;
                      }
                  }
              }
          } else {
              dep_match = false;
          }

          if (((new_curve_set.getIndependentCurves().size() > 0) && existing_curve_set.contains("independent"))) {
              for (auto& [indep_key, indep_item] : existing_curve_set["independent"].items()) {
                  auto& independents = new_curve_set.getIndependentCurves();
                  if (independents.find(indep_key) != independents.end()) {
                      size_t new_size = independents.at(indep_key).getValues().size();
                      if (largest == 0) {
                          largest = indep_item["value"].size() + static_cast<int>(new_size);
                      } else if (largest != 0 && largest != ((int)indep_item["value"].size() + static_cast<int>(new_size))) {
                         std::cerr << "Error validating independents: Record " << id << ", Curve Set " << existing_key << ", Independent " << indep_key
                                         << "'s size after append will mismatch with an earlier curve post-append." << std::endl;
                         return false;
                      }
                  } else {
                      size_t new_size = indep_item["value"].size();
                      if (largest == 0) {
                          largest = static_cast<int>(new_size);
                      } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                         std::cerr << "Error validating independents: Record " << id << ", Curve Set " << existing_key << ", Independent " << indep_key
                                         << " is not being appended to and will mismatch with curves that were." << std::endl;
                         return false;
                      }
                  }
              }

              for (auto& [new_indep_key, new_indep_item] : new_curve_set.getIndependentCurves()) {
                  auto& independents = existing_curve_set["independent"];
                  if (!independents.contains(new_indep_key)) {
                      size_t new_size = new_indep_item.getValues().size();
                      if (largest == 0) {
                          largest = static_cast<int>(new_size);
                      } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                         std::cerr << "Error validating independents: Record " << id << ", Curve Set " << existing_key << ", Independent " << new_indep_key
                                         << " will mistmatch with earlier curves in its curve_set if appended to." << std::endl;
                         return false;
                      }
                  }
              }
          } else {
              indep_match = false;
          }
          
          if (dep_match && indep_match) {
              match_found = true;
          } else if ((!(new_curve_set.getDependentCurves().size() > 0) && existing_curve_set.contains("dependent")) && (!(new_curve_set.getIndependentCurves().size() > 0) && existing_curve_set.contains("independent"))) {
              match_found = true;
          }
      }

      if (!match_found) {
          std::cerr << "No matching existing curve set found for appending." << std::endl;
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
             //--- Update Curve Sets ---
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
             
             //--- Update Data Sets ---
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
             
             //--- Update User Defined Content ---
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

bool validate_curve_sets_hdf5(const Document &newData, conduit::relay::io::IOHandle &existing_file, const conduit::Node &info) {
 std::vector<std::string> record_list;
 std::string curve_set_path;
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
                 existing_file.list_child_names("records/" + rec_name + "/curve_sets", curve_sets_list);
                 auto &new_curve_sets = record->getCurveSets();

                 for (const auto &cs_name : curve_sets_list) {
                     // Locate the new curve set information.
                     auto cs_itr = new_curve_sets.find(cs_name);
                     if (cs_itr == new_curve_sets.end())
                         continue;
                     
                     const auto &new_dependents  = cs_itr->second.getDependentCurves();
                     const auto &new_independents = cs_itr->second.getIndependentCurves();
                     
                     // Process dependent curves.
                     std::string dependent_path = cs_path + cs_name + "/dependent/";

                     bool need_dep_flag = false;
                     int size_count = -1;
                     
                     for (const auto &dep_pair : new_dependents) {
                         curve_set_path = dependent_path + dep_pair.first + "/value";

                         if(info.has_path(curve_set_path))
                         {
                             if (need_dep_flag && info[curve_set_path]["num_elements"].to_int() > 0) {
                                 std::cerr << "Error validating dependents: " << curve_set_path << " will mistmatch with earlier curves in its curve_set if appended to." << std::endl;
                                 return false;
                             } else {
                                 if (size_count == -1) {
                                     // First dependent has increased its size
                                     size_count = dep_pair.second.getValues().size() + info[curve_set_path]["num_elements"].to_int();
                                 } else if (size_count != dep_pair.second.getValues().size() + info[curve_set_path]["num_elements"].to_int()) {
                                     std::cerr << "Error validating dependents: " << curve_set_path << "'s size after append will mismatch with an earlier curve post-append." << std::endl;
                                     return false;
                                 }
                             }
                         }
                         else {
                             if (size_count == -1) {
                                 // No more dependents are allowed to increase size
                                 need_dep_flag = true;
                             } else {
                                 std::cerr << "Error validating dependents: " << curve_set_path << " is not being appended to and will mismatch with curves that were."  << std::endl;
                                 return false;
                             }
                         }
                     }
                     // Made it out of For loop, either all deps increased by same increment or none did

                     // Process independent curves.
                     std::string independent_path = cs_path + cs_name + "/independent/";

                     bool need_indep_flag = false;
                     
                     for (const auto &indep_pair : new_independents) {
                         curve_set_path = independent_path + indep_pair.first + "/value";

                         if(info.has_path(curve_set_path))
                         {
                             if (need_indep_flag && info[curve_set_path]["num_elements"].to_int() > 0) {
                                 std::cerr << "Error validating independents: " << curve_set_path << " will mistmatch with earlier curve in its curve_set if appended to." << std::endl;
                                 return false;
                             } else {
                                 if (size_count == -1) {
                                     // First independent has increased its size
                                     size_count = indep_pair.second.getValues().size() + info[curve_set_path]["num_elements"].to_int();
                                 } else if (size_count != indep_pair.second.getValues().size() + info[curve_set_path]["num_elements"].to_int()) {
                                     std::cerr << "Error validating independents: " << curve_set_path << "'s size after append will mismatch with an earlier curve post-append." << std::endl;
                                     return false;
                                 }
                             }
                         }
                         else {
                             if (size_count == -1) {
                                 // No more independents are allowed to increase size
                                 need_indep_flag = true;
                             } else {
                                 std::cerr << "Error validating independents: " << curve_set_path << " is not being appended to and will mismatch with curves that were." << std::endl;
                                 return false;
                             }
                         }
                     }
                     // Made it out of For loop, either all indeps increased by same increment or none did
                     // Since both passed, restart with the next curve set.
                 }
                 // Made it through this curve set, on to the next

             } catch (const std::exception &e) {
                 std::cerr << "Error validating curve sets: " << e.what() << std::endl;
             }
             break;
         }
         // Onto the next record
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

}  // namespace sina
}  // namespace axom
