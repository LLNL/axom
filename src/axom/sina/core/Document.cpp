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
#include "axom/json.hpp"
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
}

void protocol_warn(std::string protocol, std::string const &name) {
  size_t pos = name.rfind('.');

  if (pos != std::string::npos) {
        std::string found = name.substr(pos+1);

        if (("." + found) != protocol) {
          for (const std::string& file_type : supported_types) {
            std::string lower_case_type = file_type; 
            std::transform(lower_case_type.begin(), lower_case_type.end(), lower_case_type.begin(), [](unsigned char c) { return std::tolower(c); });

            if ((lower_case_type != protocol) && (lower_case_type == found)) {
              std::cout << "|| WARNING: INCORRECT FILE EXTENSION FOUND (FOUND: " << found << "; EXPECTED: " << protocol << "), DID YOU MEAN TO INCLUDE ONE OF " << protocol << "'s SUPPORTED TYPES? ||";
            }
          }
          std::cout << "|| WARNING: BROKEN FILE EXTENSION FOUND, SINA WILL BE APPENDING THE " << protocol << " FILE EXTENSION ||";
        } else {
          return;
        }
  } else {
    std::cout << "|| WARNING: NO FILE EXTENSION FOUND, DID YOU MEAN TO INCLUDE ONE OF " << protocol << "'s SUPPORTED TYPES?";
  }
}

void removeSlashes(const conduit::Node& originalNode, conduit::Node& modifiedNode)
{
    for (auto it = originalNode.children(); it.has_next();)
    {
        it.next();
        std::string key = it.name();
        std::string modifiedKey = key;

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
        conduit::Node modifiedRelationshipNode;
        removeSlashes(relationshipNode, modifiedRelationshipNode);
        relationshipsNode.append() = modifiedRelationshipNode;
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

  try
  {
      if (protocol == Protocol::JSON)
      {
          protocol_warn(".json", fileName);
          auto asJson = document.toJson();
          std::ofstream fout {tmpFileName};
          fout.exceptions(std::ostream::failbit | std::ostream::badbit);
          fout << asJson;
          fout.close();
      }
      else if (protocol == Protocol::HDF5)
      {
          protocol_warn(".hdf5", fileName);
          document.toHDF5(tmpFileName);
      }
      else
      {
          throw std::invalid_argument("Invalid format choice. Please enter 'json' or 'hdf5'.");
      }

      if (rename(tmpFileName.c_str(), fileName.c_str()) != 0)
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

bool validate_curve_sets_json(const DataHolder::CurveSetMap new_curve_sets, const nlohmann::json& existing_curve_sets) {
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
                            std::cerr << "Appending curve set lengths mismatch at Dependent: " << dep_key << std::endl;
                            return false;
                        }
                    } else {
                        size_t new_size = dep_item["value"].size();
                        if (largest == 0) {
                            largest = static_cast<int>(new_size);
                        } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                            std::cerr << "Appending curve set length mismatch at Existing Dependent: " << dep_key << std::endl;
                            return false;
                        }
                    }
                }

                for (auto& [new_dep_key, new_dep_item] : new_curve_set.getDependentCurves()) {
                    auto& dependents = existing_curve_set["dependent"];
                    if (dependents.contains(new_dep_key)) {
                        size_t new_size = dependents[new_dep_key]["value"].size();
                        if (largest == 0) {
                            largest = new_dep_item.getValues().size() + static_cast<int>(new_size);
                        } else if (largest != 0 && largest != ((int)new_dep_item.getValues().size() + static_cast<int>(new_size))) {
                            std::cerr << "Appending curve set lengths mismatch at Dependent: " << new_dep_key << std::endl;
                            return false;
                        }
                    } else {
                        size_t new_size = new_dep_item.getValues().size();
                        if (largest == 0) {
                            largest = static_cast<int>(new_size);
                        } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                            std::cerr << "Appending curve set length mismatch at New Dependent: " << new_dep_key << std::endl;
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
                            std::cerr << "Appending curve set lengths mismatch at Independents: " << indep_key << std::endl;
                            return false;
                        }
                    } else {
                        size_t new_size = indep_item["value"].size();
                        if (largest == 0) {
                            largest = static_cast<int>(new_size);
                        } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                            std::cerr << "Appending curve set length mismatch at Existing Independent: " << indep_key << std::endl;
                            return false;
                        }
                    }
                }

                for (auto& [new_indep_key, new_indep_item] : new_curve_set.getIndependentCurves()) {
                    auto& independents = existing_curve_set["independent"];
                    if (independents.contains(new_indep_key)) {
                        size_t new_size = independents[new_indep_key]["value"].size();
                        if (largest == 0) {
                            largest = new_indep_item.getValues().size() + static_cast<int>(new_size);
                        } else if (largest != 0 && largest != ((int)new_indep_item.getValues().size() + static_cast<int>(new_size))) {
                            std::cerr << "Appending curve set lengths mismatch at Independents: " << new_indep_key << std::endl;
                            return false;
                        }
                    } else {
                        size_t new_size = new_indep_item.getValues().size();
                        if (largest == 0) {
                            largest = static_cast<int>(new_size);
                        } else if (largest != 0 && largest != static_cast<int>(new_size)) {
                            std::cerr << "Appending curve set length mismatch at New Independent: " << new_indep_key << std::endl;
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

bool append_to_json(const std::string& jsonFilePath, Document const &newData) {
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

    for (auto& new_record : newData.getRecords()) {
        bool found = false;
        for (auto& existing_record : j["records"]) {
            if (new_record->getId().getId() == existing_record["id"]) {
                found = true;
                if ((new_record->getCurveSets().size() > 0) && existing_record.contains("curve_sets")) {
                    if (validate_curve_sets_json(new_record->getCurveSets(), existing_record["curve_sets"])) {
                        nlohmann::json& existing_curve_sets = existing_record["curve_sets"];
                        auto& new_curve_sets = new_record->getCurveSets();
                        for (auto& new_curve_set : new_curve_sets) {
                            auto& [new_cs_key, new_cs_values] = new_curve_set;
                            for (auto& new_dependent : new_cs_values.getDependentCurves()) {
                                auto& [new_dep_key, new_dep_values] = new_dependent;
                                auto& to_append = new_dep_values.getValues();
                                auto& existing_dep_values = existing_curve_sets[new_cs_key]["dependent"][new_dep_key]["value"];
                                for (const auto& value : to_append) {
                                    existing_dep_values.push_back(static_cast<double>(value));
                                }
                            }

                            for (auto& new_independent : new_cs_values.getIndependentCurves()) {
                                auto& [new_indep_key, new_indep_values] = new_independent;
                                auto& to_append = new_indep_values.getValues();
                                auto& existing_indep_values = existing_curve_sets[new_cs_key]["independent"][new_indep_key]["value"];
                                for (const auto& value : to_append) {
                                    existing_indep_values.push_back(static_cast<double>(value));
                                }
                            }
                        }
                    } else {
                        std::cerr << "There was an error validating the curve sets at id: " << existing_record["id"] << std::endl;
                        return false;
                    }
                }


                if ((new_record->getData().size() >0) && existing_record.contains("data")) {
                    nlohmann::json& existing_data_sets = existing_record["data"];
                    auto& new_data_sets = new_record->getData();
                    for (auto& new_data : new_data_sets) {
                        int data_protocol = -1;
                        auto& [new_data_key, new_data_pair] = new_data;
                        json obj = obj.parse(new_data_pair.toNode().to_json());
                        if (existing_data_sets.contains(new_data_key)) {
                            if (data_protocol == -1) {
                                std::cout << "Detected a duplicate data key, would you like to: 1 = overwrite duplicates, 2 = ignore duplicates, 3 = cancel the append";
                                std::cin >> data_protocol;
                            }

                            switch(data_protocol) {
                                case 1:
                                    existing_data_sets[new_data_key] = obj;
                                    break;
                                case 2:
                                    break;
                                case 3:
                                    std::cout << "Append Cancelled";
                                    return false;
                                default:
                                    std::cout << "Invalid Entry";
                                    return false;
                            }
                        } else {
                            existing_data_sets[new_data_key] = obj;
                        }
                    }
                }

                if ((!new_record->getUserDefinedContent().dtype().is_empty()) && existing_record.contains("user_defined")) {
                    nlohmann::json& existing_udc = existing_record["user_defined"];
                    auto& new_udc = new_record->getUserDefinedContent();
                    for (auto& udc : new_udc.children()) {
                        std::string udc_name = udc.name();
                        int udc_protocol = -1;
                        if (existing_record["user_defined"].contains(udc_name)) {
                            if (udc_protocol == -1) {
                                std::cout << "Detected a duplicate user_defined key, would you like to: 1 = overwrite duplicates, 2 = ignore duplicates, 3 = cancel the append";
                                std::cin >> udc_protocol;
                            }

                            switch(udc_protocol) {
                                case 1:
                                    existing_udc[udc_name] = udc.as_string();
                                    break;
                                case 2:
                                    break;
                                case 3:
                                    std::cout << "Append Cancelled";
                                    return false;
                                default:
                                    std::cout << "Invalid Entry";
                                    return false;
                            }
                        } else {
                            existing_udc[udc_name] = udc.as_string();
                        }
                    }
                }
            }
        }
        if (!found) {
            json obj = obj.parse(new_record->toNode().to_json());
            j["records"][j["records"].size()] = obj;
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

bool validate_curve_sets_hdf5(const DataHolder::CurveSetMap &new_curve_sets, const conduit::Node &existing_curve_sets) {
    int largest = 0;

    // Iterate over existing curve sets
    std::vector<std::string> existing_keys = existing_curve_sets.child_names();
    for (auto &existing_key : existing_keys) {
        const conduit::Node &existing_curve_set = existing_curve_sets[existing_key];
        bool match_found = true;

        if (new_curve_sets.find(existing_key) != new_curve_sets.end()) {
            bool dep_match = true;
            bool indep_match = true;
            const auto &new_curve_set = new_curve_sets.at(existing_key);

            // Validate dependent curves
            if (new_curve_set.getDependentCurves().size() > 0 && existing_curve_set.has_child("dependent")) {
                const conduit::Node &existing_dep = existing_curve_set["dependent"];
                // Check existing dependents
                for (auto dep_key : existing_dep.child_names()) {
                    const conduit::Node &dep_item = existing_dep[dep_key];
                    const auto &dependents = new_curve_set.getDependentCurves();
                    if (dependents.find(dep_key) != dependents.end()) {
                        size_t new_size = dependents.at(dep_key).getValues().size();
                        int total_size = (int)dep_item["value"].dtype().number_of_elements() + (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set lengths mismatch at Dependent: " << dep_key << std::endl;
                            return false;
                        }
                    } else {
                        size_t new_size = dep_item["value"].dtype().number_of_elements();
                        int total_size = (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set length mismatch at Existing Dependent: " << dep_key << std::endl;
                            return false;
                        }
                    }
                }

                // Check new dependents
                for (auto &new_dep : new_curve_set.getDependentCurves()) {
                    auto &new_dep_key = new_dep.first;
                    auto &new_dep_item = new_dep.second;
                    size_t new_size = new_dep_item.getValues().size();

                    if (existing_dep.has_child(new_dep_key)) {
                        // Already exists, check combined size
                        size_t old_size = existing_dep[new_dep_key]["value"].dtype().number_of_elements();
                        int total_size = (int)old_size + (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set lengths mismatch at Dependent: " << new_dep_key << std::endl;
                            return false;
                        }
                    } else {
                        // new dep curve
                        int total_size = (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set length mismatch at New Dependent: " << new_dep_key << std::endl;
                            return false;
                        }
                    }
                }
            } else {
                dep_match = false;
            }

            // Validate independent curves
            if (new_curve_set.getIndependentCurves().size() > 0 && existing_curve_set.has_child("independent")) {
                const conduit::Node &existing_indep = existing_curve_set["independent"];
                // Check existing independents
                for (auto indep_key : existing_indep.child_names()) {
                    const conduit::Node &indep_item = existing_indep[indep_key];
                    const auto &independents = new_curve_set.getIndependentCurves();
                    if (independents.find(indep_key) != independents.end()) {
                        size_t new_size = independents.at(indep_key).getValues().size();
                        int total_size = (int)indep_item["value"].dtype().number_of_elements() + (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set lengths mismatch at Independent: " << indep_key << std::endl;
                            return false;
                        }
                    } else {
                        size_t new_size = indep_item["value"].dtype().number_of_elements();
                        int total_size = (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set length mismatch at Existing Independent: " << indep_key << std::endl;
                            return false;
                        }
                    }
                }

                // Check new independents
                for (auto &new_indep : new_curve_set.getIndependentCurves()) {
                    auto &new_indep_key = new_indep.first;
                    auto &new_indep_item = new_indep.second;
                    size_t new_size = new_indep_item.getValues().size();
                    if (existing_indep.has_child(new_indep_key)) {
                        size_t old_size = existing_indep[new_indep_key]["value"].dtype().number_of_elements();
                        int total_size = (int)old_size + (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set lengths mismatch at Independent: " << new_indep_key << std::endl;
                            return false;
                        }
                    } else {
                        int total_size = (int)new_size;
                        if (largest == 0) {
                            largest = total_size;
                        } else if (largest != total_size) {
                            std::cerr << "Appending curve set length mismatch at New Independent: " << new_indep_key << std::endl;
                            return false;
                        }
                    }
                }
            } else {
                indep_match = false;
            }

            if (!dep_match && !indep_match) {
                // check if no new curves but existing sets have them
                if (!((new_curve_set.getDependentCurves().size() == 0 && existing_curve_set.has_child("dependent")) &&
                      (new_curve_set.getIndependentCurves().size() == 0 && existing_curve_set.has_child("independent")))) {
                    match_found = false;
                }
            }

        } else {
            match_found = false;
        }

        if (!match_found) {
            std::cerr << "No matching existing curve set found for appending." << std::endl;
            return false;
        }
    }

    return true;
}


bool append_to_hdf5(const std::string &hdf5FilePath, const Document &newData) {
    conduit::relay::io::IOHandle existing_file;
    conduit::Node to_load;
    conduit::Node to_set;
    // Load existing HDF5 file
    try {
        existing_file.open(hdf5FilePath);
    } catch (const std::exception &e) {
        std::cerr << "Error loading HDF5 file: " << e.what() << std::endl;
        return false;
    }

    // Check number of records
    std::vector<std::string> record_count;
    existing_file.list_child_names("records", record_count);
    int existing_record_count = record_count.size();
    int new_record_count = (int)newData.getRecords().size();
    if (existing_record_count != new_record_count) {
        std::cerr << "Mismatch in the number of records." << std::endl;
        return false;
    }

    conduit::Node opts, placeholder;
    opts["stride"] = 1;
    std::string cs_path, data_path, udc_path;

    std::vector<std::string> record_list;
    std::vector<std::string>::const_iterator rec_itr;
    existing_file.list_child_names("records", record_list);
    int extension = record_list.size();
    // TODO: Validate
    for (auto& record : newData.getRecords()) {
        bool found = false;
        for (rec_itr = record_list.begin(); rec_itr < record_list.end(); ++rec_itr) {
            std::string id_path = "records/" + *rec_itr + "/id/";
            conduit::Node id;
            existing_file.read(id_path, id);
            if (id.to_string() == "\"" + record->getId().getId() + "\"") {
                found = true;

                // DATA VALUES
                try {
                    data_path = "records/" + *rec_itr + "/data/";
                    std::vector<std::string> data_list;
                    existing_file.list_child_names("records/" + *rec_itr + "/data", data_list);
                    auto& new_data_sets = record->getData();
                    for (auto& new_data : new_data_sets) {
                        int data_protocol = -1;
                        auto& [new_data_key, new_data_pair] = new_data;
                        json obj = obj.parse(new_data_pair.toNode().to_json());
                        std::string json_str = obj.dump(4);
                        if (std::find(data_list.begin(), data_list.end(), new_data_key) != data_list.end()) {
                            if (data_protocol == -1) {
                                std::cout << "Detected a duplicate data key, would you like to: 1 = overwrite duplicates, 2 = ignore duplicates, 3 = cancel the append";
                                std::cin >> data_protocol;
                            }

                            switch(data_protocol) {
                                case 1:
                                    existing_file.remove(data_path + new_data_key);
                                    existing_file.write(new_data_pair.toNode(), data_path + new_data_key);
                                    break;
                                case 2:
                                    break;
                                case 3:
                                    std::cout << "Append Cancelled";
                                    return false;
                                default:
                                    std::cout << "Invalid Entry";
                                    return false;
                            }
                        } else {
                            existing_file.write(new_data_pair.toNode(), data_path + new_data_key);
                        }   
                    }
                } catch (const std::exception &e) {

                }


                // USER DEFINED
                try {
                    udc_path = "records/" + *rec_itr + "/user_defined/";
                    std::vector<std::string> udc_list;
                    existing_file.list_child_names("records/" + *rec_itr + "/user_defined", udc_list);
                    auto& new_udc_sets = record->getUserDefinedContent();
                    for (auto& new_udc : new_udc_sets.children()) {
                        int udc_protocol = -1;
                        std::string udc_name = new_udc.name();
                        if (std::find(udc_list.begin(), udc_list.end(), udc_name) != udc_list.end()) {
                            if (udc_protocol == -1) {
                                std::cout << "Detected a duplicate user defined entry, would you like to: 1 = overwrite duplicates, 2 = ignore duplicates, 3 = cancel the append";
                                std::cin >> udc_protocol;
                            }

                            switch(udc_protocol) {
                                case 1:
                                    existing_file.remove(udc_path + udc_name);
                                    existing_file.write(new_udc, udc_path + udc_name);
                                    break;
                                case 2:
                                    break;
                                case 3:
                                    std::cout << "Append Cancelled";
                                    return false;
                                default:
                                    std::cout << "Invalid Entry";
                                    return false;
                            }
                        } else {
                            existing_file.write(new_udc, udc_path + udc_name);
                        }   
                    }
                } catch (const std::exception &e) {

                }


                // CURVE SETS
                try {
                    cs_path = "records/" + *rec_itr + "/curve_sets/";
                    std::vector<std::string> curve_sets_list;
                    std::vector<std::string>::const_iterator cs_itr;
                    existing_file.list_child_names("records/" + *rec_itr+ "/curve_sets", curve_sets_list);
                    auto& new_curve_sets = record->getCurveSets();
                    for (cs_itr = curve_sets_list.begin(); cs_itr < curve_sets_list.end(); ++cs_itr) {
                        const auto &new_cs_values = new_curve_sets.find(*cs_itr);
                        const auto &new_dependents = new_cs_values->second.getDependentCurves();
                        const auto &new_independents = new_cs_values->second.getIndependentCurves();
                        std::string dependent_path = cs_path + *cs_itr + "/dependent/";
                        std::string independent_path = cs_path + *cs_itr + "/independent/";
                        std::vector<std::string> dependents;
                        std::vector<std::string> independents;
                        std::vector<std::string>::const_iterator dep_itr;
                        std::vector<std::string>::const_iterator indep_itr;
                        existing_file.list_child_names(dependent_path, dependents);
                        existing_file.list_child_names(independent_path, independents);
                        std::vector<double> current_array;
                        new_dependents.find(*cs_itr);
                        for (dep_itr = dependents.begin(); dep_itr < dependents.end(); ++dep_itr) {
                            std::string curve_set_path = dependent_path + *dep_itr + "/value";
                            conduit::Node* existing_dep_values = new conduit::Node();
                            existing_file.read(curve_set_path, *existing_dep_values);
                            opts["offset"] = existing_dep_values->dtype().number_of_elements();
                            delete existing_dep_values;
                            const auto &new_dep_values = new_dependents.find(*dep_itr)->second.getValues();
                            for (auto& value : new_dep_values) {
                                current_array.push_back(value);
                            }
                            to_set.set(current_array);
                            current_array.clear();
                            conduit::relay::io::hdf5_write(to_set, hdf5FilePath, curve_set_path, opts, true);
                        }

                        for (indep_itr = independents.begin(); indep_itr < independents.end(); ++indep_itr) {
                            std::string curve_set_path = independent_path + *indep_itr + "/value";
                            conduit::Node* existing_indep_values = new conduit::Node();
                            existing_file.read(curve_set_path, *existing_indep_values);
                            opts["offset"] = existing_indep_values->dtype().number_of_elements();
                            delete existing_indep_values;
                            const auto &new_indep_values = new_independents.find(*indep_itr)->second.getValues();
                            for (auto& value : new_indep_values) {
                                current_array.push_back(value);
                            }
                            to_set.set(current_array);
                            current_array.clear();
                            conduit::relay::io::hdf5_write(to_set, hdf5FilePath, curve_set_path, opts, true);
                        }
                    }
                    break;
                } catch (const std::exception &e) {

                }
            }
        }
        if (!found) {
            existing_file.write(record->toNode(), "records/" + std::to_string(extension++) + "/");
        }
    }

    conduit::relay::io::hdf5_read(hdf5FilePath, to_load);
    to_load.print();
    return false;
}


// bool append_to_hdf5_multiple(const std::string& hdf5FilePath, const json& new_records) {
//     conduit::Node data;
//     conduit::relay::io::hdf5_read(hdf5FilePath, data);

//     if (data["records"].number_of_children() != new_records["records"].size()) {
//         std::cerr << "Mismatch in the number of records." << std::endl;
//         return false;
//     }

//     for (int index = 0; index < data["records"].number_of_children(); ++index) {
//         const auto& existing_record = data["records"][index];
//         const auto& new_record = new_records["records"][index];
//         json obj = obj.parse(existing_record["curve_sets"].to_json());

//         if (!validate_curve_sets_json(obj, new_record["curve_sets"])) {
//             std::cerr << "Validation failed for curve sets at index " << index << std::endl;
//             return false;  
//         }
//     }


//     conduit::Node& records = data["records"];
//     for (int i = 0; i < records.number_of_children(); ++i) {
//         conduit::Node& record = records[i];

//         if (record.has_child("curve_sets")) {
//             conduit::Node& curve_sets = record["curve_sets"];
//             const auto& new_record = new_records["records"][i];
//             int j = 0;
//             for (const auto& new_curve_set : new_record["curve_sets"]) {
//                 conduit::Node& curve_set = curve_sets[j];
//                 if (new_curve_set.contains("dependent")) {
//                     conduit::Node& dependent = curve_set["dependent"];
//                     for (const auto& dep_item : new_curve_set["dependent"].items()) {
//                         const std::string& key = dep_item.key();
//                         conduit::Node& dep_value = dependent[key];
//                         conduit::Node& valueNode = dep_value["value"];
//                         std::vector<double> currentValues;
//                         const conduit::double_array existingArray = valueNode.as_double_array();
//                         currentValues.reserve(existingArray.number_of_elements());
//                         for (conduit::index_t i = 0; i < existingArray.number_of_elements(); ++i) {
//                             currentValues.push_back(existingArray[i]);
//                         }
//                         const auto& new_dep_value = dep_item.value()["value"];

//                         currentValues.insert(currentValues.end(), new_dep_value.begin(), new_dep_value.end());
//                         valueNode.set(currentValues);
//                     }
//                 }

//                 if (new_curve_set.contains("independent")) {
//                     conduit::Node& independent = curve_set["independent"];
//                     for (const auto& indep_item : new_curve_set["independent"].items()) {
//                         const std::string& key = indep_item.key();
//                         conduit::Node& indep_value = independent[key];
//                         conduit::Node& valueNode = indep_value["value"];
//                         std::vector<double> currentValues;
//                         const conduit::double_array existingArray = valueNode.as_double_array();
//                         currentValues.reserve(existingArray.number_of_elements());
//                         for (conduit::index_t i = 0; i < existingArray.number_of_elements(); ++i) {
//                             currentValues.push_back(existingArray[i]);
//                         }
//                         const auto& new_dep_value = indep_item.value()["value"];

//                         currentValues.insert(currentValues.end(), new_dep_value.begin(), new_dep_value.end());
//                         valueNode.set(currentValues);
//                     }
//                 }
//             }
//             j++;
//         }
//     }

//     conduit::relay::io::hdf5_write(data, hdf5FilePath);

//     return true;
// }

// bool append_to_hdf5_one_record(const std::string& hdf5FilePath, const json& new_curve_sets, size_t recordIndex) {
//     conduit::Node data;
//     conduit::relay::io::hdf5_read(hdf5FilePath, data);

//     for (int index = 0; index < data["records"].number_of_children(); ++index) {
//         const auto& existing_record = data["records"][index];
//         json obj = obj.parse(existing_record["curve_sets"].to_json());

//         if (!validate_curve_sets_json(obj, new_curve_sets)) {
//             std::cerr << "Validation failed for curve sets at index " << index << std::endl;
//             return false;  
//         }
//     }


//     conduit::Node& records = data["records"];
//     conduit::Node& record = records[recordIndex];
//     json obj = obj.parse(record["curve_sets"].to_json());

//     if (!validate_curve_sets_json(obj, new_curve_sets)) {
//         std::cerr << "Validation failed for new curve sets." << std::endl;
//         return false;  
//     }

//     if (record.has_child("curve_sets")) {
//         conduit::Node& curve_sets = record["curve_sets"];
//         int j = 0;
//         for (const auto& new_curve_set : new_curve_sets) {
//             conduit::Node& curve_set = curve_sets[j];
//             if (curve_set.has_child("dependent")) {
//                 conduit::Node& dependent = curve_set["dependent"];
//                 for (const auto& dep_item : new_curve_set["dependent"].items()) {
//                     const std::string& key = dep_item.key();
//                     conduit::Node& dep_value = dependent[key];
//                     conduit::Node& valueNode = dep_value["value"];
//                     std::vector<double> currentValues;
//                     const conduit::double_array existingArray = valueNode.as_double_array();
//                     currentValues.reserve(existingArray.number_of_elements());
//                     for (conduit::index_t i = 0; i < existingArray.number_of_elements(); ++i) {
//                         currentValues.push_back(existingArray[i]);
//                     }
//                     const auto& new_dep_value = dep_item.value()["value"];

//                     currentValues.insert(currentValues.end(), new_dep_value.begin(), new_dep_value.end());
//                     valueNode.set(currentValues);
//                 }
//             }

//                 if (curve_set.has_child("independent")) {
//                 conduit::Node& independent = curve_set["independent"];
//                 for (const auto& indep_item : new_curve_set["independent"].items()) {
//                     const std::string& key = indep_item.key();
//                     conduit::Node& indep_value = independent[key];
//                     conduit::Node& valueNode = indep_value["value"];
//                     std::vector<double> currentValues;
//                     const conduit::double_array existingArray = valueNode.as_double_array();
//                     currentValues.reserve(existingArray.number_of_elements());
//                     for (conduit::index_t i = 0; i < existingArray.number_of_elements(); ++i) {
//                         currentValues.push_back(existingArray[i]);
//                     }
//                     const auto& new_dep_value = indep_item.value()["value"];

//                     currentValues.insert(currentValues.end(), new_dep_value.begin(), new_dep_value.end());
//                     valueNode.set(currentValues);
//                 }
//             }
//         }
//         j++;
//     }

//     conduit::relay::io::hdf5_write(data, hdf5FilePath);

//     return true;
// }

}  // namespace sina
}  // namespace axom


