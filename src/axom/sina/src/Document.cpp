/// @file

#include "sina/Document.hpp"

#include <cstdio>
#include <fstream>
#include <ios>
#include <iostream>
#include <utility>
#include <sstream>
#include <stdexcept>

namespace sina {

namespace {
char const RECORDS_KEY[] = "records";
char const RELATIONSHIPS_KEY[] = "relationships";
char const SAVE_TMP_FILE_EXTENSION[] = ".sina.tmp";
}

void Document::add(std::unique_ptr<Record> record) {
    records.emplace_back(std::move(record));
}

void Document::add(Relationship relationship) {
    relationships.emplace_back(std::move(relationship));
}

conduit::Node Document::toNode() const {
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
        RecordLoader const &recordLoader) {
    if (asNode.has_child(RECORDS_KEY)) {
        conduit::Node record_nodes = asNode[RECORDS_KEY];
        if (record_nodes.dtype().is_list()) {
            auto recordIter = record_nodes.children();
            while (recordIter.has_next()){
                auto record = recordIter.next();
                add(recordLoader.load(record));
            }
    } else {
        std::ostringstream message;
        message << "The '" << RECORDS_KEY
                << "' element of a document must be an array";
        throw std::invalid_argument(message.str());
    }
    }

    if (asNode.has_child(RELATIONSHIPS_KEY)){
        conduit::Node relationship_nodes = asNode[RELATIONSHIPS_KEY];
        if (relationship_nodes.dtype().is_list()) {
          auto relationshipsIter = relationship_nodes.children();
          while (relationshipsIter.has_next()){
              auto &relationship = relationshipsIter.next();
              add(Relationship{relationship});
          }
        } else {
            std::ostringstream message;
            message << "The '" << RELATIONSHIPS_KEY
                    << "' element of a document must be an array";
            throw std::invalid_argument(message.str());
        }
    }
}

Document::Document(conduit::Node const &asNode,
        RecordLoader const &recordLoader) {
    this->createFromNode(asNode, recordLoader);
}

Document::Document(std::string const &asJson, RecordLoader const &recordLoader) {
  conduit::Node asNode;
  asNode.parse(asJson, "json");
  this->createFromNode(asNode, recordLoader);
}

std::string Document::toJson(conduit::index_t indent, conduit::index_t depth,
        const std::string &pad, const std::string &eoe) const{
  return this->toNode().to_json("json", indent, depth, pad, eoe);
}

void saveDocument(Document const &document, std::string const &fileName) {
    // It is a common use case for users to want to overwrite their files as
    // the simulation progresses. However, this operation should be atomic so
    // that if a write fails, the old file is left intact. For this reason,
    // we write to a temporary file first and then move the file. The temporary
    // file is in the same directory to ensure that it is part of the same
    // file system as the destination file so that the move operation is
    // atomic.
    std::string tmpFileName = fileName + SAVE_TMP_FILE_EXTENSION;
    auto asJson = document.toJson();
    std::ofstream fout{tmpFileName};
    fout.exceptions(std::ostream::failbit | std::ostream::badbit);
    fout << asJson;
    fout.close();

    if (rename(tmpFileName.c_str(), fileName.c_str()) != 0) {
        std::string message{"Could not save to '"};
        message += fileName;
        message += "'";
        throw std::ios::failure{message};
    }
}

Document loadDocument(std::string const &path) {
    return loadDocument(path, createRecordLoaderWithAllKnownTypes());
}

Document loadDocument(std::string const &path,
        RecordLoader const &recordLoader) {
    conduit::Node nodeFromJson;
    std::ifstream file_in{path};
    std::ostringstream file_contents;
    file_contents << file_in.rdbuf();
    file_in.close();
    nodeFromJson.parse(file_contents.str(), "json");
    return Document{nodeFromJson, recordLoader};
}

}
