/// @file

#include "sina/Record.hpp"

#include <stdexcept>
#include <utility>

#include "sina/CppBridge.hpp"
#include "sina/ConduitUtil.hpp"
#include "sina/DataHolder.hpp"
#include "sina/Run.hpp"

namespace {

char const LOCAL_ID_FIELD[] = "local_id";
char const GLOBAL_ID_FIELD[] = "id";
char const TYPE_FIELD[] = "type";
char const FILES_FIELD[] = "files";

// Used to preserve information when appending "standalone" library data
char const LIBRARY_DATA_ID_DATUM[] = "SINA_librarydata_id";
char const LIBRARY_DATA_TYPE_DATUM[] = "SINA_librarydata_type";

}

namespace sina {

Record::Record(ID id_, std::string type_) :
        DataHolder{},
        id{std::move(id_), LOCAL_ID_FIELD, GLOBAL_ID_FIELD},
        type{std::move(type_)} {}

conduit::Node Record::toNode() const {
    conduit::Node asNode = DataHolder::toNode();
    asNode[TYPE_FIELD] = type;
    id.addTo(asNode);
    // Optional fields
    if(!files.empty()){
      conduit::Node fileRef;
      for (auto &file : files) {
          auto &n = fileRef.add_child(file.getUri());
          n.set(file.toNode());
      asNode[FILES_FIELD] = fileRef;
      }
    }
    return asNode;
}

Record::Record(conduit::Node const &asNode) :
        DataHolder{asNode},
        id{asNode, LOCAL_ID_FIELD, GLOBAL_ID_FIELD},
        type{getRequiredString(TYPE_FIELD, asNode, "record")} {

    if(asNode.has_child(FILES_FIELD)){
        auto filesIter = asNode[FILES_FIELD].children();
        while(filesIter.has_next()){
            auto &namedFile = filesIter.next();
            files.insert(File(filesIter.name(), namedFile));
        }
    }
}

void Record::remove(File const &file) {
    files.erase(file);
}

void Record::add(File file) {
    files.erase(file);
    files.insert(std::move(file));
}

void Record::addRecordAsLibraryData(Record const &childRecord, std::string const &name) {
    if(!childRecord.files.empty()){
      for (auto &file : childRecord.files) {
          add(file);
      }
    }
    auto newLibData = addLibraryData(name, childRecord.toNode());
    newLibData->add(LIBRARY_DATA_ID_DATUM, Datum{childRecord.getId().getId()});
    newLibData->add(LIBRARY_DATA_TYPE_DATUM, Datum{childRecord.type});
}

void RecordLoader::addTypeLoader(std::string const &type, TypeLoader loader) {
    typeLoaders[type] = std::move(loader);
}

std::unique_ptr<Record>
RecordLoader::load(conduit::Node const &recordAsNode) const {
    auto loaderIter = typeLoaders.find(recordAsNode[TYPE_FIELD].as_string());
    if (loaderIter != typeLoaders.end()) {
        return loaderIter->second(recordAsNode);
    }
    return internal::make_unique<Record>(recordAsNode);
}

bool RecordLoader::canLoad(std::string const &type) const {
    return typeLoaders.count(type) > 0;
}

RecordLoader createRecordLoaderWithAllKnownTypes() {
    RecordLoader loader;
    addRunLoader(loader);
    return loader;
}
}
