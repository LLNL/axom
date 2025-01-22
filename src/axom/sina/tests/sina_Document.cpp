

// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <cstdio>
#include <fstream>
#include <iostream>
#include <type_traits>
#include <utility>

#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "axom/sina/core/Document.hpp"
#include "axom/sina/core/Run.hpp"
#include "axom/json.hpp"
#include "axom/sina/tests/TestRecord.hpp"
#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp"
#include "axom/sina/core/CurveSet.hpp" 
#include "axom/sina/core/Record.hpp" 
#include "conduit.hpp"
#include "conduit_relay.hpp"
#include "conduit_relay_io.hpp" 
#include "conduit_relay_io_hdf5.hpp" 

namespace axom
{
namespace sina
{
namespace testing
{
namespace
{

using ::testing::ElementsAre;
using ::testing::HasSubstr;

char const TEST_RECORD_TYPE[] = "test type";
char const EXPECTED_RECORDS_KEY[] = "records";
char const EXPECTED_RELATIONSHIPS_KEY[] = "relationships";


// Helper function to convert Conduit Node array to std::vector<double> for HDF5 assertion
std::vector<double> node_to_double_vector(const conduit::Node& node) {
    std::vector<double> result;

    if (node.dtype().is_number()) {  
        const double* intArray = node.as_double_ptr(); 
        conduit::index_t numElements = node.dtype().number_of_elements();
        for (conduit::index_t i = 0; i < numElements; ++i) {
            result.push_back(intArray[i]);
        }
    }
    return result;
}

TEST(Document, create_fromNode_empty)
{
  conduit::Node documentAsNode;
  RecordLoader loader;
  Document document {documentAsNode, loader};
  EXPECT_EQ(0u, document.getRecords().size());
  EXPECT_EQ(0u, document.getRelationships().size());
}

TEST(Document, create_fromNode_wrongRecordsType)
{
  conduit::Node recordsAsNodes;
  recordsAsNodes[EXPECTED_RECORDS_KEY] = 123;
  RecordLoader loader;
  try
  {
    Document document {recordsAsNodes, loader};
    FAIL() << "Should not have been able to parse records. Have "
           << document.getRecords().size();
  }
  catch(std::invalid_argument const &expected)
  {
    EXPECT_THAT(expected.what(), HasSubstr(EXPECTED_RECORDS_KEY));
  }
}

TEST(Document, create_fromNode_withRecords)
{
  conduit::Node recordAsNode;
  recordAsNode["type"] = "IntTestRecord";
  recordAsNode["id"] = "the ID";
  recordAsNode[TEST_RECORD_VALUE_KEY] = 123;

  conduit::Node recordsAsNodes;
  recordsAsNodes.append().set(recordAsNode);

  conduit::Node documentAsNode;
  documentAsNode[EXPECTED_RECORDS_KEY] = recordsAsNodes;

  RecordLoader loader;
  loader.addTypeLoader("IntTestRecord", [](conduit::Node const &asNode) {
    return std::make_unique<TestRecord<int>>(asNode);
  });

  Document document {documentAsNode, loader};
  auto &records = document.getRecords();
  ASSERT_EQ(1u, records.size());
  auto testRecord = dynamic_cast<TestRecord<int> const *>(records[0].get());
  ASSERT_NE(nullptr, testRecord);
  ASSERT_EQ(123, testRecord->getValue());
}

TEST(Document, create_fromNode_withRelationships)
{
  conduit::Node relationshipAsNode;
  relationshipAsNode["subject"] = "the subject";
  relationshipAsNode["object"] = "the object";
  relationshipAsNode["predicate"] = "is related to";

  conduit::Node relationshipsAsNodes;
  relationshipsAsNodes.append().set(relationshipAsNode);

  conduit::Node documentAsNode;
  documentAsNode[EXPECTED_RELATIONSHIPS_KEY] = relationshipsAsNodes;

  Document document {documentAsNode, RecordLoader {}};
  auto &relationships = document.getRelationships();
  ASSERT_EQ(1u, relationships.size());
  EXPECT_EQ("the subject", relationships[0].getSubject().getId());
  EXPECT_EQ(IDType::Global, relationships[0].getSubject().getType());
  EXPECT_EQ("the object", relationships[0].getObject().getId());
  EXPECT_EQ(IDType::Global, relationships[0].getObject().getType());
  EXPECT_EQ("is related to", relationships[0].getPredicate());
}


TEST(Document, toNode_empty)
{
  // A sina document should always have, at minimum, both records and
  // relationships as empty arrays.
  Document const document;
  conduit::Node asNode = document.toNode();
  EXPECT_TRUE(asNode[EXPECTED_RECORDS_KEY].dtype().is_list());
  EXPECT_EQ(0, asNode[EXPECTED_RECORDS_KEY].number_of_children());
  EXPECT_TRUE(asNode[EXPECTED_RELATIONSHIPS_KEY].dtype().is_list());
  EXPECT_EQ(0, asNode[EXPECTED_RELATIONSHIPS_KEY].number_of_children());
}

TEST(Document, toNode_records)
{
  Document document;
  std::string expectedIds[] = {"id 1", "id 2", "id 3"};
  std::string expectedValues[] = {"value 1", "value 2", "value 3"};

  auto numRecords = sizeof(expectedIds) / sizeof(expectedIds[0]);
  for(std::size_t i = 0; i < numRecords; ++i)
  {
    document.add(std::make_unique<TestRecord<std::string>>(expectedIds[i],
                                                           TEST_RECORD_TYPE,
                                                           expectedValues[i]));
  }

  auto asNode = document.toNode();

  auto record_nodes = asNode[EXPECTED_RECORDS_KEY];
  ASSERT_EQ(numRecords, record_nodes.number_of_children());
  for(auto i = 0; i < record_nodes.number_of_children(); ++i)
  {
    auto &actualNode = record_nodes[i];
    EXPECT_EQ(expectedIds[i], actualNode["id"].as_string());
    EXPECT_EQ(TEST_RECORD_TYPE, actualNode["type"].as_string());
    EXPECT_EQ(expectedValues[i], actualNode[TEST_RECORD_VALUE_KEY].as_string());
  }
}

TEST(Document, toNode_relationships)
{
  Document document;
  std::string expectedSubjects[] = {"subject 1", "subject 2"};
  std::string expectedObjects[] = {"object 1", "object 2"};
  std::string expectedPredicates[] = {"predicate 1", "predicate 2"};

  auto numRecords = sizeof(expectedSubjects) / sizeof(expectedSubjects[0]);
  for(unsigned long i = 0; i < numRecords; ++i)
  {
    document.add(Relationship {
      ID {expectedSubjects[i], IDType::Global},
      expectedPredicates[i],
      ID {expectedObjects[i], IDType::Global},
    });
  }

  auto asNode = document.toNode();

  auto relationship_nodes = asNode[EXPECTED_RELATIONSHIPS_KEY];
  ASSERT_EQ(numRecords, relationship_nodes.number_of_children());
  for(auto i = 0; i < relationship_nodes.number_of_children(); ++i)
  {
    auto &actualRelationship = relationship_nodes[i];
    EXPECT_EQ(expectedSubjects[i], actualRelationship["subject"].as_string());
    EXPECT_EQ(expectedObjects[i], actualRelationship["object"].as_string());
    EXPECT_EQ(expectedPredicates[i], actualRelationship["predicate"].as_string());
  }
}

/**
 * Instances of this class acquire a temporary file name when created
 * and delete the file when destructed.
 *
 * NOTE: This class uses unsafe methods and should only be used for testing
 * purposes. DO NOT move it to the main code!!!!
 */
class NamedTempFile
{
public:
  NamedTempFile();

  // As a resource-holding class, we don't want this to be copyable
  // (or movable since there is no reason to return it from a function)
  NamedTempFile(NamedTempFile const &) = delete;

  NamedTempFile(NamedTempFile &&) = delete;

  NamedTempFile &operator=(NamedTempFile const &) = delete;

  NamedTempFile &operator=(NamedTempFile &&) = delete;

  ~NamedTempFile();

  std::string const &getName() const { return fileName; }

private:
  std::string fileName;
};

NamedTempFile::NamedTempFile()
{
  std::vector<char> tmpFileName;
  tmpFileName.resize(L_tmpnam);
  // tmpnam is not the best way to do this, but it is standard and this is
  // only a test.
  if(!std::tmpnam(tmpFileName.data()))
  {
    throw std::ios::failure {"Could not get temporary file"};
  }
  fileName = tmpFileName.data();
}

NamedTempFile::~NamedTempFile() { std::remove(fileName.data()); }


TEST(Document, create_fromJson_roundtrip_json)
{
  std::string orig_json =
    "{\"records\": [{\"type\": \"test_rec\",\"id\": "
    "\"test\"}],\"relationships\": []}"; 
  axom::sina::Document myDocument =
    Document(orig_json, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(0, myDocument.getRelationships().size());
  ASSERT_EQ(1, myDocument.getRecords().size());
  EXPECT_EQ("test_rec", myDocument.getRecords()[0]->getType());
  std::string returned_json1 = myDocument.toJson(0, 0, "", "");
  EXPECT_EQ(orig_json, returned_json1);
}

TEST(Document, create_fromJson_roundtrip_hdf5)
{
  std::string orig_json =
    "{\"records\": [{\"type\": \"test_rec\",\"id\": "
    "\"test\"}],\"relationships\": []}"; 
  axom::sina::Document myDocument =
    Document(orig_json, createRecordLoaderWithAllKnownTypes());
  saveDocument(myDocument, "round_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("round_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(0, loadedDocument.getRelationships().size());
  ASSERT_EQ(1, loadedDocument.getRecords().size());
  EXPECT_EQ("test_rec", loadedDocument.getRecords()[0]->getType());
  std::string returned_json2 = loadedDocument.toJson(0, 0, "", "");
  EXPECT_EQ(orig_json, returned_json2);
}

TEST(Document, create_fromJson_full_json)
{
  std::string long_json =
    "{\"records\": [{\"type\": \"foo\",\"id\": "
    "\"test_1\",\"user_defined\":{\"name\":\"bob\"},\"files\":{\"foo/"
    "bar.png\":{\"mimetype\":\"image\"}},\"data\":{\"scalar\": {\"value\": "
    "500,\"units\": \"miles\"}}},{\"type\":\"bar\",\"id\": "
    "\"test_2\",\"data\": {\"scalar_list\": {\"value\": [1, 2, 3]}, "
    "\"string_list\": {\"value\": [\"a\",\"wonderful\",\"world\"], "
    "\"tags\":[\"observation\"]}}},{\"type\": "
    "\"run\",\"application\":\"sina_test\",\"id\": "
    "\"test_3\",\"data\":{\"scalar\": {\"value\": 12.3, \"units\": \"g/s\", "
    "\"tags\": [\"hi\"]}, \"scalar_list\": {\"value\": [1,2,3.0,4]}}}, "
    "{\"type\": \"bar\",\"id\": \"test_4\",\"data\":{\"string\": {\"value\": "
    "\"yarr\"}, \"string_list\": {\"value\": [\"y\",\"a\",\"r\"]}}, "
    "\"files\":{\"test/test.png\":{}}, "
    "\"user_defined\":{\"hello\":\"there\"}}],\"relationships\": "
    "[{\"predicate\": \"completes\",\"subject\": \"test_2\",\"object\": "
    "\"test_1\"},{\"subject\": \"test_3\", \"predicate\": \"overrides\", "
    "\"object\": \"test_4\"}]}";
  axom::sina::Document myDocument =
    Document(long_json, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(2, myDocument.getRelationships().size());
  auto &records1 = myDocument.getRecords();
  EXPECT_EQ(4, records1.size());
}

TEST(Document, create_fromJson_full_hdf5)
{
  std::string long_json =
    "{\"records\": [{\"type\": \"foo\",\"id\": "
    "\"test_1\",\"user_defined\":{\"name\":\"bob\"},\"files\":{\"foo/"
    "bar.png\":{\"mimetype\":\"image\"}},\"data\":{\"scalar\": {\"value\": "
    "500,\"units\": \"miles\"}}},{\"type\":\"bar\",\"id\": "
    "\"test_2\",\"data\": {\"scalar_list\": {\"value\": [1, 2, 3]}, "
    "\"string_list\": {\"value\": [\"a\",\"wonderful\",\"world\"], "
    "\"tags\":[\"observation\"]}}},{\"type\": "
    "\"run\",\"application\":\"sina_test\",\"id\": "
    "\"test_3\",\"data\":{\"scalar\": {\"value\": 12.3, \"units\": \"g/s\", "
    "\"tags\": [\"hi\"]}, \"scalar_list\": {\"value\": [1,2,3.0,4]}}}, "
    "{\"type\": \"bar\",\"id\": \"test_4\",\"data\":{\"string\": {\"value\": "
    "\"yarr\"}, \"string_list\": {\"value\": [\"y\",\"a\",\"r\"]}}, "
    "\"files\":{\"test/test.png\":{}}, "
    "\"user_defined\":{\"hello\":\"there\"}}],\"relationships\": "
    "[{\"predicate\": \"completes\",\"subject\": \"test_2\",\"object\": "
    "\"test_1\"},{\"subject\": \"test_3\", \"predicate\": \"overrides\", "
    "\"object\": \"test_4\"}]}";
  axom::sina::Document myDocument =
    Document(long_json, createRecordLoaderWithAllKnownTypes());
  saveDocument(myDocument, "long_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("long_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(2, loadedDocument.getRelationships().size());
  auto &records2 = loadedDocument.getRecords();
  EXPECT_EQ(4, records2.size());
}

TEST(Document, create_fromJson_value_check_json)
{
  std::string data_json =
    "{\"records\": [{\"type\": \"run\", \"application\":\"test\", \"id\": "
    "\"test_1\",\"data\":{\"int\": {\"value\": 500,\"units\": \"miles\"}, "
    "\"str/ings\": {\"value\":[\"z\", \"o\", \"o\"]}}, "
    "\"files\":{\"test/test.png\":{}}}]}";
  axom::sina::Document myDocument =
    Document(data_json, createRecordLoaderWithAllKnownTypes());
  EXPECT_EQ(0, myDocument.getRelationships().size());
  auto &records1 = myDocument.getRecords();
  EXPECT_EQ(1, records1.size());
  EXPECT_EQ(records1[0]->getType(), "run");
  auto &data1 = records1[0]->getData();
  EXPECT_EQ(data1.at("int").getScalar(), 500.0);
  std::vector<std::string> expected_string_vals = {"z", "o", "o"};
  EXPECT_EQ(data1.at("str/ings").getStringArray(), expected_string_vals);
  EXPECT_EQ(records1[0]->getFiles().count(File {"test/test.png"}), 1);
}

TEST(Document, create_fromJson_value_check_hdf5)
{
  std::string data_json =
    "{\"records\": [{\"type\": \"run\", \"application\":\"test\", \"id\": "
    "\"test_1\",\"data\":{\"int\": {\"value\": 500,\"units\": \"miles\"}, "
    "\"str/ings\": {\"value\":[\"z\", \"o\", \"o\"]}}, "
    "\"files\":{\"test/test.png\":{}}}]}";
  axom::sina::Document myDocument =
    Document(data_json, createRecordLoaderWithAllKnownTypes());
  std::vector<std::string> expected_string_vals = {"z", "o", "o"};
  saveDocument(myDocument, "data_json.hdf5", Protocol::HDF5);
  Document loadedDocument = loadDocument("data_json.hdf5", Protocol::HDF5);
  EXPECT_EQ(0, loadedDocument.getRelationships().size());
  auto &records2 = loadedDocument.getRecords();
  EXPECT_EQ(1, records2.size());
  EXPECT_EQ(records2[0]->getType(), "run");
  auto &data2 = records2[0]->getData();
  EXPECT_EQ(data2.at("int").getScalar(), 500.0);
  EXPECT_EQ(data2.at("str/ings").getStringArray(), expected_string_vals);
  EXPECT_EQ(records2[0]->getFiles().count(File {"test/test.png"}), 1);
}

TEST(Document, saveDocument_json)
{
  NamedTempFile tmpFile;

  // First, write some random stuff to the temp file to make sure it is
  // overwritten.
  {
    std::ofstream fout {tmpFile.getName()};
    fout << "Initial contents";
  }

  Document document;
  document.add(
    std::make_unique<Record>(ID {"the id", IDType::Global}, "the type"));

  saveDocument(document, tmpFile.getName() + ".json");

  conduit::Node readContents;
  {
    std::ifstream fin {tmpFile.getName() + ".json"};
    std::stringstream f_buf;
    f_buf << fin.rdbuf();
    readContents.parse(f_buf.str(), "json");
  }

  ASSERT_TRUE(readContents[EXPECTED_RECORDS_KEY].dtype().is_list());
  EXPECT_EQ(1, readContents[EXPECTED_RECORDS_KEY].number_of_children());
  auto &readRecord = readContents[EXPECTED_RECORDS_KEY][0];
  EXPECT_EQ("the id", readRecord["id"].as_string());
  EXPECT_EQ("the type", readRecord["type"].as_string());
}

TEST(Document, saveDocument_hdf5)
{
    NamedTempFile tmpFile;

    // First, write some random stuff to the temp file to make sure it is
    // overwritten.
    {
        std::ofstream fout {tmpFile.getName()};
        fout << "Initial contents";
    }

    Document document;
    document.add(
        std::make_unique<Record>(ID {"the id", IDType::Global}, "the type"));

    saveDocument(document, tmpFile.getName() + ".hdf5", Protocol::HDF5);

    conduit::Node readContents;
    conduit::relay::io::load(tmpFile.getName() + ".hdf5", "hdf5", readContents);

    ASSERT_TRUE(readContents[EXPECTED_RECORDS_KEY].dtype().is_list());
    EXPECT_EQ(1, readContents[EXPECTED_RECORDS_KEY].number_of_children());
    auto &readRecord = readContents[EXPECTED_RECORDS_KEY][0];
    EXPECT_EQ("the id", readRecord["id"].as_string());
    EXPECT_EQ("the type", readRecord["type"].as_string());
}

TEST(Document, load_specifiedRecordLoader)
{
  using RecordType = TestRecord<int>;
  auto originalRecord = std::make_unique<RecordType>("the ID", "my type", 123);
  Document originalDocument;
  originalDocument.add(std::move(originalRecord));

  NamedTempFile file;
  {
    std::ofstream fout {file.getName()};
    fout << originalDocument.toNode().to_json();
  }

  RecordLoader loader;
  loader.addTypeLoader("my type", [](conduit::Node const &asNode) {
    return std::make_unique<RecordType>(
      getRequiredString("id", asNode, "Test type"),
      getRequiredString("type", asNode, "Test type"),
      static_cast<int>(
        getRequiredField(TEST_RECORD_VALUE_KEY, asNode, "Test type").as_int64()));
  });
  Document loadedDocument = loadDocument(file.getName(), loader);
  ASSERT_EQ(1u, loadedDocument.getRecords().size());
  auto loadedRecord =
    dynamic_cast<RecordType const *>(loadedDocument.getRecords()[0].get());
  ASSERT_NE(nullptr, loadedRecord);
  EXPECT_EQ(123, loadedRecord->getValue());
}

TEST(Document, load_defaultRecordLoaders)
{
  auto originalRun =
    std::make_unique<axom::sina::Run>(ID {"the ID", IDType::Global},
                                      "the app",
                                      "1.2.3",
                                      "jdoe");
  Document originalDocument;
  originalDocument.add(std::move(originalRun));

  NamedTempFile file;
  {
    std::ofstream fout {file.getName()};
    fout << originalDocument.toNode().to_json();
  }

  Document loadedDocument = loadDocument(file.getName());
  ASSERT_EQ(1u, loadedDocument.getRecords().size());
  auto loadedRun =
    dynamic_cast<axom::sina::Run const *>(loadedDocument.getRecords()[0].get());
  EXPECT_NE(nullptr, loadedRun);
}

TEST(Document, test_append_to_json) {
    std::string jsonFilePath = "test.json";

    std::ofstream testFile(jsonFilePath);
    testFile << R"({
        "records": [
            {
                "data": {"scal1": {"value": "hello!"}},
                "type": "foo1",
                "id": "bar1",
                "user_defined":{"hello": "and"},
                "curve_sets": {
                    "1": {
                        "dependent": {
                            "0": { "value": [1, 2] },
                            "1": { "value": [10, 20] }
                        },
                        "independent": {
                            "0": { "value": [4, 5] },
                            "1": { "value": [7, 8] }
                        }
                    }
                }
            },
            {
                "type": "foo2",
                "id": "bar2",
                "curve_sets": {
                    "1": {
                        "dependent": {
                            "0": { "value": [1, 2] },
                            "1": { "value": [10, 20] }
                        },
                        "independent": {
                            "0": { "value": [4, 5] },
                            "1": { "value": [7, 8] }
                        }
                    }
                }
            }
        ]
    })";
    testFile.close();

    std::string data = "{\"records\": ["
    "   {"
    "       \"type\": \"foo1\","
    "       \"id\": \"bar1\","
    "       \"data\": {\"scal1\": {\"value\": \"goodbye!\"}, \"scal2\": {\"value\": \"goodbye!\"}},"
    "       \"user_defined\":{\"hello\": \"goodbye\"},"
    "       \"curve_sets\": {"
    "           \"1\": {"
    "               \"dependent\": {"
    "                   \"0\": {\"value\": [1]},"
    "                   \"1\": {\"value\": [10]}"
    "               },"
    "               \"independent\": {"
    "                   \"0\": {\"value\": [4]},"
    "                   \"1\": {\"value\": [7]}"
    "               }"
    "           }"
    "       }"
    "   },"
    "   {"
    "       \"type\": \"foo2\","
    "       \"id\": \"bar3\","
    "       \"curve_sets\": {"
    "           \"1\": {"
    "               \"dependent\": {"
    "                   \"0\": {\"value\": [1]},"
    "                   \"1\": {\"value\": [10]}"
    "               },"
    "               \"independent\": {"
    "                   \"0\": {\"value\": [4]},"
    "                   \"1\": {\"value\": [7]}"
    "               }"
    "           }"
    "       }"
    "   }"
    "]}";

    axom::sina::Document newData =
        Document(data, createRecordLoaderWithAllKnownTypes());

    bool result = append_to_json(jsonFilePath, newData, 1, 1);

    std::ifstream updatedFile(jsonFilePath);
    nlohmann::json j;
    updatedFile >> j;

    ASSERT_TRUE(result); 

    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["0"]["value"],
    //           nlohmann::json({1, 2, 3, 1})); 
    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["1"]["value"],
    //           nlohmann::json({10, 20, 10})); 
    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["0"]["value"],
    //           nlohmann::json({4, 5, 6, 4})); 
    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["1"]["value"],
    //           nlohmann::json({7, 8, 7})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["dependent"]["0"]["value"],
    //           nlohmann::json({1, 2, 3, 1})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["dependent"]["1"]["value"],
    //           nlohmann::json({10, 20, 10})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["independent"]["0"]["value"],
    //           nlohmann::json({4, 5, 6, 4})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["independent"]["1"]["value"],
    //           nlohmann::json({7, 8, 7}));   

    // std::remove(jsonFilePath.c_str());
}


TEST(Document, test_append_to_hdf5) {
    std::string hdf5FilePath = "test.hdf5";
    conduit::Node initialData;
    initialData["records/0/data/scal1/value"].set("hello!");
    initialData["records/0/type"].set("foo1");
    initialData["records/0/id"].set("bar1");
    initialData["records/0/user_defined/hello"].set("and");
    initialData["records/0/curve_sets/1/dependent/0/value"].set(std::vector<double>{1.0, 2.0});
    initialData["records/0/curve_sets/1/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
    initialData["records/0/curve_sets/1/independent/0/value"].set(std::vector<double>{4.0, 5.0});
    initialData["records/0/curve_sets/1/independent/1/value"].set(std::vector<double>{7.0, 8.0});
    initialData["records/1/type"].set("foo2");
    initialData["records/1/id"].set("bar2");
    initialData["records/1/curve_sets/1/dependent/0/value"].set(std::vector<double>{1.0, 2.0});
    initialData["records/1/curve_sets/1/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
    initialData["records/1/curve_sets/1/independent/0/value"].set(std::vector<double>{4.0, 5.0});
    initialData["records/1/curve_sets/1/independent/1/value"].set(std::vector<double>{7.0, 8.0});
    conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

    std::string data = "{\"records\": ["
    "   {"
    "       \"type\": \"foo1\","
    "       \"id\": \"bar1\","
    "       \"data\": {\"scal1\": {\"value\": \"goodbye!\"}, \"scal2\": {\"value\": \"goodbye!\"}},"
    "       \"user_defined\":{\"hello\": \"goodbye\"},"
    "       \"curve_sets\": {"
    "           \"1\": {"
    "               \"dependent\": {"
    "                   \"0\": {\"value\": [1]},"
    "                   \"1\": {\"value\": [10]}"
    "               },"
    "               \"independent\": {"
    "                   \"0\": {\"value\": [4]},"
    "                   \"1\": {\"value\": [7]}"
    "               }"
    "           }"
    "       }"
    "   },"
    "   {"
    "       \"type\": \"foo2\","
    "       \"id\": \"bar3\","
    "       \"curve_sets\": {"
    "           \"1\": {"
    "               \"dependent\": {"
    "                   \"0\": {\"value\": [1]},"
    "                   \"1\": {\"value\": [10]}"
    "               },"
    "               \"independent\": {"
    "                   \"0\": {\"value\": [4]},"
    "                   \"1\": {\"value\": [7]}"
    "               }"
    "           }"
    "       }"
    "   }"
    "]}";

    axom::sina::Document newData = Document(data, createRecordLoaderWithAllKnownTypes());
    bool result = append_to_hdf5(hdf5FilePath, newData, 1, 1);


    // std::ofstream testFile(hdf5FilePath);
    // testFile << R"({
    //     "records": [
    //         {
    //             "data": {"scal1": {"value": "hello!"}},
    //             "type": "foo1",
    //             "id": "bar1",
    //             "user_defined":{"hello": "and"},
    //             "curve_sets": {
    //                 "1": {
    //                     "dependent": {
    //                         "0": { "value": [1, 2] },
    //                         "1": { "value": [10, 20] }
    //                     },
    //                     "independent": {
    //                         "0": { "value": [4, 5] },
    //                         "1": { "value": [7, 8] }
    //                     }
    //                 }
    //             }
    //         },
    //         {
    //             "type": "foo2",
    //             "id": "bar2",
    //             "curve_sets": {
    //                 "1": {
    //                     "dependent": {
    //                         "0": { "value": [1, 2] },
    //                         "1": { "value": [10, 20] }
    //                     },
    //                     "independent": {
    //                         "0": { "value": [4, 5] },
    //                         "1": { "value": [7, 8] }
    //                     }
    //                 }
    //             }
    //         }
    //     ]
    // })";
    // testFile.close();

    // std::string data = "{\"records\": ["
    // "   {"
    // "       \"type\": \"foo1\","
    // "       \"id\": \"bar1\","
    // "       \"data\": {\"scal1\": {\"value\": \"goodbye!\"}, \"scal2\": {\"value\": \"goodbye!\"}},"
    // "       \"user_defined\":{\"hello\": \"goodbye\"},"
    // "       \"curve_sets\": {"
    // "           \"1\": {"
    // "               \"dependent\": {"
    // "                   \"0\": {\"value\": [1]},"
    // "                   \"1\": {\"value\": [10]}"
    // "               },"
    // "               \"independent\": {"
    // "                   \"0\": {\"value\": [4]},"
    // "                   \"1\": {\"value\": [7]}"
    // "               }"
    // "           }"
    // "       }"
    // "   },"
    // "   {"
    // "       \"type\": \"foo2\","
    // "       \"id\": \"bar3\","
    // "       \"curve_sets\": {"
    // "           \"1\": {"
    // "               \"dependent\": {"
    // "                   \"0\": {\"value\": [1]},"
    // "                   \"1\": {\"value\": [10]}"
    // "               },"
    // "               \"independent\": {"
    // "                   \"0\": {\"value\": [4]},"
    // "                   \"1\": {\"value\": [7]}"
    // "               }"
    // "           }"
    // "       }"
    // "   }"
    // "]}";

    // axom::sina::Document newData =
    //     Document(data, createRecordLoaderWithAllKnownTypes());

    // bool result = append_to_hdf5(hdf5FilePath, newData);

    ASSERT_TRUE(result); 

    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["0"]["value"],
    //           nlohmann::json({1, 2, 3, 1})); 
    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["1"]["value"],
    //           nlohmann::json({10, 20, 10})); 
    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["0"]["value"],
    //           nlohmann::json({4, 5, 6, 4})); 
    // ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["1"]["value"],
    //           nlohmann::json({7, 8, 7})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["dependent"]["0"]["value"],
    //           nlohmann::json({1, 2, 3, 1})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["dependent"]["1"]["value"],
    //           nlohmann::json({10, 20, 10})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["independent"]["0"]["value"],
    //           nlohmann::json({4, 5, 6, 4})); 
    // ASSERT_EQ(j["records"][1]["curve_sets"]["1"]["independent"]["1"]["value"],
    //           nlohmann::json({7, 8, 7}));   

    // std::remove(jsonFilePath.c_str());
}

// TEST(Document, test_append_to_json_multiple_records_invalid_error) {
//     std::string jsonFilePath = "test.json";

//     std::ofstream testFile(jsonFilePath);
//     testFile << R"({
//         "records": [
//             {
//                 "curve_sets": {
//                     "1": {
//                         "dependent": {
//                             "0": { "value": [1, 2, 3] },
//                             "1": { "value": [10, 20] }
//                         },
//                         "independent": {
//                             "0": { "value": [4, 5, 6] },
//                             "1": { "value": [7, 8] }
//                         }
//                     }
//                 }
//             },
//             {
//                 "curve_sets": {
//                     "2": {
//                         "dependent": {
//                             "0": { "value": [1, 2, 3] },
//                             "1": { "value": [10, 20] }
//                         },
//                         "independent": {
//                             "0": { "value": [4, 5, 6] },
//                             "1": { "value": [7, 8] }
//                         }
//                     }
//                 }
//             }
//         ]
//     })";
//     testFile.close();

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {"1", {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }}
//                 }}
//             },
//             {
//                 {"curve_sets", {
//                     {"2", {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"2", {{"value", {7}}}}
//                         }}
//                     }}
//                 }}
//             }
//         }}
//     };
//     bool result = append_to_json_multiple_records(jsonFilePath, newData);
//     std::ifstream updatedFile(jsonFilePath);
//     nlohmann::json j;
//     updatedFile >> j;

//     ASSERT_FALSE(result); 

//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["0"]["value"],
//               nlohmann::json({1, 2, 3})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["1"]["value"],
//               nlohmann::json({10, 20})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["0"]["value"],
//               nlohmann::json({4, 5, 6})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["1"]["value"],
//               nlohmann::json({7, 8})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["dependent"]["0"]["value"],
//               nlohmann::json({1, 2, 3})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["dependent"]["1"]["value"],
//               nlohmann::json({10, 20})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["independent"]["0"]["value"],
//               nlohmann::json({4, 5, 6})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["independent"]["1"]["value"],
//               nlohmann::json({7, 8}));     

//     std::remove(jsonFilePath.c_str());
// }

// TEST(Document, test_append_to_json_multiple_records_mismatch_error) {
//     std::string jsonFilePath = "test.json";

//     std::ofstream testFile(jsonFilePath);
//     testFile << R"({
//         "records": [
//             {
//                 "curve_sets": {
//                     "1": {
//                         "dependent": {
//                             "0": { "value": [1, 2, 3] },
//                             "1": { "value": [10, 20] }
//                         },
//                         "independent": {
//                             "0": { "value": [4, 5, 6] },
//                             "1": { "value": [7, 8] }
//                         }
//                     }
//                 }
//             },
//             {
//                 "curve_sets": {
//                     "2": {
//                         "dependent": {
//                             "0": { "value": [1, 2, 3] },
//                             "1": { "value": [10, 20] }
//                         },
//                         "independent": {
//                             "0": { "value": [4, 5, 6] },
//                             "1": { "value": [7, 8] }
//                         }
//                     }
//                 }
//             }
//         ]
//     })";
//     testFile.close();

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {"1", {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }}
//                 }}
//             }
//         }}
//     };
//     bool result = append_to_json_multiple_records(jsonFilePath, newData);

//     std::ifstream updatedFile(jsonFilePath);
//     nlohmann::json j;
//     updatedFile >> j;

//     ASSERT_FALSE(result); 

//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["0"]["value"],
//               nlohmann::json({1, 2, 3})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["1"]["value"],
//               nlohmann::json({10, 20})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["0"]["value"],
//               nlohmann::json({4, 5, 6})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["1"]["value"],
//               nlohmann::json({7, 8})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["dependent"]["0"]["value"],
//               nlohmann::json({1, 2, 3})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["dependent"]["1"]["value"],
//               nlohmann::json({10, 20})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["independent"]["0"]["value"],
//               nlohmann::json({4, 5, 6})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["independent"]["1"]["value"],
//               nlohmann::json({7, 8}));  

//     std::remove(jsonFilePath.c_str());
// }

// TEST(Document, test_append_to_json_one_record) {
//     std::string jsonFilePath = "test.json";

//     std::ofstream testFile(jsonFilePath);
//     testFile << R"({
//         "records": [
//             {
//                 "curve_sets": {
//                     "1": {
//                         "dependent": {
//                             "0": { "value": [1, 2, 3] },
//                             "1": { "value": [10, 20] }
//                         },
//                         "independent": {
//                             "0": { "value": [4, 5, 6] },
//                             "1": { "value": [7, 8] }
//                         }
//                     }
//                 }
//             },
//             {
//                 "curve_sets": {
//                     "2": {
//                         "dependent": {
//                             "0": { "value": [1, 2, 3] },
//                             "1": { "value": [10, 20] }
//                         },
//                         "independent": {
//                             "0": { "value": [4, 5, 6] },
//                             "1": { "value": [7, 8] }
//                         }
//                     }
//                 }
//             }
//         ]
//     })";
//     testFile.close();

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {"2", {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }}
//                 }}
//             }
//         }}
//     };
//     bool result = append_to_json_one_record(jsonFilePath, newData, 1);

//     std::ifstream updatedFile(jsonFilePath);
//     nlohmann::json j;
//     updatedFile >> j;

//     ASSERT_TRUE(result); 

//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["0"]["value"],
//               nlohmann::json({1, 2, 3})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["dependent"]["1"]["value"],
//               nlohmann::json({10, 20})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["0"]["value"],
//               nlohmann::json({4, 5, 6})); 
//     ASSERT_EQ(j["records"][0]["curve_sets"]["1"]["independent"]["1"]["value"],
//               nlohmann::json({7, 8})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["dependent"]["0"]["value"],
//               nlohmann::json({1, 2, 3, 1})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["dependent"]["1"]["value"],
//               nlohmann::json({10, 20, 10})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["independent"]["0"]["value"],
//               nlohmann::json({4, 5, 6, 4})); 
//     ASSERT_EQ(j["records"][1]["curve_sets"]["2"]["independent"]["1"]["value"],
//               nlohmann::json({7, 8, 7}));   
//     std::remove(jsonFilePath.c_str());
// }

// TEST(Document, test_append_to_hdf5_multiple_records) {
//     std::string hdf5FilePath = "test_append.hdf5";
//     conduit::Node initialData;
//     initialData["records/0/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/0/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/0/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/0/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});
//     initialData["records/1/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/1/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/1/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/1/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }
//                 }}
//             },
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }
//                 }}
//             }
//         }}
//     };

//     conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

//     bool result = append_to_hdf5_multiple(hdf5FilePath, newData);

//     conduit::Node updatedData;
//     conduit::relay::io::hdf5_read(hdf5FilePath, updatedData);

//     ASSERT_TRUE(result);

//     std::vector<double> expectedDependent00 = {1, 2, 3, 1};
//     std::vector<double> actualDependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent00, expectedDependent00);

//     std::vector<double> expectedDependent01 = {10, 20, 10};
//     std::vector<double> actualDependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent01, expectedDependent01);

//     std::vector<double> expectedIndependent00 = {4, 5, 6, 4};
//     std::vector<double> actualIndependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent00, expectedIndependent00);

//     std::vector<double> expectedIndependent01 = {7, 8, 7};
//     std::vector<double> actualIndependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent01, expectedIndependent01);

//     std::vector<double> expectedDependent10 = {1, 2, 3, 1};
//     std::vector<double> actualDependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent10, expectedDependent10);

//     std::vector<double> expectedDependent11 = {10, 20, 10};
//     std::vector<double> actualDependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent11, expectedDependent11);

//     std::vector<double> expectedIndependent10 = {4, 5, 6, 4};
//     std::vector<double> actualIndependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent10, expectedIndependent10);

//     std::vector<double> expectedIndependent11 = {7, 8, 7};
//     std::vector<double> actualIndependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent11, expectedIndependent11);

//     std::remove(hdf5FilePath.c_str());
// }

// TEST(Document, test_append_to_hdf5_invalid_error) {
//     std::string hdf5FilePath = "test_append.hdf5";
//     conduit::Node initialData;
//     initialData["records/0/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/0/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/0/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/0/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});
//     initialData["records/1/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/1/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/1/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/1/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }
//                 }}
//             },
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"2", {{"value", {7}}}}
//                         }}
//                     }
//                 }}
//             }
//         }}
//     };

//     conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

//     bool result = append_to_hdf5_multiple(hdf5FilePath, newData);

//     conduit::Node updatedData;
//     conduit::relay::io::hdf5_read(hdf5FilePath, updatedData);

//     ASSERT_FALSE(result);

//     std::vector<double> expectedDependent00 = {1, 2, 3};
//     std::vector<double> actualDependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent00, expectedDependent00);

//     std::vector<double> expectedDependent01 = {10, 20};
//     std::vector<double> actualDependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent01, expectedDependent01);

//     std::vector<double> expectedIndependent00 = {4, 5, 6};
//     std::vector<double> actualIndependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent00, expectedIndependent00);

//     std::vector<double> expectedIndependent01 = {7, 8};
//     std::vector<double> actualIndependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent01, expectedIndependent01);

//     std::vector<double> expectedDependent10 = {1, 2, 3};
//     std::vector<double> actualDependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent10, expectedDependent10);

//     std::vector<double> expectedDependent11 = {10, 20};
//     std::vector<double> actualDependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent11, expectedDependent11);

//     std::vector<double> expectedIndependent10 = {4, 5, 6};
//     std::vector<double> actualIndependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent10, expectedIndependent10);

//     std::vector<double> expectedIndependent11 = {7, 8};
//     std::vector<double> actualIndependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent11, expectedIndependent11);

//     std::remove(hdf5FilePath.c_str());
// }

// TEST(Document, test_append_to_hdf5_mismatch_error) {
//     std::string hdf5FilePath = "test_append.hdf5";
//     conduit::Node initialData;
//     initialData["records/0/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/0/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/0/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/0/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});
//     initialData["records/1/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/1/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/1/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/1/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                         {"independent", {
//                             {"0", {{"value", {4}}}},
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }
//                 }}
//             }
//         }}
//     };

//     conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

//     bool result = append_to_hdf5_multiple(hdf5FilePath, newData);
//     ASSERT_FALSE(result);

//     conduit::Node updatedData;
//     conduit::relay::io::hdf5_read(hdf5FilePath, updatedData);



//     std::vector<double> expectedDependent00 = {1, 2, 3};
//     std::vector<double> actualDependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent00, expectedDependent00);

//     std::vector<double> expectedDependent01 = {10, 20};
//     std::vector<double> actualDependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent01, expectedDependent01);

//     std::vector<double> expectedIndependent00 = {4, 5, 6};
//     std::vector<double> actualIndependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent00, expectedIndependent00);

//     std::vector<double> expectedIndependent01 = {7, 8};
//     std::vector<double> actualIndependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent01, expectedIndependent01);

//     std::vector<double> expectedDependent10 = {1, 2, 3};
//     std::vector<double> actualDependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent10, expectedDependent10);

//     std::vector<double> expectedDependent11 = {10, 20};
//     std::vector<double> actualDependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent11, expectedDependent11);

//     std::vector<double> expectedIndependent10 = {4, 5, 6};
//     std::vector<double> actualIndependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent10, expectedIndependent10);

//     std::vector<double> expectedIndependent11 = {7, 8};
//     std::vector<double> actualIndependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent11, expectedIndependent11);

//     std::remove(hdf5FilePath.c_str());
// }


// TEST(Document, test_append_to_hdf5_complex) {
//     std::string hdf5FilePath = "test_append.hdf5";
//     conduit::Node initialData;
//     initialData["records/0/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/0/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/0/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/0/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});
//     initialData["records/1/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/1/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/1/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/1/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});

//     nlohmann::json newData = {
//         {"records", {
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                         }},
//                         {"independent", {
//                             {"1", {{"value", {7}}}}
//                         }}
//                     }
//                 }}
//             },
//             {
//                 {"curve_sets", {
//                     {
//                         {"dependent", {
//                             {"0", {{"value", {1}}}},
//                             {"1", {{"value", {10}}}}
//                         }},
//                     }
//                 }}
//             }
//         }}
//     };

//     conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

//     bool result = append_to_hdf5_multiple(hdf5FilePath, newData);

//     conduit::Node updatedData;
//     conduit::relay::io::hdf5_read(hdf5FilePath, updatedData);

//     ASSERT_TRUE(result);

//     std::vector<double> expectedDependent00 = {1, 2, 3, 1};
//     std::vector<double> actualDependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent00, expectedDependent00);

//     std::vector<double> expectedDependent01 = {10, 20};
//     std::vector<double> actualDependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent01, expectedDependent01);

//     std::vector<double> expectedIndependent00 = {4, 5, 6};
//     std::vector<double> actualIndependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent00, expectedIndependent00);

//     std::vector<double> expectedIndependent01 = {7, 8, 7};
//     std::vector<double> actualIndependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent01, expectedIndependent01);

//     std::vector<double> expectedDependent10 = {1, 2, 3, 1};
//     std::vector<double> actualDependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent10, expectedDependent10);

//     std::vector<double> expectedDependent11 = {10, 20, 10};
//     std::vector<double> actualDependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent11, expectedDependent11);

//     std::vector<double> expectedIndependent10 = {4, 5, 6};
//     std::vector<double> actualIndependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent10, expectedIndependent10);

//     std::vector<double> expectedIndependent11 = {7, 8};
//     std::vector<double> actualIndependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent11, expectedIndependent11);

//     std::remove(hdf5FilePath.c_str());
// }

// TEST(Document, test_append_to_hdf5_one_record) {
//     std::string hdf5FilePath = "test_append.hdf5";

//     conduit::Node initialData;
//     initialData["records/0/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/0/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/0/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/0/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});
//     initialData["records/1/curve_sets/0/dependent/0/value"].set(std::vector<double>{1.0, 2.0, 3.0});
//     initialData["records/1/curve_sets/0/dependent/1/value"].set(std::vector<double>{10.0, 20.0});
//     initialData["records/1/curve_sets/0/independent/0/value"].set(std::vector<double>{4.0, 5.0, 6.0});
//     initialData["records/1/curve_sets/0/independent/1/value"].set(std::vector<double>{7.0, 8.0});

//     conduit::relay::io::hdf5_write(initialData, hdf5FilePath);

//     nlohmann::json newData = {{{"dependent", {{"0", {{"value", {1}}}}, {"1", {{"value", {10}}}}}},{"independent", {{"0", {{"value", {4}}}}, {"1", {{"value", {7}}}}}}}};

//     bool result = append_to_hdf5_one_record(hdf5FilePath, newData, 1);

//     conduit::Node updatedData;
//     conduit::relay::io::hdf5_read(hdf5FilePath, updatedData);

//     ASSERT_TRUE(result);

//     std::vector<double> expectedDependent00 = {1, 2, 3};
//     std::vector<double> actualDependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent00, expectedDependent00);

//     std::vector<double> expectedDependent01 = {10, 20};
//     std::vector<double> actualDependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent01, expectedDependent01);

//     std::vector<double> expectedIndependent00 = {4, 5, 6};
//     std::vector<double> actualIndependent00 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent00, expectedIndependent00);

//     std::vector<double> expectedIndependent01 = {7, 8};
//     std::vector<double> actualIndependent01 = node_to_double_vector(updatedData["records/0/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent01, expectedIndependent01);

//     std::vector<double> expectedDependent10 = {1, 2, 3, 1};
//     std::vector<double> actualDependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/0/value"]);
//     ASSERT_EQ(actualDependent10, expectedDependent10);

//     std::vector<double> expectedDependent11 = {10, 20,10};
//     std::vector<double> actualDependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/dependent/1/value"]);
//     ASSERT_EQ(actualDependent11, expectedDependent11);

//     std::vector<double> expectedIndependent10 = {4, 5, 6,4};
//     std::vector<double> actualIndependent10 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/0/value"]);
//     ASSERT_EQ(actualIndependent10, expectedIndependent10);

//     std::vector<double> expectedIndependent11 = {7, 8,7};
//     std::vector<double> actualIndependent11 = node_to_double_vector(updatedData["records/1/curve_sets/0/independent/1/value"]);
//     ASSERT_EQ(actualIndependent11, expectedIndependent11);

//     std::remove(hdf5FilePath.c_str());
// }

}  // namespace
}  // namespace testing
}  // namespace sina
}  // namespace axom
