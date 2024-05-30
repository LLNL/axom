#include "axom/sina.hpp"

#include <utility>
#include <memory>

namespace {

//! [create record]
void createRecord() {
    axom::sina::ID id{"some_record_id", axom::sina::IDType::Local};
    std::unique_ptr<axom::sina::Record> record{new axom::sina::Record{id, "my_record_type"}};

    // Add the record to a document
    axom::sina::Document doc;
    doc.add(std::move(record));
}
//! [create record]

//! [create run]
void createRun() {
    axom::sina::ID id{"some_run_id", axom::sina::IDType::Local};
    std::unique_ptr<axom::sina::Record> run{new axom::sina::Run{id, "My Sim Code", "1.2.3", "jdoe"}};

    // Add the record to a document
    axom::sina::Document doc;
    doc.add(std::move(run));
}
//! [create run]

//! [adding data]
void addData(axom::sina::Record &record) {
    // Add a scalar named "my_scalar" with the value 123.456
    record.add("my_scalar", axom::sina::Datum{123.456});

    // Add a string named "my_string" with the value "abc"
    record.add("my_string", axom::sina::Datum{"abc"});

    // Add a list of scalars named "my_scalar_list"
    std::vector<double> scalarList = {1.2, -3.4, 5.6};
    record.add("my_scalar_list", axom::sina::Datum{scalarList});

    // Add a list of strings named "my_string_list"
    std::vector<std::string> stringList = {"hi", "hello", "howdy"};
    record.add("my_string_list", axom::sina::Datum{stringList});
}
//! [adding data]

//! [curve sets]
void addCurveSets(axom::sina::Record &record) {
    axom::sina::CurveSet timePlots{"time_plots"};

    // Add the independent variable
    timePlots.addIndependentCurve(
        axom::sina::Curve{"time", {0.0, 0.1, 0.25, 0.3}});

    // Add some dependent variables.
    // The length of each must be the same as the length of the independent.
    timePlots.addDependentCurve(
        axom::sina::Curve{"temperature", {300.0, 310.0, 350.0, 400.0}});

    timePlots.addDependentCurve(
        axom::sina::Curve{"energy", {0.0, 10.0, 20.0, 30.0}});

    // Associate the curve sets with the record
    record.add(timePlots);
}

//! [curve sets]


//! [file add_and_remove]
void addAndRemoveFileToRecord(axom::sina::Record &run) {
    axom::sina::File my_file{"some/path.txt"};
    // Adds the file to the record's file list
    run.add(my_file);
    // Removes the file from the record's file list
    run.remove(my_file);
}

//! [file add_and_remove]

//! [relationships]
void associateRunToStudy(axom::sina::Document &doc, axom::sina::Record const &uqStudy, axom::sina::Record const &run) {
    doc.add(axom::sina::Relationship{uqStudy.getId(), "contains", run.getId()});
}
//! [relationships]


//! [library data foo]
void foo_collectData(axom::sina::DataHolder &fooData) {
    fooData.add("temperature", axom::sina::Datum{500});
    fooData.add("energy", axom::sina::Datum{1.2e10});
}
//! [library data foo]

//! [library data bar]
void bar_gatherData(axom::sina::DataHolder &barData) {
    barData.add("temperature", axom::sina::Datum{400});
    barData.add("mass", axom::sina::Datum{15});
}
//! [library data bar]

//! [library data host]
void gatherAllData(axom::sina::Record &record) {
   auto fooData = record.addLibraryData("foo");
   auto barData = record.addLibraryData("bar");

   foo_collectData(*fooData);
   bar_gatherData(*barData);

   record.add("temperature", axom::sina::Datum{450});
}
//! [library data host]

//! [io write]
void save(axom::sina::Document const &doc) {
    axom::sina::saveDocument(doc, "my_output.json");
}
//! [io write]

//! [io read]
void load() {
    axom::sina::Document doc = axom::sina::loadDocument("my_output.json");
}
//! [io read]

//! [user defined]
void addUserDefined(axom::sina::Record &record) {
    conduit::Node &userDefined = record.getUserDefinedContent();
    userDefined["var_1"] = "a";
    userDefined["var_2"] = "b";

    conduit::Node subNode;
    subNode["sub_1"] = 10;
    subNode["sub_2"] = 20;
    userDefined["sub_structure"] = subNode;
}
//! [user defined]

}

int main() {
    // Call everything to keep the compiler from complaining about unused functions
    axom::sina::Record run{axom::sina::ID{"my_record", axom::sina::IDType::Global}, "my_record_type"};
    axom::sina::Record study{axom::sina::ID{"my_run", axom::sina::IDType::Global}, "UQ study"};
    axom::sina::Document doc;
    addData(run);
    createRecord();
    createRun();
    associateRunToStudy(doc, study, run);
    gatherAllData(run);
    addCurveSets(run);
    addAndRemoveFileToRecord(run);
    addUserDefined(run);
    // TODO
    // - Add Record to doc
    // - Check output file to see if record shows up
    save(doc);
    load();
}
