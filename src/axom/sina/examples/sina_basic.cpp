#include "axom/sina.hpp"

int main (void) {
    // Create a new document
    axom::sina::Document document;
    // Create a run of "My Sim Code" version "1.2.3", which was run by "jdoe".
    // The run has an ID of "run1", which has to be unique to this file.
    axom::sina::ID runID{"run1", axom::sina::IDType::Local};
    std::unique_ptr<axom::sina::Record> run{
        new axom::sina::Run{runID, "My Sim Code", "1.2.3", "jdoe"}};
    // Add the run to the document
    document.add(std::move(run));
    // Save the document directly to a file.
    axom::sina::saveDocument(document, "MySinaData.json");
}