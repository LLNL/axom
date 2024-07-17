// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/sina.hpp"

int main(void)
{
  // Create a new document
  axom::sina::Document document;

  // Create a record of this specific study
  // This study has an ID of "study1", which has to be unique to this file
  axom::sina::ID studyID {"study1", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> study {
    new axom::sina::Record {studyID, "UQ study"}};

  // Create a run of "My Sim Code" version "1.2.3", which was run by "jdoe".
  // The run has an ID of "run1", which has to be unique to this file.
  axom::sina::ID runID {"run1", axom::sina::IDType::Local};
  std::unique_ptr<axom::sina::Record> run {
    new axom::sina::Run {runID, "My Sim Code", "1.2.3", "jdoe"}};

  // Create a relationship between the study and the run
  // Here we're saying that the study contains the run
  axom::sina::Relationship relationship {studyID, "contains", runID};

  // Add the run, study record, and relationship to the document
  document.add(std::move(run));
  document.add(std::move(study));
  document.add(relationship);

  // Save the document directly to a file.
  axom::sina::saveDocument(document, "MySinaData.json");
}