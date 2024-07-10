#include "axom/sina.hpp"

int main(void) {
    // Create 2 different files
    axom::sina::File myFile{"/path/to/file.png"};
    myFile.setMimeType("image/png");
    axom::sina::File myOtherFile{"/path/to/other/file.txt"};
    myOtherFile.setTags({"these", "are", "tags"});

    // Create a record to store the files
    axom::sina::ID myID{"my_record", axom::sina::IDType::Local};
    std::unique_ptr<axom::sina::Record> myRecord{new axom::sina::Record{myID, "my_type"}};

    // Add the files to the record
    myRecord->add(myFile);
    myRecord->add(myOtherFile);

    // Remove a file from the record
    myRecord->remove(myFile);

    std::cout << myRecord->toNode().to_json() << std::endl;
}