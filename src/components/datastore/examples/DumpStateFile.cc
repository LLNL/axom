/*******************************************************************************************************************************
 * An example algorithm for saving problem state data to a file, using the datastore for identifying and holding the data.
 *
 * This example uses a simple design, where a fresh 'restart' group is created before each dump.  It also assumes that a
 * physics package wants to write out all the data it is keeping in the datastore.
 *
 * This example also shows how a client code that is used to C++ references can convert from the pointers returned by the
 * datastore API.
 ******************************************************************************************************************************/

#include "../src/Types.hpp"
#include "../src/DataStore.hpp"
#include "../src/DataBuffer.hpp"

class PhysicsPackage
{
   public:

      PhysicsPackage(DataStoreNS::DataGroup& group):mDataGroup(group) {}

      void setup()
      {
         DataStoreNS::DataGroup& subgroup = *mDataGroup.CreateGroup("physicsB");

         DataStoreNS::DataView& dataview = *subgroup.CreateViewAndBuffer("variable1");
         DataStoreNS::DataBuffer& buffer = *dataview.GetBuffer();
         buffer.Declare(DataType::float64(100));
         buffer.Allocate();
         conduit::float64* data_ptr1 = buffer.GetNode().as_float64_ptr();
         // how do you get number of entries in conduit node? (ie, number of float64's)?
         for (size_t i=0; i < 100; ++i)
         {
            data_ptr1[i] = i;
         }

         dataview = *subgroup.CreateViewAndBuffer("variable2");
         buffer = *dataview.GetBuffer();
         buffer.Declare(DataType::float64(100));
         buffer.Allocate();
         conduit::float64* data_ptr2 = buffer.GetNode().as_float64_ptr();
         // how do you get number of entries in conduit node? (ie, number of float64's)?
         for (size_t i=0; i < 100; ++i)
         {
            data_ptr2[i] = 100-i;
         }
         
         dataview = *subgroup.CreateViewAndBuffer("dependentVariable");
         buffer = *dataview.GetBuffer();
         buffer.Declare(DataType::float64(100));
         buffer.Allocate();
         conduit::float64* data_ptr3 = buffer.GetNode().as_float64_ptr();
         // how do you get number of entries in conduit node? (ie, number of float64's)?
         for (size_t i=0; i < 100; ++i)
         {
            data_ptr3[i] = data_ptr1[i] * data_ptr2[i];
         }

      }

      void saveState(DataStoreNS::DataGroup& group)
      {
         // Since the package wants to save all it's data as-is, it can just copy over it's group (with views) to the restart group.
         // It assumes a fresh, empty restart group is provided.
         group.CopyGroup( &mDataGroup );
      }

   private:
      DataStoreNS::DataGroup& mDataGroup;

};
   
class StateFile
{
   public:
      // Iterates over everything in provided tree and adds it to file
      void save(DataStoreNS::DataGroup& group)
      {
         
         // Iterate over all groups and views to exercise needed API calls.
         // A real code would follow up by writing each item to file.
      }
      
      // Read everything from file into group.
      void restore(DataStoreNS::DataGroup& group)
      {
         // Restore state data back into group, exercising needed API calls.
         // A real code would read in each item from file first.
      }

      void close() {}

};

int main(void)
{

   // Create datastore and problem state data.
   DataStoreNS::DataStore datastore;
   DataStoreNS::DataGroup& rootGroup = *datastore.GetRoot();

   // Create a sub-tree for restart data.
   DataStoreNS::DataGroup& restartGroup = *rootGroup.CreateGroup("restart");

   // Create example physics package that will use datastore.
   PhysicsPackage physics(rootGroup);

   // Tell package to populate datastore with it's problem data.
   physics.setup();

   // Tell physics package to populate the 'restart' group with data it wants to save.
   physics.saveState( restartGroup );

   // Give 'restart' data group to another component, responsible for writing it out to file.
   StateFile file;
   file.save( restartGroup );
   file.close();

   // Clean up restart tree
   rootGroup.DestroyGroup("restart");

   return 0;
}


