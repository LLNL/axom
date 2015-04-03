/*******************************************************************************************************************************
 * An example algorithm for saving/restoring problem state data to a file, using the datastore for identifying and holding the
 * data.
 *
 ******************************************************************************************************************************/

#include "../src/Types.hpp"
#include "../src/DataStore.hpp"
#include "../src/DataBuffer.hpp"

// Setup some sample physics package classes that will use the datastore.
// This example avoids needing function registration in the datastore by having
// the physics packages implement some stock functions from an interface class.
class GenericPhysics
{
   public:
      virtual void setupProblemData(DataStoreNS::DataGroup& group) = 0;

      virtual void notifyPreStateDump(DataStoreNS::DataGroup& group) = 0;
      virtual void notifyPostStateDump(DataStoreNS::DataGroup& group) = 0;
      virtual void notifyPreStateRestore(DataStoreNS::DataGroup& group) = 0;
      virtual void notifyPostStateRestore(DataStoreNS::DataGroup& group) = 0;
};

// This algorithm performs a rather heavy weight setup and teardown of
// it's restart data.  It creates and fully tears down all the 'restart' groups and views after a state dump, and makes them
// fresh at the next start of a dump.
//
// Another version would have an application make all the groups and views, and leave them there, just clean up the data buffers
// of temp data.  This would mean the app needs to maintain the views in the restart sub-tree if any underlying buffers are
// changed though.
class PhysicsA : public GenericPhysics
{
   public:
      void setupProblemData(DataStoreNS::DataGroup& group)
      {
         DataStoreNS::DataGroup& subgroup = *group.CreateGroup("physicsA");

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
         
         dataview = *subgroup.CreateViewAndBuffer("derivedVariable");
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

      void notifyPreStateDump(DataStoreNS::DataGroup& group)
      {
         // Copy over views for persistent data to restart group
         // Create new views/buffers of temporary data.
      }

      void notifyPostStateDump(DataStoreNS::DataGroup& group)
      {
         // Remove views that were copied over
         // Clean up views/buffers of temporary data
      }

      void notifyPreStateRestore(DataStoreNS::DataGroup& group)
      {
         // TBD
      }

      void notifyPostStateRestore(DataStoreNS::DataGroup& group)
      {
         // Check data that was read into group.
         // Copy data into package if needed.
         // Re-calculated derived data if needed.
      }
};

class StateFile
{
   public:
      // Iterates over everything in provided tree and adds it to file
      void save(DataStoreNS::DataGroup& group)
      {
         // Iterate over all groups and views and save out data to file.
         // Only flesh out iteration code, don't care about actual I/O.
      }
      
      // Adds everything from file into provided tree.
      void restore(DataStoreNS::DataGroup& group)
      {
         // Iterate over all groups and views in file, add them to datastore.
         // Only flesh out code needed for iteration, don't care about actual I/O.
      }

      void close() {}

};

void saveState( std::vector< GenericPhysics*>& packages, DataStoreNS::DataGroup& group )
{
   // Create a sub-tree for restart data.
   DataStoreNS::DataGroup& restartGroup = *group.CreateGroup("restart");

   // Tell packages that we are going to perform a state dump of the restart group.  They should populate their data into
   // the group.
   for ( std::vector< GenericPhysics* >::iterator pIt = packages.begin(); pIt != packages.end(); ++pIt)
   {
      (*pIt)->notifyPreStateDump(restartGroup);
   }

   // Give group to restart component to write out data to file
   StateFile file;
   file.save( restartGroup );
   file.close();

   // Tell packages we have finished the state dump, they can perform cleanup of any temporary data.
   for ( std::vector< GenericPhysics* >::iterator pIt = packages.begin(); pIt != packages.end(); ++pIt)
   {
      (*pIt)->notifyPostStateDump(restartGroup);
   }

   // Clean up restart tree
   group.DestroyGroup("restart");
}

void restoreState( std::vector< GenericPhysics*>& packages, DataStoreNS::DataGroup& group )
{
   // Create a sub-tree for restart data.
   DataStoreNS::DataGroup& restartGroup = *group.CreateGroup("restart");

   // Tell packages that we are going to perform a state restore of the restart group.
   for ( std::vector< GenericPhysics* >::iterator pIt = packages.begin(); pIt != packages.end(); ++pIt)
   {
      (*pIt)->notifyPreStateRestore(restartGroup);
   }

   // Give group to restart component to read in data from file
   StateFile file;
   file.restore( restartGroup );
   file.close();

   // Tell packages we have restored the state data into the datastore restart group.  They should pull out the data they
   // need, copy data over to the main datastore tree, re-calculate any dependent derived data, etc.
   for ( std::vector< GenericPhysics* >::iterator pIt = packages.begin(); pIt != packages.end(); ++pIt)
   {
      (*pIt)->notifyPostStateRestore(restartGroup);
   }

   // Clean up restart tree
   group.DestroyGroup("restart");
}

int main(void)
{

   // Create datastore and problem state data.
   DataStoreNS::DataStore datastore;
   DataStoreNS::DataGroup& rootGroup = *datastore.GetRoot();

   // Create some physics packages that will use the datastore.
   std::vector< GenericPhysics* > packages;
   packages.push_back( new PhysicsA() );

   // Tell packages to populate datastore with their problem data.
   for ( std::vector< GenericPhysics* >::iterator pIt = packages.begin(); pIt != packages.end(); ++pIt)
   {
      (*pIt)->setupProblemData( rootGroup );
   }

   // Save state
   saveState( packages, rootGroup );

   // Restore state
   restoreState( packages, rootGroup );

   for ( std::vector< GenericPhysics* >::iterator pIt = packages.begin(); pIt != packages.end(); ++pIt)
   {
      delete *pIt;
   }

   return 0;
}


