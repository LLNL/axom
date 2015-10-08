#include <stdio.h>

#include "State.hpp"

#include "PolygonMeshXY.hpp"
#include "Part.hpp"

#include "slic/slic.hpp"

namespace tinyHydro {


//----------------------------------------------
State::State(PolygonMeshXY & theMesh)
    : mesh(&theMesh)
    , nParts(0)
    , maxNParts(100)
{
   printf("in State c'tor\n");
   parts = new Part[maxNParts];

   velocity = NodalVectorField( &mesh->nodes);
   position = NodalVectorField( &mesh->nodes);
}
//----------------------------------------------
State::~State(void)
{
   delete [] parts;
}

//----------------------------------------------
Part * State::getPart(int i)
{
   return &(parts[i]);
}
//----------------------------------------------
// addPart
void State::addPart(Part * partPtr)
{
   printf("in State::addPart \n");
   SLIC_ASSERT_MSG( nParts < maxNParts
              , "tried to add more than " << maxNParts <<" parts, limit is hard-wired allowed, craaashing!");

   parts[nParts] = *partPtr; // copy in part data to the list
   nParts++; // now is next unused Part ptr
}

//----------------------------------------------
// Copy operator
State::State(const State & arg)
    : mesh(arg.mesh)
    , nParts(arg.nParts)
    , maxNParts(arg.maxNParts)
{
   // printf("in State::copy \n");

   // copy over the parts
   parts = new Part[maxNParts];
   for (int i = 0; i < nParts; i++)
   {
      parts[i] = arg.parts[i];
   }

   // copy other data
   velocity = NodalVectorField(arg.velocity);
   position = NodalVectorField(arg.position);
}

//----------------------------------------------
// Assignment operator
State & State::operator=(const State & rhs)
{
   // printf("State::assignment \n");
   
   // check for self-assignment
   if (this == &rhs) return *this;
   
   // assignment can only happen when states
   // are defined on the same mesh
   SLIC_ASSERT(mesh == rhs.mesh);
   
   // free up old parts data and copy rhs's parts
   delete [] parts;
   parts = new Part[maxNParts];
   nParts = rhs.nParts;
   for (int i = 0; i < nParts; i++)
   {
      parts[i] = rhs.parts[i];
   }

   // copy other data from rhs
   velocity = NodalVectorField( rhs.velocity);
   position = NodalVectorField( rhs.position);
   
   return *this;
}
   
//----------------------------------------------
// addition in place
State & State::operator+=(const State & rhs)
{
   // printf("in State::operator+= \n");
   SLIC_ASSERT(mesh == rhs.mesh);
   SLIC_ASSERT(nParts == rhs.nParts);
   
   // part data
   for (int i = 0; i < nParts; i++)
   {
      parts[i] += rhs.parts[i];
   }
   
   // mesh-based data
   const int numNodes = mesh->numNodes();
   for (int i = 0; i < numNodes; i++)
   {
      velocity[i] += rhs.velocity[i];
      position[i] += rhs.position[i];
   }

   return *this;
}
//----------------------------------------------
// scalar multiply in place
State & State::operator*=(const double s)
{
   // printf("in State::operator*= \n");
   
   // part data
   for (int i = 0; i < nParts; i++)
   {
      parts[i] *= s;
   }
   
   // mesh-based data
   const int numNodes = mesh->numNodes();
   for (int i = 0; i < numNodes; i++)
   {
      velocity[i] *= s;
      position[i] *= s;
   }

   return *this;
}
//----------------------------------------------
//----------------------------------------------

void State::dumpState()
{
    printf("\nState has %i parts", nParts);

    printf("\n\nParts");
    for(int i=0; i< nParts; ++i)
    {
        parts[i].dumpPart();
    }

    printf("\n\nVelocities");
    for(int i=0; i< mesh->numNodes(); ++i)
    {
        VectorXY v = velocity[i];
        VectorXY p = position[i];
        printf("\n\t Node %i -- pos (%g,%g) -- vel (%g,%g) "
                , i
                , p.x, p.y
                , v.x, v.y
                );
    }

    printf("\n\n--\n\n");

}


} // end namespace tinyHydro
