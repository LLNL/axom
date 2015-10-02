#include "State.hpp"
#include "PolygonMeshXY.hpp"
#include "Part.hpp"
#include <stdio.h>

#include "slic/slic.hpp"

//----------------------------------------------
State::State(PolygonMeshXY & theMesh)
    : mesh(&theMesh)
    , nParts(0)
    , maxNParts(100)
{
   printf("in State c'tor\n");
   velocity = NodalVectorField( &mesh->nodes);
   position = NodalVectorField( &mesh->nodes);

   parts = new Part[maxNParts];
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
   if (nParts < maxNParts)
   {
      parts[nParts] = *partPtr; // copy in part data to the list
      nParts++; // now is next unused Part ptr
   }
   else
   {
      printf("tried to add more than %d parts, limit is hard-wired allowed, craaashing!\n",maxNParts);
      SLIC_ASSERT(false);
   }
}

//----------------------------------------------
// Copy operator
State::State(const State & arg):
mesh(arg.mesh),
   nParts(arg.nParts),
   maxNParts(arg.maxNParts)
{
   // printf("in State::copy \n");

   // make space for data
   velocity = NodalVectorField( &mesh->nodes);
   position = NodalVectorField( &mesh->nodes);

   parts = new Part[maxNParts];

   // copy over the parts
   for (int i = 0; i < nParts; i++)
   {
      parts[i] = arg.parts[i];
   }

   // copy over the mesh-based data
   const int nnodes = mesh->numNodes();
   for (int i = 0; i < nnodes; i++)
   {
      velocity[i] = arg.velocity[i];
      position[i] = arg.position[i];
   }
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
   
   // free up old space
   delete [] parts;
   
   // make space for data
   velocity = NodalVectorField( &mesh->nodes);
   position = NodalVectorField( &mesh->nodes);
   parts = new Part[maxNParts];

   // copy over the parts
   nParts = rhs.nParts;
   for (int i = 0; i < nParts; i++)
   {
      parts[i] = rhs.parts[i];
   }
   // copy mesh-based stuff over
   const int nn = mesh->numNodes();
   for (int i = 0; i < nn; i++)
   {
      velocity[i] = rhs.velocity[i];
      position[i] = rhs.position[i];
   }
   
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
   for (int i = 0; i < mesh->numNodes(); i++)
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
   for (int i = 0; i < mesh->numNodes(); i++)
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


