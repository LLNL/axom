// State class that holds all the material data.  In the
// multi-material world, this will mostly be a vector of Part structs
// that hold the data, plus a few functions to do mesh sums and
// averages of quantities.

#include "meshapi/Map.hpp"

#include "VectorXY.hpp"
#include "TinyHydroTypes.hpp"
#include "Part.hpp"


namespace tinyHydro {


class PolygonMeshXY;


class State
{
 public:
   State(PolygonMeshXY & theMesh);
   ~State(void);

   // Copy
   State( const State & arg);
   // Assignment
   State & operator=(const State & rhs);
   // addition in place
   State & operator+=(const State & rhs);
   // scalar multiply in place
   State & operator*=(const double s);

   // python accessors
   /* double averageRho(int i) const; */
   /* double totalE(int i) const; */
   VectorXY u(int i) const {return velocity[i];}
   Part * partBegin(void) {return parts;}
   Part * getPart(int i);
   void addPart(Part * newPart);
   void setU(int i, const VectorXY& val) { velocity[i] = val;}
   void setX(int i, const VectorXY& val) { position[i] = val;}

   PolygonMeshXY  * const mesh;

   Part * parts;
   int nParts;
   const int maxNParts;

   NodalVectorField velocity;     // MeshAPI -- Map: Nodes -> VectorXY
   NodalVectorField position;     // MeshAPI -- Map: Nodes -> VectorXY


   void dumpState();
};


} // end namespace tinyHydro
