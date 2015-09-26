// State class that holds all the material data.  In the
// multi-material world, this will mostly be a vector of Part structs
// that hold the data, plus a few functions to do mesh sums and
// averages of quantities.

#include "VectorXY.hpp"
#include "Part.hpp"
class PolygonMeshXY;

class State
{
 public:
   State(PolygonMeshXY & mesh);
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
   void setU(int i, VectorXY val) {velocity[i] = val;}
   void setX(int i, VectorXY val) {position[i] = val;}

   PolygonMeshXY  * const mesh;
   Part * parts;
   int nParts;
   const int maxNParts;
   VectorXY * velocity;
   VectorXY * position;


   void dumpState();
};
