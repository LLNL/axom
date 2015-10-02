// Hydro driver
// Fri Mar 27 14:14:28 PDT 2015

// HYDRO C: Physics-wise, this hydro is still pretty dumb: multiple
// materials, gamma-law gas.  The mesh is an arbitrary collection of
// polygons.  Architecturally, we use as few temporaries as possible.
// No sweet returning of objects when you can instead pass in an array
// and modify it in place.  No STL containers.  Explicit memory
// allocation of all the memory we ever need at problem start, and
// explicit deletes when the hydro object is destroyed.

#include "State.hpp"
#include <stdio.h>

class Hydro
{
 public:
   Hydro(State * s); // using a pointer here with PyBindGen prevents
                     // State being garbage collected by Python

   ~Hydro();

   // problem control and advance methods
   void step(double dt);
   double newDT();
   void steps(int numSteps);
   void advance(double stopTime);
   void initialize(void); // call once before taking steps
   void setBC(const char * boundary, double xVel, double yVel);
   
   // compute time derivatives from first arg, store rates in 2nd arg
   void calcDerivs(const State & state, State & dState, double dt);

   // physics methods
   void calcQ(const State & state);
   void calcForce(const State & state);
   void calcDivU(const VectorXY * u);
   void applyVelocityBC(VectorXY * u);
   double totalEnergy(const State & s) const;
   double totalEnergy(void) const;
   void calcMaxCs(void);
   
   // update all the algebraically determinable stuff, like density
   void updateConstitutives(State & state);

   // setter methods (used from Python)
   void setState(const State s) {state = s;}

   // getter methods
   State * getState(void) {return &state;}
   double getQ(int i) {return Q[i];}
   VectorXY getForce(int i) {return force[i];}
   int numBCnodes(int bc);
   int bcNode(int bc, int node);
   
   PolygonMeshXY & mesh;
   State & state;
   State halfStep; // re-used temporary for RK2
   State dState;   // re-used temporary for RK2

   double * Q; // re-used for Q                                 // MeshAPI -- Map: Zones -> scalar
   double * divu; // for dt, Q, and some PdV work schemes
   VectorXY * force; // re-used for force calc                  // MeshAPI -- Map: Nodes -> VectorXY
   VectorXY * cornerforce; // re-used for force calc            // <UNUSED>?

   VectorXY * halfStepVelocity; // re-used for work calc        // MeshAPI -- Map: Nodes -> VectorXY
   VectorXY * bcVelocity; // velocity BC values                 // MeshAPI -- Map: Domain edges -> VectorXY
                                                                // Code is missing a DomainEdge set of fixed size 4.

   double ** initialMass; // per material initial mass          // MeshAPI -- Map: Parts -> array of doubles -- relation ???
   double * totalMass; //  mass in zone                         // MeshAPI -- Map: Zones -> scalar
   double * nodeMass;  // also keep node masses                 // MeshAPI -- Map: Nodes-> scalar
   double ** pressure; // per material zone pressure            // MeshAPI -- Map: Material -> zone -> scalar
   double * totalPressure; // sum of pressures in a zone        // MeshAPI -- Map: Zones -> scalar
   double * maxCs; //

   int numBCs;                                                  // <UNUSED>?
   int numBC0Nodes, numBC1Nodes, numBC2Nodes, numBC3Nodes;      // MeshAPI -- Indirection Sets
   int * bc0Nodes; // lower boundary nodes
   int * bc1Nodes; // right boundary nodes
   int * bc2Nodes; // top boundary nodes
   int * bc3Nodes; // left boundary nodes

   int cycle;
   int dtZone;
   double time;
   double Cq, Cl; // Q coefficients
   double cfl;
   int printCycle;

   bool fuzzyEqual(double a, double b, double tol=1.0e-15)
   {
      if (fabs(a-b) < tol) return true;
      else return false;
   }
   
};


