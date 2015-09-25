#include <stdio.h>
#include "HydroC.hpp"
#include "PolygonMeshXY.hpp"
#include "Part.hpp"
#include <assert.h>
#include <string.h>
// #include "/usr/include/gperftools/profiler.h"

//----------------------------------------------
// c'tor
Hydro::Hydro(State * s):
mesh(*(s->mesh)),
   state(*s),
   halfStep(*s),
   dState(*s),
   cycle(0),
   dtZone(-1),
   time(0.0),
   Cq(1.0),
   Cl(1.0),
   cfl(0.7)
{
   printf("using HydroC\n");

   assert(state.nParts > 0);
   assert(mesh.nzones > 0);
   
   // make memory for Q calculation
   Q = new double[mesh.nzones];
   
   // make memory for Q calculation
   force = new VectorXY[mesh.nnodes];
   /* cornerforce = new VectorXY * [nParts]; */
   /* for (int p = 0; p < nParts; p++) cornerforce[p] = new VectorXY[4*state.parts[p].nzones]; */
   
   // make memory for holding zone and node masses
   const int nParts = state.nParts;
   initialMass = new double * [nParts];
   for (int p = 0; p < nParts; p++) initialMass[p] = new double[state.parts[p].nzones];
   totalMass = new double [mesh.nzones];
   totalPressure = new double [mesh.nzones];
   maxCs =  new double [mesh.nzones];
   nodeMass = new double [mesh.nnodes];
      

   // space for half-step velocity
   halfStepVelocity = new VectorXY[mesh.nnodes];

   // for dt, Q, and some PdV work schemes
   divu = new double[mesh.nzones];

   // zone pressure
   pressure = new double * [nParts];
   for (int p = 0; p < nParts; p++) pressure[p] = new double[state.parts[p].nzones];
   
   // velocity BC space.  BC's are implicity ordered as 'bottom',
   // 'right', 'top', 'left', numbered 0 through 3.
   bcVelocity = new VectorXY[4];
   // initial BC's are 0xdeadbeef, which we take as a special value
   // meaning free boundary condition.
   for (int i = 0; i < 4; i ++)
      bcVelocity[i] = VectorXY(0xdeadbeef, 0xdeadbeef);
   // find the boundary nodes if any
   numBC0Nodes = numBC1Nodes = numBC2Nodes = numBC3Nodes = 0;
   double minX = 1e99;
   double maxX = -1e99;
   double minY = 1e99;
   double maxY = -1e99;
   for (int n = 0; n < mesh.nnodes; n++)
   {
      VectorXY & pos = mesh.nodePos[n];
      if (pos.x < minX) minX = pos.x;
      if (pos.x > maxX) maxX = pos.x;
      if (pos.y < minY) minY = pos.y;
      if (pos.y > maxY) maxY = pos.y;
   }
   for (int n = 0; n < mesh.nnodes; n++)
   {
      VectorXY & pos = mesh.nodePos[n];
      if (fuzzyEqual(pos.y, minY)) numBC0Nodes++;
      if (fuzzyEqual(pos.x, maxX)) numBC1Nodes++;
      if (fuzzyEqual(pos.y, maxY)) numBC2Nodes++;
      if (fuzzyEqual(pos.x, minX)) numBC3Nodes++;
   }
   bc0Nodes = new int[numBC0Nodes];
   bc1Nodes = new int[numBC1Nodes];
   bc2Nodes = new int[numBC2Nodes];
   bc3Nodes = new int[numBC3Nodes];
   int i0 = 0; int i1 = 0; int i2 = 0; int i3 = 0;
   for (int n = 0; n < mesh.nnodes; n++)
   {
      VectorXY & pos = mesh.nodePos[n];
      if (fuzzyEqual(pos.y, minY)) {bc0Nodes[i0] = n; i0++;}
      if (fuzzyEqual(pos.x, maxX)) {bc1Nodes[i1] = n; i1++;}
      if (fuzzyEqual(pos.y, maxY)) {bc2Nodes[i2] = n; i2++;}
      if (fuzzyEqual(pos.x, minX)) {bc3Nodes[i3] = n; i3++;}
   }
}
   
//----------------------------------------------
// d'tor
Hydro::~Hydro()
{
   delete [] Q;
   delete [] force;
   /* delete [] cornerforce; */
   for (int p = 0; p < state.nParts; p++) delete [] initialMass[p];
   delete [] initialMass;
   delete [] totalMass;
   delete [] totalPressure;
   delete [] nodeMass;
   delete [] halfStepVelocity;
   delete [] divu;
   for (int p = 0; p < state.nParts; p++) delete [] pressure[p];
   delete [] pressure;
   delete [] bcVelocity;
   delete [] bc0Nodes;
   delete [] bc1Nodes;
   delete [] bc2Nodes;
   delete [] bc3Nodes;
}
//----------------------------------------------
// set hydro ready to begin taking steps
void Hydro::initialize(void)
{
   // calc new geometry data in mesh, including
   // new zone volumes
   mesh.computeNewGeometry();
   
   // get our state ready with initial mesh positions
   for (int i = 0; i < mesh.nnodes; i++)
   {
      state.position[i] = mesh.nodePos[i];
   }

   // tuck away initial zone and node masses and never forget them
   memset( (void *) totalMass, 0, sizeof(double)*mesh.nzones); // zero out before summing
   memset( (void *) nodeMass,  0, sizeof(double)*mesh.nzones); // zero out before summing
   for (int p = 0; p < state.nParts; p++)  // loop over parts
   {
      const int * zones = state.parts[p].zones;
      const double * rho = state.parts[p].density;
      for (int iz = 0; iz < state.parts[p].nzones; iz++) // loop over zones
      {
         initialMass[p][iz] = rho[iz]*mesh.zoneVolume[zones[iz]];
         totalMass[zones[iz]] += initialMass[p][iz];
         for (int in = 0; in < mesh.zNumNodes[zones[iz]]; in ++) // loop over nodes
         {
            int node = mesh.zNodes[mesh.z2firstNode[zones[iz]]+in];
            nodeMass[node] += 0.25*initialMass[p][iz];
         }
      } // end loop over zones
   } // end loop over parts
}
   
//----------------------------------------------
// computes the time derivatives, stores them in 2nd arg
void Hydro::calcDerivs(const State & s, State & dState, double dt)
{
   // printf("Hydro:: blago calcDerivs 0\n");

   const int nnodes = mesh.nnodes;
   
   // calc forces
   calcForce(s);

   // Accelerations are the delta velocities for dState, so let's
   // label that explicitly.  Re-naming makes it less confusing, at
   // least to me.
   VectorXY * accel = dState.velocity;
   
   // calc accelerations:
   for (int n = 0; n < nnodes; n++)
   {
      accel[n].x = force[n].x/nodeMass[n];
      accel[n].y = force[n].y/nodeMass[n];
   }

   // halfStepVelocity = state.u + 0.5*dt*accel
   
   // Velocities must be anchored at the half-way point in time to get
   // the conservation properties we like. Note that state.u is the
   // beginning of timestep velocity, which we use even during the
   // second partial RK2 step.
   for (int i = 0; i < mesh.nnodes; i++)
   {
      halfStepVelocity[i] = accel[i];
      halfStepVelocity[i] *= 0.5*dt;
      halfStepVelocity[i] += state.velocity[i]; 
   }
   
   // apply velocity BCs to half-step velocity to prevent
   // spurious position changes and PdV work
   applyVelocityBC(halfStepVelocity);
   
   // velocities are the delta of position for dState.
   for (int i = 0; i < mesh.nnodes; i++)
   {
      dState.position[i] = halfStepVelocity[i];
   }
   
   // calc rate of (P+Q)dV work: to be energy conservative we want
   // the velocity divergence computed on the half-step velocity
   calcDivU(halfStepVelocity);

   ////// interestingly, in the finite volume approach, the
   ////// force.dot(velocity) energy update is algebraically identical
   ////// to the simpler P*div(u) calculation.  So we do the simpler one.
   /* const int kmax = mesh.kmax; */
   /* for (int z = 0; z < mesh.nzones; z++) */
   /* { */
   /*    dState.energyPerMass[z] = -(cornerforce[4*z].dot(halfStepVelocity[LL(z)]) */
   /*                                +cornerforce[4*z+1].dot(halfStepVelocity[LR(z)]) */
   /*                                +cornerforce[4*z+2].dot(halfStepVelocity[UR(z)]) */
   /*                                +cornerforce[4*z+3].dot(halfStepVelocity[UL(z)])) / */
   /*       initialMass[z]; */
   /* } */
   
   // the energy production rate is simply -(P+Q) div(u)
   for (int p = 0; p < s.nParts; p++)
   {
      const double * rho = s.parts[p].density;
      const double * P = pressure[p];
      double * de = dState.parts[p].energyPerMass;
      const int nzones = state.parts[p].nzones;
      const int * zones = state.parts[p].zones;
      for (int z = 0; z < nzones; z++)
      {
         de[z] = -(P[z] + Q[zones[z]]) * divu[zones[z]] / rho[z];
      }
   }
   // printf("Hydro:: blago calcDerivs 9\n");

}
//----------------------------------------------
#include <algorithm>
// compute the total force that pushes the nodes
void Hydro::calcForce(const State & s)
{
   assert(s.mesh == &mesh);
   
   // compute grad of zone-centered (P+Q), return pointer to array of
   // node-centered VectorXY forces.
   std::fill(force, force+mesh.nnodes, VectorXY(0.0,0.0));

   // zero out total before summing
   for (int z = 0; z < mesh.nzones; z ++) totalPressure[z] = 0.0;

   // get pressure from EOS: gamma law gas, P = rho*e*(gamma-1)
   for (int p = 0; p < s.nParts; p++)
   {
      const int nzones = state.parts[p].nzones;
      const double gammaMinusOne = s.parts[p].gamma - 1.0;
      const double * rho = s.parts[p].density;
      const double * e = s.parts[p].energyPerMass;
      double * P = pressure[p];
      const int * zones = s.parts[p].zones;
      for (int z = 0; z < nzones; z++)
      {
         P[z] = gammaMinusOne*rho[z]*e[z];
         totalPressure[zones[z]] += P[z];
      }
   }
   
   // get Q
   calcQ(state);

   // compute grad via finite volume method
   // grad(P) = (1/Vn) sum(P*A) for P,A around node
   // for each zone, accumulate P*A to its nodes
   for (int iz = 0; iz < mesh.nzones; iz++)
   {
      for (int in = 0; in < mesh.zNumNodes[iz]; in++)
      {
         int node = mesh.zNodes[mesh.z2firstNode[iz]+in];
         int face0 = mesh.zFaces[mesh.z2firstFace[iz]+ (in + mesh.zNumNodes[iz]-1) % mesh.zNumNodes[iz]];
         int face1 = mesh.zFaces[mesh.z2firstFace[iz]+in];
         VectorXY A = mesh.faceArea[face0] + mesh.faceArea[face1];
         // multiply areas by pressure + Q (and div 2 for area factor)
         double totP = 0.5 * (totalPressure[iz] + Q[iz]);
         A *= totP;
         force[node] += A;
      }
   }
}
//----------------------------------------------
// compute the artificial viscosity
void Hydro::calcQ(const State & s)
{
   // Q = rho((L*div(u))**2 - cs*L*div(u))

   // calculate divergence of velocity
   calcDivU(s.velocity);

   // set max sound speed
   calcMaxCs();

   // compute Q for each zone
   for (int z = 0; z < mesh.nzones; z++)
   {
      if (divu[z] < 0.0)
      {
         const double vol = mesh.zoneVolume[z];
         const double L = sqrt(vol);
         const double deltaU = L*divu[z];
         double rho = totalMass[z]/vol;
         Q[z] = rho*(Cq*deltaU*deltaU - Cl*maxCs[z]*deltaU);
      }
      else
         Q[z] = 0.0;
   }
}

//----------------------------------------------
// update all the algebraically determinable stuff
void Hydro::updateConstitutives(State & s)
{
   // move mesh to new node positions
   mesh.moveNodesToPosition(s.position);

   // calc new geometry data in mesh, including
   // new zone volumes
   mesh.computeNewGeometry();
   for (int p = 0; p < s.nParts; p++)
   {
      const int * zones = s.parts[p].zones;
      const int nzones = s.parts[p].nzones;
      double * rho = s.parts[p].density;
      for (int z = 0; z < nzones; z++)
      {
         rho[z] = initialMass[p][z]/mesh.zoneVolume[zones[z]];
      }
   }

   // apply BC's to newly computed velocity
   applyVelocityBC(s.velocity);
}

//----------------------------------------------
void Hydro::step(double dt)
{
   printf("cycle %d  time = %f  dt = %f\n", cycle, time, dt);
   // printf("Hydro::step blago 1\n");

   // ----------------------------------
   // update physics variables using RK2
   // ----------------------------------
   // printf("Hydro::step blago 2\n");
   time += 0.5*dt;
   calcDerivs(state, dState, 0.5*dt); // dState = time derivatives
   // halfstep = state + dt/2 * dState
   halfStep = state;
   dState *= 0.5*dt;
   halfStep += dState;
   
   // update all the algebraically determinable stuff
   updateConstitutives(halfStep);

   // second step of our RK2 integration
   calcDerivs(halfStep, dState, dt);
   // final state = state + dt * dState(halfStep)
   dState *= dt;
   state += dState;
   
   // update all the algebraically determinable stuff
   updateConstitutives(state);

   // update control info
   time += 0.5*dt; // full time update
   cycle++;

   // printf("Hydro::step blago 9\n");
      
}
//----------------------------------------------
// calc max sound speed in each zone
void Hydro::calcMaxCs(void)
{
   // zero to start with
   memset(maxCs, 0, sizeof(double)*mesh.nzones);
   
   for (int p = 0; p < state.nParts; p++)
   {
      const int * zones = state.parts[p].zones;
      const int nzones = state.parts[p].nzones;
      const double gammaGammaMinusOne = state.parts[p].gamma*(state.parts[p].gamma - 1.0);
      const double * e = state.parts[p].energyPerMass;
      for (int z = 0; z < nzones; z++)
      {
         double cs2 = gammaGammaMinusOne*e[z];
         maxCs[zones[z]] = maxCs[zones[z]] > cs2 ? maxCs[zones[z]] : cs2;
      }
   }
   // we saved all the sqrts for the very end
   for (int z = 0; z < mesh.nzones; z++)
   {
      maxCs[z] = sqrt(maxCs[z]);
   }
}

//----------------------------------------------
double Hydro::newDT(void)
{
   double minDT = 1.0e99; // impossibly high dt
   dtZone = -1; // impossible zone index

   calcDivU(state.velocity);  // we need velocity divergence
   calcMaxCs(); // need max sound speed 
   
   for (int z = 0; z < mesh.nzones; z++)
   {
      double L2 = mesh.zoneVolume[z];
      double dt2 = L2/(maxCs[z]*maxCs[z] + L2*divu[z]*divu[z]);
      if (minDT > dt2)
      {
         minDT = dt2;
         dtZone = z;
      }
   }
   return cfl*sqrt(minDT);
}
//----------------------------------------------
// compute the divergence of velocity
void Hydro::calcDivU(const VectorXY * velocity)
{
   // div u = 1/V * sum(u dot dA)
   for (int iz = 0; iz < mesh.nzones; iz++)
   {
      divu[iz] = 0.0;
      for (int in = 0; in < mesh.zNumNodes[iz]; in++)
      {
         int node = mesh.zNodes[mesh.z2firstNode[iz]+in];
         int face0 = mesh.zFaces[mesh.z2firstFace[iz]+ (in + mesh.zNumNodes[iz]-1) % mesh.zNumNodes[iz]];
         int face1 = mesh.zFaces[mesh.z2firstFace[iz]+in];
         VectorXY A = mesh.faceArea[face0] + mesh.faceArea[face1];
         divu[iz] += velocity[node].dot(A);
      }
      divu[iz] /= 2.0*mesh.zoneVolume[iz]; // the 2 normalizes the areas
   }
}

//----------------------------------------------
void Hydro::setBC(const char * boundary, double xVel, double yVel)
{
   switch (boundary[0] )
   {
   case 'b':
      bcVelocity[0] = VectorXY(xVel, yVel);
      break;
   case 'r':
      bcVelocity[1] = VectorXY(xVel, yVel);
      break;
   case 't':
      bcVelocity[2] = VectorXY(xVel, yVel);
      break;
   case 'l':
      bcVelocity[3] = VectorXY(xVel, yVel);
      break;
   default:
      printf("ERROR, bad value specified for BC\n");
      break;
   }
}
//----------------------------------------------
void Hydro::applyVelocityBC(VectorXY * u)
{
   // bottom
   if (bcVelocity[0].x != 0xdeadbeef)
   {
      for (int in = 0; in < numBC0Nodes; in++)
         u[bc0Nodes[in]].x = bcVelocity[0].x;
   }
   if (bcVelocity[0].y != 0xdeadbeef)
   {
      for (int in = 0; in < numBC0Nodes; in++)
         u[bc0Nodes[in]].y = bcVelocity[0].y;
   }
   // right
   if (bcVelocity[1].x != 0xdeadbeef)
   {
      for (int in = 0; in < numBC1Nodes; in++)
         u[bc1Nodes[in]].x = bcVelocity[1].x;
   }
   if (bcVelocity[1].y != 0xdeadbeef)
   {
      for (int in = 0; in < numBC1Nodes; in++)
         u[bc1Nodes[in]].y = bcVelocity[1].y;
   }
   // top
   if (bcVelocity[2].x != 0xdeadbeef)
   {
      for (int in = 0; in < numBC2Nodes; in++)
         u[bc2Nodes[in]].x = bcVelocity[2].x;
   }
   if (bcVelocity[2].y != 0xdeadbeef)
   {
      for (int in = 0; in < numBC2Nodes; in++)
         u[bc2Nodes[in]].y = bcVelocity[2].y;
   }
   // left
   if (bcVelocity[3].x != 0xdeadbeef)
   {
      for (int in = 0; in < numBC3Nodes; in++)
         u[bc3Nodes[in]].x = bcVelocity[3].x;
   }
   if (bcVelocity[3].y != 0xdeadbeef)
   {
      for (int in = 0; in < numBC3Nodes; in++)
         u[bc3Nodes[in]].y = bcVelocity[3].y;
   }
}
//----------------------------------------------
void Hydro::steps(int numSteps)
{
   for (int i = 0; i < numSteps; i++)
   {
      double dt = newDT();
      step(dt);
   }
}
//----------------------------------------------
void Hydro::advance(double stopTime)
{
   /* ProfilerStart("profile.out"); */
   /* struct ProfilerState ps; */
   
   while (time < stopTime)
   {
      double dt = newDT();
      step(dt);
   }
   
   /* ProfilerGetCurrentState(&ps); */
   /* printf("profiler: enabled = %d, samples gathered = %d\n", ps.enabled, ps.samples_gathered); */
   /* ProfilerStop(); */
   
}
//----------------------------------------------
double Hydro::totalEnergy(const State & s) const
{
   double totE = 0.0;
   for (int p = 0; p < s.nParts; p++)
   {
      const int nzones = s.parts[p].nzones;
      const double * e = s.parts[p].energyPerMass;
      const double * m = initialMass[p];
      for (int z = 0; z < nzones; z++)
      {
         totE += e[z]*m[z];
      }
   }
   for (int n = 0; n < mesh.nnodes; n++)
   {
      totE += 0.5 * nodeMass[n] * s.velocity[n].dot(s.velocity[n]);
   }
   return totE;
}
//----------------------------------------------
int Hydro::numBCnodes(int bc)
{
   switch (bc)
   {
   case 0:
      return numBC0Nodes;
      break;
   case 1:
      return numBC1Nodes;
      break;
   case 2:
      return numBC2Nodes;
      break;
   case 3:
      return numBC3Nodes;
      break;
   default:
      assert(bc != 0 and bc != 1 and bc != 2 and bc != 3);
      return -1;
   }
}
//----------------------------------------------
int Hydro::bcNode(int bc, int node)
{
   switch (bc)
   {
   case 0:
      return bc0Nodes[node];
      break;
   case 1:
      return bc1Nodes[node];
      break;
   case 2:
      return bc2Nodes[node];
      break;
   case 3:
      return bc3Nodes[node];
      break;
   default:
      assert(bc != 0 and bc != 1 and bc != 2 and bc != 3);
      return -1;
   }
}
//----------------------------------------------
double Hydro::totalEnergy(void) const
{
   return totalEnergy(state);
}
//----------------------------------------------
