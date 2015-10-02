// Part class holds the material data for a single material.

#ifndef __Part__
#define __Part__

#include <vector>
#include <assert.h>

class Part
{
 public:
   Part()
     : nzones(0)
     , zones(0)
     , density(0)
     , energyPerMass(0)
     , volumeFraction(0)
     , gamma(0.)
   {}

   Part(const int * zoneList, int nzones, double userGamma = 5.0/3.0);
   Part(const std::vector<int> zoneList, double userGamma = 5.0/3.0);
   ~Part(void);

   // Copy
   Part( const Part & arg);
   // Assignment
   Part & operator=(const Part & rhs);
   // addition in place
   Part & operator+=(const Part & rhs);
   // scalar multiply in place
   Part & operator*=(const double s);

   // python accessors
   double rho(int i) const {assert (i < nzones); return density[i];}
   double e(int i) const {assert(i < nzones); return energyPerMass[i];}
   int meshZone(int i) {assert (i < nzones); return zones[i];}
   void setRho(int i, double val) {assert(i < nzones); density[i] = val;}
   void setE(int i, double val) {assert(i < nzones); energyPerMass[i] = val;}

   int nzones;                      // MeshAPI : Indirection set of zones in a part/region
   int * zones;

   double * density;                // MeshAPI : Map : Zones -> scalar
   double * energyPerMass;          // MeshAPI : Map : Zones -> scalar
   double * volumeFraction;         // MeshAPI : Map : Zones -> scalar

   double gamma;

   void dumpPart();
};

#endif
