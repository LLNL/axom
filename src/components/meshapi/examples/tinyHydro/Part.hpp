// Part class holds the material data for a single material.

#ifndef __Part__
#define __Part__

#include <vector>

#include "TinyHydroTypes.hpp"

namespace tinyHydro {

class Part
{
 public:
   Part()
     : zones(0)
     , gamma(0.)
   {}

   Part(const int * zoneList, int nzones, double userGamma = 5.0/3.0);
   Part(const std::vector<int>& zoneList, double userGamma = 5.0/3.0);
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
   double rho(int i)                const { return density[i];}
   double e(int i)                  const { return energyPerMass[i];}
   int meshZone(int i)              const { return zones[i];}
   void setRho(int i, double val)         { density[i] = val;}
   void setE(int i, double val)           { energyPerMass[i] = val;}

   int numZones() const { return zones.size(); }

   ZoneSubset zones;                // MeshAPI : Indirection set of zones in a part/region

   ZonalScalarField density;                // MeshAPI : Map : Zones -> scalar
   ZonalScalarField energyPerMass;          // MeshAPI : Map : Zones -> scalar
   ZonalScalarField volumeFraction;         // MeshAPI : Map : Zones -> scalar

   double gamma;

   void dumpPart();
};

} // end namespace tinyHydro

#endif
