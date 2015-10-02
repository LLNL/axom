// Part class holds the material data for a single material.

#ifndef __Part__
#define __Part__

#include <vector>
#include <assert.h>

#include "meshapi/IndirectionSet.hpp"
#include "meshapi/FieldRegistry.hpp"


class Part
{
 public:
    typedef asctoolkit::meshapi::VectorIndirectionSet ZoneSubset;

    typedef asctoolkit::meshapi::FieldRegistry<ZoneSubset::IndexType> SubsetRegistry;
    static SubsetRegistry setRegistry;

    typedef asctoolkit::meshapi::Map<double> ZonalScalarField;

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
   double rho(int i) const {assert (i < numZones() ); return density[i];}
   double e(int i) const {assert(i < numZones() ); return energyPerMass[i];}
   int meshZone(int i) {assert (i < numZones() ); return zones[i];}
   void setRho(int i, double val) {assert(i < numZones() ); density[i] = val;}
   void setE(int i, double val) {assert(i < numZones() ); energyPerMass[i] = val;}

   int numZones() const { return zones.size(); }

   ZoneSubset zones;                // MeshAPI : Indirection set of zones in a part/region

   ZonalScalarField density;                // MeshAPI : Map : Zones -> scalar
   ZonalScalarField energyPerMass;          // MeshAPI : Map : Zones -> scalar
   ZonalScalarField volumeFraction;         // MeshAPI : Map : Zones -> scalar

   double gamma;

   void dumpPart();
};

#endif
