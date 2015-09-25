#include "Part.hpp"
#include <stdio.h>
#include <string.h>

//----------------------------------------------
Part::Part(const int * zoneList, int nzones,
           double gamma)
    : nzones(nzones)
    , gamma(gamma)
{
   printf("in Part c'tor\n");
   density = new double[nzones];
   energyPerMass = new double[nzones];
   volumeFraction = new double[nzones];
   zones = new int[nzones];
   // copy zone IDs to the Part
   memcpy(zones, zoneList, sizeof(int)*nzones);
}
//----------------------------------------------
Part::Part(const std::vector<int> zoneList, double gamma)
    : nzones(zoneList.size())
    , gamma(gamma)
{
   const int nzones = zoneList.size();
   printf("in Part c'tor\n");
   density = new double[nzones];
   energyPerMass = new double[nzones];
   volumeFraction = new double[nzones];
   zones = new int[nzones];
   // copy zone IDs to the Part
   memcpy(zones, &(*(zoneList.begin())), sizeof(int)*nzones);
}
//----------------------------------------------
Part::~Part(void)
{
   delete [] density;
   delete [] energyPerMass;
   delete [] volumeFraction;
   delete [] zones;
}

//----------------------------------------------
// Copy
Part::Part(const Part & arg)
    : nzones(arg.nzones)
{
   // printf("in Part::copy \n");

   // copy over the gamma value
   gamma = arg.gamma;

   // make space for data
   density = new double[nzones];
   energyPerMass = new double[nzones];
   volumeFraction = new double[nzones];
   zones = new int[nzones];
   
   for (int i = 0; i < nzones; i++)
   {
      density[i] = arg.density[i];
      energyPerMass[i] = arg.energyPerMass[i];
      volumeFraction[i] = arg.volumeFraction[i];
      zones[i] = arg.zones[i];
   }
}

//----------------------------------------------
// Assignment
Part & Part::operator=(const Part & rhs)
{
   // printf("Part::assignment \n");
   
   // check for self-assignment
   if (this == &rhs) return *this;
   
   // copy over the gamma value
   gamma = rhs.gamma;
   nzones = rhs.nzones;
   
   // free up old space
   delete [] density;
   delete [] energyPerMass;
   delete [] volumeFraction;
   delete [] zones;
   
   // make space for data
   density = new double[nzones];
   energyPerMass = new double[nzones];
   volumeFraction = new double[nzones];
   zones = new int[nzones];
   
   // copy data over
   for (int i = 0; i < nzones; i++)
   {
      density[i] = rhs.density[i];
      energyPerMass[i] = rhs.energyPerMass[i];
   }
   // copy some other data using memcpy for fun
   memcpy(volumeFraction, rhs.volumeFraction, sizeof(double)*nzones);
   memcpy(zones, rhs.zones, sizeof(int)*nzones);
   
   return *this;
}
   
//----------------------------------------------
// addition in place
Part & Part::operator+=(const Part & rhs)
{
   // printf("in Part::operator+= \n");
   assert(nzones == rhs.nzones);
   
   for (int i = 0; i < nzones; i++)
   {
      density[i] += rhs.density[i];
      energyPerMass[i] += rhs.energyPerMass[i];
      volumeFraction[i] += rhs.volumeFraction[i];
   }

   return *this;
}
//----------------------------------------------
// scalar multiply in place
Part & Part::operator*=(const double s)
{
   // printf("in Part::operator*= \n");
   
   for (int i = 0; i < nzones; i++)
   {
      density[i] *= s;
      energyPerMass[i] *= s;
      volumeFraction[i] *= s;
   }

   return *this;
}
//----------------------------------------------
//----------------------------------------------

