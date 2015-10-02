#include "Part.hpp"
#include <stdio.h>
#include <string.h>

Part::SubsetRegistry Part::setRegistry;


//----------------------------------------------
Part::Part(const int * zoneList, int nzones, double gamma)
    : gamma(gamma)
{
   printf("in Part c'tor\n");

   zones = ZoneSubset(nzones);
   std::vector<int>& ptr = setRegistry.addNamelessField( &zones).data();
   memcpy(&ptr[0], zoneList, sizeof(int)*nzones);
   zones.data() = &ptr;

   density = new double[nzones];
   energyPerMass = new double[nzones];
   volumeFraction = new double[nzones];
}
//----------------------------------------------
Part::Part(const std::vector<int>& zoneList, double gamma)
    : gamma(gamma)
{
   printf("in Part c'tor\n");

   /// UGLY hack to place zoneList data in the field registry and set the data of zones to this...
//   zones = ZoneSubset::SetBuilder()
//               .size( zoneList.size() )
//               .data( &zoneList );
//   setRegistry.addNamelessField( &zones).data().swap( zoneList);

   zones = ZoneSubset( zoneList.size() );
   std::vector<int>& ptr = setRegistry.addNamelessField( &zones).data();
   memcpy(&ptr[0], &zoneList[0], sizeof(int)*zones.size());
   zones.data() = &zoneList;


   density = new double[numZones()];
   energyPerMass = new double[numZones()];
   volumeFraction = new double[numZones()];
}
//----------------------------------------------
Part::~Part(void)
{
   delete [] density;
   delete [] energyPerMass;
   delete [] volumeFraction;
}

//----------------------------------------------
// Copy
Part::Part(const Part & arg)
    :    gamma(arg.gamma)
{
   // printf("in Part::copy \n");

    zones = ZoneSubset(arg.zones.size());
    std::vector<int>& ptr = setRegistry.addNamelessField( &zones).data();
    memcpy(&ptr[0], &arg.zones[0], sizeof(int)*zones.size());
    zones.data() = &ptr;

   // make space for data
   density = new double[numZones()];
   energyPerMass = new double[numZones()];
   volumeFraction = new double[numZones()];
   

   for (int i = 0; i < numZones(); i++)
   {
      density[i] = arg.density[i];
      energyPerMass[i] = arg.energyPerMass[i];
      volumeFraction[i] = arg.volumeFraction[i];
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

   zones = ZoneSubset(rhs.zones.size());
   std::vector<int>& ptr = setRegistry.addNamelessField( &zones).data();
   memcpy(&ptr[0], &rhs.zones[0], sizeof(int)*zones.size());
   zones.data() = &ptr;

   
   // free up old space
   delete [] density;
   delete [] energyPerMass;
   delete [] volumeFraction;
   
   // make space for data
   density = new double[numZones()];
   energyPerMass = new double[numZones()];
   volumeFraction = new double[numZones()];
   
   // copy data over
   for (int i = 0; i < numZones(); i++)
   {
      density[i] = rhs.density[i];
      energyPerMass[i] = rhs.energyPerMass[i];
   }
   // copy some other data using memcpy for fun
   memcpy(volumeFraction, rhs.volumeFraction, sizeof(double)*numZones());
   
   return *this;
}
   
//----------------------------------------------
// addition in place
Part & Part::operator+=(const Part & rhs)
{
   // printf("in Part::operator+= \n");
   SLIC_ASSERT(numZones() == rhs.numZones() );
   
   for (int i = 0; i < numZones(); i++)
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
   
   for (int i = 0; i < numZones(); i++)
   {
      density[i] *= s;
      energyPerMass[i] *= s;
      volumeFraction[i] *= s;
   }

   return *this;
}
//----------------------------------------------
//----------------------------------------------


void Part::dumpPart()
{
    printf("\nPart has %i zones", numZones() );

    printf("\n\nzones");
    for(int i=0; i< numZones(); ++i)
    {
        printf("\n\t Zone %i -- idx %i"
                ,i, zones[i]
                        );
    }
    printf("\n\n--\n\n");

}
