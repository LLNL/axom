/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-741217
 *
 * All rights reserved.
 *
 * This file is part of Axom.
 *
 * For details about use and distribution, please read axom/LICENSE.
 *
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

// main for profile.  Stupid.

#include "PolygonMeshXY.hpp"
#include "HydroC.hpp"

#include<vector>
#include<cmath>

int main()  // sets up Sedov and runs it
{
   int n = 100;
   int kz = n;
   int lz = n;
   PolygonMeshXY mesh(kz+1,lz+1);

   // Create a part
   std::vector<int> zones( mesh.nzones );
   for(int i=0; i< mesh.nzones; ++i )
       zones[i] = i;

   Part p(zones);
   State s(mesh);

   double rho0 = 1.0;
   for (int i = 0; i < mesh.nzones; i++)
   {
      p.setRho(i,rho0);
      p.setE(i,0.0);
   }

   // set energy spike at the origin
   //double zonemass = rho0 * mesh.zoneVolume[0];
   //double E = 4.48333289847;
   //p.setE(0, E/(4.0*zonemass)); // puts shock at r=0.8 at t=0.3, in XY coords


   printf( "setting initial density, energy, velocity for Sedov problem\n");
   for (int i = 0; i < mesh.nnodes; i++)
   {
      VectorXY pos = mesh.getPos(i);
      double ux = pos.x - 0.5;
      double uy = pos.y - 0.5;
      s.setU(i, VectorXY(-ux,-uy));
   }
   s.addPart(&p);
   Part* pp = s.getPart(0);

   Hydro h(&s);
   h.Cl = 0.;
   h.Cq = 0.;

   double dt = 0.1;
   h.initialize();


   h.step(dt);

   while( h.time < 0.499999 )
       h.step(dt);

   printf("hydro cycle = %i hydro time  = %g\n", h.cycle, h.time);

   double rhoTheory = (1.0/(1.0-h.time));
   rhoTheory = rhoTheory*rhoTheory * rho0;

   if( std::fabs(pp->rho(0) - rhoTheory) < 1.0e-10
    && std::fabs(pp->rho(n-1) - rhoTheory) < 1.0e-10)
       printf ("\nPASS\n");
   else
       printf ("\nFAIL -- densities are not correct\n");

}
