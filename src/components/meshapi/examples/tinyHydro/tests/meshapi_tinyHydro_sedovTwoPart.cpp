
#include<vector>
#include<cmath>
#include<iostream>

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "../PolygonMeshXY.hpp"
#include "../HydroC.hpp"
#include "../myTimer.hpp"


using namespace tinyHydro;

std::ostream& operator<<(std::ostream& os, const VectorXY & vec)
{
    os <<"(" << vec.x << "," << vec.y <<")";
    return os;
}


TEST(gtest_meshapi_tinyHydro,test_sedov_2_part)
{
    std::cout <<"\n****** Sedov cylindrical square grid test." << std::endl;

    int steps = 0;
    int n = 80;
    //double dx = 1.0 / n;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> reg1, reg2;

    for(int i=0; i< mesh.numZones(); ++i )
    {
        VectorXY pos = mesh.getZonePos(i);
        if(pos.y > pos.x)
            reg1.push_back(i);
        else
            reg2.push_back(i);
    }

    State s(mesh);

    double gamma = 7.0 / 5.0;
    Part p1(reg1, gamma);
    Part p2(reg2, gamma);

    std::cout <<"**setting initial density, energy, velocity for sedov problem" << std::endl;
    double rho0 = 1.0;
    for (unsigned int i = 0; i < reg1.size(); i++)
    {
       p1.setRho(i,rho0);
       p1.setE(i, 0.);
    }
    for (unsigned int i = 0; i < reg2.size(); i++)
    {
       p2.setRho(i,rho0);
       p2.setE(i, 0.);
    }
    for (int i = 0; i < mesh.numNodes(); i++)
        s.setU(i, VectorXY(0.,0.));


    // set energy spike at the origin
    double zonemass = rho0 * mesh.zoneVol(0);
    double E = 4.48333289847;
    p2.setE(0, E/(4.0*zonemass));  // puts shock at r=0.8 at t=0.3, in XY coords

    s.addPart(&p1);
    s.addPart(&p2);


    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);

    std::cout <<"**Set boundary conditions!" << std::endl;
    h.setBC("bottom", 0xdeadbeef, 0);   // x can slide, y fixed
    h.setBC("left",   0, 0xdeadbeef);   // x fixed, y can slide
    h.setBC("top",    0, 0);            // fixed x,y
    h.setBC("right",  0, 0);            // fixed x,y


    std::cout <<"**Initialize" << std::endl;
    h.initialize();

    h.cfl = 0.1;

    double E0 = h.totalEnergy();

    timespec time1, time2; // TIMER
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);  // TIMER
    if(steps > 0)
    {
        h.steps(steps);
    }
    else
    {
        h.advance(0.3);
    }
    clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2); // TIMER
    double timeElapsed = diffSeconds(time1, time2);
    std::cout<<"Elapsed time after advancing was " << timeElapsed << " seconds." << std::endl;


    double E1 = h.totalEnergy();

    std::cout <<"\n\t Total final energy: " << E1 << std::endl;
    std::cout <<"\t (E1-E0)/E0 = : " << (E1-E0)/E0 << std::endl;

    std::cout<<"\nPASS" << std::endl;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
 int result = 0;

 ::testing::InitGoogleTest(&argc, argv);

 UnitTestLogger logger;   // create & initialize test logger,
 // finalized when exiting main scope

 result = RUN_ALL_TESTS();

 return result;
}
