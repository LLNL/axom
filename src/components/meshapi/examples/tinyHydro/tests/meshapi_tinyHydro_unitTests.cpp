
#include<vector>
#include<cmath>
#include<iostream>

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "../PolygonMeshXY.hpp"
#include "../HydroC.hpp"


using namespace tinyHydro;

std::ostream& operator<<(std::ostream& os, const VectorXY & vec)
{
    os <<"(" << vec.x << "," << vec.y <<")";
    return os;
}



TEST(gtest_meshapi_tinyHydro,test_02_density_with_prescribed_velocity)
{
    std::cout <<"\n****** move mesh with prescribed velocity, test density change." << std::endl;

    int n = 100;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> zones( mesh.numZones() );
    for(int i=0; i< mesh.numZones(); ++i )
        zones[i] = i;

    Part p(zones);
    State s(mesh);

    std::cout <<"**setting initial densities" << std::endl;
    double rho0 = 1.0;
    for (int i = 0; i < mesh.numZones(); i++)
    {
       p.setRho(i,rho0);
       p.setE(i,0.0);
    }

    // set energy spike at the origin
    //double zonemass = rho0 * mesh.zoneVolume[0];
    //double E = 4.48333289847;
    //p.setE(0, E/(4.0*zonemass)); // puts shock at r=0.8 at t=0.3, in XY coords


    std::cout <<"**setting initial velocities" << std::endl;
    for (int i = 0; i < mesh.numNodes(); i++)
    {
       VectorXY pos = mesh.getPos(i);
       double ux = pos.x - 0.5;
       double uy = pos.y - 0.5;
       s.setU(i, VectorXY(-ux,-uy));
    }

    s.addPart(&p);
    Part* pp = s.getPart(0);

    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);

    std::cout <<"**Before we begin:" << std::endl;
    h.Cl = 0.;
    h.Cq = 0.;

    double dt = 0.1;
    h.initialize();


    h.step(dt);

    std::cout <<"**Taking steps:" << std::endl;
    while( h.time < 0.499999 )
        h.step(dt);

    std::cout <<"\n"
            <<" hydro cycle = " << h.cycle
            <<" hydro time  = " << h.time
            << std::endl;

    double rhoTheory = (1.0/(1.0-h.time));
    rhoTheory = rhoTheory*rhoTheory * rho0;

    double tol = 1.0e-10;
    SLIC_ASSERT_MSG(
               std::fabs(pp->rho(0) - rhoTheory) < tol
            && std::fabs(pp->rho(n-1) - rhoTheory) < tol
            , "FAIL -- densities are not correct\n");

    std::cout<<"\nPASS" << std::endl;
}


TEST(gtest_meshapi_tinyHydro,test_03_gradAndForce)
{
    std::cout <<"\n****** test gradient operator and force calculation." << std::endl;

    typedef std::vector<int> IntVec;
    IntVec vals;
    vals.push_back(0);
    vals.push_back(1);
    vals.push_back(10);
    vals.push_back(11);
    vals.push_back(12);
    vals.push_back(13);

    int n = 10;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> zones( mesh.numZones() );
    for(int i=0; i< mesh.numZones(); ++i )
        zones[i] = i;

    Part p(zones);
    State s(mesh);

    std::cout <<"**setting initial density, energy, to make a linear pressure profile" << std::endl;
    double rho0 = 1.0;
    for (int i = 0; i < mesh.numZones(); i++)
    {
       p.setRho(i,rho0);
       VectorXY zPos = mesh.getZonePos(i);
       p.setE(i, zPos.x);
    }

    s.addPart(&p);
    double dt = 0.1;

    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);
    h.setBC("bottom", 0xdeadbeef, 0); // fix bottom nodes in y


    std::cout <<"**Before we begin:" << std::endl;
    h.initialize();

    for(IntVec::iterator it = vals.begin(); it != vals.end(); ++it)
    {
        std::cout <<"density[ " << *it << "] = " <<  p.rho( *it)
                  << " energy[ " << *it << "] = " << p.e(*it)
                  << std::endl;
    }
    std::cout <<"e0, e1, e2, e3, e4 = "
              << p.e(0) << " "
              << p.e(1) << " "
              << p.e(2) << " "
              << p.e(3) << " "
              << p.e(4) << " "
              << std::endl
              ;

    std::cout <<"**Computing forces:" << std::endl;

    h.calcForce(s);

    for(IntVec::iterator it = vals.begin(); it != vals.end(); ++it)
    {
        VectorXY f = h.getForce(*it);
        std::cout<<"\tForce on node " << *it << " = " << f <<"\n";
    }


    double tol = 1e-12;
    VectorXY f12 = h.getForce(12);
    VectorXY f13 = h.getForce(13);
    SLIC_ASSERT_MSG(    std::fabs(f12.x - f13.x) < tol
                     && std::fabs(f12.y) < tol
                     && std::fabs(f13.y) < tol
        , "force calculation FAILS -- f12-f13: " << f12-f13 );

    std::cout <<"**Testing acceleration:" << std::endl;

    h.Cq = 0.0; // turn off Q
    h.Cl = 0.0; // turn off Q
    std::cout <<"\ttaking a step.." << std::endl;
    h.step(dt);


    for(IntVec::iterator it = vals.begin(); it != vals.end(); ++it)
    {
        VectorXY f = s.u(*it);
        std::cout<<"\tNew velocity on node " << *it << ": " << f <<")\n";
    }

    VectorXY u0 = s.u(0);
    VectorXY u11 = s.u(11);
    VectorXY u12 = s.u(12);
    SLIC_ASSERT_MSG(    std::fabs(u0.x - u12.x) < tol
                     && std::fabs(u11.y) < tol
                     && std::fabs(u12.y) < tol
        , "acceleration/velocity calculation FAILS");


    std::cout<<"\nPASS" << std::endl;

    if(false) { std::cout<< f12 << f13 << u0 << u11 << u12; } // tmp workaround for compiler error about unused variables

}



TEST(gtest_meshapi_tinyHydro,test_04_BC)
{
    std::cout <<"\n****** testing boundary conditions." << std::endl;

    int n = 10;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> zones( mesh.numZones() );
    for(int i=0; i< mesh.numZones(); ++i )
        zones[i] = i;

    Part p(zones);
    State s(mesh);

    std::cout <<"**setting initial density, energy, to make a linear pressure profile" << std::endl;
    double rho0 = 1.0;
    for (int i = 0; i < mesh.numZones(); i++)
    {
       p.setRho(i,rho0);
       VectorXY zPos = mesh.getZonePos(i);
       p.setE(i, zPos.x);
    }

    s.addPart(&p);
    double dt = 0.1;

    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);


    std::cout <<"**setting boundary conditions" << std::endl;
    h.setBC("bottom", 0xdeadbeef, 0); // fix bottom nodes in y
    std::cout <<"\tbottom free in x, y is fixed" << std::endl;

    h.setBC("top", 0xdeadbeef,0);
    std::cout <<"\ttop free in x, y is fixed"<< std::endl;

    //h.setBC("left", 0,0); # leave unset to test free condition default
    std::cout <<"\tleft free in x, free in y" << std::endl;

    h.setBC("right", -1,0);
    std::cout <<"\tright -1 in x, y is fixed" << std::endl;

    h.initialize();


    std::cout <<"**stepping" << std::endl;
    h.step(dt);

    VectorXY u0 = s.u(0);
    VectorXY u11 = s.u(11);
    VectorXY u21 = s.u(21);
    VectorXY u120 = s.u(120);

    VectorXY n0 = mesh.getPos(0);
    VectorXY n11 = mesh.getPos(11);
    VectorXY n21 = mesh.getPos(21);
    VectorXY n120 = mesh.getPos(120);


    std::cout <<" Bottom left -- u: "   << u0   <<" -- n: " << n0 << std::endl;
    std::cout <<" Left side -- u: "     << u11  <<" -- n: " << n11 << std::endl;
    std::cout <<" Right side -- u: "    << u21  <<" -- n: " << n21 << std::endl;
    std::cout <<" Top right -- u: "     << u120 <<" -- n: " << n120 << std::endl;


    double tol = 1e-21;

    SLIC_ASSERT_MSG(    std::fabs(u0.y) < tol
                     && std::fabs(u11.x) > tol
                     && std::fabs(u21.x + 1.0) < tol
                     && std::fabs(u120.x + 1.0) < tol
                     && std::fabs(u120.y ) < tol
        , "BC test FAILS ");

    std::cout<<"\nPASS" << std::endl;
}



TEST(gtest_meshapi_tinyHydro,test_05_newDT_Noh)
{
    std::cout <<"\n****** testing newDT method -- Noh." << std::endl;

    int n = 10;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> zones( mesh.numZones() );
    for(int i=0; i< mesh.numZones(); ++i )
        zones[i] = i;

    Part p(zones);
    State s(mesh);

    std::cout <<"**setting initial density, energy, velocity for Noh problem" << std::endl;
    double rho0 = 1.0;
    for (int i = 0; i < mesh.numZones(); i++)
    {
       p.setRho(i,rho0);
       //VectorXY zPos = mesh.getZonePos(i);
       p.setE(i, 0.);
    }
    for (int i = 0; i < mesh.numNodes(); i++)
    {
        s.setU(i,VectorXY(-1.,0.));
    }

    s.addPart(&p);


    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);

    for (int i = 0; i < h.numBCnodes(3); i++)   //set left nodes (left=3) to zero velocity
    {
        h.getState()->setU(h.bcNode(3,i),VectorXY(0,0));
    }

    std::cout <<"**Set boundary conditions!" << std::endl;
    h.setBC("bottom", 0xdeadbeef, 0);
    h.setBC("top",    0xdeadbeef, 0);
    h.setBC("left",            0, 0xdeadbeef);
    h.setBC("right",          -1, 0xdeadbeef);


    std::cout <<"**Printing positions of boundary condition nodes:" << std::endl;
    for(int j=0; j< 4; ++j)
    {
        switch (j)
        {
        case 0: std::cout <<"\n  bottom boundary nodes\n"; break;
        case 1: std::cout <<"\n  right boundary nodes\n"; break;
        case 2: std::cout <<"\n  top boundary nodes\n"; break;
        case 3: std::cout <<"\n  left boundary nodes\n"; break;
        }
        for (int i = 0; i < h.numBCnodes(j); ++i)   //set left nodes (left=3) to zero velocity
        {
            int nodeNum = h.bcNode(j,i);
            std::cout << "\t" << nodeNum << " -- " <<   s.u(nodeNum) ;

        }
    }
    std::cout <<std::endl;



    h.initialize();

    std::cout <<"**Get DT!" << std::endl;
    double dt = h.newDT();
    std::cout <<"\tnewDT = " << dt << std::endl;

    double tol=1.0e-16;
    double expDT = 0.1*h.cfl;
    SLIC_ASSERT_MSG( std::fabs(dt - expDT) < tol
        , " newDT calculation FAILS -- expected dt = " << expDT << " but got " << dt << " instead, leaving " << dt - expDT);

    std::cout<<"\nPASS" << std::endl;
}



TEST(gtest_meshapi_tinyHydro,test_05_newDT_Sedov)
{
    std::cout <<"\n****** testing newDT method -- Sedov." << std::endl;

    int n = 10;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> zones( mesh.numZones() );
    for(int i=0; i< mesh.numZones(); ++i )
        zones[i] = i;

    Part p(zones);
    State s(mesh);

    std::cout <<"**setting initial density, energy, velocity for Sedov problem" << std::endl;
    double rho0 = 1.0;
    for (int i = 0; i < mesh.numZones(); i++)
    {
       p.setRho(i,rho0);
       //kz = i%n;
       //lz = i / n;
       p.setE(i, 0.);
    }

    // set energy spike at the origin
    double zonemass = rho0 * mesh.zoneVol(0);
    double E = 4.48333289847;
    p.setE(0, E/(4.0*zonemass)); // puts shock at r=0.8 at t=0.3, in XY coords


    s.addPart(&p);


    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);

    for (int i = 0; i < mesh.numNodes(); i++)
        s.setU(i,VectorXY(0,0));

    std::cout <<"**set boundary conditions" << std::endl;
    h.setBC("bottom", 0xdeadbeef, 0);
    h.setBC("top",    0xdeadbeef, 0xdeadbeef);
    h.setBC("left",            0, 0xdeadbeef);
    h.setBC("right",  0xdeadbeef, 0xdeadbeef);

    h.initialize();

    std::cout <<"**Get DT!" << std::endl;
    double dt = h.newDT();
    std::cout <<"\tnewDT = " << dt << std::endl;

    double L = 0.1;
    double cfl = 0.7;
    double cs = sqrt(10*E/(4*zonemass*9));
    double theoryDT = cfl * L / cs;

    double tol=1.0e-16;

    SLIC_ASSERT_MSG( std::fabs(dt - theoryDT) < tol
        , " newDT calculation FAILS -- expected dt = " << theoryDT << " but code got " << dt << ". Diff:" <<  dt - theoryDT);

    std::cout<<"\nPASS" << std::endl;

}




TEST(gtest_meshapi_tinyHydro,test_06_PdV_work)
{
    std::cout <<"\n****** testing PdV work." << std::endl;

    int steps = 100;

    int n = 1;
    int kz = n;
    int lz = n;
    PolygonMeshXY mesh(kz+1,lz+1);

    // Create a part
    std::vector<int> zones( mesh.numZones() );
    for(int i=0; i< mesh.numZones(); ++i )
        zones[i] = i;

    Part p(zones);
    State s(mesh);

    std::cout <<"**setting initial density, energy, velocity for 1x1 compression problem" << std::endl;
    double rho0 = 1.0;
    double e0 = 1.0;
    double u0 = -0.1;

    for (int i = 0; i < mesh.numZones(); i++)
    {
       p.setRho(i,rho0);
       p.setE(i, e0);
    }

    s.addPart(&p);

    // right nodes compress at speed u0
    s.setU(1,VectorXY(u0,0));
    s.setU(3,VectorXY(u0,0));
    // set left nodes to zero velocity
    s.setU(0, VectorXY(0,0));
    s.setU(2, VectorXY(0,0));



    std::cout <<"**making hydro" << std::endl;
    Hydro h(&s);

    std::cout <<"**Set boundary conditions!" << std::endl;
    h.setBC("bottom", 0xdeadbeef, 0);
    h.setBC("top",    0xdeadbeef, 0);
    h.setBC("left",            0, 0xdeadbeef);
    h.setBC("right",          u0, 0xdeadbeef);


    std::cout <<"**Turn off Q!" << std::endl;
    h.Cq = 0.;
    h.Cl = 0.;

    std::cout <<"**Initialize" << std::endl;
    h.initialize();

    double dt = 0.01;

    std::cout <<"**take " << steps << " steps with dt = " << dt << std::endl;
    for (int i = 0; i < mesh.numZones(); i++)
        h.step(dt);

    std::cout <<"**done stepping. " << std::endl;

    double tol=1.0e-6;
    double theoryRho = 1.0/(1.0 + h.time*u0);
    double rhoCode = h.getState()->getPart(0)->rho(0);
    SLIC_ASSERT_MSG( std::fabs(theoryRho - rhoCode) < tol
        , "density calculation FAILS -- rhoCode = " << rhoCode << " but should be " << theoryRho );


    double theoryE = std::pow(theoryRho/rho0, 2.0/3.0); // e = e0*(rho/rho0)**(gamma-1)
    double eCode = h.getState()->getPart(0)->e(0);

    std::cout <<"\t h.getState().getPart(0).e(0) " << eCode << std::endl;
    std::cout <<"\t theory e " << theoryE << std::endl;

    SLIC_ASSERT_MSG( std::fabs(theoryE - eCode) < tol
                    , "PdV calculation FAILS");



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
