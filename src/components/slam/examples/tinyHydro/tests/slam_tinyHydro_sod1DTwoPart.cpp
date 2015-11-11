/*
 * Copyright (c) 2015, Lawrence Livermore National Security, LLC.
 * Produced at the Lawrence Livermore National Laboratory.
 *
 * All rights reserved.
 *
 * This source code cannot be distributed without permission and further
 * review from Lawrence Livermore National Laboratory.
 */



#include <vector>
#include <cmath>
#include <iostream>

#include "gtest/gtest.h"

#include "slic/slic.hpp"

#include "../PolygonMeshXY.hpp"
#include "../HydroC.hpp"
#include "../myTimer.hpp"

using namespace tinyHydro;

std::ostream& operator<<(std::ostream& os, const VectorXY & vec)
{
  os << "(" << vec.x << "," << vec.y << ")";
  return os;
}


TEST(gtest_slam_tinyHydro,test_sod1D_2_part)
{
  std::cout << "\n****** Sod 1D test, two parts." << std::endl;


  double goalTime = 0.171;   // puts shock right on 0.80
  double gamma = 7.0 / 5.0;
  double steps = 0;

  // zoning
  int n = 200;
  double dx = 1.0 / static_cast<double>(n);
  int kz = n;
  int lz = 1;

  // initial parameters
  double Pp = 1.0;   // left pressure
  double dp = 1.0;   // left density
  double ep = Pp / (dp * (gamma - 1)); //left specific energy

  //double c0 = sqrt(gamma*Pp/dp); // left sound speed

  double Pm = 1.0 / 10.0; // right pressure
  double dm = 1.0 / 8.0; // right density
  double em = Pm / (dm * (gamma - 1)); // right specific energy

  PolygonMeshXY mesh(kz + 1,lz + 1);

  for(int i = mesh.numNodes() / 2; i< mesh.numNodes(); ++i )
  {
    VectorXY pos = mesh.getPos(i);
    mesh.setPos(i, VectorXY(pos.x, dx));
  }

  // Create a part
  std::vector<int> reg1, reg2;
  int halfZones = mesh.numZones() / 2;
  for(int i = 0; i< halfZones; ++i )
  {
    reg1.push_back(i);
    reg2.push_back(i + halfZones);
  }

  SLIC_ASSERT(  reg1.size() == reg2.size());
  SLIC_ASSERT(  reg1.size() + reg2.size() == static_cast<size_t>( mesh.numZones() ));

  Part p1(reg1, gamma);
  Part p2(reg2, gamma);
  State s(mesh);

  std::cout << "**setting initial density, energy, velocity for sod problem" << std::endl;
  for (unsigned int i = 0; i < reg1.size(); i++)
  {
    p1.setRho(i,dp);
    p1.setE(i, ep);
  }
  for (unsigned int i = 0; i < reg2.size(); i++)
  {
    p2.setRho(i,dm);
    p2.setE(i, em);
  }
  for (int i = 0; i < mesh.numNodes(); i++)
    s.setU(i, VectorXY(0.,0.));

  s.addPart(&p1);
  s.addPart(&p2);

  std::cout << "**making hydro" << std::endl;
  Hydro h(&s);


  std::cout << "**Set boundary conditions!" << std::endl;
  h.setBC("bottom", 0xdeadbeef, 0);
  h.setBC("top",    0xdeadbeef, 0);
  h.setBC("left",   0,          0);
  h.setBC("right",  0,          0);


  std::cout << "**Initialize" << std::endl;
  h.initialize();

  double dt = h.newDT();
  std::cout << "  newDT is " << dt << std::endl;
  h.cfl = 0.2;


  //bool testQ = false;
  //bool testP = false;
  bool special = false;
  double E0 = 0.;   // starting energy

  Timer timer;
  timer.start();

  if (special)
  {
    dt = 0.2;
    h.Cq = 0.0;
    h.Cl = 0.0;
    E0 = h.totalEnergy();
    std::cout << "  total starting energy = " << E0 << std::endl;
    std::cout << "**stepping..." << std::endl;
    h.step(dt);
  }
  else if(steps > 0)
  {
    E0 = h.totalEnergy();
    std::cout << "  total starting energy = " << E0 << std::endl;
    std::cout << "**stepping..." << std::endl;
    h.steps(steps);
  }
  else
  {
    E0 = h.totalEnergy();
    std::cout << "  total starting energy = " << E0 << std::endl;
    std::cout << "**advancing..." << std::endl;
    h.advance(goalTime);     // puts shock at 0.8 (p.57, Triangular Hydro Schemes notebook)
  }
  double E1 = h.totalEnergy();

  timer.stop();
  std::cout << "Elapsed time after advancing was " << timer.getElapsedTime() << " seconds." << std::endl;

  std::cout << "\n\t Total final energy: " << E1 << std::endl;
  std::cout << "\t (E1-E0)/E0 = : " << (E1 - E0) / E0 << std::endl;

  std::cout << "\nPASS" << std::endl;
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
#include "slic/UnitTestLogger.hpp"
using asctoolkit::slic::UnitTestLogger;

int main(int argc, char * argv[])
{
  int result = 0;

  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,

  // finalized when exiting main scope

  result = RUN_ALL_TESTS();

  return result;
}
