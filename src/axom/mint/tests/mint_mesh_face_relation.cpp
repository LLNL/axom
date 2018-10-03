/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2017, Lawrence Livermore National Security, LLC.
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
#include "axom/mint/config.hpp"                      /* for mint defintions */
#include "axom/mint/mesh/UnstructuredMesh.hpp"       /* for mint::UnstructuredMesh */
#include "axom/mint/mesh/internal/MeshHelpers.hpp"   /* for mint::initFaces */

#include "axom/core/utilities/Utilities.hpp"         /* for utilities::max */

#include "axom/slic/core/UnitTestLogger.hpp"         /* for UnitTestLogger */
#include "axom/slic/interface/slic.hpp"              /* for slic macros */

#include "gtest/gtest.h"                             /* for TEST and EXPECT_* macros */

#include <sstream>
#include <vector>

namespace axom
{
namespace mint
{

namespace internal
{

struct MeshFaceTest
{
  MeshFaceTest(): name(""), mesh(nullptr), initShouldSucceed(false),
                  totalFaceCount(-1)
  { }

  MeshFaceTest(std::string thename, Mesh * themesh, bool itsInitShouldSucceed,
               IndexType theTotalFaceCount,
               std::vector< IndexType > theCellFaceCount,
               std::vector< IndexType > theCellNeighbors):
    name(thename), mesh(themesh), initShouldSucceed(itsInitShouldSucceed),
    totalFaceCount(theTotalFaceCount)
  {
    cellFaceCount = theCellFaceCount;
    cellNeighbors = theCellNeighbors;
  }

  ~MeshFaceTest()
  {
    delete mesh;
  }

  std::string name;
  Mesh * mesh;
  bool initShouldSucceed;
  IndexType totalFaceCount;
  std::vector< IndexType > cellFaceCount;
  std::vector< IndexType > cellNeighbors;
};

} /* end namespace internal */

std::vector<internal::MeshFaceTest *> generateFaceTestCases()
{
  std::vector<internal::MeshFaceTest *> tests;

  // Each test mesh in "Face relation tests.svg" is instantiated
  // and put into tests.
  constexpr int TWO_D = 2;
  constexpr int THREE_D = 3;
  constexpr bool EXPECT_INIT_SUCCESS = true;
  constexpr bool EXPECT_INIT_FAILURE = false;

  {
    // 1. tri  ===============================================================
    UnstructuredMesh< SINGLE_SHAPE > * tri =
      new UnstructuredMesh< SINGLE_SHAPE >(TWO_D, TRIANGLE);
    double    trinodes[] = { 0, 0, 1, 0, 0, 1 } ;
    IndexType tricells[] = { 0, 1, 2 } ;
  
    tri->appendNodes(trinodes, 3);
    tri->appendCells(tricells, 1);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces, face count,
                    ("tri", tri, EXPECT_INIT_SUCCESS, 3,
                     // per-cell face count, neighbor cells for each cell.
                     { 3 },  { -1, -1, -1 }));
  }

  {
    // 2. two tris  ==========================================================
    UnstructuredMesh< SINGLE_SHAPE > * twotris =
      new UnstructuredMesh< SINGLE_SHAPE >(TWO_D, TRIANGLE);
    double    twotrisnodes[] = { 0, 0, 1, 0, 0, 1, 0.8, 1.2 } ;
    IndexType twotriscells[] = { 0, 1, 2,
                                 1, 3, 2 } ;
  
    twotris->appendNodes(twotrisnodes, 4);
    twotris->appendCells(twotriscells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("two tris", twotris, EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     5, { 3, 3 },
                     // neighbor cells for each cell.
                     { -1,  1, -1,
                       -1, -1,  0 }));
  }

  {
    // 3. three quads and a tri  =============================================
    UnstructuredMesh< MIXED_SHAPE > * thrqtri =
      new UnstructuredMesh< MIXED_SHAPE >(TWO_D);
    double    thrqtrixs[] = { -1, -1, -1,  0, 0, 0, 1, 1 };
    double    thrqtriys[] = { -1,  0,  1, -1, 0, 1, 0, 1 };
    IndexType thrqtricells[] = { 0, 3, 4, 1, // lower-left quad
                                 1, 4, 5, 2, // upper-left quad
                                 3, 6, 4,    // lower-right tri
                                 4, 6, 7, 5  } ; // upper-right quad
    // since thrqtri is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType thrqtrioffs[] = { 0, 4, 8, 11, 15 };
    CellType  thrqtritype[] = { QUAD, QUAD, TRIANGLE, QUAD };

    thrqtri->appendNodes(thrqtrixs, thrqtriys, 8);
    thrqtri->appendCells(thrqtricells, 4, thrqtrioffs, thrqtritype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("three quads and a tri", thrqtri,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     11, { 4, 4, 3, 4 },
                     // neighbors at each face of the cells
                     { -1,  2,  1, -1,
                        0,  3, -1, -1,
                       -1,  3,  0,
                        2, -1, -1,  1 }));
  }

  {
    // 4. four tris, with a hole  ============================================
    UnstructuredMesh< SINGLE_SHAPE > * fourtris =
      new UnstructuredMesh< SINGLE_SHAPE >(TWO_D, TRIANGLE);
    double    fourtrisxs[] = {  -1,  0, 0, -.2,  .2, 1 };
    double    fourtrisys[] = { -.1, -1, 0,   1, -.2, 0 };
    IndexType fourtriscells[] = { 0, 1, 2,
                                  0, 2, 3,
                                  1, 5, 4,
                                  2, 5, 3 };
  
    fourtris->appendNodes(fourtrisxs, fourtrisys, 6);
    fourtris->appendCells(fourtriscells, 4);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("four tris, with a hole", fourtris,
                     // should init faces, total face count,
                     EXPECT_INIT_SUCCESS, 10,
                     // per-cell face count,
                     { 3, 3, 3, 3 },
                     // neighbors at each face of each cell.
                     { -1, -1,  1,
                        0,  3, -1,
                       -1, -1, -1,
                       -1, -1,  1 }));
  }

  {
    // 5. one 3D tri  ========================================================
    UnstructuredMesh< SINGLE_SHAPE > * threeDtri =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TRIANGLE);
    double    threeDtrinodes[] = { -1, 0, 0, 0, 1, 0.5, 1.2, -.2, 3 };
    IndexType threeDtricells[] = { 0, 1, 2 };
  
    threeDtri->appendNodes(threeDtrinodes, 3);
    threeDtri->appendCells(threeDtricells, 1);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("one 3D tri", threeDtri,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     3, { 3 },
                     // neighbors at each face of each cell.
                     { -1, -1, -1 }));
  }

  {
    // 6. a tri not far from a quad  =========================================
    UnstructuredMesh< MIXED_SHAPE > * qandtri =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    double    qandtrixs[] = { -1, 1,   0, 0.3, 0.9,   2, 1.4 };
    double    qandtriys[] = {  0, 0, 0.8, 0.9, 0.2, 0.5, 1.5 };
    double    qandtrizs[] = {  1, 0,   0,   0,   0,   0,   0 };
    IndexType qandtricells[] = { 0, 1, 2,       // left: the triangle
                                 3, 4, 5, 6 } ; // right: the quad
    // since qandtri is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType qandtrioffs[] = { 0, 3, 7 };
    CellType  qandtritype[] = { TRIANGLE, QUAD };

    qandtri->appendNodes(qandtrixs, qandtriys, qandtrizs, 7);
    qandtri->appendCells(qandtricells, 2, qandtrioffs, qandtritype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("tri separated from quad", qandtri,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     7, { 3, 4 },
                     // neighbors at each face of the cells
                     { -1, -1, -1,
                       -1, -1, -1, -1 }));
  }

  {
    // 7. tet from four tris  ================================================
    UnstructuredMesh< SINGLE_SHAPE > * tettris =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TRIANGLE);
    double    tettrisxs[] = { 0,   1,   1, 1 };
    double    tettrisys[] = { 0,   0,   1, 1 };
    double    tettriszs[] = { 0, -.1, 0.2, 1 };
    IndexType tettriscells[] = { 0, 2, 1, 0, 1, 3, 1, 2, 3, 2, 0, 3 };
  
    tettris->appendNodes(tettrisxs, tettrisys, tettriszs, 4);
    tettris->appendCells(tettriscells, 4);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("tet from four tris", tettris,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     6, { 3, 3, 3, 3 },
                     // neighbors at each face of each cell.
                     { 3, 2, 1,
                       0, 2, 3,
                       0, 3, 1,
                       0, 1, 2 }));
  }

  {
    // 8. hex from six quads  ================================================
    UnstructuredMesh< SINGLE_SHAPE > * hexquads =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, QUAD);
    double    hexquadsxs[] = { 0, 1, 1, 0, 0, 1, 1, 0 };
    double    hexquadsys[] = { 0, 0, 1, 1, 0, 0, 1, 1 };
    double    hexquadszs[] = { 0, 0, 0, 0, 1, 1, 1, 1 };
    IndexType hexquadscells[] = { 0, 3, 2, 1,
                                  0, 1, 5, 4,
                                  1, 2, 6, 5,
                                  2, 3, 7, 6,
                                  3, 0, 4, 7,
                                  4, 5, 6, 7 };
  
    hexquads->appendNodes(hexquadsxs, hexquadsys, hexquadszs, 8);
    hexquads->appendCells(hexquadscells, 6);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("hex from six quads", hexquads,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     12, { 4, 4, 4, 4, 4, 4 },
                     // neighbors at each face of each cell.
                     { 4, 3, 2, 1,
                       0, 2, 5, 4,
                       0, 3, 5, 1,
                       0, 4, 5, 2,
                       0, 1, 5, 3,
                       1, 2, 3, 4 }));
  }

  {
    // 9. pyramid from a quad and four tris  =================================
    UnstructuredMesh< MIXED_SHAPE > * pyr =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    double    pyrxs[] = { -1, 0, 0, 1, 0 };
    double    pyrys[] = { 0, -1, 1, 0, 0 };
    double    pyrzs[] = { 0,  0, 0, 0, 1 };
    IndexType pyrcells[] = { 0, 1, 3, 2,       // quad for the base
                             0, 1, 4,          // four tris for the sides
                             1, 3, 4,
                             3, 2, 4,
                             2, 0, 4 } ;
    // since pyr is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType pyroffs[] = { 0, 4, 7, 10, 13, 16 };
    CellType  pyrtype[] = { QUAD, TRIANGLE, TRIANGLE, TRIANGLE, TRIANGLE };

    pyr->appendNodes(pyrxs, pyrys, pyrzs, 5);
    pyr->appendCells(pyrcells, 5, pyroffs, pyrtype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("hollow pyramid", pyr, EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     8, { 4, 3, 3, 3, 3 },
                     // neighbors at each face of the cells
                     { 1, 2, 3, 4,
                       0, 2, 4,
                       0, 3, 1,
                       0, 4, 2,
                       0, 1, 3 }));
  }

  {
    // 10. two tris back to back, forming a closed surface  ==================
    UnstructuredMesh< SINGLE_SHAPE > * b2btris =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TRIANGLE);
    double    b2btrisnodes[] = { 0,   0,   0,
                                 2, -.3, -.1,
                                 1,   1,   1 };
    IndexType b2btriscells[] = { 0, 1, 2,
                                 0, 2, 1 };
  
    b2btris->appendNodes(b2btrisnodes, 3);
    b2btris->appendCells(b2btriscells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("back-to-back tris", b2btris,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     3, { 3, 3 },
                     // neighbors at each face of each cell.
                     { 1, 1, 1,
                       0, 0, 0 }));
  }

  {
    // 11. three quads (corner of a box)  ====================================
    UnstructuredMesh< MIXED_SHAPE > * threeq =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    // for variety, you can use a MIXED_SHAPE for a homogeneous mesh---
    // it just means a little more typing.
    double    threeqxs[] = { 0, 1, 1, 0, 1.4, 1.4, 0.4 };
    double    threeqys[] = { 0, 0, 1, 1, 0.4, 1.4, 1.4 };
    double    threeqzs[] = { 0, 0, 0, 0, 0.4, 0.4, 0.4 };
    IndexType threeqcells[] = { 0, 1, 2, 3,
                                1, 4, 5, 2,
                                2, 5, 6, 3 };
    // since threeq is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType threeqoffs[] = { 0, 4, 8, 12 };
    CellType  threeqtype[] = { QUAD, QUAD, QUAD };

    threeq->appendNodes(threeqxs, threeqys, threeqzs, 7);
    threeq->appendCells(threeqcells, 3, threeqoffs, threeqtype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("box corner", threeq, EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     9, { 4, 4, 4 },
                     // neighbors at each face of the cells
                     { -1,  1,  2, -1,
                       -1, -1,  2,  0,
                        1, -1, -1,  0 }));
  }

  {
    // 12. two quads, two tris  ==============================================
    UnstructuredMesh< MIXED_SHAPE > * twoqtwot =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    double    twoqtwotxs[] = { 0, 1, 1, 0, 1.4, 1.4, 0.4 };
    double    twoqtwotys[] = { 0, 0, 1, 1, 0.4, 1.4, 1.4 };
    double    twoqtwotzs[] = { 0, 0, 0, 0, 0.4, 0.4, 0.4 };
    IndexType twoqtwotcells[] = { 0, 1, 2,
                                  1, 4, 5, 2,
                                  2, 5, 6, 3,
                                  3, 0, 2 };
    // since twoqtwot is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType twoqtwotoffs[] = { 0, 3, 7, 11, 14 };
    CellType  twoqtwottype[] = { TRIANGLE, QUAD, QUAD, TRIANGLE };

    twoqtwot->appendNodes(twoqtwotxs, twoqtwotys, twoqtwotzs, 7);
    twoqtwot->appendCells(twoqtwotcells, 4, twoqtwotoffs, twoqtwottype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("box corner (tris and quads)", twoqtwot,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     10, { 3, 4, 4, 3 },
                     // neighbors at each face of the cells
                     { -1,  1,  3,
                       -1, -1,  2, 0,
                        1, -1, -1, 3,
                       -1,  0,  2 }));
  }

  {
    // 13. two quads, two tris forming a prism open at both ends  ============
    UnstructuredMesh< MIXED_SHAPE > * oprism =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    double    oprismxs[] = { 0, 1, 2, 0, 1, 2 };
    double    oprismys[] = { 0, -1, 1, 0, -1, 1 };
    double    oprismzs[] = { 0, 0, 0, 1, 1, 1 };
    IndexType oprismcells[] = { 0, 1, 4, 3,
                                1, 2, 4,
                                2, 5, 4,
                                5, 2, 0, 3 };
    // since oprism is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType oprismoffs[] = { 0, 4, 7, 10, 14 };
    CellType  oprismtype[] = { QUAD, TRIANGLE, TRIANGLE, QUAD };

    oprism->appendNodes(oprismxs, oprismys, oprismzs, 6);
    oprism->appendCells(oprismcells, 4, oprismoffs, oprismtype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("open prism (quads, tris)", oprism,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     10, { 4, 3, 3, 4 },
                     // neighbors at each face of the cells
                     { -1,  1, -1,  3,
                       -1,  2,  0,
                        3, -1,  1,
                        2, -1,  0, -1}));
  }

  {
    // 14. cracked tet  ======================================================
    UnstructuredMesh< SINGLE_SHAPE > * crackedtet =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TRIANGLE);
    double    crackedtetxs[] = { 0,   1,   1, 1, 0.9 };
    double    crackedtetys[] = { 0,   0,   1, 1, 0.9 };
    double    crackedtetzs[] = { 0, -.1, 0.2, 1, 1 };
    IndexType crackedtetcells[] = { 0, 2, 1, 0, 1, 4, 1, 2, 3, 2, 0, 3 };
  
    crackedtet->appendNodes(crackedtetxs, crackedtetys, crackedtetzs, 5);
    crackedtet->appendCells(crackedtetcells, 4);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("cracked tet", crackedtet,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     8, { 3, 3, 3, 3 },
                     // neighbors at each face of each cell.
                     { 3,  2,  1,
                       0, -1, -1,
                       0,  3, -1,
                       0, -1,  2 }));
  }

  {
    // 15. cracked pyramid  ==================================================
    UnstructuredMesh< MIXED_SHAPE > * crackedpyr =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    double    crackedpyrxs[] = { -1, 0, 0, 1, 0,  0.2 };
    double    crackedpyrys[] = { 0, -1, 1, 0, 0, -0.2 };
    double    crackedpyrzs[] = { 0,  0, 0, 0, 1,    1 };
    IndexType crackedpyrcells[] = { 0, 1, 3, 2,     // quad for the base
                                    0, 1, 4,        // four tris for the sides
                                    1, 3, 5,
                                    3, 2, 4,
                                    2, 0, 4 } ;
    // since crackedpyr is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType crackedpyroffs[] = { 0, 4, 7, 10, 13, 16 };
    CellType  crackedpyrtype[] = { QUAD, TRIANGLE, TRIANGLE, TRIANGLE,
                                   TRIANGLE };

    crackedpyr->appendNodes(crackedpyrxs, crackedpyrys, crackedpyrzs, 6);
    crackedpyr->appendCells(crackedpyrcells, 5, crackedpyroffs, crackedpyrtype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("cracked pyramid", crackedpyr,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     10, { 4, 3, 3, 3, 3 },
                     // neighbors at each face of the cells
                     { 1, 2, 3, 4,
                       0, -1,  4,
                       0, -1, -1,
                       0,  4, -1,
                       0,  1,  3 }));
  }

  {
    // 16. not a manifold  ===================================================
    UnstructuredMesh< SINGLE_SHAPE > * notmanf =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TRIANGLE);
    double    notmanfxs[] = {  -1, 0, 0, 0.2, 1 };
    double    notmanfys[] = {   0, 0, 0, -.6, 0 };
    double    notmanfzs[] = { 0.6, 0, 1, 0.4, 0.4 };
    IndexType notmanfcells[] = { 0, 1, 2, 
                                 3, 1, 2,
                                 4, 1, 2 };
  
    notmanf->appendNodes(notmanfxs, notmanfys, notmanfzs, 5);
    notmanf->appendCells(notmanfcells, 3);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("not a manifold", notmanf,
                     // should NOT init faces,
                     EXPECT_INIT_FAILURE,
                     // total face count, per-cell face count,
                     7, { 3, 3, 3 },
                     // neighbors at each face of each cell
                     // can't be properly expressed.
                     { -1, -1, -1,
                       -1, -1, -1,
                       -1, -1, -1 }));
  }

  {
    // 17. egregiously not a manifold  =======================================
    UnstructuredMesh< MIXED_SHAPE > * egreg =
      new UnstructuredMesh< MIXED_SHAPE >(THREE_D);
    double    egregxs[] = { 0, 1, 1, 0,   0, 1, 0,    1,    1 };
    double    egregys[] = { 0, 0, 1, 1, 0.8, 0, 0, -0.4, -0.9 };
    double    egregzs[] = { 0, 0, 0, 0, 0.4, 1, 1,  1.2,  0.5 };
    IndexType egregcells[] = { 0, 1, 2, 3,
                               0, 1, 4, 
                               0, 1, 5, 6,
                               0, 1, 7,
                               0, 1, 8 } ;
    // since egreg is a mixed topology mesh, we need to give offsets and
    // shapes to properly define the cells.
    IndexType egregoffs[] = { 0, 4, 7, 11, 14, 17 };
    CellType  egregtype[] = { QUAD, TRIANGLE, QUAD, TRIANGLE, TRIANGLE };

    egreg->appendNodes(egregxs, egregys, egregzs, 9);
    egreg->appendCells(egregcells, 5, egregoffs, egregtype);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("egregiously not a manifold", egreg,
                     //  should NOT init faces,
                     EXPECT_INIT_FAILURE,
                     // total face count, per-cell face count,
                     10, { 4, 3, 3, 3, 3 },
                     // neighbors at each face of each cell
                     // can't be properly expressed.
                     { -1, -1, -1, -1,
                       -1, -1, -1,
                       -1, -1, -1, -1,
                       -1, -1, -1,
                       -1, -1, -1 }));
  }

  {
    // 18. 3D tet  ===========================================================
    UnstructuredMesh< SINGLE_SHAPE > * tet =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TET);
    double    tetxs[] = { 0,   1,   1, 1 };
    double    tetys[] = { 0,   0,   1, 1 };
    double    tetzs[] = { 0, -.1, 0.2, 1 };
    IndexType tetcells[] = { 0, 1, 2, 3 };
  
    tet->appendNodes(tetxs, tetys, tetzs, 4);
    tet->appendCells(tetcells, 1);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("3D tet", tet, EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     4, { 4 },
                     // neighbors at each face of each cell.
                     { -1, -1, -1, -1 }));
  }

  {
    // 19. two hexs  =========================================================
    UnstructuredMesh< SINGLE_SHAPE > * hexs =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, HEX);
    double    hexsxs[] = { 0,   1, 2, 0,   1, 3, 0,   1, 2, 0,   1, 3 };
    double    hexsys[] = { 0,   0, 0, 1, 0.8, 1, 0, 0.2, 0, 1,   1, 1 };
    double    hexszs[] = { 0, 0.2, 0, 0,   0, 0, 1,   1, 1, 1, 0.8, 1 };
    IndexType hexscells[] = { 0, 1, 4, 3, 6, 7, 10, 9,
                              1, 2, 5, 4, 7, 8, 11, 10 };
  
    hexs->appendNodes(hexsxs, hexsys, hexszs, 12);
    hexs->appendCells(hexscells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("two hexs", hexs, EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     11, { 6, 6 },
                     // neighbors at each face of each cell.
                     { -1, -1,  1, -1, -1, -1,
                       -1, -1, -1, -1,  0, -1 }));
  }

  {
    // 20. three hexs  =======================================================
    UnstructuredMesh< SINGLE_SHAPE > * hexs3 =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, HEX);
    double    hexs3xs[] = { 0, 0,  1, 1, 1, 2, 2, 0, 0,  1, 1, 1, 2, 2 };
    double    hexs3ys[] = { 0, 0,  0, 0, 0, 0, 0, 1, 1,  1, 1, 1, 1, 1 };
    double    hexs3zs[] = { 0, 2, -1, 1, 3, 0, 2, 0, 2, -1, 1, 3, 0, 2 };
    IndexType hexs3cells[] = { 0, 3, 10,  7, 1, 4, 11,  8,
                               3, 5, 12, 10, 4, 6, 13, 11,
                               0, 2,  9,  7, 3, 5, 12, 10 };
  
    hexs3->appendNodes(hexs3xs, hexs3ys, hexs3zs, 14);
    hexs3->appendCells(hexs3cells, 3);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer, should init faces,
                    ("three hexs", hexs3, EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     15, { 6, 6, 6 },
                     // neighbors at each face of each cell.
                     { 2, -1,  1, -1, -1, -1,
                       2, -1, -1, -1,  0, -1,
                      -1, -1, -1, -1,  0,  1 }));
  }

  {
    // 21. two coincident tets, one inside-out, forming a closed manifold  ===
    UnstructuredMesh< SINGLE_SHAPE > * mtet =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TET);
    double    mtetxs[] = { 0,   1,   1, 1 };
    double    mtetys[] = { 0,   0,   1, 1 };
    double    mtetzs[] = { 0, -.1, 0.2, 1 };
    IndexType mtetcells[] = { 0, 1, 2, 3,
                              0, 2, 1, 3 };
  
    mtet->appendNodes(mtetxs, mtetys, mtetzs, 4);
    mtet->appendCells(mtetcells, 2);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("two \"back-to-back\" tets", mtet,
                     // should init faces,
                     EXPECT_INIT_SUCCESS,
                     // total face count, per-cell face count,
                     4, { 4, 4 },
                     // neighbors at each face of each cell.
                     { 1, 1, 1, 1,
                       0, 0, 0, 0 }));
  }

  {
    // 22. bad tet mesh (not a manifold)  ====================================
    UnstructuredMesh< SINGLE_SHAPE > * badtets =
      new UnstructuredMesh< SINGLE_SHAPE >(THREE_D, TET);
    double    badtetsxs[] = { 0,   1,   1, 1,   0, 0.3 };
    double    badtetsys[] = { 0,   0,   1, 1,   1, 1.2 };
    double    badtetszs[] = { 0, -.1, 0.2, 1, 0.5, 0.8 };
    IndexType badtetscells[] = { 0, 1, 2, 3,
                                 0, 3, 2, 4,
                                 0, 3, 2, 5 };
  
    badtets->appendNodes(badtetsxs, badtetsys, badtetszs, 6);
    badtets->appendCells(badtetscells, 3);
    tests.push_back(new internal::MeshFaceTest
                    // test name, test pointer,
                    ("non-manifold tet mesh", badtets,
                     // should NOT init faces,
                     EXPECT_INIT_FAILURE,
                     // total face count, per-cell face count,
                     10, { 4, 4, 4 },
                     // neighbors at each face of each cell can't be properly
                     // expressed.
                     { -1, -1, -1, -1,
                       -1, -1, -1, -1,
                       -1, -1, -1, -1 }));
  }

  return tests;
}

bool verifyNeighbors(IndexType facecount,
                     IndexType * testnbrs,
                     IndexType * nbrs)
{
  std::map<IndexType, int> testnbrset, nbrset;

  for (IndexType i = 0; i < facecount; ++i)
  {
    if (testnbrset.count(testnbrs[i]) < 1) { testnbrset[testnbrs[i]] = 1; }
    else { testnbrset[testnbrs[i]] = testnbrset[testnbrs[i]] + 1; }

    if (nbrset.count(nbrs[i]) < 1) { nbrset[nbrs[i]] = 1; }
    else { nbrset[nbrs[i]] = nbrset[nbrs[i]] + 1; }
  }

  for (auto nbrcount : nbrset)
  {
    if (testnbrset.count(nbrcount.first) != 1) { return false; }
    if (testnbrset[nbrcount.first] != nbrcount.second) { return false; }
    testnbrset.erase(nbrcount.first);
  }
  for (auto testnbrcount : testnbrset)
  {
    if (testnbrcount.first != -1) { return false; }
  }

  return true;
}

void runMeshFaceTest(internal::MeshFaceTest * t)
{
  int facecount = -1;
  IndexType * f2c = nullptr;
  IndexType * c2f = nullptr;
  IndexType * c2n = nullptr;
  IndexType * c2foffsets = nullptr;

  bool initresult =
    internal::initFaces(t->mesh, facecount, f2c, c2f, c2n, c2foffsets);

  if (!initresult && !t->initShouldSucceed)
  {
    SCOPED_TRACE(t->name);
    SUCCEED();
  }
  else if (!initresult && t->initShouldSucceed)
  {
    FAIL() << "test mesh \"" << t->name << 
      "\" call to initFaces() failed but should have succeeded.";
  }
  else if (initresult && !(t->initShouldSucceed))
  {
    delete [] f2c;
    delete [] c2f;
    delete [] c2n;
    delete [] c2foffsets;

    FAIL() << "test mesh \"" << t->name << 
      "\" call to initFaces() succeeded but should have failed.";
  }
  else
  {
    // initFaces() succeeded where it should have!
    SCOPED_TRACE(t->name);

    Mesh * m = t->mesh;

    // do we have enough faces on each of the cells?
    IndexType cellcount = m->getNumberOfCells();
    for (IndexType c = 0; c < cellcount; ++c)
    {
      IndexType thisCellFaceCount = c2foffsets[c + 1] - c2foffsets[c];
      EXPECT_EQ(thisCellFaceCount, t->cellFaceCount[c]) <<
        "thisCellFaceCount " << thisCellFaceCount <<
        " differs from t->cellFaceCount[c] " << t->cellFaceCount[c] <<
        " with c == " << c;
    }
    EXPECT_EQ(facecount, t->totalFaceCount);

    // does each cell have the neighbors we expect?
    for (IndexType c = 0; c < cellcount; ++c)
    {
      std::stringstream errmesg;

      IndexType thisCellFaceCount = c2foffsets[c + 1] - c2foffsets[c];
      bool neighborsMatched =
        verifyNeighbors(thisCellFaceCount, &c2n[c2foffsets[c]],
                       &(t->cellNeighbors[c2foffsets[c]]));
      if (!neighborsMatched)
      {
        errmesg << "Cell " << c << " expected neighbors ";
        IndexType facecount = t->cellFaceCount[c];
        for (int i = 0; i < facecount; ++i) 
        {
          errmesg << t->cellNeighbors[c2foffsets[c] + i] << " ";
        }
        errmesg << std::endl << "found neighbors ";
        for (int i = 0; i < thisCellFaceCount; ++i)
        {
          errmesg << c2n[c2foffsets[c] + i] << " ";
        }
      }
      EXPECT_TRUE(neighborsMatched) << errmesg.str();
    }

    delete [] f2c;
    delete [] c2f;
    delete [] c2n;
    delete [] c2foffsets;
  }
}

} /* end namespace mint */
} /* end namespace axom */

/* Here are tests for the mint face relation method.
 * 
 * In general, the steps are
 * 1. generate a test mesh with some pre-written answers
 *    1. the mesh itself
 *    2. the number of faces for each cell and the mesh's total number of
 *       faces
 *    3. the mesh's cell-to-cell neighbor relation
 * 2. run axom::mint::internal::initFaces
 * 3. from the face-cell and cell-face relations, generate the cell-cell
 *    neighbor relation
 * 4. make sure the entire mesh and each of its cells has the right number
 *    of faces
 * 5. make sure each cell has the right neighbors
 * 6. clean up.
 *
 * The function std::vector<internal::MeshFaceTest> generateFaceTestCases()
 * will produce mesh face test cases corresponding to the figures in the file
 * mint_mesh_face_relation.svg along with the test answers (step 1 above).
 *
 * The function void runMeshFaceTest(internal::MeshFaceTest *) will call
 * internal::initFaces() on the test mesh, generate the neighbor relation,
 * and verify the face counts and neighbor relation as scoped tests
 * (steps 2--6 above).
 *
 * The driver for all the mesh face relation tests, calling
 * generateTestCase() and runMeshFaceTest(), is the TEST function below.
 */

TEST( mint_mesh_face_relation, correct_construction )
{
  std::vector<axom::mint::internal::MeshFaceTest *> tests =
    axom::mint::generateFaceTestCases();

  for (axom::mint::internal::MeshFaceTest * t : tests)
  {
    axom::mint::runMeshFaceTest(t);
    delete t;
  }
}

//------------------------------------------------------------------------------
#include "axom/slic/core/UnitTestLogger.hpp"
using axom::slic::UnitTestLogger;

int main(int argc, char* argv[])
{
  ::testing::InitGoogleTest(&argc, argv);

  UnitTestLogger logger;  // create & initialize test logger,
  // finalized when exiting main scope

  return RUN_ALL_TESTS();
}
