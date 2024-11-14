// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "opencascade/BRepAdaptor_Curve.hxx"
#include "opencascade/BRepBuilderAPI_NurbsConvert.hxx"
#include "opencascade/BRepTools.hxx"
#include "opencascade/BRep_Tool.hxx"
#include "opencascade/GeomAbs_CurveType.hxx"
#include "opencascade/Geom_BSplineCurve.hxx"
#include "opencascade/Geom_Curve.hxx"
#include "opencascade/STEPControl_Reader.hxx"
#include "opencascade/TopAbs.hxx"
#include "opencascade/TopExp_Explorer.hxx"
#include "opencascade/TopoDS_Edge.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include "opencascade/TopoDS.hxx"
#include "opencascade/Geom2d_Curve.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/Geom2d_Curve.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/BRep_Tool.hxx"

#include <iostream>

TopoDS_Shape processStepFile(const std::string& filename)
{
  STEPControl_Reader reader;
  IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());

  if(status != IFSelect_RetDone)
  {
    std::cerr << "Error: Cannot read the file." << std::endl;
    return TopoDS_Shape();
  }

  Standard_Integer numRoots = reader.NbRootsForTransfer();
  reader.TransferRoots();
  TopoDS_Shape shape = reader.OneShape();

  if(shape.IsNull())
  {
    std::cerr << "Error: No shape found in the file." << std::endl;
    return TopoDS_Shape();
  }

  SLIC_INFO(axom::fmt::format("Successfully read the STEP file with {} roots",
                              numRoots));

  return shape;
}

void printShapeStats(const TopoDS_Shape& shape)
{
  if(shape.IsNull())
  {
    std::cerr << "Error: The shape is invalid or empty." << std::endl;
    return;
  }

  int numPatches = 0;
  int numTrimmingCurves = 0;

  for(TopExp_Explorer exp(shape, TopAbs_FACE); exp.More(); exp.Next())
  {
    numPatches++;
  }

  for(TopExp_Explorer exp(shape, TopAbs_EDGE); exp.More(); exp.Next())
  {
    numTrimmingCurves++;
  }
  std::map<int, int> patchTrimmingCurves;
  int patchIndex = 0;

  for(TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More();
      faceExp.Next(), ++patchIndex)
  {
    int trimmingCurvesCount = 0;
    for(TopExp_Explorer edgeExp(faceExp.Current(), TopAbs_EDGE); edgeExp.More();
        edgeExp.Next())
    {
      trimmingCurvesCount++;
    }
    patchTrimmingCurves[patchIndex] = trimmingCurvesCount;
  }

  axom::fmt::memory_buffer out;
  axom::fmt::format_to(std::back_inserter(out), "Shape statistics:\n");
  axom::fmt::format_to(std::back_inserter(out),
                       " - Number of patches: {}\n",
                       numPatches);
  axom::fmt::format_to(std::back_inserter(out),
                       " - Number of trimming curves: {}\n",
                       numTrimmingCurves);
  axom::fmt::format_to(std::back_inserter(out), "---- \n");
  axom::fmt::format_to(std::back_inserter(out),
                       " - Trimming curves per patch:\n");
  for(const auto& entry : patchTrimmingCurves)
  {
    axom::fmt::format_to(std::back_inserter(out),
                         "\t Patch {}: {} trimming curves\n",
                         entry.first,
                         entry.second);
  }
  SLIC_INFO(axom::fmt::to_string(out));
}

std::map<int, axom::Array<axom::primal::NURBSCurve<double, 2>>>
convertTrimmingCurvesToNurbs(const TopoDS_Shape& shape)
{
  using NCurve = axom::primal::NURBSCurve<double, 2>;
  using PointType = axom::primal::Point<double, 2>;
  using BBox = axom::primal::BoundingBox<double, 2>;

  std::map<int, axom::Array<NCurve>> patchTrimmingCurves;
  int patchIndex = 0;

  for(TopExp_Explorer faceExp(shape, TopAbs_FACE); faceExp.More();
      faceExp.Next(), ++patchIndex)
  {
    // Get span of this patch in u and v directions
    BBox patchBbox;
    BBox expandedPatchBbox;
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());
      Handle(Geom_Surface) surface = BRep_Tool::Surface(face);

      Standard_Real u1, u2, v1, v2;
      surface->Bounds(u1, u2, v1, v2);
      patchBbox.addPoint(PointType {u1, v1});
      patchBbox.addPoint(PointType {u2, v2});

      expandedPatchBbox = patchBbox;
      expandedPatchBbox.scale(1. + 1e-3);

      SLIC_INFO(
        axom::fmt::format("[Patch {}]: U span [{}, {}], V span [{}, {}], BBox "
                          "{}; expanded BBox {}",
                          patchIndex,
                          u1,
                          u2,
                          v1,
                          v2,
                          patchBbox,
                          expandedPatchBbox));
    }

    axom::Array<axom::primal::NURBSCurve<double, 2>> curves;

    for(TopExp_Explorer edgeExp(faceExp.Current(), TopAbs_EDGE); edgeExp.More();
        edgeExp.Next())
    {
      TopoDS_Edge edge = TopoDS::Edge(edgeExp.Current());

      BRepAdaptor_Curve curveAdaptor(edge);
      GeomAbs_CurveType curveType = curveAdaptor.GetType();

      if(curveType == GeomAbs_BSplineCurve || curveType == GeomAbs_BezierCurve)
      {
        Standard_Real first, last;
        // Handle(Geom_Curve) geomCurve = BRep_Tool::Curve(edge, first, last);
        // Handle(Geom_BSplineCurve) bsplineCurve =
        //   Handle(Geom_BSplineCurve)::DownCast(geomCurve);

        //Handle(Geom_Surface) surface = BRep_Tool::Surface(TopoDS::Face(faceExp.Current()));

        Handle(Geom2d_Curve) parametricCurve =
          BRep_Tool::CurveOnSurface(edge,
                                    TopoDS::Face(faceExp.Current()),
                                    first,
                                    last);
        Handle(Geom2d_BSplineCurve) bsplineCurve =
          Handle(Geom2d_BSplineCurve)::DownCast(parametricCurve);
        if(!parametricCurve.IsNull() && !bsplineCurve.IsNull())
        {
          // Extract the control points in parametric space
          TColgp_Array1OfPnt2d paraPoints(1, bsplineCurve->NbPoles());
          bsplineCurve->Poles(paraPoints);

          axom::Array<PointType> controlPoints;
          for(Standard_Integer i = paraPoints.Lower(); i <= paraPoints.Upper();
              ++i)
          {
            gp_Pnt2d paraPt = paraPoints(i);
            auto pt = PointType {paraPt.X(), paraPt.Y()};
            SLIC_DEBUG_IF(
              !expandedPatchBbox.contains(pt),
              axom::fmt::format(
                "Distance of {} to {} is {}",
                pt,
                patchBbox,
                sqrt(axom::primal::squared_distance(pt, patchBbox))));
            controlPoints.emplace_back(pt);
          }

          SLIC_INFO(
            axom::fmt::format("[Patch {} Curve {}] Control Points: [{}]",
                              patchIndex,
                              curves.size(),
                              axom::fmt::join(controlPoints, ", ")));

          // Extract the weights for the control points
          TColStd_Array1OfReal weights(1, bsplineCurve->NbPoles());
          bsplineCurve->Weights(weights);

          bool isRational = false;
          axom::Array<double> weightsVector(weights.Length());
          for(int i = 1; i <= weights.Length(); ++i)
          {
            weightsVector[i - 1] = weights(i);
            if(!isRational && i > 2 && weightsVector[i - 1] != weightsVector[0])
            {
              isRational = true;
            }
          }
          SLIC_INFO(axom::fmt::format(
            "[Patch {} Curve {}] Weights: [{}]; spline {} rational",
            patchIndex,
            curves.size(),
            axom::fmt::join(weightsVector, ", "),
            isRational ? "is" : "is not"));

          // Extract the knots and their multiplicities
          TColStd_Array1OfReal knots(1, bsplineCurve->NbKnots());
          bsplineCurve->Knots(knots);

          TColStd_Array1OfInteger multiplicities(1, bsplineCurve->NbKnots());
          bsplineCurve->Multiplicities(multiplicities);

          SLIC_ASSERT(knots.Length() == multiplicities.Length());
          const int curveDegree = bsplineCurve->Degree();
          const bool isClamped = (multiplicities.First() == curveDegree + 1) &&
            (multiplicities.Last() == curveDegree + 1);

          axom::Array<double> knotVector;
          axom::Array<double> debug_multipliciesVector;  // delete me!
          for(int i = 1; i <= knots.Length(); ++i)
          {
            debug_multipliciesVector.push_back(multiplicities(i));
            for(int j = 0; j < multiplicities(i); ++j)
            {
              knotVector.push_back(knots(i));
            }
          }
          SLIC_INFO(axom::fmt::format("[Patch {} Curve {}] Knots: [{}]",
                                      patchIndex,
                                      curves.size(),
                                      axom::fmt::join(knotVector, ", ")));

          SLIC_INFO(axom::fmt::format("[Patch {} Curve {}] Degree: {}",
                                      patchIndex,
                                      curves.size(),
                                      curveDegree));

          SLIC_INFO(
            axom::fmt::format("[Patch {} Curve {}] knot multiplicities: {}, "
                              "degree: {}, is clamped: {}",
                              patchIndex,
                              curves.size(),
                              axom::fmt::join(debug_multipliciesVector, ", "),
                              curveDegree,
                              isClamped));

          if(!isClamped)
          {
            SLIC_WARNING(axom::fmt::format(
              "[Patch {} Curve {}] skipping curve -- Axom only currently "
              "supports clamped trimming curves",
              patchIndex,
              curves.size()));
            SLIC_INFO("---");
            continue;
          }
          else
          {
            SLIC_INFO("---");
          }

          if(isRational)
          {
            NCurve nurbs {controlPoints, weightsVector, knotVector};
            SLIC_ASSERT(nurbs.isValidNURBS());
            SLIC_ASSERT(nurbs.getDegree() == curveDegree);

            curves.emplace_back(nurbs);
          }
          else
          {
            NCurve nurbs {controlPoints, knotVector};
            SLIC_ASSERT(nurbs.isValidNURBS());
            SLIC_ASSERT(nurbs.getDegree() == curveDegree);

            curves.emplace_back(nurbs);
          }
          AXOM_UNUSED_VAR(curveDegree);
        }
      }
    }

    patchTrimmingCurves[patchIndex] = curves;
  }

  return patchTrimmingCurves;
}

TopoDS_Shape convertToNURBS(const TopoDS_Shape& shape)
{
  BRepBuilderAPI_NurbsConvert converter(shape);
  TopoDS_Shape nurbsShape = converter.Shape();

  if(nurbsShape.IsNull())
  {
    std::cerr << "Error: Conversion to NURBS failed." << std::endl;
  }
  else
  {
    std::cout << "Successfully converted to NURBS." << std::endl;
  }

  return nurbsShape;
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  axom::CLI::App app {"Quest Step File Example"};

  std::string filename;
  app.add_option("-f,--file", filename, "Input file")->required();

  CLI11_PARSE(app, argc, argv);

  std::cout << "Processing file: " << filename << std::endl;

  TopoDS_Shape shape = processStepFile(filename);
  if(shape.IsNull())
  {
    std::cerr << "Error: The shape is invalid or empty." << std::endl;
    return 1;
  }

  // convert to NURBS and print some stats
  auto nurbs_shape = convertToNURBS(shape);
  printShapeStats(nurbs_shape);

  const auto nurbsCurvesMap = convertTrimmingCurvesToNurbs(nurbs_shape);

  for(const auto& entry : nurbsCurvesMap)
  {
    std::cout << "Patch " << entry.first << " has " << entry.second.size()
              << " NURBS curves." << std::endl;
  }

  return 0;
}