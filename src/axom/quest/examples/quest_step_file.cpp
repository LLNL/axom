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
#include "opencascade/TopoDS_Wire.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/Poly_Triangulation.hxx"
#include "opencascade/Poly_Array1OfTriangle.hxx"
#include "opencascade/TColgp_Array1OfPnt.hxx"
#include "opencascade/TopLoc_Location.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/Precision.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/BRep_Tool.hxx"
#include "opencascade/Geom2dConvert.hxx"
#include "opencascade/TopExp.hxx"
#include "opencascade/TopTools_IndexedMapOfShape.hxx"
#include "opencascade/Geom2dAPI_ExtremaCurveCurve.hxx"
#include "opencascade/BRepAdaptor_Surface.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/BRepBndLib.hxx"
#include "opencascade/Bnd_Box.hxx"
#include "opencascade/BRepAdaptor_Surface.hxx"

#include <iostream>

struct PatchData
{
  int patchIndex {-1};
  axom::Array<double> uKnots;
  axom::Array<double> vKnots;
  axom::primal::BoundingBox<double, 2> parametricBBox;
  axom::primal::BoundingBox<double, 3> physicalBBox;
  axom::Array<axom::primal::NURBSCurve<double, 2>> trimmingCurves;
};

class StepFileProcessor
{
public:
  using PatchDataMap = std::map<int, PatchData>;

  enum class LoadStatus
  {
    UNINITIALIZED = 0,
    SUCEESS = 1 << 0,
    FAILED_TO_READ = 1 << 1,
    FAILED_NO_SHAPES = 1 << 2,
    FAILED_TO_CONVERT = 1 << 3,
    FAILED = FAILED_TO_READ | FAILED_TO_CONVERT
  };

  static constexpr int CurveDim = 2;
  using NCurve = axom::primal::NURBSCurve<double, CurveDim>;
  using PointType = axom::primal::Point<double, CurveDim>;
  using PointType2D = axom::primal::Point<double, 2>;
  using PointType3D = axom::primal::Point<double, 3>;
  using BBox2D = axom::primal::BoundingBox<double, 2>;
  using BBox3D = axom::primal::BoundingBox<double, 3>;
  using NCurveArray = axom::Array<NCurve>;
  using PatchToTrimmingCurvesMap = std::map<int, NCurveArray>;

private:
  BBox2D faceBoundingBox(const TopoDS_Face& face) const
  {
    BBox2D bbox;

    Handle(Geom_Surface) surface = BRep_Tool::Surface(face);

    // TODO: Do we need to handle closed and/or periodicity for U or V?
    Standard_Real u1, u2, v1, v2;
    surface->Bounds(u1, u2, v1, v2);
    bbox.addPoint(PointType {u1, v1});
    bbox.addPoint(PointType {u2, v2});

    return bbox;
  }

  class PatchProcessor
  {
  public:
    PatchProcessor(const TopoDS_Face& face) : m_face(face)
    {
      extractKnotVectors();
      extractBoundingBoxes();
    }

    const axom::Array<double>& getUKnots() const { return m_uKnots; }
    const axom::Array<double>& getVKnots() const { return m_vKnots; }
    const BBox2D& getParametricBBox() const { return m_parametricBBox; }
    const BBox3D& getPhysicalBBox() const { return m_physicalBBox; }

  private:
    void extractKnotVectors()
    {
      Handle(Geom_Surface) surface = BRep_Tool::Surface(m_face);
      if(surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface)))
      {
        Handle(Geom_BSplineSurface) bsplineSurface =
          Handle(Geom_BSplineSurface)::DownCast(surface);
        TColStd_Array1OfReal uKnots(1, bsplineSurface->NbUKnots());
        TColStd_Array1OfReal vKnots(1, bsplineSurface->NbVKnots());
        bsplineSurface->UKnots(uKnots);
        bsplineSurface->VKnots(vKnots);

        m_uKnots.resize(uKnots.Length());
        m_vKnots.resize(vKnots.Length());
        for(int i = 1; i <= uKnots.Length(); ++i)
        {
          m_uKnots[i - 1] = uKnots(i);
        }
        for(int i = 1; i <= vKnots.Length(); ++i)
        {
          m_vKnots[i - 1] = vKnots(i);
        }
      }
    }

    void extractBoundingBoxes()
    {
      Handle(Geom_Surface) surface = BRep_Tool::Surface(m_face);
      Standard_Real u1, u2, v1, v2;
      surface->Bounds(u1, u2, v1, v2);
      m_parametricBBox.addPoint(PointType2D {u1, v1});
      m_parametricBBox.addPoint(PointType2D {u2, v2});

      BRepAdaptor_Surface surfaceAdaptor(m_face);
      Bnd_Box bndBox;
      BRepBndLib::Add(m_face, bndBox);
      Standard_Real xMin, yMin, zMin, xMax, yMax, zMax;
      bndBox.Get(xMin, yMin, zMin, xMax, yMax, zMax);
      m_physicalBBox.addPoint(PointType3D {xMin, yMin, zMin});
      m_physicalBBox.addPoint(PointType3D {xMax, yMax, zMax});
    }

  private:
    TopoDS_Face m_face;
    axom::Array<double> m_uKnots;
    axom::Array<double> m_vKnots;
    BBox2D m_parametricBBox;
    BBox3D m_physicalBBox;
  };

  /// Helper class in support of extracting necessary information from trimming curves
  class CurveProcessor
  {
  public:
    CurveProcessor(const Handle(Geom2d_BSplineCurve) & curve) : m_curve(curve)
    { }

    const Handle(Geom2d_BSplineCurve) & getCurve() const { return m_curve; }

    // note: should only be called after extractKnots()
    bool isClamped() const { return m_isClamped; }

    void convertPeriodicToClamped()
    {
      if(!m_curve->IsPeriodic())
      {
        return;
      }

      // Copy original poles, knots, and multiplicities
      axom::Array<PointType> orig_poles(m_curve->NbPoles());
      axom::Array<double> orig_knots(m_curve->NbKnots());
      axom::Array<int> orig_mults(m_curve->NbKnots());
      int orig_knots_total = 0;
      {
        TColgp_Array1OfPnt2d originalPoles(1, m_curve->NbPoles());
        m_curve->Poles(originalPoles);
        for(int i = 1; i <= originalPoles.Length(); ++i)
        {
          orig_poles[i - 1] =
            PointType {originalPoles(i).X(), originalPoles(i).Y()};
        }

        TColStd_Array1OfReal originalKnots(1, m_curve->NbKnots());
        TColStd_Array1OfInteger originalMultiplicities(1, m_curve->NbKnots());
        m_curve->Knots(originalKnots);
        m_curve->Multiplicities(originalMultiplicities);
        for(int i = 1; i <= originalKnots.Length(); ++i)
        {
          orig_knots[i - 1] = originalKnots(i);
          orig_mults[i - 1] = originalMultiplicities(i);
          orig_knots_total += originalMultiplicities(i);
        }
      }

      // Get the value of the curve at the first and last knot
      PointType orig_first, orig_last;
      {
        gp_Pnt2d firstKnotPoint = m_curve->Value(orig_knots.front());
        orig_first = PointType {firstKnotPoint.X(), firstKnotPoint.Y()};

        gp_Pnt2d lastKnotPoint = m_curve->Value(orig_knots.back());
        orig_last = PointType {lastKnotPoint.X(), lastKnotPoint.Y()};
      }

      // Set the curve to not periodic
      m_curve->SetNotPeriodic();

      // Copy new poles, knots, and multiplicities
      axom::Array<PointType> new_poles(m_curve->NbPoles());
      axom::Array<double> new_knots(m_curve->NbKnots());
      axom::Array<int> new_mults(m_curve->NbKnots());
      int new_knots_total = 0;
      {
        TColgp_Array1OfPnt2d newPoles(1, m_curve->NbPoles());
        m_curve->Poles(newPoles);
        for(int i = 1; i <= newPoles.Length(); ++i)
        {
          new_poles[i - 1] = PointType {newPoles(i).X(), newPoles(i).Y()};
        }

        TColStd_Array1OfReal newKnots(1, m_curve->NbKnots());
        TColStd_Array1OfInteger newMultiplicities(1, m_curve->NbKnots());
        m_curve->Knots(newKnots);
        m_curve->Multiplicities(newMultiplicities);
        for(int i = 1; i <= newKnots.Length(); ++i)
        {
          new_knots[i - 1] = newKnots(i);
          new_mults[i - 1] = newMultiplicities(i);
          new_knots_total += newMultiplicities(i);
        }
      }

      // Create axom arrays for updated knots, multiplicities, and poles
      auto updated_poles = new_poles;
      auto updated_knots = orig_knots;

      axom::Array<int> updated_mults(new_mults.size() - 2);
      for(int i = 0; i < updated_mults.size(); ++i)
      {
        updated_mults[i] = new_mults[i + 1];
      }
      updated_mults.front() += new_mults.front();
      updated_mults.back() += new_mults.back();
      SLIC_ASSERT(updated_mults.size() == updated_knots.size());

      // Copy updated multiplicities, knots, and poles into OpenCascade arrays
      TColStd_Array1OfReal updatedKnots(1, updated_knots.size());
      TColStd_Array1OfInteger updatedMultiplicities(1, updated_mults.size());
      TColgp_Array1OfPnt2d updatedPoles(1, updated_poles.size());

      for(int i = 1; i <= updated_knots.size(); ++i)
      {
        updatedKnots.SetValue(i, updated_knots[i - 1]);
      }

      for(int i = 1; i <= updated_mults.size(); ++i)
      {
        updatedMultiplicities.SetValue(i, updated_mults[i - 1]);
      }

      for(int i = 1; i <= updated_poles.size(); ++i)
      {
        updatedPoles.SetValue(
          i,
          gp_Pnt2d(updated_poles[i - 1][0], updated_poles[i - 1][1]));
      }

      SLIC_INFO(
        axom::fmt::format("\nOrig knots ({}): [{}]"
                          "\nNew knots ({}): [{}]",
                          orig_knots.size(),
                          axom::fmt::join(orig_knots, ", "),
                          new_knots.size(),
                          axom::fmt::join(new_knots, ", ")));

      SLIC_INFO(
        axom::fmt::format("\nOrig multiplicities ({}): [{}]"
                          "\nNew multiplicities ({}): [{}]"
                          "\nUpdated multiplicities ({}): [{}]",
                          orig_mults.size(),
                          axom::fmt::join(orig_mults, ", "),
                          new_mults.size(),
                          axom::fmt::join(new_mults, ", "),
                          updated_mults.size(),
                          axom::fmt::join(updated_mults, ", ")));

      SLIC_INFO(
        axom::fmt::format("\nOrig poles ({}): [{}]"
                          "\nNew poles ({}): [{}]"
                          "\nUpdated poles ({}): [{}]",
                          orig_poles.size(),
                          axom::fmt::join(orig_poles, ", "),
                          new_poles.size(),
                          axom::fmt::join(new_poles, ", "),
                          updated_poles.size(),
                          axom::fmt::join(updated_poles, ", ")));

      const int degree = m_curve->Degree();
      SLIC_INFO(
        axom::fmt::format("Degree (d): {}"
                          "\nOrig curve total knots (m): {} and poles (n): {}"
                          "\nNew curve total knots (m): {} and poles (n): {}",
                          degree,
                          orig_knots_total,
                          orig_poles.size(),
                          new_knots_total,
                          new_poles.size()));

      // Save the points as axom PointType
      PointType new_first, new_last;
      {
        gp_Pnt2d firstKnotPoint = m_curve->Value(new_knots.front());
        new_first = PointType {firstKnotPoint.X(), firstKnotPoint.Y()};

        gp_Pnt2d lastKnotPoint = m_curve->Value(new_knots.back());
        new_last = PointType {lastKnotPoint.X(), lastKnotPoint.Y()};
      }

      axom::slic::flushStreams();

      // Copy updated weights into OpenCascade array
      TColStd_Array1OfReal updatedWeights(1, updated_poles.size());
      //if (m_curve->IsRational())
      {
        m_curve->Weights(updatedWeights);
      }

      // try
      // {
      Handle(Geom2d_BSplineCurve) clamped_curve =
        new Geom2d_BSplineCurve(updatedPoles,
                                updatedWeights,
                                updatedKnots,
                                updatedMultiplicities,
                                degree);

      // Save the points as axom PointType
      PointType updated_first, updated_last;
      {
        gp_Pnt2d firstKnotPoint = clamped_curve->Value(updated_knots.front());
        updated_first = PointType {firstKnotPoint.X(), firstKnotPoint.Y()};

        gp_Pnt2d lastKnotPoint = clamped_curve->Value(updated_knots.back());
        updated_last = PointType {lastKnotPoint.X(), lastKnotPoint.Y()};
      }

      SLIC_INFO(axom::fmt::format(
        "Original, new, and updated curve at first knot: {} -> {} -> {}",
        orig_first,
        new_first,
        updated_first));
      SLIC_INFO(axom::fmt::format(
        "Original, new, and updated curve at last knot: {} -> {} -> {}",
        orig_last,
        new_last,
        updated_last));

      // } catch (Standard_Failure& e) {
      //     std::cerr << "***  Error: " << e.GetMessageString() << std::endl;
      //     std::cout << "**** Error: " << e.GetMessageString() << std::endl;
      //     throw std::runtime_error(e.GetMessageString());
      // }

      // Use Geom2dAPI_ExtremaCurveCurve to find the extrema between the original and clamped curves
      {
        using RangeType = axom::primal::BoundingBox<double, 1>;
        using RangePoint = axom::primal::Point<double, 1>;

        Geom2dAPI_ExtremaCurveCurve extrema(m_curve,
                                            clamped_curve,
                                            m_curve->FirstParameter(),
                                            m_curve->LastParameter(),
                                            clamped_curve->FirstParameter(),
                                            clamped_curve->LastParameter());

        RangeType range;
        for(int i = 1; i <= extrema.NbExtrema(); ++i)
        {
          Standard_Real U1, U2;
          extrema.Parameters(i, U1, U2);
          gp_Pnt2d P1 = m_curve->Value(U1);
          gp_Pnt2d P2 = clamped_curve->Value(U2);
          SLIC_INFO(axom::fmt::format(
            "Extrema {}: Distance = {}, Location on original curve = ({}, {}), "
            "Location on clamped curve = ({}, {})",
            i,
            extrema.Distance(i),
            P1.X(),
            P1.Y(),
            P2.X(),
            P2.Y()));
          range.addPoint(RangePoint {extrema.Distance(i)});
        }

        SLIC_INFO(axom::fmt::format(
          "Distance between original and clamped curve using extrema: {}",
          range));
      }
      // Swap clamped_curve and m_curve
      m_curve = clamped_curve;
      m_isClamped = true;
    }

    // note: bounding boxes are used as debug checks that control points
    // lie within the parameter space of the parent patch
    axom::Array<PointType> extractControlPoints(const BBox2D& bbox,
                                                const BBox2D& expandedBBox) const
    {
      axom::Array<PointType> controlPoints;
      TColgp_Array1OfPnt2d paraPoints(1, m_curve->NbPoles());
      m_curve->Poles(paraPoints);

      for(Standard_Integer i = paraPoints.Lower(); i <= paraPoints.Upper(); ++i)
      {
        gp_Pnt2d paraPt = paraPoints(i);
        auto pt = PointType2D {paraPt.X(), paraPt.Y()};
        SLIC_DEBUG_IF(
          !expandedBBox.contains(pt),
          axom::fmt::format("Distance of {} to {} is {}",
                            pt,
                            bbox,
                            sqrt(axom::primal::squared_distance(pt, bbox))));
        controlPoints.emplace_back(pt);
      }

      AXOM_UNUSED_VAR(bbox);
      AXOM_UNUSED_VAR(expandedBBox);

      return controlPoints;
    }

    axom::Array<double> extractWeights() const
    {
      axom::Array<double> weights;
      if(m_curve->IsRational())
      {
        TColStd_Array1OfReal curveWeights(1, m_curve->NbPoles());
        m_curve->Weights(curveWeights);
        weights.resize(curveWeights.Length());
        for(int i = 1; i <= curveWeights.Length(); ++i)
        {
          weights[i - 1] = curveWeights(i);
        }
      }
      return weights;
    }

    axom::Array<double> extractKnots()
    {
      axom::Array<double> knots;

      TColStd_Array1OfReal curveKnots(1, m_curve->NbKnots());
      m_curve->Knots(curveKnots);

      TColStd_Array1OfInteger multiplicities(1, m_curve->NbKnots());
      m_curve->Multiplicities(multiplicities);

      SLIC_ASSERT(curveKnots.Length() == multiplicities.Length());

      const int curveDegree = m_curve->Degree();
      m_isClamped = (multiplicities.First() == curveDegree + 1) &&
        (multiplicities.Last() == curveDegree + 1);

      for(int i = 1; i <= curveKnots.Length(); ++i)
      {
        for(int j = 0; j < multiplicities(i); ++j)
        {
          knots.push_back(curveKnots(i));
        }
      }

      return knots;
    }

  private:
    Handle(Geom2d_BSplineCurve) m_curve;
    bool m_isClamped {false};
  };

public:
  StepFileProcessor() = delete;

  StepFileProcessor(const std::string& filename, bool verbose = false)
    : m_verbose(verbose)
  {
    m_shape = loadStepFile(filename);
  }

  void setVerbosity(bool verbose) { m_verbose = verbose; }

  const TopoDS_Shape& getShape() const { return m_shape; }

  bool isLoaded() const { return m_loadStatus == LoadStatus::SUCEESS; }

  void printShapeStats() const
  {
    if(m_shape.IsNull())
    {
      std::cerr << "Error: The shape is invalid or empty." << std::endl;
      return;
    }

    int numPatches = 0;
    int numTrimmingCurves = 0;

    for(TopExp_Explorer exp(m_shape, TopAbs_FACE); exp.More(); exp.Next())
    {
      numPatches++;
    }

    for(TopExp_Explorer exp(m_shape, TopAbs_EDGE); exp.More(); exp.Next())
    {
      numTrimmingCurves++;
    }
    std::map<int, int> patchTrimmingCurves;
    int patchIndex = 0;

    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
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

  void extractPatches()
  {
    int patchIndex = 0;

    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      const TopoDS_Face& face = TopoDS::Face(faceExp.Current());

      PatchProcessor patchProcessor(face);

      PatchData& patchData = m_patchData[patchIndex];
      patchData.patchIndex = patchIndex;
      patchData.uKnots = patchProcessor.getUKnots();
      patchData.vKnots = patchProcessor.getVKnots();
      patchData.parametricBBox = patchProcessor.getParametricBBox();
      patchData.physicalBBox = patchProcessor.getPhysicalBBox();
    }
  }

  void extractTrimmingCurves()
  {
    int patchIndex = 0;

    std::map<GeomAbs_CurveType, std::string> curveTypeMap = {
      {GeomAbs_Line, "Line"},
      {GeomAbs_Circle, "Circle"},
      {GeomAbs_Ellipse, "Ellipse"},
      {GeomAbs_Hyperbola, "Hyperbola"},
      {GeomAbs_Parabola, "Parabola"},
      {GeomAbs_BezierCurve, "Bezier Curve"},
      {GeomAbs_BSplineCurve, "BSpline Curve"},
      {GeomAbs_OffsetCurve, "Offset Curve"},
      {GeomAbs_OtherCurve, "Other Curve"}};

    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      const TopoDS_Face& face = TopoDS::Face(faceExp.Current());

      PatchData& patchData = m_patchData[patchIndex];

      // Get span of this patch in u and v directions
      BBox2D patchBbox = faceBoundingBox(face);
      auto expandedPatchBbox = patchBbox;
      expandedPatchBbox.scale(1. + 1e-3);

      SLIC_INFO_IF(
        m_verbose,
        axom::fmt::format(
          "[Patch {}]: BBox in parametric space: {}; expanded BBox {}",
          patchIndex,
          patchBbox,
          expandedPatchBbox));

      NCurveArray curves;
      int wireIndex = 0;
      for(TopExp_Explorer wireExp(faceExp.Current(), TopAbs_WIRE); wireExp.More();
          wireExp.Next(), ++wireIndex)
      {
        const TopoDS_Wire& wire = TopoDS::Wire(wireExp.Current());

        int edgeIndex = 0;
        for(TopExp_Explorer edgeExp(wire, TopAbs_EDGE); edgeExp.More();
            edgeExp.Next(), ++edgeIndex)
        {
          const TopoDS_Edge& edge = TopoDS::Edge(edgeExp.Current());

          TopAbs_Orientation orientation = edge.Orientation();
          const bool isReversed = orientation == TopAbs_REVERSED;

          BRepAdaptor_Curve curveAdaptor(edge);
          GeomAbs_CurveType curveType = curveAdaptor.GetType();

          std::string curveTypeStr = curveTypeMap[curveType];
          SLIC_INFO_IF(
            m_verbose,
            axom::fmt::format("[Patch {} Wire {} Edge {} Curve {}] Processing "
                              "edge with curve type: {}",
                              patchIndex,
                              wireIndex,
                              edgeIndex,
                              curves.size(),
                              curveTypeStr));

          Standard_Real first, last;
          Handle(Geom2d_Curve) parametricCurve =
            BRep_Tool::CurveOnSurface(edge,
                                      TopoDS::Face(faceExp.Current()),
                                      first,
                                      last);
          Handle(Geom2d_BSplineCurve) bsplineCurve =
            Geom2dConvert::CurveToBSplineCurve(parametricCurve);

          if(!parametricCurve.IsNull() && !bsplineCurve.IsNull())
          {
            CurveProcessor curveProcessor(bsplineCurve);
            const int curveDegree = bsplineCurve->Degree();

            curveProcessor.convertPeriodicToClamped();
            SLIC_ASSERT(!curveProcessor.getCurve()->IsPeriodic());

            auto controlPoints =
              curveProcessor.extractControlPoints(patchBbox, expandedPatchBbox);
            SLIC_INFO_IF(
              m_verbose,
              axom::fmt::format(
                "[Patch {} Wire {} Edge {} Curve {}] Control Points: [{}]",
                patchIndex,
                wireIndex,
                edgeIndex,
                curves.size(),
                axom::fmt::join(controlPoints, ", ")));

            // Check if the B-spline curve is rational
            // and extract the weights for the control points
            const bool isRational = bsplineCurve->IsRational();
            axom::Array<double> weightsVector = curveProcessor.extractWeights();

            SLIC_INFO_IF(
              m_verbose && isRational,
              axom::fmt::format("[Patch {} Wire {} Edge {} Curve {}] Weights: "
                                "[{}]",
                                patchIndex,
                                wireIndex,
                                edgeIndex,
                                curves.size(),
                                axom::fmt::join(weightsVector, ", ")));

            // Extract the knots and their multiplicities
            axom::Array<double> knotVector = curveProcessor.extractKnots();
            SLIC_INFO_IF(
              m_verbose,
              axom::fmt::format(
                "[Patch {} Wire {} Edge {} Curve {}] Degree: {}, Knots: [{}]",
                patchIndex,
                wireIndex,
                edgeIndex,
                curves.size(),
                curveDegree,
                axom::fmt::join(knotVector, ", ")));

            // Note: This should never trigger since we've enforced this above
            if(!curveProcessor.isClamped())
            {
              SLIC_WARNING(
                axom::fmt::format("[Patch {} Wire {} Edge {} Curve {}] "
                                  "skipping curve -- Axom only currently "
                                  "supports clamped trimming curves",
                                  patchIndex,
                                  wireIndex,
                                  edgeIndex,
                                  curves.size()));
              SLIC_INFO_IF(m_verbose, "---");
              continue;
            }
            else
            {
              SLIC_INFO_IF(m_verbose, "---");
            }

            // Generate a NURBSCurve and add it to the list
            isRational ? curves.emplace_back(
                           NCurve {controlPoints, weightsVector, knotVector})
                       : curves.emplace_back(NCurve {controlPoints, knotVector});

            // Ensure consistency of the curves w.r.t. the patch
            if(isReversed)
            {
              curves.back().reverseOrientation();
            }
            SLIC_ASSERT(curves.back().isValidNURBS());
            SLIC_ASSERT(curves.back().getDegree() == curveDegree);

            AXOM_UNUSED_VAR(curveDegree);
          }
        }
      }

      patchData.trimmingCurves = curves;
    }
  }

  const PatchDataMap& getPatchDataMap() const { return m_patchData; }

private:
  TopoDS_Shape loadStepFile(const std::string& filename)
  {
    STEPControl_Reader reader;

    IFSelect_ReturnStatus status = reader.ReadFile(filename.c_str());
    if(status != IFSelect_RetDone)
    {
      m_loadStatus = LoadStatus::FAILED_TO_READ;
      std::cerr << "Error: Cannot read the file." << std::endl;
      return TopoDS_Shape();
    }

    Standard_Integer numRoots = reader.NbRootsForTransfer();
    reader.TransferRoots();
    TopoDS_Shape shape = reader.OneShape();
    if(shape.IsNull())
    {
      m_loadStatus = LoadStatus::FAILED_NO_SHAPES;
      std::cerr << "Error: No shape found in the file." << std::endl;
      return TopoDS_Shape();
    }

    // convert to NURBS
    // TODO: Check if this also converts trimming curves
    BRepBuilderAPI_NurbsConvert converter(shape);
    TopoDS_Shape nurbsShape = converter.Shape();

    if(nurbsShape.IsNull())
    {
      m_loadStatus = LoadStatus::FAILED_TO_CONVERT;
      std::cerr << "Error: Conversion to NURBS failed." << std::endl;
      return TopoDS_Shape();
    }

    m_loadStatus = LoadStatus::SUCEESS;
    SLIC_INFO_IF(
      m_verbose,
      axom::fmt::format("Successfully read the STEP file with {} roots",
                        numRoots));

    // Initialize m_faceMap with faces from nurbsShape
    TopExp::MapShapes(nurbsShape, TopAbs_FACE, m_faceMap);

    return nurbsShape;
  }

private:
  TopoDS_Shape m_shape;
  bool m_verbose {false};
  LoadStatus m_loadStatus {LoadStatus::UNINITIALIZED};
  TopTools_IndexedMapOfShape m_faceMap;

  PatchDataMap m_patchData;
};

std::string nurbsCurveToSVGPath(const axom::primal::NURBSCurve<double, 2>& curve)
{
  using PointType = axom::primal::Point<double, 2>;

  const int degree = curve.getDegree();
  const auto& knotVector = curve.getKnots();
  const bool isRational = curve.isRational();

  axom::fmt::memory_buffer svgPath;
  axom::fmt::format_to(std::back_inserter(svgPath),
                       "<path class=\"{} degree-{}\" d=\"",
                       isRational ? "rational" : "non-rational",
                       degree);

  if(curve.isRational() || degree > 3)
  {
    const int numSamples = 100;
    const double tMin = knotVector[0];
    const double tMax = knotVector[knotVector.getNumKnots() - 1];

    for(int i = 0; i <= numSamples; ++i)
    {
      const double t =
        axom::utilities::lerp(tMin, tMax, static_cast<double>(i) / numSamples);

      PointType pt = curve.evaluate(t);
      if(i == 0)
      {
        axom::fmt::format_to(std::back_inserter(svgPath), "M {} {} ", pt[0], pt[1]);
      }
      else
      {
        axom::fmt::format_to(std::back_inserter(svgPath), "L {} {} ", pt[0], pt[1]);
      }
    }
  }
  else
  {
    auto bezierCurves = curve.extractBezier();
    for(const auto& bezier : bezierCurves)
    {
      const auto& bezierControlPoints = bezier.getControlPoints();
      if(degree == 2)
      {
        for(int i = 0; i < bezierControlPoints.size(); ++i)
        {
          const PointType& pt = bezierControlPoints[i];
          if(i == 0)
          {
            axom::fmt::format_to(std::back_inserter(svgPath),
                                 "M {} {} ",
                                 pt[0],
                                 pt[1]);
          }
          else if(i == 2)
          {
            axom::fmt::format_to(std::back_inserter(svgPath),
                                 "Q {} {} {} {} ",
                                 bezierControlPoints[1][0],
                                 bezierControlPoints[1][1],
                                 pt[0],
                                 pt[1]);
          }
        }
      }
      else if(degree == 3)
      {
        for(int i = 0; i < bezierControlPoints.size(); ++i)
        {
          const PointType& pt = bezierControlPoints[i];
          if(i == 0)
          {
            axom::fmt::format_to(std::back_inserter(svgPath),
                                 "M {} {} ",
                                 pt[0],
                                 pt[1]);
          }
          else if(i == 3)
          {
            axom::fmt::format_to(std::back_inserter(svgPath),
                                 "C {} {} {} {} {} {} ",
                                 bezierControlPoints[1][0],
                                 bezierControlPoints[1][1],
                                 bezierControlPoints[2][0],
                                 bezierControlPoints[2][1],
                                 pt[0],
                                 pt[1]);
          }
        }
      }
      else
      {
        for(int i = 0; i < bezierControlPoints.size(); ++i)
        {
          const PointType& pt = bezierControlPoints[i];
          if(i == 0)
          {
            axom::fmt::format_to(std::back_inserter(svgPath),
                                 "M {} {} ",
                                 pt[0],
                                 pt[1]);
          }
          else
          {
            axom::fmt::format_to(std::back_inserter(svgPath),
                                 "L {} {} ",
                                 pt[0],
                                 pt[1]);
          }
        }
      }
    }
  }

  axom::fmt::format_to(std::back_inserter(svgPath), "\"  />");

  return axom::fmt::to_string(svgPath);
}

void generateSVGForPatch(int patchIndex, const PatchData& patchData)
{
  using PointType = axom::primal::Point<double, 2>;

  const auto& parametricBBox = patchData.parametricBBox;

  SLIC_INFO(axom::fmt::format("Parametric BBox for patch {}: {}",
                              patchIndex,
                              parametricBBox));

  const auto& curves = patchData.trimmingCurves;
  axom::fmt::memory_buffer svgContent;

  // Create a new bounding box by scaling and translating the parametricBBox by a factor of two, keeping the center the same
  PointType center = parametricBBox.getCentroid();
  PointType min = parametricBBox.getMin();
  PointType max = parametricBBox.getMax();

  PointType newMin = center + 2.0 * (min - center);
  PointType newMax = center + 2.0 * (max - center);

  axom::primal::BoundingBox<double, 2> scaledParametricBBox;
  scaledParametricBBox.addPoint(newMin);
  scaledParametricBBox.addPoint(newMax);

  SLIC_INFO(
    axom::fmt::format("Scaled and translated parametric BBox for patch {}: {}",
                      patchIndex,
                      scaledParametricBBox));
  axom::fmt::format_to(
    std::back_inserter(svgContent),
    "<svg xmlns=\"http://www.w3.org/2000/svg\" "
    "version=\"1.1\" viewBox=\"{} {} {} {}\">\n",
    scaledParametricBBox.getMin()[0],
    scaledParametricBBox.getMin()[1],
    scaledParametricBBox.getMax()[0] - scaledParametricBBox.getMin()[0],
    scaledParametricBBox.getMax()[1] - scaledParametricBBox.getMin()[1]);

  axom::fmt::format_to(std::back_inserter(svgContent), R"raw(
  <style>
    path {{ fill:none; stroke:black; stroke-width:.1; marker-end:url(#arrow);paint-order:fill stroke markers;stroke-linejoin:round;stroke-linecap:round; }}
    rect {{ fill: white; stroke: gray; stroke-width: 0.1; }}
  </style>
  )raw");
  // add a marker for the arrow's head to indicate the orientation
  axom::fmt::format_to(std::back_inserter(svgContent), R"raw(
  <defs>
    <marker id="arrow" style="overflow:visible" orient="auto-start-reverse"
        refX="0" refY="0"
        markerWidth="3.3239999" markerHeight="3.8427744"
        viewBox="0 0 5.3244081 6.1553851">
      <path
          transform="scale(0.5)"
          style="fill:context-stroke;fill-rule:evenodd;stroke:none"
          d="M 5.77,0 -2.88,5 V -5 Z" />
    </marker>
  </defs>
  )raw");

  axom::fmt::format_to(
    std::back_inserter(svgContent),
    "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" />\n",
    parametricBBox.getMin()[0],
    parametricBBox.getMin()[1],
    parametricBBox.getMax()[0] - parametricBBox.getMin()[0],
    parametricBBox.getMax()[1] - parametricBBox.getMin()[1]);

  for(const auto& curve : curves)
  {
    std::string pathData = nurbsCurveToSVGPath(curve);
    axom::fmt::format_to(std::back_inserter(svgContent), "{}\n", pathData);
  }

  axom::fmt::format_to(std::back_inserter(svgContent), "</svg>");

  std::string svgFilename = axom::fmt::format("patch_{}.svg", patchIndex);
  std::ofstream svgFile(svgFilename);
  if(svgFile.is_open())
  {
    svgFile << axom::fmt::to_string(svgContent);
    svgFile.close();
    std::cout << "SVG file generated: " << svgFilename << std::endl;
  }
  else
  {
    std::cerr << "Error: Unable to open file " << svgFilename << " for writing."
              << std::endl;
  }
}

void triangulateTrimmedPatchAndWriteSTL(const TopoDS_Shape& shape, int patchIndex)
{
  // Get the patch (face) associated with patchIndex
  TopExp_Explorer faceExp(shape, TopAbs_FACE);
  for(int i = 0; i < patchIndex && faceExp.More(); ++i)
  {
    faceExp.Next();
  }
  if(!faceExp.More())
  {
    SLIC_WARNING(
      axom::fmt::format("Error: Patch index {} is out of range.", patchIndex));
    return;
  }
  TopoDS_Face face = TopoDS::Face(faceExp.Current());

  // Create a triangulation of this patch and write its triangles in the STL format
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, loc);

  if(triangulation.IsNull())
  {
    SLIC_WARNING(
      axom::fmt::format("Error: STL could not be generated for patch {}",
                        patchIndex));
    return;
  }

  axom::fmt::memory_buffer stlContent;
  axom::fmt::format_to(std::back_inserter(stlContent),
                       "solid patch_{}\n",
                       patchIndex);

  const int numTriangles = triangulation->NbTriangles();

  for(int i = 1; i <= numTriangles; ++i)
  {
    Poly_Triangle triangle = triangulation->Triangle(i);
    int n1, n2, n3;
    triangle.Get(n1, n2, n3);

    gp_Pnt p1 = triangulation->Node(n1);
    gp_Pnt p2 = triangulation->Node(n2);
    gp_Pnt p3 = triangulation->Node(n3);

    gp_Vec v1(p1, p2);
    gp_Vec v2(p1, p3);
    gp_Vec normal = v1.Crossed(v2);
    normal.Normalize();

    axom::fmt::format_to(std::back_inserter(stlContent),
                         "facet normal {} {} {}\n",
                         normal.X(),
                         normal.Y(),
                         normal.Z());
    axom::fmt::format_to(std::back_inserter(stlContent), "  outer loop\n");
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "    vertex {} {} {}\n",
                         p1.X(),
                         p1.Y(),
                         p1.Z());
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "    vertex {} {} {}\n",
                         p2.X(),
                         p2.Y(),
                         p2.Z());
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "    vertex {} {} {}\n",
                         p3.X(),
                         p3.Y(),
                         p3.Z());
    axom::fmt::format_to(std::back_inserter(stlContent), "  endloop\n");
    axom::fmt::format_to(std::back_inserter(stlContent), "endfacet\n");
  }

  axom::fmt::format_to(std::back_inserter(stlContent),
                       "endsolid patch_{}\n",
                       patchIndex);

  // write the STL file
  std::string stlFilename = axom::fmt::format("patch_{}.stl", patchIndex);
  {
    std::ofstream stlFile(stlFilename);
    if(!stlFile.is_open())
    {
      std::cerr << "Error: Unable to open file " << stlFilename
                << " for writing." << std::endl;
      return;
    }
    stlFile << axom::fmt::to_string(stlContent);
    stlFile.close();
  }
  std::cout << "STL file generated: " << stlFilename << std::endl;
}

void triangulateUntrimmedPatchAndWriteSTL(const TopoDS_Shape& shape,
                                          int patchIndex)
{
  // Get the patch (face) associated with patchIndex
  TopExp_Explorer faceExp(shape, TopAbs_FACE);
  for(int i = 0; i < patchIndex && faceExp.More(); ++i)
  {
    faceExp.Next();
  }
  if(!faceExp.More())
  {
    SLIC_WARNING(
      axom::fmt::format("Error: Patch index {} is out of range.", patchIndex));
    return;
  }
  TopoDS_Face face = TopoDS::Face(faceExp.Current());

  // Get the underlying surface of the face
  Handle(Geom_Surface) surface = BRep_Tool::Surface(face);

  // Optionally, you can create a rectangular trimmed surface if needed
  // Here, we assume the surface is already suitable for creating a new face
  Standard_Real u1, u2, v1, v2;
  surface->Bounds(u1, u2, v1, v2);
  Handle(Geom_RectangularTrimmedSurface) untrimmedSurface =
    new Geom_RectangularTrimmedSurface(surface, u1, u2, v1, v2);

  // Create a new face from the untrimmed surface
  TopoDS_Face newFace =
    BRepBuilderAPI_MakeFace(untrimmedSurface, Precision::Confusion());

  // Mesh the new face
  Standard_Real deflection = 0.1;
  Standard_Real angle = 0.5;
  BRepMesh_IncrementalMesh mesh(newFace, deflection, Standard_False, angle);

  // Now you can access the triangulation of the new face
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation =
    BRep_Tool::Triangulation(newFace, loc);

  if(triangulation.IsNull())
  {
    SLIC_WARNING(
      axom::fmt::format("Error: STL could not be generated for patch {}",
                        patchIndex));
    return;
  }

  axom::fmt::memory_buffer stlContent;
  axom::fmt::format_to(std::back_inserter(stlContent),
                       "solid patch_{}\n",
                       patchIndex);

  const int numTriangles = triangulation->NbTriangles();

  for(int i = 1; i <= numTriangles; ++i)
  {
    Poly_Triangle triangle = triangulation->Triangle(i);
    int n1, n2, n3;
    triangle.Get(n1, n2, n3);

    gp_Pnt p1 = triangulation->Node(n1);
    gp_Pnt p2 = triangulation->Node(n2);
    gp_Pnt p3 = triangulation->Node(n3);

    gp_Vec v1(p1, p2);
    gp_Vec v2(p1, p3);
    gp_Vec normal = v1.Crossed(v2);
    normal.Normalize();

    axom::fmt::format_to(std::back_inserter(stlContent),
                         "facet normal {} {} {}\n",
                         normal.X(),
                         normal.Y(),
                         normal.Z());
    axom::fmt::format_to(std::back_inserter(stlContent), "  outer loop\n");
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "    vertex {} {} {}\n",
                         p1.X(),
                         p1.Y(),
                         p1.Z());
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "    vertex {} {} {}\n",
                         p2.X(),
                         p2.Y(),
                         p2.Z());
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "    vertex {} {} {}\n",
                         p3.X(),
                         p3.Y(),
                         p3.Z());
    axom::fmt::format_to(std::back_inserter(stlContent), "  endloop\n");
    axom::fmt::format_to(std::back_inserter(stlContent), "endfacet\n");
  }

  axom::fmt::format_to(std::back_inserter(stlContent),
                       "endsolid patch_{}\n",
                       patchIndex);

  // write the STL file
  std::string stlFilename =
    axom::fmt::format("patch_untrimmed_{}.stl", patchIndex);
  {
    std::ofstream stlFile(stlFilename);
    if(!stlFile.is_open())
    {
      std::cerr << "Error: Unable to open file " << stlFilename
                << " for writing." << std::endl;
      return;
    }
    stlFile << axom::fmt::to_string(stlContent);
    stlFile.close();
  }
  std::cout << "STL file generated: " << stlFilename << std::endl;
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  axom::CLI::App app {"Quest Step File Example"};

  std::string filename;
  app.add_option("-f,--file", filename, "Input file")->required();

  bool verbosity {false};
  app.add_flag("-v,--verbose", verbosity, "Enable verbose output");

  CLI11_PARSE(app, argc, argv);

  std::cout << "Processing file: " << filename << std::endl;

  StepFileProcessor stepProcessor(filename, verbosity);
  if(!stepProcessor.isLoaded())
  {
    std::cerr << "Error: The shape is invalid or empty." << std::endl;
    return 1;
  }

  stepProcessor.printShapeStats();

  stepProcessor.extractPatches();
  stepProcessor.extractTrimmingCurves();

  for(const auto& entry : stepProcessor.getPatchDataMap())
  {
    std::cout << "Patch " << entry.first << " has "
              << entry.second.trimmingCurves.size() << " NURBS curves."
              << std::endl;
  }

  std::string cwd = axom::utilities::filesystem::getCWD();
  std::cout << "Current working directory: " << cwd << std::endl;

  for(const auto& entry : stepProcessor.getPatchDataMap())
  {
    generateSVGForPatch(entry.first, entry.second);
  }

  auto& nurbs_shape = stepProcessor.getShape();
  Standard_Real deflection = 0.1;
  Standard_Real angle = 0.5;
  BRepMesh_IncrementalMesh mesh(nurbs_shape, deflection, Standard_False, angle);

  if(!mesh.IsDone())
  {
    throw std::runtime_error("Mesh generation failed.");
  }

  for(const auto& entry : stepProcessor.getPatchDataMap())
  {
    triangulateTrimmedPatchAndWriteSTL(nurbs_shape, entry.first);
    triangulateUntrimmedPatchAndWriteSTL(nurbs_shape, entry.first);
  }

  return 0;
}