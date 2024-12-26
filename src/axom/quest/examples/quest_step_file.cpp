// Copyright (c) 2017-2024, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "axom/config.hpp"
#include "axom/core.hpp"
#include "axom/slic.hpp"
#include "axom/primal.hpp"
#include "axom/mint.hpp"

#include "axom/CLI11.hpp"
#include "axom/fmt.hpp"

#include "opencascade/BRepAdaptor_Curve.hxx"
#include "opencascade/BRepAdaptor_Surface.hxx"
#include "opencascade/BRepBndLib.hxx"
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_NurbsConvert.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/BRepTools.hxx"
#include "opencascade/BRep_Tool.hxx"
#include "opencascade/Bnd_Box.hxx"
#include "opencascade/Geom2dAPI_ExtremaCurveCurve.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/Geom2d_Curve.hxx"
#include "opencascade/Geom2dConvert.hxx"
#include "opencascade/GeomAbs_CurveType.hxx"
#include "opencascade/Geom_BSplineCurve.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_Curve.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/Poly_Array1OfTriangle.hxx"
#include "opencascade/Poly_Triangulation.hxx"
#include "opencascade/Precision.hxx"
#include "opencascade/STEPControl_Reader.hxx"
#include "opencascade/TColgp_Array1OfPnt.hxx"
#include "opencascade/TColgp_Array2OfPnt.hxx"
#include "opencascade/TopAbs.hxx"
#include "opencascade/TopExp.hxx"
#include "opencascade/TopExp_Explorer.hxx"
#include "opencascade/TopLoc_Location.hxx"
#include "opencascade/TopTools_IndexedMapOfShape.hxx"
#include "opencascade/TopoDS.hxx"
#include "opencascade/TopoDS_Edge.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include "opencascade/TopoDS_Wire.hxx"

#include <iostream>

struct PatchData
{
  int patchIndex {-1};
  bool wasOriginallyPeriodic_u {false};
  bool wasOriginallyPeriodic_v {false};
  axom::primal::NURBSPatch<double, 3> nurbsPatch;
  axom::primal::BoundingBox<double, 2> parametricBBox;
  axom::primal::BoundingBox<double, 3> physicalBBox;
  axom::Array<axom::primal::NURBSCurve<double, 2>> trimmingCurves;
};

/**
 * Class to read in STEP files representing trimmed NURBS meshes using OpenCASCADE 
 * and convert the patches and trimming curves to axom's NURBSPatch and NURBSCurves.
 * 
 * Since axom's primitives do not support periodic knots, we must convert the OpenCASCADE
 * analogues to the open/clamped representation, where appropriate.
 */
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
  static constexpr int SpaceDim = 3;
  using NCurve = axom::primal::NURBSCurve<double, CurveDim>;
  using NPatch = axom::primal::NURBSPatch<double, SpaceDim>;
  using PointType = axom::primal::Point<double, CurveDim>;
  using PointType2D = axom::primal::Point<double, CurveDim>;
  using PointType3D = axom::primal::Point<double, SpaceDim>;
  using BBox2D = axom::primal::BoundingBox<double, CurveDim>;
  using BBox3D = axom::primal::BoundingBox<double, SpaceDim>;
  using NCurveArray = axom::Array<NCurve>;
  using PatchToTrimmingCurvesMap = std::map<int, NCurveArray>;

private:
  /// Returns a bounding box convering the patch's knot spans in 2D parametric space
  BBox2D faceBoundingBox(const TopoDS_Face& face) const
  {
    BBox2D bbox;

    Handle(Geom_Surface) surface = BRep_Tool::Surface(face);

    Standard_Real u1, u2, v1, v2;
    surface->Bounds(u1, u2, v1, v2);
    bbox.addPoint(PointType {u1, v1});
    bbox.addPoint(PointType {u2, v2});

    return bbox;
  }

  /// Helper class to convert faces to valid NURBSPatch instances
  /// The constructor converts the surface to a clamped (non-periodic) representation, if necessary
  class PatchProcessor
  {
  public:
    PatchProcessor() = delete;

    PatchProcessor(const Handle(Geom_BSplineSurface) & surface,
                   bool verbose = false)
      : m_surface(surface)
      , m_verbose(verbose)
    {
      m_inputSurfaceWasPeriodic_u = m_surface->IsUPeriodic();
      m_inputSurfaceWasPeriodic_v = m_surface->IsVPeriodic();

      ensureClamped();
    }

    const Handle(Geom_BSplineSurface) & getSurface() const { return m_surface; }

    /// Returns a representation of the surface as an axom::primal::NURBSPatch
    NPatch nurbsPatch() const
    {
      // Check if the surface is periodic in u or v
      const bool isUPeriodic = m_surface->IsUPeriodic();
      const bool isVPeriodic = m_surface->IsVPeriodic();
      SLIC_ERROR_IF(isUPeriodic || isVPeriodic,
                    "Axom's NURBSPatch only supports non-periodic patches");

      // Extract weights, if the surface is rational
      const bool isRational =
        m_surface->IsURational() || m_surface->IsVRational();

      // Create the NURBSPatch from control points, weights, and knots
      return isRational ? NPatch(extractControlPoints(),
                                 extractWeights(),
                                 extractCombinedKnots_u(),
                                 extractCombinedKnots_v())
                        : NPatch(extractControlPoints(),
                                 extractCombinedKnots_u(),
                                 extractCombinedKnots_v());
    }

    bool patchWasOriginallyPeriodic_u() const
    {
      return m_inputSurfaceWasPeriodic_u;
    }
    bool patchWasOriginallyPeriodic_v() const
    {
      return m_inputSurfaceWasPeriodic_v;
    }

    /**
     * Utility function to compare the represented curve to \a otherSurface
     * at \a numSamples uniformly sampled points in parameter space
     * Returns true when the sum of distances is less than the provided tolerance
     */
    bool compareToSurface(Handle(Geom_BSplineSurface) otherSurface,
                          int numSamples,
                          double sq_tol = 1e-8) const
    {
      if(numSamples <= 1)
      {
        return true;
      }

      auto knot_vals_u = extractKnotValues_u();
      auto knot_vals_v = extractKnotValues_v();

      axom::Array<double> params_u(numSamples);
      axom::Array<double> params_v(numSamples);

      auto evaluateSurface = [](auto surface, double u, double v) {
        gp_Pnt point;
        surface->D0(u, v, point);
        return PointType3D {point.X(), point.Y(), point.Z()};
      };

      double squared_sum = 0.;
      for(auto u : params_u)
      {
        for(auto v : params_v)
        {
          auto my_val = evaluateSurface(m_surface, u, v);
          auto other_val = evaluateSurface(otherSurface, u, v);

          const double sq_dist =
            axom::primal::squared_distance(my_val, other_val);
          squared_sum += sq_dist;
          SLIC_WARNING_IF(
            m_verbose && sq_dist > sq_tol,
            axom::fmt::format("Distance between surfaces at evaluated param "
                              "({},{}) exceeded tolerance {}.\n"
                              "Point on my surface {}; Point on other surface "
                              "{}; Squared distance {}  (running sum {})",
                              u,
                              v,
                              sq_tol,
                              my_val,
                              other_val,
                              sq_dist,
                              squared_sum));
        }
      }

      return squared_sum <= sq_tol;
    }

    void printSurfaceStats() const
    {
      const bool isUPeriodic = m_surface->IsUPeriodic();
      const bool isVPeriodic = m_surface->IsVPeriodic();

      SLIC_INFO(axom::fmt::format("Patch is periodic in u: {}", isUPeriodic));
      SLIC_INFO(axom::fmt::format("Patch is periodic in v: {}", isVPeriodic));

      const auto patch_control_points = extractControlPoints();
      SLIC_INFO(axom::fmt::format("Patch control points ({} x {}): {}",
                                  patch_control_points.shape()[0],
                                  patch_control_points.shape()[1],
                                  patch_control_points));

      const bool isRational =
        m_surface->IsURational() || m_surface->IsVRational();

      const auto patch_weights =
        isRational ? extractWeights() : axom::Array<double, 2>();
      if(isRational)
      {
        SLIC_INFO(axom::fmt::format("Patch weights ({} x {}): {}",
                                    patch_weights.shape()[0],
                                    patch_weights.shape()[1],
                                    patch_weights));
      }
      else
      {
        SLIC_INFO("Patch is polynomial (uniform weights)");
      }
    }

  private:
    /// converts the surface from periodic knots to clamped knots, when necessary
    void ensureClamped()
    {
      const bool isUPeriodic = m_surface->IsUPeriodic();
      const bool isVPeriodic = m_surface->IsVPeriodic();
      if(!isUPeriodic && !isVPeriodic)
      {
        return;  // nothing to do, return
      }

      const bool isRational =
        m_surface->IsURational() || m_surface->IsVRational();

      const int uDegree = m_surface->UDegree();
      const int vDegree = m_surface->VDegree();
      SLIC_INFO_IF(m_verbose, axom::fmt::format("  --> U degree: {}", uDegree));
      SLIC_INFO_IF(m_verbose, axom::fmt::format("  --> V degree: {}", vDegree));

      if(isUPeriodic)
      {
        m_surface->SetUNotPeriodic();

        // Modify knots and mults to ensure proper clamping
        // We remove the first and last knot and increase the multiplicity of the second and second-to-last knots
        auto mod_knots_u = extractKnotValues_u();
        auto mod_mults_u = extractKnotMultiplicities_u();
        SLIC_ASSERT(mod_knots_u.size() >= 4);
        SLIC_ASSERT(mod_mults_u.size() == mod_knots_u.size());

        const int last_idx = mod_knots_u.size() - 1;
        mod_mults_u[1] += mod_mults_u[0];
        mod_mults_u[last_idx - 1] += mod_mults_u[last_idx];
        mod_knots_u.erase(mod_knots_u.end() - 1);
        mod_knots_u.erase(mod_knots_u.begin());
        mod_mults_u.erase(mod_mults_u.end() - 1);
        mod_mults_u.erase(mod_mults_u.begin());

        try
        {
          // Create a new BSpline surface from the modified knots and control points
          m_surface = createBSplineSurfaceFromAxomArrays(
            extractControlPoints(),
            isRational ? extractWeights() : axom::Array<double, 2>(),
            mod_knots_u,
            mod_mults_u,
            extractKnotValues_v(),
            extractKnotMultiplicities_v(),
            false,
            isVPeriodic,
            uDegree,
            vDegree);
        }
        catch(const Standard_Failure& e)
        {
          std::cerr << "Error: " << e.GetMessageString() << std::endl;
          throw;
        }
      }

      if(isVPeriodic)
      {
        m_surface->SetVNotPeriodic();

        // Modify knots and mults to ensure proper clamping
        // We remove the first and last knot and increase the multiplicity of the second and second-to-last knots
        auto mod_knots_v = extractKnotValues_v();
        auto mod_mults_v = extractKnotMultiplicities_v();
        SLIC_ASSERT(mod_knots_v.size() >= 4);
        SLIC_ASSERT(mod_mults_v.size() == mod_knots_v.size());

        const int last_idx = mod_knots_v.size() - 1;
        mod_mults_v[1] += mod_mults_v[0];
        mod_mults_v[last_idx - 1] += mod_mults_v[last_idx];
        mod_knots_v.erase(mod_knots_v.end() - 1);
        mod_knots_v.erase(mod_knots_v.begin());
        mod_mults_v.erase(mod_mults_v.end() - 1);
        mod_mults_v.erase(mod_mults_v.begin());

        try
        {
          // Create a new BSpline surface from the modified knots and control points
          m_surface = createBSplineSurfaceFromAxomArrays(
            extractControlPoints(),
            isRational ? extractWeights() : axom::Array<double, 2>(),
            extractKnotValues_u(),
            extractKnotMultiplicities_u(),
            mod_knots_v,
            mod_mults_v,
            false,
            false,
            uDegree,
            vDegree);
        }
        catch(const Standard_Failure& e)
        {
          std::cerr << "Error: " << e.GetMessageString() << std::endl;
          throw;
        }
      }
    }

    /// generates a new BSpline surface from data in axom::Arrays
    Handle(Geom_BSplineSurface) createBSplineSurfaceFromAxomArrays(
      const axom::Array<PointType3D, 2>& controlPoints,
      const axom::Array<double, 2>& weights,  // will be empty if surface is non-rational
      const axom::Array<double>& uKnots,
      const axom::Array<int>& uMults,
      const axom::Array<double>& vKnots,
      const axom::Array<int>& vMults,
      bool isUPeriodic,
      bool isVPeriodic,
      int uDegree,
      int vDegree)
    {
      // Convert control points to OpenCascade array
      TColgp_Array2OfPnt poles(1,
                               controlPoints.shape()[0],
                               1,
                               controlPoints.shape()[1]);
      for(int i = 0; i < controlPoints.shape()[0]; ++i)
      {
        for(int j = 0; j < controlPoints.shape()[1]; ++j)
        {
          const PointType3D& pt = controlPoints(i, j);
          poles.SetValue(i + 1, j + 1, gp_Pnt(pt[0], pt[1], pt[2]));
        }
      }

      // Convert uKnots and uMults to OpenCascade arrays
      TColStd_Array1OfReal occUKnots(1, uKnots.size());
      TColStd_Array1OfInteger occUMults(1, uMults.size());
      for(int i = 0; i < uKnots.size(); ++i)
      {
        occUKnots.SetValue(i + 1, uKnots[i]);
        occUMults.SetValue(i + 1, uMults[i]);
      }

      // Convert vKnots and vMults to OpenCascade arrays
      TColStd_Array1OfReal occVKnots(1, vKnots.size());
      TColStd_Array1OfInteger occVMults(1, vMults.size());
      for(int i = 0; i < vKnots.size(); ++i)
      {
        occVKnots.SetValue(i + 1, vKnots[i]);
        occVMults.SetValue(i + 1, vMults[i]);
      }

      // Create and return the BSpline surface
      Handle(Geom_BSplineSurface) bsplineSurface;
      if(!weights.empty())
      {
        TColStd_Array2OfReal occWeights(1,
                                        controlPoints.shape()[0],
                                        1,
                                        controlPoints.shape()[1]);
        for(int i = 0; i < controlPoints.shape()[0]; ++i)
        {
          for(int j = 0; j < controlPoints.shape()[1]; ++j)
          {
            occWeights.SetValue(i + 1, j + 1, weights(i, j));
          }
        }
        bsplineSurface = new Geom_BSplineSurface(poles,
                                                 occWeights,
                                                 occUKnots,
                                                 occVKnots,
                                                 occUMults,
                                                 occVMults,
                                                 uDegree,
                                                 vDegree,
                                                 isUPeriodic,
                                                 isVPeriodic);
      }
      else
      {
        bsplineSurface = new Geom_BSplineSurface(poles,
                                                 occUKnots,
                                                 occVKnots,
                                                 occUMults,
                                                 occVMults,
                                                 uDegree,
                                                 vDegree,
                                                 isUPeriodic,
                                                 isVPeriodic);
      }

      return bsplineSurface;
    }

    /// extracts control points (poles) from the patch as a 2D axom::Array
    axom::Array<PointType3D, 2> extractControlPoints() const
    {
      axom::Array<PointType3D, 2> patch_control_points;

      TColgp_Array2OfPnt poles(1, m_surface->NbUPoles(), 1, m_surface->NbVPoles());
      m_surface->Poles(poles);

      patch_control_points.resize(axom::ArrayOptions::Uninitialized {},
                                  poles.ColLength(),
                                  poles.RowLength());
      for(int i = poles.LowerRow(); i <= poles.UpperRow(); ++i)
      {
        for(int j = poles.LowerCol(); j <= poles.UpperCol(); ++j)
        {
          gp_Pnt pole = poles(i, j);
          patch_control_points(i - 1, j - 1) =
            PointType3D {pole.X(), pole.Y(), pole.Z()};
        }
      }

      return patch_control_points;
    }

    /// extracts weights from the patch as a 2D axom::Array,
    /// each weight corresponds to a control point
    axom::Array<double, 2> extractWeights() const
    {
      axom::Array<double, 2> patch_weights;

      TColStd_Array2OfReal weights(1,
                                   m_surface->NbUPoles(),
                                   1,
                                   m_surface->NbVPoles());
      m_surface->Weights(weights);

      patch_weights.resize(axom::ArrayOptions::Uninitialized {},
                           weights.ColLength(),
                           weights.RowLength());
      for(int i = weights.LowerRow(); i <= weights.UpperRow(); ++i)
      {
        for(int j = weights.LowerCol(); j <= weights.UpperCol(); ++j)
        {
          patch_weights(i - 1, j - 1) = weights(i, j);
        }
      }

      return patch_weights;
    }

    /// extracts the u knot vector from the patch, accounting for multiplicities
    axom::Array<double> extractCombinedKnots_u() const
    {
      const auto vals = extractKnotValues_u();
      const auto mults = extractKnotMultiplicities_u();
      const int total_knots = std::accumulate(mults.begin(), mults.end(), 0);

      axom::Array<double> knots_u(0, total_knots);
      for(int i = 0; i < vals.size(); ++i)
      {
        knots_u.insert(knots_u.end(), mults[i], vals[i]);
      }
      SLIC_ASSERT(knots_u.size() == total_knots);

      return knots_u;
    }

    /// extracts the v knot vector from the patch, accounting for multiplicities
    axom::Array<double> extractCombinedKnots_v() const
    {
      const auto vals = extractKnotValues_v();
      const auto mults = extractKnotMultiplicities_v();
      const int total_knots = std::accumulate(mults.begin(), mults.end(), 0);

      axom::Array<double> knots_v(0, total_knots);
      for(int i = 0; i < vals.size(); ++i)
      {
        knots_v.insert(knots_v.end(), mults[i], vals[i]);
      }
      SLIC_ASSERT(knots_v.size() == total_knots);

      return knots_v;
    }

    /// converts u knot values to axom Array w/o accounting for multiplicity
    /// /sa extractKnotMultiplicities_u
    axom::Array<double> extractKnotValues_u() const
    {
      const int num_knots = m_surface->NbUKnots();
      axom::Array<double> uKnotsArray(0, num_knots);

      TColStd_Array1OfReal uKnots(1, num_knots);
      m_surface->UKnots(uKnots);

      for(int i = uKnots.Lower(); i <= uKnots.Upper(); ++i)
      {
        uKnotsArray.push_back(uKnots(i));
      }

      return uKnotsArray;
    }

    /// converts u knot multiplicities to axom Array
    /// /sa extractKnotValues_u
    axom::Array<int> extractKnotMultiplicities_u() const
    {
      const int num_knots = m_surface->NbUKnots();
      axom::Array<int> uMultsArray(0, num_knots);

      TColStd_Array1OfInteger uMults(1, num_knots);
      m_surface->UMultiplicities(uMults);

      for(int i = uMults.Lower(); i <= uMults.Upper(); ++i)
      {
        uMultsArray.push_back(uMults(i));
      }

      return uMultsArray;
    }

    /// converts v knot values to axom Array w/o accounting for multiplicity
    /// /sa extractKnotMultiplicities_v
    axom::Array<double> extractKnotValues_v() const
    {
      const int num_knots = m_surface->NbVKnots();
      axom::Array<double> vKnotsArray(0, num_knots);

      TColStd_Array1OfReal vKnots(1, num_knots);
      m_surface->VKnots(vKnots);

      for(int i = vKnots.Lower(); i <= vKnots.Upper(); ++i)
      {
        vKnotsArray.push_back(vKnots(i));
      }

      return vKnotsArray;
    }

    /// converts v knot multiplicities to axom Array
    /// /sa extractKnotValues_v
    axom::Array<int> extractKnotMultiplicities_v() const
    {
      const int num_knots = m_surface->NbVKnots();
      axom::Array<int> vMultsArray(0, num_knots);

      TColStd_Array1OfInteger vMults(1, num_knots);
      m_surface->VMultiplicities(vMults);

      for(int i = vMults.Lower(); i <= vMults.Upper(); ++i)
      {
        vMultsArray.push_back(vMults(i));
      }

      return vMultsArray;
    }

  private:
    Handle(Geom_BSplineSurface) m_surface;

    bool m_verbose {false};
    bool m_inputSurfaceWasPeriodic_u {false};
    bool m_inputSurfaceWasPeriodic_v {false};
  };

  /// Helper class in support of extracting necessary information from trimming curves
  /// The constructor converts the curve to a clamped (non-periodic) representation, if necessary
  class CurveProcessor
  {
  public:
    CurveProcessor() = delete;

    CurveProcessor(const Handle(Geom2d_BSplineCurve) & curve, bool verbose = false)
      : m_curve(curve)
      , m_verbose(verbose)
    {
      ensureClamped();
    }

    const Handle(Geom2d_BSplineCurve) & getCurve() const { return m_curve; }

    NCurve nurbsCurve() const
    {
      const bool isPeriodic = m_curve->IsPeriodic();
      SLIC_ERROR_IF(isPeriodic,
                    "Axom's NURBSCurve only supports non-periodic curves");

      return m_curve->IsRational()
        ? NCurve(extractControlPoints(), extractWeights(), extractCombinedKnots())
        : NCurve(extractControlPoints(), extractCombinedKnots());
    }

    /// Function to check that the current curve is reasonably close to another curve
    bool compareToCurve(Handle(Geom2d_BSplineCurve) & otherCurve,
                        int numSamples,
                        double sq_tol = 1e-8) const
    {
      if(numSamples <= 1)
      {
        return true;
      }

      auto knot_vals = extractKnotValues();

      axom::Array<double> params(numSamples);
      axom::numerics::linspace(knot_vals.front(),
                               knot_vals.back(),
                               params.data(),
                               numSamples);

      auto evaluateCurve = [](auto curve, double t) {
        const gp_Pnt2d knot_point = curve->Value(t);
        return PointType {knot_point.X(), knot_point.Y()};
      };

      double squared_sum = 0.;
      for(auto val : params)
      {
        auto my_val = evaluateCurve(m_curve, val);
        auto other_val = evaluateCurve(otherCurve, val);

        const double sq_dist = axom::primal::squared_distance(my_val, other_val);
        squared_sum += sq_dist;
        SLIC_WARNING_IF(
          m_verbose && sq_dist > sq_tol,
          axom::fmt::format("Distance between curves at evaluated param {} "
                            "exceeded tolerance {}.\n"
                            "Point on my curve {}; Point on other curve {}; "
                            "Squared distance {}  (running sum {})",
                            val,
                            sq_tol,
                            my_val,
                            other_val,
                            sq_dist,
                            squared_sum));
      }

      return squared_sum <= sq_tol;
    }

  private:
    /// converts the curve from periodic knots to clamped knots, when necessary
    void ensureClamped()
    {
      if(!m_curve->IsPeriodic())
      {
        return;
      }

      const int degree = m_curve->Degree();

      const bool isRational = m_curve->IsRational();

      // Set the curve to not periodic; this can change the poles, knots and multiplicities
      m_curve->SetNotPeriodic();

      // Modify knots and mults to ensure proper clamping
      // We remove the first and last knot and increase the multiplicity of the second and second-to-last knots
      auto mod_knots = extractKnotValues();
      auto mod_mults = extractKnotMultiplicities();
      SLIC_ASSERT(mod_knots.size() >= 4);
      SLIC_ASSERT(mod_mults.size() == mod_knots.size());

      const int last_idx = mod_knots.size() - 1;
      mod_mults[1] += mod_mults[0];
      mod_mults[last_idx - 1] += mod_mults[last_idx];
      mod_knots.erase(mod_knots.end() - 1);
      mod_knots.erase(mod_knots.begin());
      mod_mults.erase(mod_mults.end() - 1);
      mod_mults.erase(mod_mults.begin());

      try
      {
        // Create a new BSpline surface from the modified knots and control points
        m_curve = createBSplineCurveFromAxomArrays(
          extractControlPoints(),
          isRational ? extractWeights() : axom::Array<double>(),
          mod_knots,
          mod_mults,
          degree);
      }
      catch(const Standard_Failure& e)
      {
        std::cerr << "Error: " << e.GetMessageString() << std::endl;
        throw;
      }
    }

    /// generates a new BSpline surface from data in axom::Arrays
    Handle(Geom2d_BSplineCurve) createBSplineCurveFromAxomArrays(
      const axom::Array<PointType>& controlPoints,
      const axom::Array<double>& weights,  // will be empty if curve is non-rational
      const axom::Array<double>& knots,
      const axom::Array<int>& mults,
      int degree)
    {
      // Convert control points to OpenCascade array
      TColgp_Array1OfPnt2d occPoles(1, controlPoints.size());
      for(int i = 0; i < controlPoints.size(); ++i)
      {
        occPoles.SetValue(i + 1,
                          gp_Pnt2d(controlPoints[i][0], controlPoints[i][1]));
      }

      // Convert knots and mults to OpenCascade arrays
      TColStd_Array1OfReal occKnots(1, knots.size());
      TColStd_Array1OfInteger occMults(1, mults.size());
      for(int i = 0; i < knots.size(); ++i)
      {
        occKnots.SetValue(i + 1, knots[i]);
        occMults.SetValue(i + 1, mults[i]);
      }

      Handle(Geom2d_BSplineCurve) clamped_curve;

      // Copy updated weights into OpenCascade array
      if(!weights.empty())
      {
        TColStd_Array1OfReal occWeights(1, weights.size());
        for(int i = 0; i < weights.size(); ++i)
        {
          occWeights.SetValue(i + 1, weights[i]);
        }

        clamped_curve =
          new Geom2d_BSplineCurve(occPoles, occWeights, occKnots, occMults, degree);
      }
      else
      {
        clamped_curve =
          new Geom2d_BSplineCurve(occPoles, occKnots, occMults, degree);
      }

      return clamped_curve;
    }

    // note: bounding boxes are used as debug checks that control points
    // lie within the parameter space of the parent patch
    axom::Array<PointType> extractControlPoints() const
    {
      axom::Array<PointType> controlPoints;

      TColgp_Array1OfPnt2d paraPoints(1, m_curve->NbPoles());
      m_curve->Poles(paraPoints);

      for(Standard_Integer i = paraPoints.Lower(); i <= paraPoints.Upper(); ++i)
      {
        gp_Pnt2d paraPt = paraPoints(i);
        controlPoints.emplace_back(PointType2D {paraPt.X(), paraPt.Y()});
      }

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

    axom::Array<double> extractCombinedKnots() const
    {
      const auto vals = extractKnotValues();
      const auto mults = extractKnotMultiplicities();
      const int total_knots = std::accumulate(mults.begin(), mults.end(), 0);

      axom::Array<double> knots(0, total_knots);
      for(int i = 0; i < vals.size(); ++i)
      {
        knots.insert(knots.end(), mults[i], vals[i]);
      }
      SLIC_ASSERT(knots.size() == total_knots);

      return knots;
    }

    /// converts knot values to axom Array w/o accounting for multiplicity
    /// /sa extractKnotMultiplicities
    axom::Array<double> extractKnotValues() const
    {
      const int num_knots = m_curve->NbKnots();
      axom::Array<double> knots(0, num_knots);

      TColStd_Array1OfReal occ_knots(1, num_knots);
      m_curve->Knots(occ_knots);
      for(int i = 1; i <= occ_knots.Length(); ++i)
      {
        knots.push_back(occ_knots(i));
      }

      return knots;
    }

    /// converts knot multiplicities to axom Array
    /// /sa extractKnotValues
    axom::Array<int> extractKnotMultiplicities() const
    {
      const int num_knots = m_curve->NbKnots();
      axom::Array<int> mults(0, num_knots);

      TColStd_Array1OfInteger occ_mults(1, num_knots);
      m_curve->Multiplicities(occ_mults);
      for(int i = 1; i <= occ_mults.Length(); ++i)
      {
        mults.push_back(occ_mults(i));
      }

      return mults;
    }

  private:
    Handle(Geom2d_BSplineCurve) m_curve;
    bool m_verbose {false};
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

  /// Extracts data from the faces of the mesh and converts each patch to a NURBSPatch
  void extractPatches()
  {
    int patchIndex = 0;

    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      const TopoDS_Face& face = TopoDS::Face(faceExp.Current());

      Handle(Geom_Surface) surface = BRep_Tool::Surface(face);
      if(surface->IsKind(STANDARD_TYPE(Geom_BSplineSurface)))
      {
        Handle(Geom_BSplineSurface) bsplineSurface =
          Handle(Geom_BSplineSurface)::DownCast(surface);

        PatchProcessor patchProcessor(bsplineSurface);

        SLIC_INFO(axom::fmt::format("---"));
        SLIC_INFO("*** Processing patch " << patchIndex);
        SLIC_INFO(axom::fmt::format("---"));

        PatchData& patchData = m_patchData[patchIndex];
        patchData.patchIndex = patchIndex;
        patchData.nurbsPatch = patchProcessor.nurbsPatch();
        patchData.wasOriginallyPeriodic_u =
          patchProcessor.patchWasOriginallyPeriodic_u();
        patchData.wasOriginallyPeriodic_v =
          patchProcessor.patchWasOriginallyPeriodic_v();
        patchData.parametricBBox = faceBoundingBox(face);
        patchData.physicalBBox = patchData.nurbsPatch.boundingBox();

        if(patchData.wasOriginallyPeriodic_u || patchData.wasOriginallyPeriodic_v)
        {
          Handle(Geom_BSplineSurface) origSurface =
            Handle(Geom_BSplineSurface)::DownCast(surface);
          const bool withinThreshold =
            patchProcessor.compareToSurface(origSurface, 25);

          SLIC_WARNING_IF(!withinThreshold,
                          axom::fmt::format("[Patch {}] Patch geometry was not "
                                            "within threshold after clamping.",
                                            patchIndex,
                                            patchData.nurbsPatch));
        }
      }
    }
  }

  /// Extracts data from the trimming curves of each patch and converts the curves to a NURBSCurve representation
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
      PatchData& patchData = m_patchData[patchIndex];
      NCurveArray& curves = patchData.trimmingCurves;

      // Get span of this patch in u and v directions
      BBox2D patchBbox = patchData.parametricBBox;
      auto expandedPatchBbox = patchBbox;
      expandedPatchBbox.scale(1. + 1e-3);

      SLIC_INFO_IF(
        m_verbose,
        axom::fmt::format(
          "[Patch {}]: BBox in parametric space: {}; expanded BBox {}",
          patchIndex,
          patchBbox,
          expandedPatchBbox));

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
            const bool originalCurvePeriodic = bsplineCurve->IsPeriodic();

            CurveProcessor curveProcessor(bsplineCurve, m_verbose);
            curves.emplace_back(curveProcessor.nurbsCurve());

            if(isReversed)  // Ensure consistency of curve w.r.t. patch
            {
              curves.back().reverseOrientation();
            }
            SLIC_ASSERT(curves.back().isValidNURBS());
            SLIC_ASSERT(curves.back().getDegree() == bsplineCurve->Degree());

            SLIC_INFO_IF(
              m_verbose,
              axom::fmt::format(
                "[Patch {} Wire {} Edge {} Curve {}] Added curve: {}",
                patchIndex,
                wireIndex,
                edgeIndex,
                curves.size(),
                curves.back()));

            // Check to ensure that curve did not change geometrically after making non-periodic
            if(originalCurvePeriodic)
            {
              Handle(Geom2d_BSplineCurve) origCurve =
                Geom2dConvert::CurveToBSplineCurve(parametricCurve);
              const bool withinThreshold =
                curveProcessor.compareToCurve(origCurve, 25);
              SLIC_WARNING_IF(
                !withinThreshold,
                axom::fmt::format(
                  "[Patch {} Wire {} Edge {} Curve {}] Trimming curve was not "
                  "within threshold after clamping.",
                  patchIndex,
                  wireIndex,
                  edgeIndex,
                  curves.size(),
                  curves.back()));
            }
          }
        }
      }
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

/**
 * Utility function to represent a NURBSCurve as an SVG path
 *
 * Since SVG only represents polynomial Bezier splines up to order 3,
 * this function distretizes rational curves and linear curves with order above three
 * to a polyline representation
 */
std::string nurbsCurveToSVGPath(const axom::primal::NURBSCurve<double, 2>& curve)
{
  using PointType = axom::primal::Point<double, 2>;

  const int degree = curve.getDegree();
  const auto& knotVector = curve.getKnots();
  const bool isRational = curve.isRational();

  axom::fmt::memory_buffer svgPath;
  axom::fmt::format_to(std::back_inserter(svgPath),
                       "  <path class='{} degree-{}' d='",
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

  // add the closing tags for the path
  axom::fmt::format_to(std::back_inserter(svgPath), "' />");

  return axom::fmt::to_string(svgPath);
}

/**
 * Utility function to generate an SVG file over the parametric space of a trimmed NURBS patch
 *
 * Uses a <rect> for the bounding box in parameter space; 
 * adds a <line> for each knot vector in u- and v-
 * and a <path> for each oriented trimming curve
 */
void generateSVGForPatch(int patchIndex, const PatchData& patchData)
{
  const auto& parametricBBox = patchData.parametricBBox;

  SLIC_INFO(axom::fmt::format("Parametric BBox for patch {}: {}",
                              patchIndex,
                              parametricBBox));

  const auto& curves = patchData.trimmingCurves;
  axom::fmt::memory_buffer svgContent;

  // Create a new bounding box by scaling and translating the parametricBBox
  auto scaledParametricBBox = parametricBBox;
  scaledParametricBBox.scale(1.25);

  SLIC_INFO(
    axom::fmt::format("Scaled and translated parametric BBox for patch {}: {}",
                      patchIndex,
                      scaledParametricBBox));

  // add the SVG header
  axom::fmt::format_to(
    std::back_inserter(svgContent),
    "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' \n"
    "     width='{}in' height='{}in' \n"
    "     viewBox='{} {} {} {}' >\n",
    scaledParametricBBox.range()[0],
    scaledParametricBBox.range()[1],
    scaledParametricBBox.getMin()[0],
    scaledParametricBBox.getMin()[1],
    scaledParametricBBox.range()[0],
    scaledParametricBBox.range()[1]);

  // add some CSS styles
  axom::fmt::format_to(std::back_inserter(svgContent), R"raw(
  <style>
    path {{ fill:none; stroke:black; stroke-width:.03; marker-end:url(#arrow); paint-order:fill stroke markers; stroke-linejoin:round; stroke-linecap:round; }}
    rect {{ fill: white; stroke: gray; stroke-width: 0.05; }}
    .u-line {{ fill: none; stroke: gray; stroke-width: 0.01; }}
    .v-line {{ fill: none; stroke: gray; stroke-width: 0.01; }}
  </style>
  )raw");

  // add a marker for the arrow's head to indicate the orientation
  axom::fmt::format_to(std::back_inserter(svgContent), R"raw(
  <defs>
    <marker id='arrow' style='overflow:visible' orient='auto-start-reverse'
        refX='0' refY='0'
        markerWidth='3.3239999' markerHeight='3.8427744'
        viewBox='0 0 5.3244081 6.1553851'>
      <path
          transform='scale(0.8)'
          style='fill:context-stroke;fill-rule:evenodd;stroke:none'
          d='M 5.77,0 L -2.88,4.5 L -1.44,0 L -2.88,-4.5 Z' />
    </marker>
  </defs>
  )raw");

  // add a rectangle for the parametric bounding box and a comment for its bounding boxes
  axom::fmt::format_to(
    std::back_inserter(svgContent),
    "  <!-- Bounding box of ({},{})-degree patch in parametric space: {}; \n"
    "       BBox in physical space: {} -->\n",
    patchData.nurbsPatch.getDegree_u(),
    patchData.nurbsPatch.getDegree_v(),
    patchData.parametricBBox,
    patchData.physicalBBox);

  axom::fmt::format_to(std::back_inserter(svgContent),
                       "  <rect x='{}' y='{}' width='{}' height='{}' />\n",
                       parametricBBox.getMin()[0],
                       parametricBBox.getMin()[1],
                       parametricBBox.range()[0],
                       parametricBBox.range()[1]);

  // add lines for the u- and v- knots
  axom::fmt::format_to(std::back_inserter(svgContent),
                       "  <!-- Lines for u- and v- knots -->\n");

  auto unique_knots_and_multiplicities =
    [](const axom::Array<double>& knots_vector) {
      axom::Array<std::pair<double, int>> uniqueCounts;
      if(knots_vector.size() == 0) return uniqueCounts;

      double currentValue = knots_vector[0];
      int count = 1;

      for(int i = 1; i < knots_vector.size(); ++i)
      {
        if(knots_vector[i] == currentValue)
        {
          ++count;
        }
        else
        {
          uniqueCounts.emplace_back(currentValue, count);
          currentValue = knots_vector[i];
          count = 1;
        }
      }
      uniqueCounts.emplace_back(currentValue, count);

      return uniqueCounts;
    };

  for(const auto& u :
      unique_knots_and_multiplicities(patchData.nurbsPatch.getKnotsArray_u()))
  {
    axom::fmt::format_to(
      std::back_inserter(svgContent),
      "  <line class='u-line mult-{}' x1='{}' y1='{}' x2='{}' y2='{}' />\n",
      u.second,
      u.first,
      parametricBBox.getMin()[1],
      u.first,
      parametricBBox.getMax()[1]);
  }

  for(const auto& v :
      unique_knots_and_multiplicities(patchData.nurbsPatch.getKnotsArray_v()))
  {
    axom::fmt::format_to(
      std::back_inserter(svgContent),
      "  <line class='v-line mult-{}' x1='{}' y1='{}' x2='{}' y2='{}' />\n",
      v.second,
      parametricBBox.getMin()[0],
      v.first,
      parametricBBox.getMax()[0],
      v.first);
  }

  // add a path for each trimming curve
  // add lines for the u- and v- knots
  axom::fmt::format_to(std::back_inserter(svgContent),
                       "  <!-- Paths for patch trimming curves -->\n");
  for(const auto& curve : curves)
  {
    std::string pathData = nurbsCurveToSVGPath(curve);
    axom::fmt::format_to(std::back_inserter(svgContent), "{}\n", pathData);
  }

  // close the image and write to disk
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

class PatchTriangulator
{
public:
  PatchTriangulator() = delete;

  PatchTriangulator(const TopoDS_Shape& shape,
                    double deflection,
                    double angularDeflection)
    : m_shape(shape)
    , m_deflection(deflection)
    , m_angularDeflection(angularDeflection)
  {
    BRepMesh_IncrementalMesh mesh(m_shape,
                                  m_deflection,
                                  Standard_False,
                                  m_angularDeflection);
    if(!mesh.IsDone())
    {
      throw std::runtime_error("Mesh generation failed.");
    }
  }

  /// Utility function to triangulate a trimmed patch and write it to disk as as STL mesh
  /// Uses OpenCASCADE functionality to triangulate the patch
  void triangulateTrimmedPatches()
  {
    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());

      // Create a triangulation of this patch and write its triangles in the STL format
      TopLoc_Location loc;
      Handle(Poly_Triangulation) triangulation =
        BRep_Tool::Triangulation(face, loc);

      if(triangulation.IsNull())
      {
        SLIC_WARNING(
          axom::fmt::format("Error: STL could not be generated for patch {}",
                            patchIndex));
        break;
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
        axom::fmt::format_to(std::back_inserter(stlContent),
                             triangleAsSTLString(triangulation->Node(n1),
                                                 triangulation->Node(n2),
                                                 triangulation->Node(n3)));
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
          break;
        }
        stlFile << axom::fmt::to_string(stlContent);
      }
      std::cout << "STL file generated: " << stlFilename << std::endl;
    }
  }

  void triangulateUntrimmedPatches()
  {
    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
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
      BRepMesh_IncrementalMesh mesh(newFace,
                                    m_deflection,
                                    Standard_False,
                                    m_angularDeflection);

      // Now you can access the triangulation of the new face
      TopLoc_Location loc;
      Handle(Poly_Triangulation) triangulation =
        BRep_Tool::Triangulation(newFace, loc);

      if(triangulation.IsNull())
      {
        SLIC_WARNING(
          axom::fmt::format("Error: STL could not be generated for patch {}",
                            patchIndex));
        break;
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

        axom::fmt::format_to(std::back_inserter(stlContent),
                             triangleAsSTLString(triangulation->Node(n1),
                                                 triangulation->Node(n2),
                                                 triangulation->Node(n3)));
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
          break;
        }
        stlFile << axom::fmt::to_string(stlContent);
        stlFile.close();
      }
      std::cout << "STL file generated: " << stlFilename << std::endl;
    }
  }

  void triangulateFullMesh()
  {
    // Create an unstructured mesh with 3D vertices and triangular cells
    axom::mint::UnstructuredMesh<axom::mint::SINGLE_SHAPE> mesh(
      3,
      axom::mint::TRIANGLE);

    std::vector<int> patch_id;

    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      TopoDS_Face face = TopoDS::Face(faceExp.Current());

      // Create a triangulation of this patch
      TopLoc_Location loc;
      Handle(Poly_Triangulation) triangulation =
        BRep_Tool::Triangulation(face, loc);

      if(triangulation.IsNull())
      {
        SLIC_WARNING(axom::fmt::format(
          "Error: Triangulation could not be generated for patch {}",
          patchIndex));
        continue;
      }

      const int numTriangles = triangulation->NbTriangles();
      for(int i = 1; i <= numTriangles; ++i)
      {
        Poly_Triangle triangle = triangulation->Triangle(i);
        int n1, n2, n3;
        triangle.Get(n1, n2, n3);

        gp_Pnt p1 = triangulation->Node(n1);
        gp_Pnt p2 = triangulation->Node(n2);
        gp_Pnt p3 = triangulation->Node(n3);

        axom::IndexType v1 = mesh.appendNode(p1.X(), p1.Y(), p1.Z());
        axom::IndexType v2 = mesh.appendNode(p2.X(), p2.Y(), p2.Z());
        axom::IndexType v3 = mesh.appendNode(p3.X(), p3.Y(), p3.Z());

        axom::IndexType cell[3] = {v1, v2, v3};
        mesh.appendCell(cell);
        patch_id.push_back(patchIndex);
      }
    }

    // Add a field to store the patch index for each cell
    auto* patchIndexField =
      mesh.createField<int>("patch_index", axom::mint::CELL_CENTERED);

    for(axom::IndexType i = 0; i < mesh.getNumberOfCells(); ++i)
    {
      patchIndexField[i] = patch_id[i];
    }

    axom::mint::write_vtk(&mesh, "triangulated_mesh.vtk");
  }

private:
  // Format the triangle (represented as three points) as an STL triangle
  std::string triangleAsSTLString(gp_Pnt p1, gp_Pnt p2, gp_Pnt p3)
  {
    gp_Vec v1(p1, p2);
    gp_Vec v2(p1, p3);
    gp_Vec normal = v1.Crossed(v2);
    //normal.Normalize(); // note: this can throw an error

    axom::fmt::memory_buffer stlContent;
    axom::fmt::format_to(std::back_inserter(stlContent),
                         "facet normal {} {} {}\n",
                         normal.X(),
                         normal.Y(),
                         normal.Z());
    axom::fmt::format_to(std::back_inserter(stlContent),  //
                         "  outer loop\n");
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
    axom::fmt::format_to(std::back_inserter(stlContent),  //
                         "  endloop\n");
    axom::fmt::format_to(std::back_inserter(stlContent),  //
                         "endfacet\n");

    return axom::fmt::to_string(stlContent);
  }

private:
  TopoDS_Shape m_shape;
  double m_deflection;
  double m_angularDeflection;
};

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger(axom::slic::message::Info);

  axom::CLI::App app {"Quest Step File Example"};

  std::string filename;
  app.add_option("-f,--file", filename, "Input file")->required();

  bool verbosity {false};
  app.add_flag("-v,--verbose", verbosity)
    ->description("Enable verbose output")
    ->capture_default_str();

  double deflection {.1};
  app.add_option("--deflection", deflection)
    ->description(
      "Max distance between actual geometry and triangulated geometry")
    ->capture_default_str();

  double angular_deflection {0.5};
  app.add_option("--angular_deflection", angular_deflection)
    ->description(
      "Angular deflection between adjacent normals when triangulating surfaces")
    ->capture_default_str();

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
  PatchTriangulator patchTriangulator(nurbs_shape, deflection, angular_deflection);
  patchTriangulator.triangulateTrimmedPatches();
  patchTriangulator.triangulateFullMesh();
  patchTriangulator.triangulateUntrimmedPatches();

  return 0;
}
