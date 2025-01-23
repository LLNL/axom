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
#include "opencascade/BRepBuilderAPI_MakeFace.hxx"
#include "opencascade/BRepBuilderAPI_NurbsConvert.hxx"
#include "opencascade/BRepMesh_IncrementalMesh.hxx"
#include "opencascade/BRepTools.hxx"
#include "opencascade/BRep_Tool.hxx"
#include "opencascade/Geom2d_BSplineCurve.hxx"
#include "opencascade/Geom2d_Curve.hxx"
#include "opencascade/Geom2dConvert.hxx"
#include "opencascade/GeomAbs_CurveType.hxx"
#include "opencascade/Geom_BSplineSurface.hxx"
#include "opencascade/Geom_RectangularTrimmedSurface.hxx"
#include "opencascade/Geom_Surface.hxx"
#include "opencascade/Interface_Static.hxx"
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
#include "opencascade/TopoDS_Face.hxx"
#include "opencascade/TopoDS_Shape.hxx"
#include "opencascade/TopoDS_Wire.hxx"
#include <iostream>

bool USE_SUBSET = false;
const int NUM_SUBSET = 1;
int RELEVANT_INDICES[NUM_SUBSET] = {11};
/**
 * /file quest_jacob_step_file.cpp
 * /brief Example that loads in a STEP file and converts the surface patches and curves to Axom's NURBS representations
 *
 * This specific version contains some of Jacob's visualization and experiment code.
 *  If you find this file in a real PR, please berate Jacob.
 * This example reads in STEP files representing trimmed NURBS meshes using OpenCASCADE, 
 * converts the patches and trimming curves to Axom's NURBSPatch and NURBSCurve primitives, 
 * and generates various outputs including SVG and STL files.
 *
 * /note This example requires Axom to be configured with OpenCASCADE enabled.
 */

/// Struct to hold data associated with each surface patch of the mesh
struct PatchData
{
  int patchIndex {-1};
  bool wasOriginallyPeriodic_u {false};
  bool wasOriginallyPeriodic_v {false};
  axom::primal::NURBSPatch<double, 3> nurbsPatch;
  axom::primal::NURBSPatchData<double> nurbsPatchData;
  axom::primal::BoundingBox<double, 2> parametricBBox;
  axom::primal::BoundingBox<double, 3> physicalBBox;
  // axom::Array<axom::primal::NURBSCurve<double, 2>> trimmingCurves;
  axom::Array<bool> trimmingCurves_originallyPeriodic;
};

/**
 * Class to read in a STEP file representing trimmed NURBS meshes using OpenCASCADE 
 * and convert the patches and trimming curves to axom's NURBSPatch and NURBSCurve primitives.
 * 
 * Implementation note: Since axom's primitives do not support periodic knots, 
 * we must convert the OpenCASCADE analogues to a open/clamped representation, when necessary.
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
  // using PatchToTrimmingCurvesMap = std::map<int, NCurveArray>;

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

  /// Helper class to convert faces of CAD mesh to valid NURBSPatch instances
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
     * Utility function to compare the represented surface to \a otherSurface
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
      axom::Array<double> params_u(numSamples);
      axom::numerics::linspace(knot_vals_u.front(),
                               knot_vals_u.back(),
                               params_u.data(),
                               numSamples);

      auto knot_vals_v = extractKnotValues_v();
      axom::Array<double> params_v(numSamples);
      axom::numerics::linspace(knot_vals_v.front(),
                               knot_vals_v.back(),
                               params_v.data(),
                               numSamples);

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

    /// Logs some information about the patch
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

      if(isRational)
      {
        const auto patch_weights = extractWeights();
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

  /// Helper class to convert trimming curves to valid NURBSCurve instances
  /// The constructor converts the curve to a clamped (non-periodic) representation, if necessary
  class CurveProcessor
  {
  public:
    CurveProcessor() = delete;

    CurveProcessor(const Handle(Geom2d_BSplineCurve) & curve, bool verbose = false)
      : m_curve(curve)
      , m_verbose(verbose)
    {
      m_inputSurfaceWasPeriodic = m_curve->IsPeriodic();

      ensureClamped();
    }

    /// Returns a representation of the trimming curve as an axom::primal::NURBSCurve
    NCurve nurbsCurve() const
    {
      const bool isPeriodic = m_curve->IsPeriodic();
      SLIC_ERROR_IF(isPeriodic,
                    "Axom's NURBSCurve only supports non-periodic curves");

      return m_curve->IsRational()
        ? NCurve(extractControlPoints(), extractWeights(), extractCombinedKnots())
        : NCurve(extractControlPoints(), extractCombinedKnots());
    }

    bool curveWasOriginallyPeriodic() const
    {
      return m_inputSurfaceWasPeriodic;
    }

    /**
     * Utility function to compare the represented curve to \a otherCurve
     * at \a numSamples uniformly sampled points in parameter space
     * Returns true when the sum of distances is less than the provided tolerance
     */
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
    bool m_inputSurfaceWasPeriodic {false};
  };

public:
  StepFileProcessor() = delete;

  StepFileProcessor(const std::string& filename, bool verbose = false)
    : m_verbose(verbose)
  {
    m_shape = loadStepFile(filename);
  }

  void setVerbosity(bool verbose) { m_verbose = verbose; }

  bool isLoaded() const { return m_loadStatus == LoadStatus::SUCEESS; }

  const TopoDS_Shape& getShape() const { return m_shape; }

  int getNumberOfPatches() const { return m_patchData.size(); }

  /// Logs some information about the loaded mesh
  void printMeshInfo() const
  {
    // Helper struct for simple stats over a collection of integers
    struct AccumStatistics
    {
      int min;
      int max;
      double mean;
      double stddev;
    };

    // Lambda to generate AccumStatistics for a list of integers
    auto computeStatistics = [](const std::vector<int>& data) -> AccumStatistics {
      AccumStatistics stats;
      stats.min = *std::min_element(data.begin(), data.end());
      stats.max = *std::max_element(data.begin(), data.end());

      const double sum = std::accumulate(data.begin(), data.end(), 0.0);
      const double sumSquared =
        std::accumulate(data.begin(), data.end(), 0.0, [](double a, double b) {
          return a + b * b;
        });
      stats.mean = sum / data.size();
      stats.stddev =
        std::sqrt(sumSquared / data.size() - stats.mean * stats.mean);
      return stats;
    };

    axom::fmt::memory_buffer out;

    axom::fmt::format_to(std::back_inserter(out),
                         "Details about loaded mesh:\n");

    // summarize the number of patches and trimming curves
    {
      int totalTrimmingCurves = 0;
      for(const auto& kv : m_patchData)
      {
        totalTrimmingCurves += kv.second.nurbsPatch.getNumTrimmingCurves();
      }

      axom::fmt::format_to(
        std::back_inserter(out),
        " - Mesh has {} patches with a total of {} trimming curves\n",
        m_patchData.size(),
        totalTrimmingCurves);
    }

    // compute and print the bounding box of the mesh in physical space
    {
      BBox3D meshBBox;
      for(const auto& kv : m_patchData)
      {
        meshBBox.addBox(kv.second.physicalBBox);
        // axom::fmt::format_to(
        //   std::back_inserter(out),
        //   " - Bounding box of patch {} in physical space: {}\n",
        //   kv.second.patchIndex,
        //   kv.second.physicalBBox);
      }

      axom::fmt::format_to(
        std::back_inserter(out),
        " - Bounding box of the mesh in physical space (in {}): {}\n",
        m_fileUnits,
        meshBBox);
    }

    // compute a histogram of the patch degrees w/ some additional info
    {
      struct counts
      {
        int total {0};
        int rational {0};
        int periodic_u {0};
        int periodic_v {0};
      };

      std::map<std::pair<int, int>, counts> patchDegrees;
      for(const auto& kv : m_patchData)
      {
        const auto& patch = kv.second.nurbsPatch;
        auto& c =
          patchDegrees[std::make_pair(patch.getDegree_u(), patch.getDegree_v())];
        c.total++;
        if(patch.isRational())
        {
          c.rational++;
        }
        if(kv.second.wasOriginallyPeriodic_u)
        {
          c.periodic_u++;
        }
        if(kv.second.wasOriginallyPeriodic_v)
        {
          c.periodic_v++;
        }
      }

      axom::fmt::format_to(std::back_inserter(out),
                           " - Patch degree histogram:\n");
      for(const auto& entry : patchDegrees)
      {
        axom::fmt::format_to(
          std::back_inserter(out),
          "   - Degree (u={}, v={}): {} patches ({} rational{}{})\n",
          entry.first.first,   // degree u
          entry.first.second,  // degree v
          entry.second.total,
          entry.second.rational,
          entry.second.periodic_u > 0
            ? axom::fmt::format("; {} originally periodic in u",
                                entry.second.periodic_u)
            : "",
          entry.second.periodic_v > 0
            ? axom::fmt::format("; {} originally periodic in v",
                                entry.second.periodic_v)
            : "");
      }
    }

    // Compute statistics on the number of spans per patch
    {
      std::vector<int> uSpansPerPatch;
      std::vector<int> vSpansPerPatch;
      for(const auto& kv : m_patchData)
      {
        const auto& patch = kv.second.nurbsPatch;
        uSpansPerPatch.push_back(patch.getKnots_u().getNumKnotSpans());
        vSpansPerPatch.push_back(patch.getKnots_v().getNumKnotSpans());
      }

      AccumStatistics uSpansStats = computeStatistics(uSpansPerPatch);
      AccumStatistics vSpansStats = computeStatistics(vSpansPerPatch);

      axom::fmt::format_to(std::back_inserter(out),
                           " - Number of u-spans per patch:  min: {}, max: {}, "
                           "mean: {:.2f} (stdev: {:.2f})\n",
                           uSpansStats.min,
                           uSpansStats.max,
                           uSpansStats.mean,
                           uSpansStats.stddev);

      axom::fmt::format_to(std::back_inserter(out),
                           " - Number of v-spans per patch:  min: {}, max: {}, "
                           "mean: {:.2f} (stdev: {:.2f})\n",
                           vSpansStats.min,
                           vSpansStats.max,
                           vSpansStats.mean,
                           vSpansStats.stddev);
    }

    // Compute statistics on the number of trimming curves per patch
    {
      std::vector<int> trimmingCurvesPerPatch;
      for(const auto& kv : m_patchData)
      {
        trimmingCurvesPerPatch.push_back(
          kv.second.nurbsPatch.getNumTrimmingCurves());
      }

      AccumStatistics trimmingCurvesStats =
        computeStatistics(trimmingCurvesPerPatch);

      axom::fmt::format_to(std::back_inserter(out),
                           " - Number of trimming curves per patch:  min: {}, "
                           "max: {}, mean: {:.2f} (stdev: {:.2f})\n",
                           trimmingCurvesStats.min,
                           trimmingCurvesStats.max,
                           trimmingCurvesStats.mean,
                           trimmingCurvesStats.stddev);
    }
    // Compute statistics on the degrees of the trimming curves in the mesh
    {
      struct counts
      {
        int total {0};
        int rational {0};
        int periodic {0};
      };

      std::map<int, counts> curveDegreeCounts;
      std::vector<int> curveDegreeList;
      for(const auto& kv : m_patchData)
      {
        const auto& curves = kv.second.nurbsPatch.getTrimmingCurves();
        for(axom::IndexType i = 0; i < curves.size(); i++)
        {
          const int degree = curves[i].getDegree();
          auto& c = curveDegreeCounts[degree];
          c.total++;
          if(curves[i].isRational())
          {
            c.rational++;
          }
          //   if(kv.second.trimmingCurves_originallyPeriodic[i])
          //   {
          //     c.periodic++;
          //   }

          curveDegreeList.push_back(degree);
        }
      }

      AccumStatistics curveDegreeStats = computeStatistics(curveDegreeList);

      // Output the results for curve orders
      axom::fmt::format_to(std::back_inserter(out),
                           " - Mesh trimming curve degree histogram:\n");
      for(const auto& entry : curveDegreeCounts)
      {
        axom::fmt::format_to(std::back_inserter(out),
                             "   - Degree {}: {} curves ({} rational{})\n",
                             entry.first,  // degree
                             entry.second.total,
                             entry.second.rational,
                             entry.second.periodic > 0
                               ? axom::fmt::format("; {} originally periodic",
                                                   entry.second.periodic)
                               : "");
      }

      axom::fmt::format_to(
        std::back_inserter(out),
        "   - Average trimming curve order: {:.2f} (stdev: {:.2f})\n",
        curveDegreeStats.mean,
        curveDegreeStats.stddev);
    }

    // Compute statistics on the number of spans in the trimming curves
    {
      std::vector<int> spansPerCurve;
      for(const auto& kv : m_patchData)
      {
        const auto& curves = kv.second.nurbsPatch.getTrimmingCurves();
        for(const auto& curve : curves)
        {
          spansPerCurve.push_back(curve.getKnots().getNumKnotSpans());
        }
      }

      AccumStatistics spansStats = computeStatistics(spansPerCurve);

      axom::fmt::format_to(std::back_inserter(out),
                           " - Number of spans per trimming curve:  min: {}, "
                           "max: {}, mean: {:.2f} (stdev: {:.2f})\n",
                           spansStats.min,
                           spansStats.max,
                           spansStats.mean,
                           spansStats.stddev);
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

        SLIC_INFO_IF(m_verbose, axom::fmt::format("---"));
        SLIC_INFO_IF(m_verbose, "*** Processing patch " << patchIndex);
        SLIC_INFO_IF(m_verbose, axom::fmt::format("---"));

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
      NCurveArray& curves = patchData.nurbsPatch.getTrimmingCurves();

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
            if(!curveProcessor.nurbsCurve().isValidNURBS())
            {
              std::cout << "OH NO CURVE IS INVALID!!!" << std::endl;
            }

            patchData.trimmingCurves_originallyPeriodic.push_back(
              curveProcessor.curveWasOriginallyPeriodic());

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

      // patchData.nurbsPatch.normalize();
      // if(patchIndex == 40)
      // {
      //   patchData.nurbsPatch.makeUntrimmed();
      //   patchData.nurbsPatch.makeSimpleTrimmed();
      // }

      if(patchData.nurbsPatch.isTrimmed())
      {
        auto min_u = patchData.nurbsPatch.getMinKnot_u();
        auto max_u = patchData.nurbsPatch.getMaxKnot_u();

        auto min_v = patchData.nurbsPatch.getMinKnot_v();
        auto max_v = patchData.nurbsPatch.getMaxKnot_v();

        while(true)
        {
          double test_u = axom::utilities::random_real(min_u, max_u);
          double test_v = axom::utilities::random_real(min_v, max_v);

          axom::primal::Point2D p {test_u, test_v};
          auto gwn = std::lround(patchData.nurbsPatch.parameterGWN(p[0], p[1]));

          if(gwn < 0)
          {
            patchData.nurbsPatch.flipNormals();
            break;
          }
          if(gwn > 0)
          {
            break;
          }
          else if(gwn == 0)
          {
            // patchData.nurbsPatch.printTrimmingCurves(
            // "C:\\Users\\Fireh\\Code\\winding_number_code\\figures_2d\\trimming_"
            // "curve_examples\\unsure_" + std::to_string(patchIndex) + ".txt");
            std::cout << "Unsure if it's the right orientation! Trying again..."
                      << patchIndex << std::endl;
          }
        }
      }
    }
  }

  void precomputePatchData()
  {
    for(auto& kv : m_patchData)
    {
      kv.second.nurbsPatchData =
        axom::primal::NURBSPatchData<double>(kv.first, kv.second.nurbsPatch);
      kv.second.nurbsPatchData.curve_quadrature_maps.resize(
        kv.second.nurbsPatch.getNumTrimmingCurves());
    }
  }

  const PatchDataMap& getPatchDataMap() const { return m_patchData; }
  PatchDataMap& getMutablePatchDataMap() { return m_patchData; }

  std::string getFileUnits() const { return m_fileUnits; }

private:
  /// Returns the canoninical representation of a unit string (e.g. "centimeter" -> "cm")
  std::string getCanonicalUnit(const std::string& unit) const
  {
    // we'll convert all units to lower case
    auto toLower = [](std::string str) {
      std::transform(str.begin(), str.end(), str.begin(), ::tolower);
      return str;
    };

    // start with imperial units
    std::map<std::string, std::string> unitCanonicalMap = {{"inch", "in"},
                                                           {"inches", "in"},
                                                           {"in", "in"},
                                                           {"foot", "ft"},
                                                           {"feet", "ft"},
                                                           {"ft", "ft"},
                                                           {"mile", "mi"},
                                                           {"miles", "mi"},
                                                           {"mi", "mi"}};

    // now add the SI units w/ several suffixes
    // we're going to reverse this for the map to canonical units
    std::map<std::string, std::string> prefixes = {
      {"am", "atto"},
      {"fm", "femto"},
      {"pm", "pico"},
      {"nm", "nano"},
      {"um", "micro"},
      {"mm", "milli"},
      {"cm", "centi"},
      {"dm", "deci"},
      {"m", ""},
      {"dam", "deca"},
      {"hm", "hecto"},
      {"km", "kilo"},
    };

    for(const auto& kv : prefixes)
    {
      const std::string& canonical = kv.first;
      const std::string& prefix = kv.second;
      unitCanonicalMap[canonical] = canonical;
      for(const std::string& suffix : {"meter", "meters", "metre", "metres"})
      {
        unitCanonicalMap[prefix + suffix] = canonical;
      }
    }

    return unitCanonicalMap[toLower(unit)];
  }

  /**
   *  Returns the conversion factor from an input unit to an output unit
   * 
   * \note Converts the units to their canonical form
   * \sa getCanonicalUnit
   */
  double getConversionFactor(const std::string& fileUnits,
                             const std::string& defaultUnits = "mm") const
  {
    std::map<std::string, double> unitConversionMap = {{"am", 1e-15},
                                                       {"fm", 1e-12},
                                                       {"pm", 1e-9},
                                                       {"nm", 1e-6},
                                                       {"um", 1e-3},
                                                       {"mm", 1.0},
                                                       {"cm", 10.0},
                                                       {"dm", 100.0},
                                                       {"m", 1e3},
                                                       {"dam", 1e4},
                                                       {"hm", 1e5},
                                                       {"km", 1e6},
                                                       {"in", 25.4},
                                                       {"ft", 304.8},
                                                       {"mi", 1609344.0}};

    const double fileUnitFactor = unitConversionMap[getCanonicalUnit(fileUnits)];
    const double defaultUnitFactor =
      unitConversionMap[getCanonicalUnit(defaultUnits)];

    return fileUnitFactor / defaultUnitFactor;
  };

  /// Loads the step file \a filename from disk
  /// Uses the units from \a filename
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

    // adjust the units, as needed
    TColStd_SequenceOfAsciiString anUnitLengthNames;
    TColStd_SequenceOfAsciiString anUnitAngleNames;
    TColStd_SequenceOfAsciiString anUnitSolidAngleNames;
    reader.FileUnits(anUnitLengthNames, anUnitAngleNames, anUnitSolidAngleNames);
    if(anUnitLengthNames.Size() > 0)
    {
      m_fileUnits = getCanonicalUnit(anUnitLengthNames(1).ToCString());
      std::string defaultUnit = Interface_Static::CVal("xstep.cascade.unit");
      const double lengthUnit = getConversionFactor(m_fileUnits, defaultUnit);
      reader.SetSystemLengthUnit(lengthUnit);
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

    // Convert to NURBS
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

  std::string m_fileUnits {"mm"};

  PatchDataMap m_patchData;
};

/**
 * Class that processes NURBS patches and generates to  an SVG file over the parametric space of a trimmed NURBS patch
 *
 * Uses a <rect> for the bounding box in parameter space; 
 * adds a <line> for each knot vector in u- and v-
 * and a <path> for each oriented trimming curve
 */
class PatchParametricSpaceProcessor
{
public:
  PatchParametricSpaceProcessor() { }

  void setOutputDirectory(const std::string& dir) { m_outputDirectory = dir; }
  void setUnits(const std::string& units) { m_units = units; }
  void setVerbosity(bool verbosityFlag) { m_verbose = verbosityFlag; }
  void setNumFillZeros(int num)
  {
    if(num >= 0)
    {
      m_numFillZeros = num;
    }
  }

  void generateSVGForPatch(int patchIndex, const PatchData& patchData)
  {
    const auto& parametricBBox = patchData.parametricBBox;

    SLIC_INFO_IF(m_verbose,
                 axom::fmt::format("Parametric BBox for patch {}: {}",
                                   patchIndex,
                                   parametricBBox));

    const auto& curves = patchData.nurbsPatch.getTrimmingCurves();
    axom::fmt::memory_buffer svgContent;

    // Create a new bounding box by scaling and translating the parametricBBox
    auto scaledParametricBBox = parametricBBox;
    scaledParametricBBox.scale(1.25);

    SLIC_INFO_IF(m_verbose,
                 axom::fmt::format(
                   "Scaled and translated parametric BBox for patch {}: {}",
                   patchIndex,
                   scaledParametricBBox));

    // add the SVG header
    axom::fmt::format_to(
      std::back_inserter(svgContent),
      "<svg xmlns='http://www.w3.org/2000/svg' version='1.1' \n"
      "     width='{0}{2}' height='{1}{2}' \n"
      "     viewBox='{3} {4} {0} {1}' >\n",
      scaledParametricBBox.range()[0],
      scaledParametricBBox.range()[1],
      m_units,
      scaledParametricBBox.getMin()[0],
      scaledParametricBBox.getMin()[1]);

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

    std::string svgFilename = axom::utilities::filesystem::joinPath(
      m_outputDirectory,
      axom::fmt::format("patch_{:0{}}.svg", patchIndex, m_numFillZeros));
    std::ofstream svgFile(svgFilename);
    if(svgFile.is_open())
    {
      svgFile << axom::fmt::to_string(svgContent);
      svgFile.close();
      SLIC_INFO_IF(m_verbose, "SVG file generated: " << svgFilename);
    }
    else
    {
      std::cerr << "Error: Unable to open file " << svgFilename
                << " for writing." << std::endl;
    }
  }

private:
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

private:
  std::string m_outputDirectory {"."};
  std::string m_units {"mm"};
  bool m_verbose {false};
  int m_numFillZeros {0};
};

/**
 * Class to assist with triangulating STEP files
 * 
 * This class uses OpenCASCADE's triangulation functionality to generate the triangle meshes
 * Supported triangulations:
 *  - triangulateTrimmedPatch: Separately triangulate each trimmed patch, and output as an STL file ('patch_NNN.stl')
 *  - triangulateUntrimmedPatch: Separately triangulate each patch, ignoring the trimming curves, and output as an STL file ('patch_untrimmed_NNN.stl')
 *  - triangulateFullMesh: Output a triangulation of the entire (trimmed) mesh as a VTK file ('triangulated_mesh.vtk').
 *    The output mesh has a field 'patch_index' that associates each triangle with the index of its original patch in the input mesh
 */
class PatchTriangulator
{
public:
  PatchTriangulator() = delete;

  PatchTriangulator(const TopoDS_Shape& shape,
                    double deflection,
                    double angularDeflection,
                    bool verbose)
    : m_shape(shape)
    , m_deflection(deflection)
    , m_angularDeflection(angularDeflection)
    , m_verbose(verbose)
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

  void setOutputDirectory(const std::string& dir) { m_outputDirectory = dir; }
  void setNumFillZeros(int num)
  {
    if(num >= 0)
    {
      m_numFillZeros = num;
    }
  }

  /// Utility function to triangulate each trimmed patch and write it to disk as as STL mesh
  void triangulateTrimmedPatches()
  {
    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      if(USE_SUBSET)
      {
        bool found = false;
        for(int i = 0; i < NUM_SUBSET; i++)
        {
          if(patchIndex == RELEVANT_INDICES[i])
          {
            found = true;
            break;
          }
        }
        if(!found) continue;
      }

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
      std::string stlFilename = axom::utilities::filesystem::joinPath(
        m_outputDirectory,
        axom::fmt::format("patch_{:0{}}.stl", patchIndex, m_numFillZeros));
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
      SLIC_INFO_IF(m_verbose, "STL file generated: " << stlFilename);
    }
  }

  /// Utility function to triangulate each patch, ignoring the trimming curves, and write it to disk as as STL mesh
  void triangulateUntrimmedPatches()
  {
    int patchIndex = 0;
    for(TopExp_Explorer faceExp(m_shape, TopAbs_FACE); faceExp.More();
        faceExp.Next(), ++patchIndex)
    {
      if(USE_SUBSET)
      {
        bool found = false;
        for(int i = 0; i < NUM_SUBSET; i++)
        {
          if(patchIndex == RELEVANT_INDICES[i])
          {
            found = true;
            break;
          }
        }
        if(!found) continue;
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
      std::string stlFilename = axom::utilities::filesystem::joinPath(
        m_outputDirectory,
        axom::fmt::format("patch_untrimmed_{:0{}}.stl",
                          patchIndex,
                          m_numFillZeros));
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
      SLIC_INFO_IF(m_verbose, "STL file generated: " << stlFilename);
    }
  }

  /// Triangulates the entire mesh as a single VTK file
  /// The association to the original patches is tracked via the patch_index field
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
      if(USE_SUBSET)
      {
        bool found = false;
        for(int i = 0; i < NUM_SUBSET; i++)
        {
          if(patchIndex == RELEVANT_INDICES[i])
          {
            found = true;
            break;
          }
        }
        if(!found) continue;
      }

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

    const std::string filename =
      axom::utilities::filesystem::joinPath(m_outputDirectory,
                                            "triangulated_mesh.vtk");
    axom::mint::write_vtk(&mesh, filename);
    SLIC_INFO_IF(m_verbose,
                 "VTK triangle mesh of entire model generated: " << filename);
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
  bool m_verbose {false};
  std::string m_outputDirectory {"."};
  int m_numFillZeros {0};
};

StepFileProcessor import_step_file(std::string prefix,
                                   std::string filename,
                                   bool getTriangulation = false,
                                   double deflection = 0.1,
                                   double angular_deflection = 0.5,
                                   bool export_vtk = true)
{
  AXOM_ANNOTATE_SCOPE("Load step file");

  const std::string stepFilePath =
    axom::utilities::filesystem::joinPath(prefix, filename + ".step");
  StepFileProcessor stepProcessor(stepFilePath, false);
  if(!stepProcessor.isLoaded())
  {
    std::cerr << "Error: The shape is invalid or empty." << std::endl;
    return stepProcessor;
  }

  stepProcessor.extractPatches();
  stepProcessor.extractTrimmingCurves();
  stepProcessor.precomputePatchData();
  stepProcessor.printMeshInfo();

  const int numPatches = stepProcessor.getPatchDataMap().size();
  const int numFillZeros = static_cast<int>(std::log10(numPatches)) + 1;

  // Generate outputs
  std::string output_dir = prefix + filename + "_output";

  // Ensure output directory exists
  if(!axom::utilities::filesystem::pathExists(output_dir))
  {
    axom::utilities::filesystem::makeDirsForPath(output_dir);
  }

  {
    AXOM_ANNOTATE_SCOPE("Generate SVGs");
    PatchParametricSpaceProcessor patchProcessor;
    patchProcessor.setUnits(stepProcessor.getFileUnits());
    patchProcessor.setVerbosity(false);
    patchProcessor.setOutputDirectory(output_dir);
    patchProcessor.setNumFillZeros(numFillZeros);
    SLIC_INFO(
      axom::fmt::format("Generating SVG meshes for patches and their trimming "
                        "curves in '{}' directory",
                        output_dir));
    for(const auto& entry : stepProcessor.getPatchDataMap())
    {
      patchProcessor.generateSVGForPatch(entry.first, entry.second);
    }
  }
  if(getTriangulation)
  {
    AXOM_ANNOTATE_BEGIN("Triangulate STEP");

    AXOM_UNUSED_VAR(export_vtk);

    auto& nurbs_shape = stepProcessor.getShape();
    PatchTriangulator patchTriangulator(nurbs_shape,
                                        deflection,
                                        angular_deflection,
                                        false);
    patchTriangulator.setOutputDirectory(output_dir);
    patchTriangulator.setNumFillZeros(numFillZeros);
    SLIC_INFO(
      axom::fmt::format("Generating triangles meshes for trimmed and untrimmed "
                        "patches in '{}' directory",
                        output_dir));
    patchTriangulator.triangulateTrimmedPatches();
    patchTriangulator.triangulateFullMesh();
    patchTriangulator.triangulateUntrimmedPatches();
  }

  return stepProcessor;
};

void nut_3d_example()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_gwn_"
    "example\\";

  std::string filename = "nut";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      //auto new_query =
      //  axom::primal::Point<double, 3> {-0.000705868, 0.0389087, -0.00646476};
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  timer.start();
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "nut_slice.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0, 0.0, 0.79375 / 1000},
    axom::primal::Vector<double, 3> {-0.10624472993169809,
                                     0.35459739950848745,
                                     -0.928963261718976},
    1.5 * the_range,
    1.5 * the_range,
    100,
    100);
  timer.stop();

  //auto elapsed_time = timer.elapsedTimeInSec();
}

void nut_2d_example()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_gwn_"
    "example\\";

  std::string filename = "nut";

  axom::Array<axom::primal::BezierCurve<double, 2>> curves;
  convert_from_svg(prefix + filename + ".svg", curves);
  std::ofstream curve_out(prefix + filename + "_curves.txt");
  for(auto& curve : curves)
  {
    curve_out << curve << std::endl;
  }

  std::ofstream wn_out(prefix + filename + "_wn.csv");

  auto bbox = curves_bbox(curves, 1.2, true);

  simple_grid_test(curves, bbox, 300, 300, wn_out);
}

void graphical_abstract_watertight()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
    "abstract\\";

  std::string filename = "machine_part_rhino";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      //auto new_query =
      //  axom::primal::Point<double, 3> {-0.000705868, 0.0389087, -0.00646476};
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  timer.start();
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_watertight_slice1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0021127422476921076, 0.0404, 0.0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    the_range,
    the_range,
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_watertight_slice2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {-0.01048739455123153, 0.0404, 0.0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    the_range,
    the_range,
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_watertight_slice3.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0, 0.0404, 0.0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    the_range,
    the_range,
    100,
    100);
  timer.stop();

  //auto elapsed_time = timer.elapsedTimeInSec();
}

void graphical_abstract_exploded()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
    "abstract\\";

  std::string filename = "machine_part_exploded";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      //auto new_query =
      //  axom::primal::Point<double, 3> {-0.000705868, 0.0389087, -0.00646476};
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  timer.start();
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_explodedt_slice1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0021127422476921076, -0.005, 0.0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    the_range,
    the_range,
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_explodedt_slice2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {-0.01048739455123153, -0.005, 0.0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    the_range,
    the_range,
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_explodedt_slice3.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0, -0.005, 0.0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    the_range,
    the_range,
    100,
    100);
  timer.stop();

  // auto elapsed_time = timer.elapsedTimeInSec();
}

void graphical_abstract_hollow()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\graphical_"
    "abstract\\";

  std::string filename = "machine_part_hollow";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      //auto new_query =
      //  axom::primal::Point<double, 3> {-0.000705868, 0.0389087, -0.00646476};
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  timer.start();
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_hollow_slice1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0021127422476921076, 0.0404, 0.0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    the_range,
    the_range,
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_hollow_slice2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {-0.01048739455123153, 0.0404, 0.0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    the_range,
    the_range,
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_hollow_slice3.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.0, 0.0404, 0.0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    the_range,
    the_range,
    100,
    100);
  timer.stop();

  //auto elapsed_time = timer.elapsedTimeInSec();
}

// void discretized_curve_gwn()
// {
// std::string prefix =
// "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\wn_comparison\\";

// std::string filename = "floppy_curve";

// axom::Array<axom::primal::BezierCurve<double, 2>> curves;
// convert_from_svg(prefix + filename + ".svg", curves);
// std::ofstream curve_out(prefix + filename + "_curves.txt");
// for(auto& curve : curves)
// {
// curve_out << curve << std::endl;
// }

// std::ofstream wn_out(prefix + filename + "_wn.csv");

// axom::prieam curve_out2(prefix + filename + "_curves.txt");
// for(auto& curve : curves2)
// {
// curve_out2 << curve << std::endl;
// }

// std::ofstream wn_out2(prefix + filename + "_wn.csv");

// simple_grid_test(curves2, bbox, 300, 300, wn_out2);mal::BoundingBox<double, 2> bbox;
// bbox.addPoint(axom::primal::Point<double, 2> {78.0, 74.0});
// bbox.addPoint(axom::primal::Point<double, 2> {40.0, 105.0});
// bbox.expand(1.3);

// simple_grid_test(curves, bbox, 300, 300, wn_out);

// filename = "discretized_curve";

// axom::Array<axom::primal::BezierCurve<double, 2>> curves2;
// convert_from_svg(prefix + filename + ".svg", curves2);
// std::ofstr
// }

void discretized_surface_gwn()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\wn_comparison\\";

  std::string filename = "floppy_surface";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-16;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();

  // Use this for u
  // u = Vector<T, 3>( {-0.6536123100031038, -0.10738450531201763, -0.7491725543766935});
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_comparison_slice.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.4226723313331604 / 1000,
                                    -0.08878263831138611 / 1000,
                                    0.0012638476388358486},
    axom::primal::Vector<double, 3> {-0.13465792265548568,
                                     0.9906402898209166,
                                     -0.022339651958791562},
    0.9 * the_range,
    0.9 * the_range,
    100,
    100);

  // ---------------------------------------------------------

  std::string filename2 = "discretized_mesh";

  auto stepProcessor2 = import_step_file(prefix, filename2);

  //auto wn_field2 = [&stepProcessor2, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field2 = [&stepProcessor2,
                    &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor2.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox2;
  for(const auto& kv : stepProcessor2.getPatchDataMap())
  {
    meshBBox2.addBox(kv.second.physicalBBox);
  }

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename2 + "_comparison_slice.vtk",
    wn_field2,
    axom::primal::Point<double, 3> {0.4226723313331604 / 1000,
                                    -0.08878263831138611 / 1000,
                                    0.0012638476388358486},
    axom::primal::Vector<double, 3> {-0.13465792265548568,
                                     0.9906402898209166,
                                     -0.022339651958791562},
    0.9 * the_range,
    0.9 * the_range,
    100,
    100);
}

void quadrature_is_bad()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\quadrature_bad\\";

  std::string filename = "random_surface";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field_quad = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field_quad =
    [&stepProcessor](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      //int integrated_trimming_curves;
      double the_val =
        axom::primal::detail::surface_winding_number(query,
                                                     kv.second.nurbsPatch,
                                                     30);
      wn += the_val;
    }

    return wn;
  };

  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  //auto the_range = 0.5 * meshBBox.range().norm();

  // Use this for u
  // u = Vector<T, 3>( {-0.6536123100031038, -0.10738450531201763, -0.7491725543766935});
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_quadrature_slice_1.vtk",
    wn_field_quad,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    200,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_quadrature_slice_3.vtk",
    wn_field_quad,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[2] - meshBBox.getMin()[2],
    200,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_casting_slice_1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    200,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_casting_slice_3.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[2] - meshBBox.getMin()[2],
    200,
    100);
}

void closure_intuition_3d()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\stokes_"
    "intuition\\";

  std::string filename = "intuition_3d";
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field_boundary = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field_boundary =
    [&stepProcessor](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      //int integrated_trimming_curves;
      double the_val = axom::primal::detail::stokes_winding_number_cached(
        query,
        kv.second.nurbsPatchData,
        axom::primal::detail::DiscontinuityAxis::z,
        15,
        quad_tol);
      wn += the_val;
    }

    return wn;
  };

  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  //auto the_range = 0.5 * meshBBox.range().norm();

  // Use this for u
  // u = Vector<T, 3>( {-0.6536123100031038, -0.10738450531201763, -0.7491725543766935});
  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_gwn_slice_1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {0.18489071073484686, 0.9824752787680089, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    2 * (meshBBox.getMax()[1] - meshBBox.getMin()[1]),
    100,
    200);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_boundary_slice_1.vtk",
    wn_field_boundary,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {0.18489071073484686, 0.9824752787680089, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    2 * (meshBBox.getMax()[1] - meshBBox.getMin()[1]),
    100,
    200);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_gwn_slice_2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, -4.5},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_boundary_slice_2.vtk",
    wn_field_boundary,
    axom::primal::Point<double, 3> {0, 0, -4.5},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    100,
    100);
}

void closure_intuition_2d()
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\stokes_"
    "intuition\\";

  std::string filename = "intuition_2d";

  axom::Array<axom::primal::BezierCurve<double, 2>> curves;
  axom::Array<axom::primal::BezierCurve<double, 2>> curves2;
  convert_from_svg(prefix + filename + ".svg", curves);
  std::ofstream curve_out(prefix + filename + "_curves.txt");
  std::ofstream curve_out2(prefix + filename + "_curves2.txt");
  for(int i = 1; i < curves.size(); i++)
  {
    curve_out << curves[i] << std::endl;
  }
  curves2.push_back(curves[0]);
  curve_out2 << curves[0] << std::endl;

  curves.erase(curves.begin());

  std::ofstream wn_out(prefix + filename + "_wn.csv");
  std::ofstream wn_out2(prefix + filename + "_wn2.csv");

  auto bbox = curves_bbox(curves, 1.2, true);

  simple_grid_test(curves, bbox, 300, 300, wn_out);
  simple_grid_test(curves2, bbox, 300, 300, wn_out2);
}

void plot_trimming_curves(const std::string& prefix)
{
  std::string filename = "trimmed_surface";
  auto stepProcessor = import_step_file(prefix, filename);

  axom::primal::NURBSPatch<double, 3> the_patch;
  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    std::cout << kv.first << std::endl;
    meshBBox.addBox(kv.second.physicalBBox);
    the_patch = kv.second.nurbsPatch;
  }

  // auto the_patch = stepProcessor.getPatchDataMap()[0].nurbsPatch;
  axom::primal::NURBSPatch<double, 3> n1, n2;
  the_patch.diskSplit(0.4, 0.2, 0.08, n1, n2);

  using axom::utilities::filesystem::joinPath;
  the_patch.printTrimmingCurves(joinPath(prefix, "original.txt"));
  n1.printTrimmingCurves(joinPath(prefix, "disk.txt"));
  n2.printTrimmingCurves(joinPath(prefix, "remaining.txt"));
}

// Surface curvature evaluation at parameter (3.86209, 3.73582):
// Surface curvature evaluation at parameter (1.29495, 1.83515):
// Surface curvature evaluation at parameter (2.58151, 0.527985):

void coincident_point()
{
  axom::primal::Point<double, 2> centers[3] = {
    axom::primal::Point<double, 2> {3.86209, 3.73582},
    axom::primal::Point<double, 2> {1.29495, 1.83515},
    axom::primal::Point<double, 2> {2.58151, 0.52798}};

  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_"
    "code\\siggraph25\\intersections\\";
  std::string filename = "squashed_surface";

  auto stepProcessor = import_step_file(prefix, filename);

  axom::primal::NURBSPatch<double, 3> the_patch;
  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    // std::cout << kv.first << std::endl;
    meshBBox.addBox(kv.second.physicalBBox);
    the_patch = kv.second.nurbsPatch;
  }

  const int N = 10;
  double radii[N] =
    {0.125, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001};

  std::ofstream radii_out(prefix + "radii.txt");
  for(int i = 0; i < N; ++i)
  {
    radii_out << radii[i] << std::endl;
  }

  for(int i = 0; i < 3; ++i)
  {
    axom::primal::Vector<double, 3> direction =  //axom::primal::Vector<double, 3> {0.0, 0.0, 1.0};
      -the_patch.normal(centers[i][0], centers[i][1]);

    axom::primal::Point<double, 3> point =
      the_patch.evaluate(centers[i][0], centers[i][1]);

    std::cout << point << std::endl;

    std::ofstream out(prefix + "coincident_point_" + std::to_string(i + 1) +
                      ".txt");

    for(int j = 0; j < N; ++j)
    {
      // Remove a disk from the surface
      axom::primal::NURBSPatch<double, 3> disk, n1;
      the_patch.diskSplit(centers[i][0], centers[i][1], radii[j], n1, disk);

      the_patch.printTrimmingCurves(prefix + "original.txt");
      disk.printTrimmingCurves(prefix + "disk.txt");
      n1.printTrimmingCurves(prefix + "remaining.txt");

      // Requires setting onSurface to true
      out << simple_coincident_wn(point, disk, direction, 0) << std::endl;
    }
  }
}

void nearby_point()
{
  axom::primal::Point<double, 2> centers[3] = {
    axom::primal::Point<double, 2> {3.86209, 3.73582},
    axom::primal::Point<double, 2> {1.29495, 1.83515},
    axom::primal::Point<double, 2> {2.58151, 0.52798}};

  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_"
    "code\\siggraph25\\intersections\\";
  std::string filename = "squashed_surface";

  auto stepProcessor = import_step_file(prefix, filename);

  axom::primal::NURBSPatch<double, 3> the_patch;
  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    // std::cout << kv.first << std::endl;
    meshBBox.addBox(kv.second.physicalBBox);
    the_patch = kv.second.nurbsPatch;
  }

  const int N = 10;
  double radii[N] =
    {0.125, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001};

  std::ofstream radii_out(prefix + "radii.txt");
  for(int i = 0; i < N; ++i)
  {
    radii_out << radii[i] << std::endl;
  }

  for(int i = 0; i < 3; ++i)
  {
    axom::primal::Vector<double, 3> direction =  //axom::primal::Vector<double, 3> {0.0, 0.0, 1.0};
      -the_patch.normal(centers[i][0], centers[i][1]);

    // axom::primal::Vector<double, 3> other_direction =
    //   the_patch.du(centers[i][0], centers[i][1]);

    the_patch.evaluate(centers[i][0], centers[i][1]);

    // axom::primal::Point<double, 3> other_point {point[0] + 0.1 * direction[0],
    // point[1] + 0.1 * direction[1],
    // point[2] + 0.1 * direction[2]};

    std::ofstream out(prefix + "nearby_point_" + std::to_string(i + 1) + ".txt");

    for(int j = 0; j < N; ++j)
    {
      axom::primal::Point<double, 3> random_point {
        axom::utilities::random_real(meshBBox.getMin()[0], meshBBox.getMax()[0]),
        axom::utilities::random_real(meshBBox.getMin()[1], meshBBox.getMax()[1]),
        axom::utilities::random_real(meshBBox.getMin()[2], meshBBox.getMax()[2])};

      // Remove a disk from the surface
      axom::primal::NURBSPatch<double, 3> disk, n1;
      the_patch.diskSplit(centers[i][0], centers[i][1], radii[j], n1, disk);

      // the_patch.printTrimmingCurves(prefix + "original.txt");
      // disk.printTrimmingCurves(prefix + "disk.txt");
      // n1.printTrimmingCurves(prefix + "remaining.txt");

      // Requires setting onSurface to true
      out << simple_coincident_wn(random_point, disk, direction, 0) << std::endl;
    }
  }
}

void spring_example(bool process_only = false)
{
  // ABC example 86
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "spring";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, -0.0019539435993821794, 0},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[2] - meshBBox.getMin()[2],
    50,
    50);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 25.39554751828993 / 1000},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    50,
    50);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_3.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, -25.23807955173339 / 1000},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    50,
    50);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_4.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[2] - meshBBox.getMin()[2],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    50,
    50);
}

void generic_timing_test(std::string prefix,
                         std::string filename,
                         axom::primal::BoundingBox<double, 3> bbox,
                         const axom::primal::Point<int, 3>& resolution)
{
  auto stepProcessor = import_step_file(prefix, filename);

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  int xSteps = resolution[0];
  int ySteps = resolution[1];
  int zSteps = resolution[2];

  double dx =
    (xSteps > 1) ? (bbox.getMax()[0] - bbox.getMin()[0]) / (xSteps - 1) : 0.0;
  double dy =
    (ySteps > 1) ? (bbox.getMax()[1] - bbox.getMin()[1]) / (ySteps - 1) : 0.0;
  double dz =
    (zSteps > 1) ? (bbox.getMax()[2] - bbox.getMin()[2]) / (zSteps - 1) : 0.0;

  int case_patch_totals[4] = {0, 0, 0, 0};
  int case_curves_totals[4] = {0, 0, 0, 0};
  double case_time_totals[4] = {0.0, 0.0, 0.0, 0.0};

  // Write VTK header
  // clang-format off
  axom::fmt::memory_buffer vtk;
  axom::fmt::format_to(std::back_inserter(vtk), "# vtk DataFile Version 3.0\n");
  axom::fmt::format_to(std::back_inserter(vtk), "Scalar field data\n");
  axom::fmt::format_to(std::back_inserter(vtk), "ASCII\n");
  axom::fmt::format_to(std::back_inserter(vtk), "DATASET STRUCTURED_POINTS\n");
  axom::fmt::format_to(std::back_inserter(vtk), "DIMENSIONS {} {} {}\n", xSteps, ySteps, zSteps);
  axom::fmt::format_to(std::back_inserter(vtk), "ORIGIN {} {} {}\n", bbox.getMin()[0], bbox.getMin()[1], bbox.getMin()[2]);
  axom::fmt::format_to(std::back_inserter(vtk), "SPACING {} {} {}\n", dx, dy, dz);
  axom::fmt::format_to(std::back_inserter(vtk), "POINT_DATA {}\n", xSteps * ySteps * zSteps);
  axom::fmt::format_to(std::back_inserter(vtk), "SCALARS scalars float\n");
  axom::fmt::format_to(std::back_inserter(vtk), "LOOKUP_TABLE default\n");
  // clang-format on

  const bool output_query_file = false;
  axom::fmt::memory_buffer res;
  if(output_query_file)
  {
    // write the header for the results table
    axom::fmt::format_to(std::back_inserter(res),
                         "x, y, z, patch_index, integrated_trimming_curves, "
                         "winding_number, case_code, elapsed_time\n");
  }

  AXOM_ANNOTATE_BEGIN("GWN query");
  for(int k = 0; k < zSteps; ++k)
  {
    axom::primal::printLoadingBar(k, zSteps);

    for(int j = 0; j < ySteps; ++j)
    {
      for(int i = 0; i < xSteps; ++i)
      {
        double x = bbox.getMin()[0] + i * dx;
        double y = bbox.getMin()[1] + j * dy;
        double z = bbox.getMin()[2] + k * dz;

        auto query = axom::primal::Point<double, 3> {x, y, z};
        int case_code = -1;
        //int integrated_trimming_curves = 0;
        axom::utilities::Timer timer(false);

        double wn = 0.0;
        for(const auto& kv : stepProcessor.getPatchDataMap())
        {
          int integrated_trimming_curves = 0;
          timer.start();
          double the_val =
            axom::primal::winding_number_casting(query,
                                                 kv.second.nurbsPatchData,
                                                 case_code,
                                                 integrated_trimming_curves,
                                                 edge_tol,
                                                 quad_tol,
                                                 EPS);
          timer.stop();

          case_patch_totals[case_code]++;
          case_curves_totals[case_code] += integrated_trimming_curves;
          case_time_totals[case_code] += timer.elapsedTimeInSec();
          wn += the_val;

          if(output_query_file)
          {
            axom::fmt::format_to(
              std::back_inserter(res),
              "{:.15f}, {:.15f}, {:.15f}, {}, {}, {:.15f}, {}, {:.15f}\n",
              x,
              y,
              z,
              kv.first,
              integrated_trimming_curves,
              the_val,
              case_code,
              timer.elapsedTimeInSec());
          }
        }
        axom::fmt::format_to(std::back_inserter(vtk), "{:15f}\n", wn);
      }
    }
  }
  AXOM_ANNOTATE_END("GWN query");

  using axom::utilities::filesystem::joinPath;
  if(output_query_file)
  {
    std::ofstream results(joinPath(prefix, filename + "_timing_results.csv"));
    results << axom::fmt::to_string(res);
  }
  {
    std::ofstream gwn_field(joinPath(prefix, filename + "_gwn_field.vtk"));
    gwn_field << axom::fmt::to_string(vtk);
  }

  {
    std::ofstream summary(joinPath(prefix, filename + "_timing_summary.csv"));
    summary << "Case, Total patches, Total curves, Time Totals (s)\n";
    summary << axom::fmt::format("Case 0: Outside AABB: {}, {}, {:.15f}\n",
                                 case_patch_totals[0],
                                 case_curves_totals[0],
                                 case_time_totals[0]);

    summary << axom::fmt::format("Case 1: Outside OBB: {}, {}, {:.15f}\n",
                                 case_patch_totals[1],
                                 case_curves_totals[1],
                                 case_time_totals[1]);

    summary << axom::fmt::format("Case 2: Casting Necessary: {}, {}, {:.15f}\n",
                                 case_patch_totals[2],
                                 case_curves_totals[2],
                                 case_time_totals[2]);

    summary << axom::fmt::format(
      "Case 3: Trimming Curve Subdivision: {}, {}, {:.15f}\n",
      case_patch_totals[3],
      case_curves_totals[3],
      case_time_totals[3]);
  }
}

void query_timing_test(const std::string& in_prefix,
                       const std::string& filename,
                       const std::string& out_prefix,
                       const std::vector<axom::primal::Point<double, 3>>& query_points,
                       const std::string& out_suffix = "")
{
  auto stepProcessor = import_step_file(in_prefix, filename);

  constexpr double quad_tol = 1e-12;
  constexpr double EPS = 1e-12;
  constexpr double edge_tol = 1e-12;

  int case_patch_totals[4] = {0, 0, 0, 0};
  int case_curves_totals[4] = {0, 0, 0, 0};
  double case_time_totals[4] = {0.0, 0.0, 0.0, 0.0};

  AXOM_ANNOTATE_BEGIN("GWN query");

  const int num_query_pts = query_points.size();
  const int progress_idx = [num_query_pts](int div) {
    return num_query_pts > div
      ? static_cast<int>(num_query_pts / static_cast<double>(div))
      : num_query_pts;
  }(20);
  axom::Array<double> gwn(num_query_pts);

  for(int idx = 0; idx < num_query_pts; ++idx)
  {
    if(idx % progress_idx == 0)
    {
      axom::primal::printLoadingBar(idx, num_query_pts);
    }

    const auto& query = query_points[idx];
    int case_code = -1;
    axom::utilities::Timer timer(false);

    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves = 0;
      timer.start();
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      timer.stop();

      case_patch_totals[case_code]++;
      case_curves_totals[case_code] += integrated_trimming_curves;
      case_time_totals[case_code] += timer.elapsedTimeInSec();
      wn += the_val;
    }
    gwn[idx] = wn;
  }
  AXOM_ANNOTATE_END("GWN query");

  using axom::utilities::filesystem::joinPath;
  {
    // Write query points and GWN field to CSV
    axom::fmt::memory_buffer csv;
    axom::fmt::format_to(
      std::back_inserter(csv),
      "x, y, z, radius, winding_number, expected_inside, rounded\n");

    int outside_count = 0;
    for(int idx = 0; idx < num_query_pts; ++idx)
    {
      const auto& query = query_points[idx];
      const double radius = std::sqrt(
        query[0] * query[0] + query[1] * query[1] + query[2] * query[2]);
      const bool inside = radius > 1.0 ? false : true;
      const int rounded = std::lround(gwn[idx]);
      axom::fmt::format_to(
        std::back_inserter(csv),
        "{:.15f}, {:.15f}, {:.15f}, {:.15f}, {:.15f}, {}, {}\n",
        query[0],
        query[1],
        query[2],
        radius,
        gwn[idx],
        inside,
        rounded);

      if(rounded == 0)
      {
        ++outside_count;
      }
    }

    std::string csvFilename =
      joinPath(out_prefix,
               axom::fmt::format("{}_gwn_field{}.csv", filename, out_suffix));
    std::ofstream csvFile(csvFilename);
    csvFile << axom::fmt::to_string(csv);
    SLIC_INFO(axom::fmt::format("CSV file generated: '{}'", csvFilename));
    const int inside_count = num_query_pts - outside_count;
    SLIC_INFO(axom::fmt::format(
      axom::utilities::locale(),
      "Out of {:L} query points, {:L} were inside ({:.2f}%) and {:L} were "
      "outside ({:.2f}%)",
      num_query_pts,
      inside_count,
      (static_cast<double>(inside_count) / num_query_pts) * 100.0,
      outside_count,
      (static_cast<double>(outside_count) / num_query_pts) * 100.0));
  }

  {
    std::string out_file = joinPath(
      out_prefix,
      axom::fmt::format("{}_timing_summary{}.csv", filename, out_suffix));

    SLIC_INFO(axom::fmt::format("Writing results to '{}'", out_file));
    std::ofstream summary(out_file);
    summary << "Case, Total patches, Total curves, Time Totals (s)\n";
    summary << axom::fmt::format("Case 0: Outside AABB: {}, {}, {:.15f}\n",
                                 case_patch_totals[0],
                                 case_curves_totals[0],
                                 case_time_totals[0]);

    summary << axom::fmt::format("Case 1: Outside OBB: {}, {}, {:.15f}\n",
                                 case_patch_totals[1],
                                 case_curves_totals[1],
                                 case_time_totals[1]);

    summary << axom::fmt::format("Case 2: Casting Necessary: {}, {}, {:.15f}\n",
                                 case_patch_totals[2],
                                 case_curves_totals[2],
                                 case_time_totals[2]);

    summary << axom::fmt::format(
      "Case 3: Trimming Curve Subdivision: {}, {}, {:.15f}\n",
      case_patch_totals[3],
      case_curves_totals[3],
      case_time_totals[3]);
  }
}

void van_example(bool process_only = false)
{
  // ABC example 4192 i think
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "van";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 1.2728508154400588},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    50,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 1.2728508154400588},
    axom::primal::Vector<double, 3> {0, 0, 1},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[2] - meshBBox.getMin()[2],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    50,
    100);
}

void spring_example_volume(bool process_only = false)
{
  // ABC example 4192 i think
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "spring_subset";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  axom::primal::BoundingBox<double, 3> smaller_box;
  smaller_box.addPoint(
    axom::primal::Point<double, 3> {-0.0385 + 0.005, 0.0 + 0.005, 0.025 + 0.005});
  smaller_box.addPoint(
    axom::primal::Point<double, 3> {-0.0385 - 0.005, 0.0 - 0.005, 0.025 - 0.005});

  axom::primal::exportScalarFieldToVTK<double>(prefix + filename + "_volume.vtk",
                                               wn_field,
                                               smaller_box,
                                               100,
                                               100,
                                               100);
}

void bobbin_example_volume(bool process_only = false)
{
  // ABC example 4192 i think
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "bobbin";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      int integrated_trimming_curves;
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  axom::primal::BoundingBox<double, 3> smaller_box;
  smaller_box.addPoint(
    axom::primal::Point<double, 3> {0.047 + 0.01, -0.005 + 0.01, 0.05 + 0.01});
  smaller_box.addPoint(
    axom::primal::Point<double, 3> {0.047 - 0.01, -0.005 - 0.01, 0.05 - 0.01});

  axom::primal::exportScalarFieldToVTK<double>(prefix + filename + "_volume.vtk",
                                               wn_field,
                                               smaller_box,
                                               100,
                                               100,
                                               100);
}
//0.2741457196484367

void bobbin_example(bool process_only = false)
{
  // ABC example 4192 i think
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "bobbin";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      // 14

      int integrated_trimming_curves;  // 0.062689486818427134, 0, -0.00067335004127314774
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0.02539999969303608},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 0, 1},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[2] - meshBBox.getMin()[2],
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0.04906700113117939},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    100,
    100);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_3.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, -0.0018166745596870217},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    meshBBox.getMax()[0] - meshBBox.getMin()[0],
    meshBBox.getMax()[1] - meshBBox.getMin()[1],
    100,
    100);
}

void complex_gear_example(bool process_only = false)
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "complex_gear";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      // 14

      int integrated_trimming_curves;  // 0.062689486818427134, 0, -0.00067335004127314774
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_slice_1.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0, 0, 0.0005342458753727963},
    axom::primal::Vector<double, 3> {1, 0, 0},
    axom::primal::Vector<double, 3> {0, 1, 0},
    0.2 * (meshBBox.getMax()[0] - meshBBox.getMin()[0]),
    0.2 * (meshBBox.getMax()[1] - meshBBox.getMin()[1]),
    100,
    100);
}

void complex_gear_subset_example(bool process_only = false)
{
  std::string prefix =
    "C:\\Users\\Fireh\\Code\\winding_number_code\\siggraph25\\simple_results\\";

  std::string filename = "complex_gear";
  auto stepProcessor = import_step_file(prefix, filename);

  if(process_only)
  {
    return;
  }

  constexpr double quad_tol = 1e-5;
  constexpr double EPS = 1e-10;
  constexpr double edge_tol = 1e-6;

  // (!bBox, !oBox, casting, noCache)
  int case_code = -1;
  //auto wn_field = [&stepProcessor, &edge_tol, &quad_tol, &EPS, &case_code](
  auto wn_field = [&stepProcessor,
                   &case_code](axom::primal::Point<double, 3> query) -> double {
    double wn = 0.0;
    for(const auto& kv : stepProcessor.getPatchDataMap())
    {
      // 14

      int integrated_trimming_curves;  // 0.062689486818427134, 0, -0.00067335004127314774
      double the_val =
        axom::primal::winding_number_casting(query,
                                             kv.second.nurbsPatchData,
                                             case_code,
                                             integrated_trimming_curves,
                                             edge_tol,
                                             quad_tol,
                                             EPS);
      wn += the_val;
    }

    return wn;
  };

  axom::primal::BoundingBox<double, 3> meshBBox;
  for(const auto& kv : stepProcessor.getPatchDataMap())
  {
    meshBBox.addBox(kv.second.physicalBBox);
  }

  auto the_range = 0.5 * meshBBox.range().norm();
  meshBBox.expand(0.1 * the_range);

  //axom::primal::Point<double, 3> origin = meshBBox.getCentroid();
  axom::primal::Vector<double, 3> normal = {1.0, 1.0, 1.0};

  axom::utilities::Timer timer(false);

  // axom::primal::exportSliceScalarFieldToVTK<double>(
  //   prefix + filename + "_subset_slice_1.vtk",
  //   wn_field,
  //   axom::primal::Point<double, 3> {0.044658458307408416, 0.015908217648718987, 0.0076132148484121575},
  //   axom::primal::Vector<double, 3> {0.9815644542845647, 0.051996087231145134, 0.183922888727031},
  //   0.02,
  //   0.02,
  //   200,
  //   200);

  axom::primal::exportSliceScalarFieldToVTK<double>(
    prefix + filename + "_subset_slice_2.vtk",
    wn_field,
    axom::primal::Point<double, 3> {0.043949259805650946,
                                    0.01569962565063724,
                                    0.0113693041165032},
    axom::primal::Vector<double, 3> {0.05241730314620468,
                                     -0.7681876302012245,
                                     -0.6380753804502299},
    0.02 * 0.075,
    0.02 * 0.075,
    200,
    200);
}

/**
 * Generates a given number of sample points near the surface of a unit sphere
 * 
 * @param numSamples The number of sample points to generate.
 * @param EPS Allowed distance to sphere (inside the sphere)
 * @return A vector of sample points near unit sphere
 */
std::vector<axom::primal::Point<double, 3>> generateSamplePointsOnSphere(
  int numSamples,
  double EPS = 1e-3)
{
  std::vector<axom::primal::Point<double, 3>> samplePoints;
  samplePoints.reserve(numSamples);

  const double eps_lower = 1. - 1.1 * EPS;
  const double eps_upper = 1. - 0.9 * EPS;
  SLIC_INFO(
    axom::fmt::format("Generating {} sample points near unit sphere with EPS: "
                      "{} and range: [{}, {}]",
                      numSamples,
                      EPS,
                      eps_lower,
                      eps_upper));

  for(int i = 0; i < numSamples; ++i)
  {
    const double phi = axom::utilities::random_real(0., 2. * M_PI);
    const double cosTheta = axom::utilities::random_real(-1., 1.);
    const double sinTheta = std::sqrt(1. - cosTheta * cosTheta);

    const double mag = axom::utilities::random_real(eps_lower, eps_upper);

    samplePoints.emplace_back(
      axom::primal::Point<double, 3> {mag * sinTheta * std::cos(phi),
                                      mag * sinTheta * std::sin(phi),
                                      mag * cosTheta});
  }

  return samplePoints;
}

int main(int argc, char** argv)
{
  axom::slic::SimpleLogger logger;

  axom::CLI::App app {"Winding number experiments"};

  app.set_config("--config");

  std::string in_prefix;
  app.add_option("--input-directory", in_prefix, "Prefix path for input file")
    ->required();

  std::string out_prefix {"winding_results"};
  app.add_option("--output-prefix", out_prefix, "Prefix path for output files")
    ->capture_default_str();

  std::string filename;
  app.add_option("--filename", filename, "Filename string")->required();

  std::vector<double> min_bbox(3), max_bbox(3);
  app.add_option("--min-bbox", min_bbox, "BoundingBox min values")
    ->required()
    ->expected(3);
  app.add_option("--max-bbox", max_bbox, "BoundingBox max values")
    ->required()
    ->expected(3);

  std::vector<int> resolution {50, 50, 50};
  app.add_option("--resolution", resolution, "Resolution values")
    ->capture_default_str()
    ->expected(3);

  int num_query_pts {1000};
  app.add_option("-q,--num-query-pts", num_query_pts, "Number of query points")
    ->capture_default_str();

  double distance_to_surface {1e-3};
  app
    .add_option("--distance-to-surface",
                distance_to_surface,
                "Distance to surface")
    ->capture_default_str();

  std::string annotationMode {"none"};
#ifdef AXOM_USE_CALIPER
  app.add_option("--caliper", annotationMode)
    ->description(
      "caliper annotation mode. Valid options include 'none' and 'report'. "
      "Use 'help' to see full list.")
    ->capture_default_str()
    ->check(axom::utilities::ValidCaliperMode);
#endif

  double scale_factor {1.01};
  app.add_option("--scale_factor", scale_factor, "BoundingBox scale factor")
    ->capture_default_str();

  CLI11_PARSE(app, argc, argv);

  axom::primal::BoundingBox<double, 3> unscaled_bbox;
  unscaled_bbox.addPoint(
    axom::primal::Point<double, 3> {min_bbox[0], min_bbox[1], min_bbox[2]});
  unscaled_bbox.addPoint(
    axom::primal::Point<double, 3> {max_bbox[0], max_bbox[1], max_bbox[2]});
  axom::primal::BoundingBox<double, 3> bbox = unscaled_bbox;
  bbox.scale(scale_factor);

  axom::primal::Point<int, 3> res {resolution[0], resolution[1], resolution[2]};

  SLIC_INFO(axom::fmt::format(R"raw(Input parameters:
   Input directory prefix: {}
   Output directory prefix: {}
   Filename: {}
   BBox scale factor: {}
   Unscaled BBox: {}
   Scaled BBox: {}
   Resolution: {}
   Num query pts: {}
   Distance from surface: {}
)raw",
                              in_prefix,
                              out_prefix,
                              filename,
                              scale_factor,
                              unscaled_bbox,
                              bbox,
                              res,
                              num_query_pts,
                              distance_to_surface));

  axom::utilities::raii::AnnotationsWrapper annotations_raii_wrapper(
    annotationMode);

  //generic_timing_test(prefix, filename, bbox, res);
  if(!axom::utilities::filesystem::pathExists(out_prefix))
  {
    axom::utilities::filesystem::makeDirsForPath(out_prefix);
  }

  auto query_pts =
    generateSamplePointsOnSphere(num_query_pts, distance_to_surface);
  query_timing_test(in_prefix,
                    filename,
                    out_prefix,
                    query_pts,
                    axom::fmt::format("_within_{}", distance_to_surface));

  // nut_3d_example();
  // nut_2d_example();
  // discretized_curve_gwn();
  // discretized_surface_gwn();
  // quadrature_is_bad();
  // closure_intuition_2d();
  // coincident_point();
  // nearby_point();
  // spring_timing_example();
  // van_example(true);
  // bobbin_example(true);
  // complex_gear_example(true);
  // complex_gear_subset_example();
  // bobbin_example_volume();
  // van_example_volume();
  // Overlap of 9995 and 4196
  // 9979 at varying levels zooming

  return 0;
}
