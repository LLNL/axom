#ifndef SINA_CURVESET_HPP
#define SINA_CURVESET_HPP

/**
 * @file
 *
 * Contains the definition of the CurveSet class. 
 */

#include <string>
#include <unordered_map>

#include "conduit.hpp"

#include "sina/Curve.hpp"

namespace sina {

/**
 * A CurveSet represents an entry in a record's "curve_set".
 *
 * A CurveSet consist of a set of independent and dependent curves. Each curve
 * is a list of numbers along with optional units and tags.
 *
 * \see Record
 * \see Curve
 */
class CurveSet {
public:
    using CurveMap = std::unordered_map<std::string, Curve>;

    /**
     * Create a CurveSet with the given name
     *
     * @param name the name of the CurveSet
     */
    explicit CurveSet(std::string name);

    /**
     * Create a CurveSet from the given Conduit node.
     *
     * @param name the name of the CurveSet
     * @param node the Conduit node representing the CurveSet
     */
    CurveSet(std::string name, conduit::Node const &node);

    /**
     * Get the name of the this CurveSet.
     *
     * @return the curve set's name
     */
    std::string const & getName() const {
        return name;
    }

    /**
     * Add an independent curve.
     *
     * @param curve the curve to add
     */
    void addIndependentCurve(Curve curve);

    /**
     * Add a dependent curve.
     *
     * @param curve the curve to add
     */
    void addDependentCurve(Curve curve);

    /**
     * Get a map of all the independent curves.
     *
     * @return a map of all the independent curves
     */
    CurveMap const &getIndependentCurves() const {
        return independentCurves;
    }

    /**
     * Get a map of all the dependent curves.
     *
     * @return a map of all the dependent curves
     */
    CurveMap const &getDependentCurves() const {
        return dependentCurves;
    }

    /**
     * Convert his CurveSet to a Conduit node.
     *
     * @return the Node representation of this CurveSet
     */
    conduit::Node toNode() const;

private:
    std::string name;
    CurveMap independentCurves;
    CurveMap dependentCurves;
};

}

#endif //SINA_CURVESET_HPP
