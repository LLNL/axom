// Copyright (c) 2017-2020, Lawrence Livermore National Security, LLC and
// other Axom Project Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef AXOM_KLEE_GEOMETRYOPERATOR_HPP
#define AXOM_KLEE_GEOMETRYOPERATOR_HPP

#include <memory>
#include <vector>

#include "axom/core/numerics/Matrix.hpp"
#include "axom/klee/Dimensions.hpp"
#include "axom/primal/geometry/Point.hpp"
#include "axom/primal/geometry/Vector.hpp"

namespace axom { namespace klee {

class GeometryOperatorVisitor;

/**
 * A GeometryOperator describes an operation to perform on the Geometry
 * of a Shape.
 *
 * There is a subclass of GeometryOperator for each operator defined in
 * the format specification. You can figure out which one you're working
 * with via the accept() method. This follows the standard visitor pattern.
 */
class GeometryOperator {
public:
    virtual ~GeometryOperator() = default;

    /**
     * Get the number of dimensions that the geometry must be in before
     * this operator is applied.
     * \return always 2 or 3
     */
    virtual Dimensions startDims() const = 0;

    /**
     * Get the number of dimensions that the geometry will be in after
     * this operator is applied.
     * \return always 2 or 3
     */
    virtual Dimensions endDims() const = 0;

    /**
     * Accept the given visitor. The appropriate visit() method will
     * be called on the visitor based on the run-time type of this operator.
     *
     * \param visitor the visitor to accept.
     */
    virtual void accept(GeometryOperatorVisitor &visitor) const = 0;
};

/**
 * A MatrixOperator is a type of GeometryOperator whose operation can be
 * expressed as a 4x4 affine transformation matrix.
 */
class MatrixOperator : public virtual GeometryOperator {
public:
    /**
     * Convert this operator to its matrix representation.
     *
     * \return a 4x4 affine transformation matrix
     */
    virtual numerics::Matrix<double> toMatrix() const = 0;
};

/**
 * A CompositeOperator is a GeometryOperator which consists of a list of
 * other operators.
 */
class CompositeOperator : public GeometryOperator {
public:
    using OpPtr = std::shared_ptr<const GeometryOperator>;

    Dimensions startDims() const override;

    Dimensions endDims() const override;

    void accept(GeometryOperatorVisitor &visitor) const override;

    /**
     * Add the given operator to the end of the list of operators in this
     * composite.
     *
     * \param op the operator to add
     */
    void addOperator(const OpPtr &op);

    /**
     * Get a list of all the operators. They should be applied in order.
     *
     * \return the list of operators
     */
    const std::vector<OpPtr> &getOperators() const {
        return m_operators;
    }

private:
    std::vector<OpPtr> m_operators;
};

/**
 * A ConstantDimensionOperator is an abstract GeometryOperator which can
 * be used by subclasses that operate in a single dimension. This is the
 * case for all transformations (e.g. rotations, translations, etc).
 */
class ConstantDimensionOperator : public virtual GeometryOperator {
public:
    /**
     * Create a ConstantDimensionOperator which works in the specified
     * number of dimensions.
     *
     * \param dims the number of dimensions
     */
    explicit ConstantDimensionOperator(Dimensions dims) : m_dims{dims} {}

    Dimensions startDims() const override {
        return m_dims;
    }

    Dimensions endDims() const override {
        return m_dims;
    }

private:
    Dimensions m_dims;
};

/**
 * A Translation is a GeometryOperator which translates points.
 */
class Translation : public ConstantDimensionOperator, public MatrixOperator {
public:
    /**
     * Create a Translation.
     *
     * \param offset the amount by which to offset points
     * \param dims the number of dimensions the translation was meant to
     * operate in. If 2, the 3rd entry in the offset should be zero, but this
     * is not checked.
     */
    Translation(const primal::Vector3D &offset, Dimensions dims);

    /**
     * Get the amount by which to offset points.
     *
     * \return a vector by which points should be offset
     */
    const primal::Vector3D &getOffset() const {
        return m_offset;
    }

    numerics::Matrix<double> toMatrix() const override;

    void accept(GeometryOperatorVisitor &visitor) const override;

private:
    primal::Vector3D m_offset;
};

/**
 * A Rotation is a GeometryOperator which rotates points about a given
 * axis.
 */
class Rotation : public ConstantDimensionOperator, public MatrixOperator {
public:
    /**
     * Create a Rotation.
     *
     * \param angle the angle, in degrees, by which to rotate. Rotations are
     * counter-clockwise.
     * \param center the center of rotation
     * \param axis the axis about which to rotate points
     * \param dims the number of dimensions the translation was meant to
     * operate in. If 2, the axis should be [0, 0, 1], but this is not checked.
     */
    Rotation(double angle, const primal::Point3D &center,
            const primal::Vector3D &axis, Dimensions dims);

    /**
     * Get the angle of rotation.
     *
     * \return the amount by which to rotate in degrees.
     */
    double getAngle() const {
        return m_angle;
    }

    /**
     * Get the center of rotation.
     *
     * \return the point about which to rotate in 2D, and in 3D, the point
     * which defines the axis of rotation along with getAxis().
     */
    const primal::Point3D &getCenter() const {
        return m_center;
    }

    /**
     * The direction of the axis of rotation.
     *
     * \return the vector, which when combined with the center, defines the
     * axis of rotation.
     */
    const primal::Vector3D &getAxis() const {
        return m_axis;
    }

    numerics::Matrix<double> toMatrix() const override;

    void accept(GeometryOperatorVisitor &visitor) const override;

private:
    double m_angle;
    primal::Point3D m_center;
    primal::Vector3D m_axis;
};

/**
 * A GeometryOperator for scaling shapes.
 */
class Scale : public ConstantDimensionOperator, public MatrixOperator {
public:
    /**
     * Create a new Scale operator.
     *
     * \param xFactor the amount by which to scale in the x direction
     * \param yFactor the amount by which to scale in the y direction
     * \param zFactor the amount by which to scale in the z direction
     * \param dims the number of dimensions. If 3, zFactor should be 1.0,
     * but this is not enforced.
     */
    Scale(double xFactor, double yFactor, double zFactor, Dimensions dims);

    /**
     * Get the scale factor in the x direction.
     *
     * \return the x scale factor
     */
    double getXFactor() const {
        return m_xFactor;
    }

    /**
     * Get the scale factor in the y direction.
     *
     * \return the y scale factor
     */
    double getYFactor() const {
        return m_yFactor;
    }

    /**
     * Get the scale factor in the z direction.
     *
     * \return the z scale factor
     */
    double getZFactor() const {
        return m_zFactor;
    }

    numerics::Matrix<double> toMatrix() const override;

    void accept(GeometryOperatorVisitor &visitor) const override;

private:
    double m_xFactor;
    double m_yFactor;
    double m_zFactor;
};

/**
 * An ArbitraryMatrixOperator is an operator which a user specified via
 * an arbitrary affine matrix.
 */
class ArbitraryMatrixOperator : public ConstantDimensionOperator,
                                public MatrixOperator {
public:
    /**
     * Create an ArbitraryMatrixOperator from the given transformation matrix.
     *
     * \param transformation the matrix which defines the transformation
     * \param dims the number of dimensions the matrix operates on. This is
     * not enforced.
     */
    ArbitraryMatrixOperator(const numerics::Matrix<double> &transformation,
            Dimensions dims);

    numerics::Matrix<double> toMatrix() const override;

    void accept(GeometryOperatorVisitor &visitor) const override;

private:
    numerics::Matrix<double> m_transformation;
};

/**
 * A SliceOperator takes a 3D shape and makes it 2D by defining a cut plane.
 * This plane is augmented with an origin and an "up" direction to orient
 * the slice to the x-y plane.
 */
class SliceOperator : public MatrixOperator {
public:
    /**
     * Create a new Slice.
     *
     * \param origin the origin of the coordinate system
     * \param normal a vector normal to the slice plane. Cannot be a zero
     * vector.
     * \param up the direction of the positive Y axis. Must be normal to
     * the "normal" vector.
     */
    SliceOperator(const primal::Point3D &origin, const primal::Vector3D &normal,
            const primal::Vector3D &up);

    /**
     * Get the origin of the coordinate system.
     *
     * \return the system's origin
     */
    const primal::Point3D &getOrigin() const {
        return m_origin;
    }

    /**
     * Get a vector normal to the slice plane.
     *
     * \return a vector normal to the slice plane
     */
    const primal::Vector3D &getNormal() const {
        return m_normal;
    }

    /**
     * Get a vector in the direction of the positive Y axis.
     *
     * \return the direction of the positive Y axis
     */
    const primal::Vector3D &getUp() const {
        return m_up;
    }

    Dimensions startDims() const override;
    Dimensions endDims() const override;

    numerics::Matrix<double> toMatrix() const override;

    void accept(GeometryOperatorVisitor &visitor) const override;

private:
    numerics::Matrix<double> createRotation() const;
    numerics::Matrix<double> createTranslationToOrigin() const;
    primal::Vector3D calculateRightVector() const;

    primal::Point3D m_origin;
    primal::Vector3D m_normal;
    primal::Vector3D m_up;
};

/**
 * A GeometryOperatorVisitor defines visitor for interacting with
 * instances of GeometryOperator. It defines a "visit()" method for each
 * type of operator that a user can specify in the input file, as well
 * as one for CompositeOperator to work with a whole list of operators.
 *
 */
class GeometryOperatorVisitor {
public:
    virtual ~GeometryOperatorVisitor() = default;

    virtual void visit(const Translation &translation) = 0;

    virtual void visit(const Rotation &rotation) = 0;

    virtual void visit(const Scale &scale) = 0;

    virtual void visit(const ArbitraryMatrixOperator &op) = 0;

    virtual void visit(const CompositeOperator &composite) = 0;

    virtual void visit(const SliceOperator &slice) = 0;
};

}}

#endif //AXOM_KLEE_GEOMETRYOPERATOR_HPP
