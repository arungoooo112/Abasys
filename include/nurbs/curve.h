/**
 * The Curve class represents a non-uniform polynomial B-spline curve, while the
 * RationalCurve class represents a non-uniform rational B-spline (NURBS) curve.
 *
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef TINYNURBS_CURVE_H
#define TINYNURBS_CURVE_H

#include <exception>
#include <stdexcept>
#include <vector>

#include "util/util.h"
#include "util/coord.h"
#include "refine.h"

namespace tinynurbs
{

// Forward declaration
template <typename T> struct RationalCurve;

/**
Struct for holding a polynomial B-spline curve
@tparam T Data type of control points and knots (float or double)
*/
template <typename T> struct Curve
{
    unsigned int degree;
    std::vector<T> knots;
    std::vector<vec3<T>> control_points;

    Curve() = default;
    Curve(const RationalCurve<T> &crv) : Curve(crv.degree, crv.knots, crv.control_points) {}
    Curve(unsigned int degree, const std::vector<T> &knots,
          const std::vector<vec3<T>> &control_points)
        : degree(degree), knots(knots), control_points(control_points)
    {
    }
};

template <typename T> struct HCurve
{
    unsigned int degree;
    std::vector<T> knots;
    std::vector<hvec3<T>> control_points;

    HCurve() = default;
    HCurve(unsigned int degree, const std::vector<T> &knots,
          const std::vector<hvec3<T>> &control_points)
        : degree(degree), knots(knots), control_points(control_points)
    {
    }
    HCurve(const RationalCurve<T> &crv) : degree(crv.degree), knots(crv.knots) {
        for (int i = 0; i < crv.control_points.size(); i++) {
            control_points.emplace_back(crv.control_points[i].x * crv.weights[i],
                                        crv.control_points[i].y * crv.weights[i],
                                        crv.control_points[i].z * crv.weights[i], crv.weights[i]);
        }
    }

    HCurve& operator=(const RationalCurve<T> &crv) {return static_cast<HCurve<T>>(crv);}


};

/**
Struct for holding a rational B-spline curve
@tparam T Data type of control points and knots (float or double)
*/
template <typename T> struct RationalCurve
{
    unsigned int degree;
    std::vector<T> knots;
    std::vector<vec3<T>> control_points;
    std::vector<T> weights;

    RationalCurve() = default;
    RationalCurve(const Curve<T> &crv)
        : RationalCurve(crv, std::vector<T>(crv.control_points.size(), 1.0))
    {
    }
    RationalCurve(const Curve<T> &crv, const std::vector<T> &weights)
        : RationalCurve(crv.degree, crv.knots, crv.control_points, weights)
    {
    }
    RationalCurve(unsigned int degree, const std::vector<T> &knots,
                  const std::vector<vec3<T>> &control_points, const std::vector<T> weights)
        : degree(degree), knots(knots), control_points(control_points), weights(weights)
    {
    }
    void HRefine(const std::vector<T>& X) {
        vector<T> Ubar;
        vector<hvec3<T>> Qw;

        auto hpoints = cartesianToHomogenous(control_points, weights);
        tinynurbs::RefineKnotVectCurve(degree, knots, hpoints, X, Ubar, Qw);
        knots = Ubar;
        tinynurbs::homogenousToCartesian(Qw, control_points, weights);
    }



};

// Typedefs for ease of use
typedef Curve<float> Curve3f;
typedef Curve<double> Curve3d;
typedef HCurve<float> HCurve3f;
typedef HCurve<double> HCurve3d;
typedef RationalCurve<float> RCurve3f;
typedef RationalCurve<double> RCurve3d;

} // namespace tinynurbs

#endif // TINYNURBS_CURVE_H
