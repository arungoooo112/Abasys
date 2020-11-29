/**
 * Functionality for checking validity and properties of NURBS curves and
 * surfaces
 *
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef TINYNURBS_CHECK_H
#define TINYNURBS_CHECK_H

#include <vector>
#include "curve.h"
#include "surface.h"
#include "detail/check.h"

namespace tinynurbs
{

/////////////////////////////////////////////////////////////////////

/**
 * Returns the mulitplicity of the knot at index
 * @tparam Type of knot values
 * @param[in] knots Knot vector
 * @param[in] index Index of knot of interest
 * @return Multiplicity (>= 1)
 */
template <typename T> int knotMultiplicity(const std::vector<T> &knots, unsigned int index)
{
    T u = knots[index];
    int mult = 0;
    int beg = 0;
    for (; beg < knots.size() && abs(knots[beg], u) > 1e-6; beg++);
    for (; beg < knots.size() && abs(knots[beg], u) <= 1e-6; beg++, mult++);
    return mult;
}

/**
 * Returns whether the curve is valid
 * @tparam T Type of control point coordinates, knot values
 * @param[in] crv Curve object
 * @return Whether valid
 */
template <typename T> bool isValidCurve(const Curve<T> &crv)
{
    return funs::isValidCurve<T>(crv.degree, crv.knots, crv.control_points);
}

/**
 * Returns whether the curve is valid
 * @tparam T Type of control point coordinates, knot values
 * @param[in] crv RationalCurve object
 * @return Whether valid
 */
template <typename T> bool isValidCurve(const RationalCurve<T> &crv)
{
    return funs::isValidCurve(crv.degree, crv.knots, crv.control_points, crv.weights);
}

/**
 * Returns whether the surface is valid
 * @tparam T Type of control point coordinates, knot values
 * @param srf Surface object
 * @return Whether valid
 */
template <typename T> bool isValidSurface(const Surface<T> &srf)
{
    return funs::isValidSurface(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v,
                                    srf.control_points);
}

/**
 * Returns whether the rational surface is valid
 * @tparam T Type of control point coordinates, knot values
 * @param[in] srf RationalSurface object
 * @return Whether valid
 */
template <typename T> bool isValidSurface(const RationalSurface<T> &srf)
{
    return funs::isValidSurface(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v,
                                    srf.control_points, srf.weights);
}

// /**
//  * Checks whether the curve is closed
//  * @param[in] crv Curve object
//  * @return  Whether closed
//  */
// template <typename T> bool curveIsClosed(const Curve<T> &crv)
// {
//     return isArray1Closed(crv.degree, crv.control_points) &&
//            isKnotVectorClosed(crv.degree, crv.knots);
// }

// /**
//  * Checks whether the rational curve is closed
//  * @param[in] crv RationalCurve object
//  * @return  Whether closed
//  */
// template <typename T> bool curveIsClosed(const RationalCurve<T> &crv)
// {
//     return isArray1Closed(crv.degree, crv.control_points) &&
//            isArray1Closed(crv.degree, crv.weights) &&
//            isKnotVectorClosed(crv.degree, crv.knots);
// }

// /**
//  * Checks whether the surface is closed along u-direction
//  * @param[in] srf Surface object
//  * @return  Whether closed along u-direction
//  */
// template <typename T> bool surfaceIsClosedU(const Surface<T> &srf)
// {
//     return isArray2ClosedU(srf.degree_u, srf.control_points) &&
//            isKnotVectorClosed(srf.degree_u, srf.knots_u);
// }

// /**
//  * Checks whether the surface is closed along v-direction
//  * @param[in] srf Surface object
//  * @return  Whether closed along v-direction
//  */
// template <typename T> bool surfaceIsClosedV(const Surface<T> &srf)
// {
//     return isArray2ClosedV(srf.degree_v, srf.control_points) &&
//            isKnotVectorClosed(srf.degree_v, srf.knots_v);
// }

// /**
//  * Checks whether the rational surface is closed along u-direction
//  * @param[in] srf RationalSurface object
//  * @return  Whether closed along u-direction
//  */
// template <typename T> bool surfaceIsClosedU(const RationalSurface<T> &srf)
// {
//     return isArray2ClosedU(srf.degree_u, srf.control_points) &&
//            isKnotVectorClosed(srf.degree_u, srf.knots_u) &&
//            isArray2ClosedU(srf.degree_u, srf.weights);
// }

// /**
//  * Checks whether the rational surface is closed along v-direction
//  * @param[in] srf RationalSurface object
//  * @return  Whether closed along v-direction
//  */
// template <typename T> bool surfaceIsClosedV(const RationalSurface<T> &srf)
// {
//     return isArray2ClosedV(srf.degree_v, srf.control_points) &&
//            isKnotVectorClosed(srf.degree_v, srf.knots_v) &&
//            isArray2ClosedV(srf.degree_v, srf.weights);
// }

} // namespace tinynurbs

#endif // TINYNURBS_CHECK_H
