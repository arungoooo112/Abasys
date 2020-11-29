/**
 * Core functionality for evaluating points, derivatives and related
 * quantities on NURBS curves and surfaces.
 *
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef TINYNURBS_EVALUATE_H
#define TINYNURBS_EVALUATE_H

#include <tuple>
#include <vector>
#include "basis.h"
#include "curve.h"
#include "surface.h"
#include "./coord.h"
#include "./array2.h"
#include "util/util.h"
#include "detail/evaluate.h"

namespace tinynurbs
{


/////////////////////////////////////////////////////////////////////

/**
Evaluate point on a nonrational NURBS curve
@param[in] crv Curve object
@param[in] u Parameter to evaluate the curve at.
@return point Resulting point on the curve at parameter u.
*/
template <typename T> vec3<T> curvePoint(const Curve<T> &crv, T u)
{
    return funs::curvePoint(crv.degree, crv.knots, crv.control_points, u);
}

/**
 * Evaluate point on a rational NURBS curve
 * @param[in] crv RationalCurve object
 * @param[in] u Parameter to evaluate the curve at.
 * @return point Resulting point on the curve.
 */
template <typename T> vec3<T> curvePoint(const RationalCurve<T> &crv, T u)
{

    typedef hvec3<T> tvecnp1;

    // Compute homogenous coordinates of control points
    std::vector<tvecnp1> Cw;
    Cw.reserve(crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); i++)
    {
        Cw.push_back(tvecnp1(cartesianToHomogenous(crv.control_points[i], crv.weights[i])));
    }

    // Compute point using homogenous coordinates
    tvecnp1 pointw = funs::curvePoint(crv.degree, crv.knots, Cw, u);

    // Convert back to cartesian coordinates
    return homogenousToCartesian(pointw);
}

/**
 * Evaluate derivatives of a non-rational NURBS curve
 * @param[in] crv Curve object
 * @param[in] num_ders Number of times to derivate.
 * @param[in] u Parameter to evaluate the derivatives at.
 * @return curve_ders Derivatives of the curve at u.
 * E.g. curve_ders[n] is the nth derivative at u, where 0 <= n <= num_ders.
 */
template <typename T>
std::vector<vec3<T>> curveDerivatives(const Curve<T> &crv, int num_ders, T u)
{
    return funs::curveDerivatives(crv.degree, crv.knots, crv.control_points, num_ders, u);
}

/**
 * Evaluate derivatives of a rational NURBS curve
 * @param[in] u Parameter to evaluate the derivatives at.
 * @param[in] knots Knot vector of the curve.
 * @param[in] control_points Control points of the curve.
 * @param[in] weights Weights corresponding to each control point.
 * @param[in] num_ders Number of times to differentiate.
 * @param[inout] curve_ders Derivatives of the curve at u.
 * E.g. curve_ders[n] is the nth derivative at u, where n is between 0 and
 * num_ders-1.
 */
template <typename T>
std::vector<vec3<T>> curveDerivatives(const RationalCurve<T> &crv, int num_ders, T u)
{

    typedef vec3<T> tvecn;
    typedef hvec3<T> tvecnp1;

    std::vector<tvecn> curve_ders;
    curve_ders.reserve(num_ders + 1);

    // Compute homogenous coordinates of control points
    std::vector<tvecnp1> Cw;
    Cw.reserve(crv.control_points.size());
    for (int i = 0; i < crv.control_points.size(); i++)
    {
        Cw.push_back(cartesianToHomogenous(crv.control_points[i], crv.weights[i]));
    }

    // Derivatives of Cw
    std::vector<tvecnp1> Cwders =
        funs::curveDerivatives(crv.degree, crv.knots, Cw, num_ders, u);

    // Split Cwders into coordinates and weights
    std::vector<tvecn> Aders;
    std::vector<T> wders;
    for (const auto &val : Cwders)
    {
        Aders.push_back(truncateHomogenous(val));
        wders.push_back(val.w);
    }

    // Compute rational derivatives
    for (int k = 0; k <= num_ders; k++)
    {
        tvecn v = Aders[k];
        for (int i = 1; i <= k; i++)
        {
            v -= static_cast<T>(binomial(k, i)) * wders[i] * curve_ders[k - i];
        }
        curve_ders.push_back(v / wders[0]);
    }
    return curve_ders;
}

/**
 * Evaluate the tangent of a B-spline curve
 * @param[in] crv Curve object
 * @return Unit tangent of the curve at u.
 */
template <typename T> vec3<T> curveTangent(const Curve<T> &crv, T u)
{
    std::vector<vec3<T>> ders = curveDerivatives(crv, 1, u);
    vec3<T> du = ders[1];
    T du_len = glm::length(du);
    if (!close(du_len, T(0)))
    {
        du /= du_len;
    }
    return du;
}

/**
 * Evaluate the tangent of a rational B-spline curve
 * @param[in] crv RationalCurve object
 * @return Unit tangent of the curve at u.
 */
template <typename T> vec3<T> curveTangent(const RationalCurve<T> &crv, T u)
{
    std::vector<vec3<T>> ders = curveDerivatives(crv, 1, u);
    vec3<T> du = ders[1];
    T du_len = glm::length(du);
    if (!close(du_len, T(0)))
    {
        du /= du_len;
    }
    return du;
}

/**
 * Evaluate point on a nonrational NURBS surface
 * @param[in] srf Surface object
 * @param[in] u Parameter to evaluate the surface at.
 * @param[in] v Parameter to evaluate the surface at.
 * @return Resulting point on the surface at (u, v).
 */
template <typename T> vec3<T> surfacePoint(const Surface<T> &srf, T u, T v)
{
    return funs::surfacePoint(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v,
                                  srf.control_points, u, v);
}

/**
 * Evaluate point on a non-rational NURBS surface
 * @param[in] srf RationalSurface object
 * @param[in] u Parameter to evaluate the surface at.
 * @param[in] v Parameter to evaluate the surface at.
 * @return Resulting point on the surface at (u, v).
 */
template <typename T> vec3<T> surfacePoint(const RationalSurface<T> &srf, T u, T v)
{

    typedef hvec3<T> tvecnp1;

    // Compute homogenous coordinates of control points
    array2<tvecnp1> Cw;
    Cw.resize(srf.control_points.rows(), srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); i++)
    {
        for (int j = 0; j < srf.control_points.cols(); j++)
        {
            Cw(i, j) =
                tvecnp1(cartesianToHomogenous(srf.control_points(i, j), srf.weights(i, j)));
        }
    }

    // Compute point using homogenous coordinates
    tvecnp1 pointw =
        funs::surfacePoint(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, Cw, u, v);

    // Convert back to cartesian coordinates
    return homogenousToCartesian(pointw);
}

/**
 * Evaluate derivatives on a non-rational NURBS surface
 * @param[in] degree_u Degree of the given surface in u-direction.
 * @param[in] degree_v Degree of the given surface in v-direction.
 * @param[in] knots_u Knot vector of the surface in u-direction.
 * @param[in] knots_v Knot vector of the surface in v-direction.
 * @param[in] control_points Control points of the surface in a 2D array.
 * @param[in] num_ders Number of times to differentiate
 * @param[in] u Parameter to evaluate the surface at.
 * @param[in] v Parameter to evaluate the surface at.
 * @return surf_ders Derivatives of the surface at (u, v).
 */
template <typename T>
array2<vec3<T>> surfaceDerivatives(const Surface<T> &srf, int num_ders, T u, T v)
{
    return funs::surfaceDerivatives(srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v,
                                        srf.control_points, num_ders, u, v);
}

/**
 * Evaluate derivatives on a rational NURBS surface
 * @param[in] srf RationalSurface object
 * @param[in] u Parameter to evaluate the surface at.
 * @param[in] v Parameter to evaluate the surface at.
 * @param[in] num_ders Number of times to differentiate
 * @return Derivatives on the surface at parameter (u, v).
 */
template <typename T>
array2<vec3<T>> surfaceDerivatives(const RationalSurface<T> &srf, int num_ders, T u, T v)
{

    using namespace std;
    using namespace glm;

    typedef vec<3, T> tvecn;
    typedef vec<4, T> tvecnp1;

    array2<tvecnp1> homo_cp;
    homo_cp.resize(srf.control_points.rows(), srf.control_points.cols());
    for (int i = 0; i < srf.control_points.rows(); ++i)
    {
        for (int j = 0; j < srf.control_points.cols(); ++j)
        {
            homo_cp(i, j) =
                cartesianToHomogenous(srf.control_points(i, j), srf.weights(i, j));
        }
    }

    array2<tvecnp1> homo_ders = funs::surfaceDerivatives(
        srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, homo_cp, num_ders, u, v);

    array2<tvecn> Aders;
    Aders.resize(num_ders + 1, num_ders + 1);
    for (int i = 0; i < homo_ders.rows(); ++i)
    {
        for (int j = 0; j < homo_ders.cols(); ++j)
        {
            Aders(i, j) = truncateHomogenous(homo_ders(i, j));
        }
    }

    array2<tvecn> surf_ders(num_ders + 1, num_ders + 1);
    for (int k = 0; k < num_ders + 1; ++k)
    {
        for (int l = 0; l < num_ders - k + 1; ++l)
        {
            auto der = Aders(k, l);

            for (int j = 1; j < l + 1; ++j)
            {
                der -= (T)binomial(l, j) * homo_ders(0, j).w * surf_ders(k, l - j);
            }

            for (int i = 1; i < k + 1; ++i)
            {
                der -= (T)binomial(k, i) * homo_ders(i, 0).w * surf_ders(k - i, l);

                tvecn tmp((T)0.0);
                for (int j = 1; j < l + 1; ++j)
                {
                    tmp -= (T)binomial(l, j) * homo_ders(i, j).w * surf_ders(k - 1, l - j);
                }

                der -= (T)binomial(k, i) * tmp;
            }

            der *= 1 / homo_ders(0, 0).w;
            surf_ders(k, l) = der;
        }
    }
    return surf_ders;
}

/**
 * Evaluate the two orthogonal tangents of a non-rational surface at the given
 * parameters
 * @param[in] srf Surface object
 * @param u Parameter in the u-direction
 * @param v Parameter in the v-direction
 * @return Tuple with unit tangents along u- and v-directions
 */
template <typename T>
std::tuple<vec3<T>, vec3<T>> surfaceTangent(const Surface<T> &srf, T u, T v)
{
    array2<vec3<T>> ptder = surfaceDerivatives(srf, 1, u, v);
    vec3<T> du = ptder(1, 0);
    vec3<T> dv = ptder(0, 1);
    T du_len = glm::length(ptder(1, 0));
    T dv_len = glm::length(ptder(0, 1));
    if (!close(du_len, T(0)))
    {
        du /= du_len;
    }
    if (!close(dv_len, T(0)))
    {
        dv /= dv_len;
    }
    return std::make_tuple(std::move(du), std::move(dv));
}

/**
 * Evaluate the two orthogonal tangents of a rational surface at the given
 * parameters
 * @param[in] srf Rational Surface object
 * @param u Parameter in the u-direction
 * @param v Parameter in the v-direction
 * @return Tuple with unit tangents along u- and v-directions
 */
template <typename T>
std::tuple<vec3<T>, vec3<T>> surfaceTangent(const RationalSurface<T> &srf, T u, T v)
{
    array2<vec3<T>> ptder = surfaceDerivatives(srf, 1, u, v);
    vec3<T> du = ptder(1, 0);
    vec3<T> dv = ptder(0, 1);
    T du_len = glm::length(ptder(1, 0));
    T dv_len = glm::length(ptder(0, 1));
    if (!close(du_len, T(0)))
    {
        du /= du_len;
    }
    if (!close(dv_len, T(0)))
    {
        dv /= dv_len;
    }
    return std::make_tuple(std::move(du), std::move(dv));
}

/**
 * Evaluate the normal a non-rational surface at the given parameters
 * @param[in] srf Surface object
 * @param u Parameter in the u-direction
 * @param v Parameter in the v-direction
 * @param[inout] normal Unit normal at of the surface at (u, v)
 */
template <typename T> vec3<T> surfaceNormal(const Surface<T> &srf, T u, T v)
{
    array2<vec3<T>> ptder = surfaceDerivatives(srf, 1, u, v);
    vec3<T> n = glm::cross(ptder(0, 1), ptder(1, 0));
    T n_len = glm::length(n);
    if (!close(n_len, T(0)))
    {
        n /= n_len;
    }
    return n;
}

/**
 * Evaluate the normal of a rational surface at the given parameters
 * @param[in] srf Rational Surface object
 * @param u Parameter in the u-direction
 * @param v Parameter in the v-direction
 * @return Unit normal at of the surface at (u, v)
 */
template <typename T> vec3<T> surfaceNormal(const RationalSurface<T> &srf, T u, T v)
{
    array2<vec3<T>> ptder = surfaceDerivatives(srf, 1, u, v);
    vec3<T> n = glm::cross(ptder(0, 1), ptder(1, 0));
    T n_len = glm::length(n);
    if (!close(n_len, T(0)))
    {
        n /= n_len;
    }
    return n;
}

} // namespace tinynurbs

#endif // TINYNURBS_EVALUATE_H
