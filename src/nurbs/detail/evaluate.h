/**
 * Core functionality for evaluating points, derivatives and related
 * quantities on NURBS curves and surfaces.
 *
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef FUNS_EVALUATE_H
#define FUNS_EVALUATE_H

#include <tuple>
#include <vector>
#include "../array2.h"
#include "../coord.h"
#include "basis.h"
#include "curve.h"
#include "surface.h"
#include "util/util.h"


namespace funs
{

/////////////////////////////////////////////////////////////////////

/**
 * Evaluate point on a nonrational NURBS curve
 * @param[in] degree Degree of the given curve.
 * @param[in] knots Knot vector of the curve.
 * @param[in] control_points Control points of the curve.
 * @param[in] u Parameter to evaluate the curve at.
 * @return point Resulting point on the curve at parameter u.
 */
template <int dim, typename T>
glm::vec<dim, T> curvePoint(unsigned int degree, const std::vector<T> &knots,
                            const std::vector<glm::vec<dim, T>> &control_points, T u)
{
    // Initialize result to 0s
    glm::vec<dim, T> point(T(0));

    // Find span and corresponding non-zero basis functions
    int span = findSpan(degree, knots, u);
    std::vector<T> N = bsplineBasis(degree, span, knots, u);

    // Compute point
    for (unsigned int j = 0; j <= degree; j++)
    {
        point += static_cast<T>(N[j]) * control_points[span - degree + j];
    }
    return point;
}

/**
 * Evaluate derivatives of a non-rational NURBS curve
 * @param[in] degree Degree of the curve
 * @param[in] knots Knot vector of the curve.
 * @param[in] control_points Control points of the curve.
 * @param[in] num_ders Number of times to derivate.
 * @param[in] u Parameter to evaluate the derivatives at.
 * @return curve_ders Derivatives of the curve at u.
 * E.g. curve_ders[n] is the nth derivative at u, where 0 <= n <= num_ders.
 */
template <int dim, typename T>
std::vector<glm::vec<dim, T>> curveDerivatives(unsigned int degree, const std::vector<T> &knots,
                                               const std::vector<glm::vec<dim, T>> &control_points,
                                               int num_ders, T u)
{

    typedef glm::vec<dim, T> tvecn;
    using std::vector;

    std::vector<glm::vec<dim, T>> curve_ders;
    curve_ders.resize(num_ders + 1);

    // Assign higher order derivatives to zero
    for (int k = degree + 1; k <= num_ders; k++)
    {
        curve_ders[k] = tvecn(0.0);
    }

    // Find the span and corresponding non-zero basis functions & derivatives
    int span = findSpan(degree, knots, u);
    array2<T> ders = bsplineDerBasis<T>(degree, span, knots, u, num_ders);

    // Compute first num_ders derivatives
    int du = num_ders < degree ? num_ders : degree;
    for (int k = 0; k <= du; k++)
    {
        curve_ders[k] = tvecn(0.0);
        for (int j = 0; j <= degree; j++)
        {
            curve_ders[k] += static_cast<T>(ders(k, j)) * control_points[span - degree + j];
        }
    }
    return curve_ders;
}

/**
 * Evaluate point on a nonrational NURBS surface
 * @param[in] degree_u Degree of the given surface in u-direction.
 * @param[in] degree_v Degree of the given surface in v-direction.
 * @param[in] knots_u Knot vector of the surface in u-direction.
 * @param[in] knots_v Knot vector of the surface in v-direction.
 * @param[in] control_points Control points of the surface in a 2d array.
 * @param[in] u Parameter to evaluate the surface at.
 * @param[in] v Parameter to evaluate the surface at.
 * @return point Resulting point on the surface at (u, v).
 */
template <int dim, typename T>
glm::vec<dim, T> surfacePoint(unsigned int degree_u, unsigned int degree_v,
                              const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                              const array2<glm::vec<dim, T>> &control_points, T u, T v)
{

    // Initialize result to 0s
    glm::vec<dim, T> point(T(0.0));

    // Find span and non-zero basis functions
    int span_u = findSpan(degree_u, knots_u, u);
    int span_v = findSpan(degree_v, knots_v, v);
    std::vector<T> Nu = bsplineBasis(degree_u, span_u, knots_u, u);
    std::vector<T> Nv = bsplineBasis(degree_v, span_v, knots_v, v);

    for (int l = 0; l <= degree_v; l++)
    {
        glm::vec<dim, T> temp(0.0);
        for (int k = 0; k <= degree_u; k++)
        {
            temp += static_cast<T>(Nu[k]) *
                    control_points(span_u - degree_u + k, span_v - degree_v + l);
        }

        point += static_cast<T>(Nv[l]) * temp;
    }
    return point;
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
 * @param[out] surf_ders Derivatives of the surface at (u, v).
 */
template <int dim, typename T>
array2<glm::vec<dim, T>>
surfaceDerivatives(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
                   const std::vector<T> &knots_v, const array2<glm::vec<dim, T>> &control_points,
                   unsigned int num_ders, T u, T v)
{

    array2<glm::vec<dim, T>> surf_ders(num_ders + 1, num_ders + 1, glm::vec<dim, T>(0.0));

    // Set higher order derivatives to 0
    for (int k = degree_u + 1; k <= num_ders; k++)
    {
        for (int l = degree_v + 1; l <= num_ders; l++)
        {
            surf_ders(k, l) = glm::vec<dim, T>(0.0);
        }
    }

    // Find span and basis function derivatives
    int span_u = findSpan(degree_u, knots_u, u);
    int span_v = findSpan(degree_v, knots_v, v);
    array2<T> ders_u = bsplineDerBasis(degree_u, span_u, knots_u, u, num_ders);
    array2<T> ders_v = bsplineDerBasis(degree_v, span_v, knots_v, v, num_ders);

    // Number of non-zero derivatives is <= degree
    unsigned int du = std::min(num_ders, degree_u);
    unsigned int dv = std::min(num_ders, degree_v);

    std::vector<glm::vec<dim, T>> temp;
    temp.resize(degree_v + 1);
    // Compute derivatives
    for (int k = 0; k <= du; k++)
    {
        for (int s = 0; s <= degree_v; s++)
        {
            temp[s] = glm::vec<dim, T>(0.0);
            for (int r = 0; r <= degree_u; r++)
            {
                temp[s] += static_cast<T>(ders_u(k, r)) *
                           control_points(span_u - degree_u + r, span_v - degree_v + s);
            }
        }

        int dd = std::min(num_ders - k, dv);

        for (int l = 0; l <= dd; l++)
        {
            for (int s = 0; s <= degree_v; s++)
            {
                surf_ders(k, l) += ders_v(l, s) * temp[s];
            }
        }
    }
    return surf_ders;
}


} // namespace funs

#endif // 
