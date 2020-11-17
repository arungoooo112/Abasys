/**
 * Functionality for checking validity and properties of NURBS curves and
 * surfaces
 *
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef BASES_CHECK_H
#define BASES_CHECK_H

#include <algorithm> // issort 
#include <vector>
#include "util/coord.h"
#include "util/array2.h"


namespace tinynurbs
{

/////////////////////////////////////////////////////////////////////

/**
 * Checks if the relation between degree, number of knots, and
 * number of control points is valid
 * @param[in] degree Degree
 * @param[in] num_knots Number of knot values
 * @param[in] num_ctrl_pts Number of control points
 * @return Whether the relationship is valid
 */
inline bool isValidRelation(int degree, int num_knots, int num_ctrl_pts)
{
    return (num_knots - degree - 1) == num_ctrl_pts;
}

/**
 * isValidKnotVector returns whether the knots are in valid with degree
 * @tparam Type of knot values
 * @param[in] knot vector
 * @param[in] degree
 * @return Whether valid
 */
template <typename T> 
bool isValidKnotVector(const std::vector<T> &knots, int degree)
{
    if (knots.size() < (degree + 1) * 2) return false;
    int i = 0, j = knots.size() - 1;
    while (degree--)
    {
        if (abs(knots[i] - knots[i+1]) > 1e-6) return false;
        if (abs(knots[j] - knots[j-1]) > 1e-6) return false;
        i++; j--;
    }
    if (abs(knots[i] - knots[i+1]) <= 1e-6) return false;
    if (abs(knots[j] - knots[j-1]) <= 1e-6) return false;

    return std::is_sorted(knots.begin(), knots.end());
}

/**
 * Returns whether the curve is valid
 * @tparam T Type of control point coordinates, knot values
 * @param[in] degree Degree of curve
 * @param[in] knots Knot vector of curve
 * @param[in] control_points Control points of curve
 * @return Whether valid
 */
template <typename T>
bool isValidCurve(int degree, const std::vector<T> &knots,
                  const std::vector<vec3<T>> &control_points)
{
    if (degree < 1 || degree > 9)
    {
        return false;
    }
    if (!isValidRelation(degree, knots.size(), control_points.size()))
    {
        return false;
    }
    if (!isValidKnotVector(knots, degree))
    {
        return false;
    }
    return true;
}

/**
 * Returns whether the curve is valid
 * @tparam T Type of control point coordinates, knot values and weights
 * @param[in] degree Degree of curve
 * @param[in] knots Knot vector of curve
 * @param[in] control_points Control points of curve
 * @return Whether valid
 */
template <typename T>
bool isValidCurve(unsigned int degree, const std::vector<T> &knots,
                  const std::vector<vec3<T>> &control_points, const std::vector<T> &weights)
{
    if (weights.size() != control_points.size())
    {
        return false;
    }
    if (!isValidCurve(degree, knots, control_points))
    {
        return false;
    }

    return true;
}

/**
 * Returns whether the surface is valid
 * @tparam T Type of control point coordinates, knot values
 * @param[in] degree_u Degree of surface along u-direction
 * @param[in] degree_v Degree of surface along v-direction
 * @param[in] knots_u Knot vector of surface along u-direction
 * @param[in] knots_v Knot vector of surface along v-direction
 * @param[in] control_points Control points grid of surface
 * @return Whether valid
 */
template <typename T>
bool isValidSurface(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
                    const std::vector<T> &knots_v, const array2<vec3<T>> &control_points)
{
    if (degree_u < 1 || degree_u > 9 || degree_v < 1 || degree_v > 9)
    {
        return false;
    }
    if (!isValidRelation(degree_u, knots_u.size(), control_points.rows()) ||
        !isValidRelation(degree_v, knots_v.size(), control_points.cols()))
    {
        return false;
    }
    if (!isValidKnotVector(knots_u, degree_u) || !isValidKnotVector(knots_v, degree_v))
    {
        return false;
    }
    return true;
}

/**
 * Returns whether the rational surface is valid
 * @tparam T Type of control point coordinates, knot values
 * @param[in] degree_u Degree of surface along u-direction
 * @param[in] degree_v Degree of surface along v-direction
 * @param[in] knots_u Knot vector of surface along u-direction
 * @param[in] knots_v Knot vector of surface along v-direction
 * @param[in] control_points Control points grid of surface
 * @param[in] weights Weights corresponding to control point grid of surface
 * @return Whether valid
 */
template <typename T>
bool isValidSurface(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
                    const std::vector<T> &knots_v, const array2<vec3<T>> &control_points,
                    const array2<T> &weights)
{
    if (control_points.rows() != weights.rows() || control_points.cols() != weights.cols())
    {
        return false;
    }
    if (!isValidSurface(degree_u, degree_v, knots_u, knots_v, control_points))
    {
        return false;
    }

    return true;
}

// /**
//  * Returns whether the given knot vector is closed by checking the
//  * periodicity of knot vectors near the start and end
//  * @param[in] degree Degree of curve/surface
//  * @param[in] knots Knot vector of curve/surface
//  * @return Whether knot vector is closed
//  */
// template <typename T> bool isKnotVectorClosed(int degree, const std::vector<T> &knots)
// {
//     for (int i = 0; i < degree - 1; ++i)
//     {
//         int j = knots.size() - degree + i;
//         if (std::abs((knots[i + 1] - knots[i]) - (knots[j + 1] - knots[j])) >
//             std::numeric_limits<T>::epsilon())
//         {
//             return false;
//         }
//     }
//     return true;
// }

// /**
//  * Returns whether the given knot vector is closed by checking the
//  * periodicity of knot vectors near the start and end
//  * @param[in] degree Degree of curve/surface
//  * @param[in] vec Array of any control points/weights
//  * @return Whether knot vector is closed
//  */
// template <typename T> bool isArray1Closed(int degree, const std::vector<T> &vec)
// {
//     for (int i = 0; i < degree; ++i)
//     {
//         int j = vec.size() - degree + i;
//         if (glm::length(vec[i] - vec[j]) > 1e-5)
//         {
//             return false;
//         }
//     }
//     return true;
// }

// /**
//  * Returns whether the 2D array is closed along the u-direction
//  * i.e., along rows.
//  * @param[in] degree_u Degree along u-direction
//  * @param[in] arr 2D array of control points / weights
//  * @return Whether closed along u-direction
//  */

// template <typename T> bool isArray2ClosedU(unsigned int degree_u, const array2<T> &arr)
// {
//     for (int i = 0; i < degree_u; ++i)
//     {
//         for (int j = 0; j < arr.cols(); ++j)
//         {
//             int k = arr.cols() - degree_u + i;
//             if (glm::length(arr(i, j) - arr(k, j)) > 1e-5)
//             {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

// /**
//  * Returns whether the 2D array is closed along the v-direction
//  * i.e., along columns.
//  * @param[in] degree_v Degree along v-direction
//  * @param[in] arr 2D array of control points / weights
//  * @return Whether closed along v-direction
//  */
// template <typename T> bool isArray2ClosedV(unsigned int degree_v, const array2<T> &arr)
// {
//     for (int i = 0; i < arr.rows(); ++i)
//     {
//         for (int j = 0; j < degree_v; j++)
//         {
//             int k = arr.rows() - degree_v + i;
//             if (glm::length(arr(i, j) - arr(i, k)) > 1e-5)
//             {
//                 return false;
//             }
//         }
//     }
//     return true;
// }

} 
#endif 