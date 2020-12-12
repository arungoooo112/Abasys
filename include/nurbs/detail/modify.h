/**
 * Functions for modifying NURBS curves and surfaces.
 * 
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef TINYNURBS_MODIFY_H
#define TINYNURBS_MODIFY_H

#include <tuple>
#include <vector>
#include "util/coord.h"
#include "check.h"

#include "util/util.h"


namespace funs
{

/////////////////////////////////////////////////////////////////////


/**
 * Insert knots in the curve
 * @param[in] deg Degree of the curve
 * @param[in] knots Knot vector of the curve
 * @param[in] cp Control points of the curve
 * @param[in] u Parameter to insert knot(s) at
 * @param[in] r Number of times to insert knot
 * @param[out] new_knots Updated knot vector
 * @param[out] new_cp Updated control points
 */
template <typename T>
void curveKnotInsert(unsigned int deg, const std::vector<T> &knots,
                     const std::vector<hvec3<T>> &cp, T u, unsigned int r,
                     std::vector<T> &new_knots, std::vector<hvec3<T>> &new_cp)
{
    int k = findSpan(deg, knots, u);
    unsigned int s = knotMultiplicity(knots, k);
    if (s == deg)
    {
        return;
    }
    if ((r + s) > deg)
    {
        r = deg - s;
    }

    // Insert new knots between span and (span + 1)
    new_knots.resize(knots.size() + r);
    for (int i = 0; i < k + 1; ++i)
    {
        new_knots[i] = knots[i];
    }
    for (unsigned int i = 1; i < r + 1; ++i)
    {
        new_knots[k + i] = u;
    }
    for (int i = k + 1; i < knots.size(); ++i)
    {
        new_knots[i + r] = knots[i];
    }
    // Copy unaffected control points
    new_cp.resize(cp.size() + r);
    for (int i = 0; i < k - deg + 1; ++i)
    {
        new_cp[i] = cp[i];
    }
    for (int i = k - s; i < cp.size(); ++i)
    {
        new_cp[i + r] = cp[i];
    }
    // Copy affected control points
    std::vector<hvec3<T>> tmp;
    tmp.resize(deg - s + 1);
    for (int i = 0; i < deg - s + 1; ++i)
    {
        tmp[i] = cp[k - deg + i];
    }
    // Modify affected control points
    for (int j = 1; j <= r; ++j)
    {
        int L = k - deg + j;
        for (int i = 0; i < deg - j - s + 1; ++i)
        {
            T a = (u - knots[L + i]) / (knots[i + k + 1] - knots[L + i]);
            tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
        }
        new_cp[L] = tmp[0];
        new_cp[k + r - j - s] = tmp[deg - j - s];
    }
    int L = k - deg + r;
    for (int i = L + 1; i < k - s; ++i)
    {
        new_cp[i] = tmp[i - L];
    }
}

/**
 * Insert knots in the surface along one direction
 * @param[in] degree Degree of the surface along which to insert knot
 * @param[in] knots Knot vector
 * @param[in] cp 2D array of control points
 * @param[in] knot Knot value to insert
 * @param[in] r Number of times to insert
 * @param[in] along_u Whether inserting along u-direction
 * @param[out] new_knots Updated knot vector
 * @param[out] new_cp Updated control points
 */
template <typename T>
void surfaceKnotInsert(unsigned int degree, const std::vector<T> &knots,
                       const array2<hvec3<T>> &cp, T knot, unsigned int r, bool along_u,
                       std::vector<T> &new_knots, array2<hvec3<T>> &new_cp)
{
    int span = findSpan(degree, knots, knot);
    unsigned int s = knotMultiplicity(knots, span);
    if (s == degree)
    {
        return;
    }
    if ((r + s) > degree)
    {
        r = degree - s;
    }

    // Create a new knot vector
    new_knots.resize(knots.size() + r);
    for (int i = 0; i <= span; ++i)
    {
        new_knots[i] = knots[i];
    }
    for (int i = 1; i <= r; ++i)
    {
        new_knots[span + i] = knot;
    }
    for (int i = span + 1; i < knots.size(); ++i)
    {
        new_knots[i + r] = knots[i];
    }
    // Compute alpha
    array2<T> alpha(degree - s, r + 1, T(0));
    for (int j = 1; j <= r; ++j)
    {
        int L = span - degree + j;
        for (int i = 0; i <= degree - j - s; ++i)
        {
            alpha(i, j) = (knot - knots[L + i]) / (knots[i + span + 1] - knots[L + i]);
        }
    }

    // Create a temporary container for affected control points per row/column
    std::vector<hvec3<T>> tmp(degree + 1);

    if (along_u)
    {
        // Create new control points with additional rows
        new_cp.resize(cp.rows() + r, cp.cols());

        // Update control points
        // Each row is a u-isocurve, each col is a v-isocurve
        for (int col = 0; col < cp.cols(); ++col)
        {
            // Copy unaffected control points
            for (int i = 0; i <= span - degree; ++i)
            {
                new_cp(i, col) = cp(i, col);
            }
            for (int i = span - s; i < cp.rows(); ++i)
            {
                new_cp(i + r, col) = cp(i, col);
            }
            // Copy affected control points to temp array
            for (int i = 0; i < degree - s + 1; ++i)
            {
                tmp[i] = cp(span - degree + i, col);
            }
            // Insert knot
            for (int j = 1; j <= r; ++j)
            {
                int L = span - degree + j;
                for (int i = 0; i <= degree - j - s; ++i)
                {
                    T a = alpha(i, j);
                    tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                }
                new_cp(L, col) = tmp[0];
                new_cp(span + r - j - s, col) = tmp[degree - j - s];
            }
            int L = span - degree + r;
            for (int i = L + 1; i < span - s; ++i)
            {
                new_cp(i, col) = tmp[i - L];
            }
        }
    }
    else
    {
        // Create new control points with additional columns
        new_cp.resize(cp.rows(), cp.cols() + r);

        // Update control points
        // Each row is a u-isocurve, each col is a v-isocurve
        for (int row = 0; row < cp.rows(); ++row)
        {
            // Copy unaffected control points
            for (int i = 0; i <= span - degree; ++i)
            {
                new_cp(row, i) = cp(row, i);
            }
            for (int i = span - s; i < cp.cols(); ++i)
            {
                new_cp(row, i + r) = cp(row, i);
            }
            // Copy affected control points to temp array
            for (int i = 0; i < degree - s + 1; ++i)
            {
                tmp[i] = cp(row, span - degree + i);
            }
            // Insert knot
            for (int j = 1; j <= r; ++j)
            {
                int L = span - degree + j;
                for (int i = 0; i <= degree - j - s; ++i)
                {
                    T a = alpha(i, j);
                    tmp[i] = (1 - a) * tmp[i] + a * tmp[i + 1];
                }
                new_cp(row, L) = tmp[0];
                new_cp(row, span + r - j - s) = tmp[degree - j - s];
            }
            int L = span - degree + r;
            for (int i = L + 1; i < span - s; ++i)
            {
                new_cp(row, i) = tmp[i - L];
            }
        }
    }
}

/**
 * Split the curve into two
 * @param[in] degree Degree of curve
 * @param[in] knots Knot vector
 * @param[in] control_points Array of control points
 * @param[in] u Parameter to split curve
 * @param[out] left_knots Knots of the left part of the curve
 * @param[out] left_control_points Control points of the left part of the curve
 * @param[out] right_knots Knots of the right part of the curve
 * @param[out] right_control_points Control points of the right part of the curve
 */
template <typename T>
void curveSplit(unsigned int degree, const std::vector<T> &knots,
                const std::vector<hvec3<T>> &control_points, T u,
                std::vector<T> &left_knots, std::vector<hvec3<T>> &left_control_points,
                std::vector<T> &right_knots, std::vector<hvec3<T>> &right_control_points)
{
    std::vector<T> tmp_knots;
    std::vector<hvec3<T>> tmp_cp;

    int span = findSpan(degree, knots, u);
    int r = degree - knotMultiplicity(knots, span);

    internal::curveKnotInsert(degree, knots, control_points, u, r, tmp_knots, tmp_cp);

    left_knots.clear();
    right_knots.clear();
    left_control_points.clear();
    right_control_points.clear();

    int span_l = findSpan(degree, tmp_knots, u) + 1;
    for (int i = 0; i < span_l; ++i)
    {
        left_knots.push_back(tmp_knots[i]);
    }
    left_knots.push_back(u);

    for (int i = 0; i < degree + 1; ++i)
    {
        right_knots.push_back(u);
    }
    for (int i = span_l; i < tmp_knots.size(); ++i)
    {
        right_knots.push_back(tmp_knots[i]);
    }

    int ks = span - degree + 1;
    for (int i = 0; i < ks + r; ++i)
    {
        left_control_points.push_back(tmp_cp[i]);
    }
    for (int i = ks + r - 1; i < tmp_cp.size(); ++i)
    {
        right_control_points.push_back(tmp_cp[i]);
    }
}

/**
 * Split the surface into two along given parameter direction
 * @param[in] degree Degree of surface along given direction
 * @param[in] knots Knot vector of surface along given direction
 * @param[in] control_points Array of control points
 * @param[in] param Parameter to split curve
 * @param[in] along_u Whether the direction to split along is the u-direction
 * @param[out] left_knots Knots of the left part of the curve
 * @param[out] left_control_points Control points of the left part of the curve
 * @param[out] right_knots Knots of the right part of the curve
 * @param[out] right_control_points Control points of the right part of the curve
 */
template <typename T>
void surfaceSplit(unsigned int degree, const std::vector<T> &knots,
                  const array2<hvec3<T>> &control_points, T param, bool along_u,
                  std::vector<T> &left_knots, array2<hvec3<T>> &left_control_points,
                  std::vector<T> &right_knots, array2<hvec3<T>> &right_control_points)
{
    std::vector<T> tmp_knots;
    array2<hvec3<T>> tmp_cp;

    int span = findSpan(degree, knots, param);
    unsigned int r = degree - knotMultiplicity(knots, span);
    internal::surfaceKnotInsert(degree, knots, control_points, param, r, along_u, tmp_knots,
                                tmp_cp);

    left_knots.clear();
    right_knots.clear();

    int span_l = findSpan(degree, tmp_knots, param) + 1;
    for (int i = 0; i < span_l; ++i)
    {
        left_knots.push_back(tmp_knots[i]);
    }
    left_knots.push_back(param);

    for (int i = 0; i < degree + 1; ++i)
    {
        right_knots.push_back(param);
    }
    for (int i = span_l; i < tmp_knots.size(); ++i)
    {
        right_knots.push_back(tmp_knots[i]);
    }

    int ks = span - degree + 1;
    if (along_u)
    {
        size_t ii = 0;
        left_control_points.resize(ks + r, tmp_cp.cols());
        for (int i = 0; i < ks + r; ++i)
        {
            for (int j = 0; j < tmp_cp.cols(); ++j)
            {
                left_control_points[ii++] = tmp_cp(i, j);
            }
        }
        ii = 0;
        right_control_points.resize(tmp_cp.rows() - ks - r + 1, tmp_cp.cols());
        for (int i = ks + r - 1; i < tmp_cp.rows(); ++i)
        {
            for (int j = 0; j < tmp_cp.cols(); ++j)
            {
                right_control_points[ii++] = tmp_cp(i, j);
            }
        }
    }
    else
    {
        size_t ii = 0;
        left_control_points.resize(tmp_cp.rows(), ks + r);
        for (int i = 0; i < tmp_cp.rows(); ++i)
        {
            for (int j = 0; j < ks + r; ++j)
            {
                left_control_points[ii++] = tmp_cp(i, j);
            }
        }
        ii = 0;
        right_control_points.resize(tmp_cp.rows(), tmp_cp.cols() - ks - r + 1);
        for (int i = 0; i < tmp_cp.rows(); ++i)
        {
            for (int j = ks + r - 1; j < tmp_cp.cols(); ++j)
            {
                right_control_points[ii++] = tmp_cp(i, j);
            }
        }
    }
}

} // namespace internal

#endif 