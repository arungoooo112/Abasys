/**
 * The Surface and RationalSurface classes represent non-rational and rational
 * NURBS surfaces, respectively.
  *
 * Use of this source code is governed by a BSD-style license that can be found in
 * the LICENSE file.
 */

#ifndef TINYNURBS_SURFACE_H
#define TINYNURBS_SURFACE_H

#include <stdexcept>
#include <vector>

#include "util/util.h"
#include "util/array2.h"
#include "util/coord.h"
#include "refine.h"

namespace tinynurbs
{

// Forward declaration
template <typename T> struct RationalSurface;

/**
Struct for representing a non-rational NURBS surface
\tparam T Data type of control points and weights (float or double)
*/
template <typename T> struct Surface
{
    unsigned int degree_u, degree_v;
    vector<T> knots_u, knots_v;
    array2<vec3<T>> control_points;

    Surface() = default;
    Surface(const RationalSurface<T> &srf)
        : degree_u(srf.degree_u), degree_v(srf.degree_v), knots_u(srf.knots_u),
          knots_v(srf.knots_v), control_points(srf.control_points)
    {
    }
    Surface(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
            const std::vector<T> &knots_v, array2<vec3<T>> control_points)
        : degree_u(degree_u), degree_v(degree_v), knots_u(knots_u), knots_v(knots_v),
          control_points(control_points)
    {
    }
};

/**
Struct for representing a rational NURBS surface
\tparam T Data type of control points and weights (float or double)
*/
template <typename T> struct RationalSurface
{
    int degree_u, degree_v;
    vector<T> knots_u, knots_v;
    array2<vec3<T>> control_points;
    array2<T> weights;

    RationalSurface() = default;
    RationalSurface(const Surface<T> &srf, const array2<T> &weights)
        : degree_u(srf.degree_u), degree_v(srf.degree_v), knots_u(srf.knots_u),
          knots_v(srf.knots_v), control_points(srf.control_points), weights(weights)
    {
    }
    RationalSurface(const Surface<T> &srf)
        : RationalSurface(srf, array2<T>(srf.control_points.rows(), srf.control_points.cols(), T(1.0)))
    {
    }
    RationalSurface(unsigned int degree_u, unsigned int degree_v, const std::vector<T> &knots_u,
                    const std::vector<T> &knots_v, const array2<vec3<T>> &control_points,
                    const array2<T> &weights)
        : degree_u(degree_u), degree_v(degree_v), knots_u(knots_u), knots_v(knots_v),
          control_points(control_points), weights(weights)
    {
    }

public: 
    // n1, n2 每个单元再细化后的个数
   void HRefine(int n1, int n2);
};

template <typename T> 
void RationalSurface<T>::HRefine(int n1, int n2) 
{
    array2<hvec3<T>> hpoints = cartesianToHomogenous(control_points, weights); 
    vector<T> new_knotsu;
    vector<T> new_knotsv;
    array2<hvec3<T>> new_hpoints;
    auto X = getinsertKnotVect(knots_u, n1 - 1);
    auto Y = getinsertKnotVect(knots_v, n2 - 1);
    tinynurbs::RefineKnotVectSurface(degree_u, degree_v, knots_u, knots_v, hpoints,
                            X, Y, new_knotsu, new_knotsv, new_hpoints);
    knots_u = new_knotsu;
    knots_v = new_knotsv;
    tinynurbs::homogenousToCartesian(new_hpoints, control_points, weights);
}

// template <typename T> struct HSurface
// {
//     int degree_u, degree_v;
//     vector<T> knots_u, knots_v;
//     array2<hvec3<T>> control_points;

//     HSurface() = default;
//     HSurface(const RationalSurface<T>& RSurface) 
//         : degree_u(RSurface.degree_u), degree_v(RSurface.degree_v), knots_u(RSurface.knots_u), knots_v(RSurface.knots_v) {
//             control_points.resize(RSurface.control_points.rows(), RSurface.control_points.cols(), hvec3<T>(0));
//             for (int i = 0; i < RSurface.control_points.rows(); i++) {
//                 for (int j = 0; j < RSurface.control_points.cols(); j++) {
//                     control_points(i, j).x = RSurface.control_points(i, j).x * RSurface.weights(i, j);
//                     control_points(i, j).y = RSurface.control_points(i, j).y * RSurface.weights(i, j);
//                     control_points(i, j).z = RSurface.control_points(i, j).z * RSurface.weights(i, j);
//                     control_points(i, j).w = RSurface.weights(i, j);
//                 }
//             }
//         }

//     HSurface(int degree_u, int degree_v, const vector<T> &knots_u,
//             const vector<T> &knots_v, array2<hvec3<T>> control_points)
//         : degree_u(degree_u), degree_v(degree_v), knots_u(knots_u), knots_v(knots_v),
//           control_points(control_points) {}

//     // void hRefine(int n1, int n2) {
//     //     vector<T> new_knotsu;
//     //     vector<T> new_knotsv;
//     //     array2<hvec3<T>> new_CtrlPoints;
//     //     RefineKnotVectSurface(degree_u, degree_v, knots_u, knots_v, control_points,
//     //                           getinsertKnotVect(knots_u, n1 - 1), getinsertKnotVect(knots_v, n2 - 1),
//     //                           new_knotsu,
//     //                           new_knotsv, new_CtrlPoints);
//     //     knots_u = new_knotsu;
//     //     knots_v = new_knotsv;
//     //     control_points = new_CtrlPoints;
//     // }
// };

// Typedefs for ease of use
// typedef HSurface<float>  HSurface3f;
// typedef HSurface<double> HSurface3d;
typedef Surface<float>  Surface3f;
typedef Surface<double> Surface3d;
typedef RationalSurface<float>  RSurface3f;
typedef RationalSurface<double> RSurface3d;
typedef RationalSurface<float>  RationalSurface3f;
typedef RationalSurface<double> RationalSurface3d;

} // namespace tinynurbs

#endif // TINYNURBS_SURFACE_H
