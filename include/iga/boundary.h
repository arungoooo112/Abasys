#ifndef TINYNURBS_BOUNDARY_H
#define TINYNURBS_BOUNDARY_H

#include <vector>
#include <stdexcept>
#include "nurbs/surface.h"
#include "nurbs/curve.h"
#include "util/coord.h"

namespace abab{

enum Boundary {U0, U1, V0, V1};

template<typename T> tinynurbs::RationalCurve<T> getBoundaryU0(const tinynurbs::RationalSurface<T>& surf)
{
    unsigned int degree = surf.degree_v;
    std::vector<T> knots = surf.knots_v;
    std::vector<vec3<T>> points;
    std::vector<T> weights;
    for (int c = 0; c < surf.control_points.cols(); c++) {
        points.push_back(surf.control_points(0, c));
        weights.push_back(surf.weights(0, c));
    }
    return tinynurbs::RationalCurve<T>(degree, knots, points, weights);
}

template<typename T> tinynurbs::RationalCurve<T> getBoundaryU1(const tinynurbs::RationalSurface<T> &surf)
{
    unsigned int degree = surf.degree_v;
    std::vector<T> knots = surf.knots_v;
    std::vector<vec3<T>> points;
    std::vector<T> weights;
    int rn = surf.control_points.rows() - 1;
    for (int c = 0; c < surf.control_points.cols(); c++)
    {
        points.push_back(surf.control_points(rn, c));
        weights.push_back(surf.weights(rn, c));
    }
    return tinynurbs::RationalCurve<T>(degree, knots, points, weights);
}


template<typename T>
tinynurbs::RationalCurve<T> getBoundaryV0(const tinynurbs::RationalSurface<T>& surf)
{
    unsigned int degree = surf.degree_u;
    std::vector<T> knots = surf.knots_u;
    std::vector<vec3<T>> points;
    std::vector<T> weights;
    for (int r = 0; r < surf.control_points.rows(); r++) {
        points.push_back(surf.control_points(r, 0));
        weights.push_back(surf.weights(r, 0));
    }
    return tinynurbs::RationalCurve<T>(degree, knots, points, weights);
}


template<typename T>
tinynurbs::RationalCurve<T> getBoundaryV1(const tinynurbs::RationalSurface<T>& surf)
{
    unsigned int degree = surf.degree_u;
    std::vector<T> knots = surf.knots_u;
    std::vector<vec3<T>> points;
    std::vector<T> weights;
    int cn = surf.control_points.cols() - 1;
    for (int r = 0; r < surf.control_points.rows(); r++) {
        points.push_back(surf.control_points(r, cn));
        weights.push_back(surf.weights(r, cn));
    }
    return tinynurbs::RationalCurve<T>(degree, knots, points, weights);
}

template<typename T>
tinynurbs::RationalCurve<T> getBoundary(const tinynurbs::RationalSurface<T>& surf, Boundary b)
{
    switch (b)
    {
    case U0:
        return getBoundaryU0(surf);
        break;
    case U1:
        return getBoundaryU1(surf);
        break;
    case V0:
        return getBoundaryV0(surf);
        break;
    case V1:
        return getBoundaryV1(surf);
        break;
    default:
        throw std::runtime_error("invalid boundry in surface");
    }
}



} //namespace abab


#endif //TINYNURBS_BOUNDARY_H