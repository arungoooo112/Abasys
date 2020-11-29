#ifndef ABASYS_IMPOSEBC_H
#define ABASYS_IMPOSEBC_H

#include <vector>
#include <Eigen/Dense>

#include "nurbs/surface.h"
#include "nurbs/curve.h"
#include "boundary.h"
#include "./coord.h"
#include "./array2.h"
#include "assembly.h"
#include "guass.h"

using std::vector;
using Eigen::VectorX;

template<typename T, typename Function>
VectorX<T> calcElementeEquivalentVals2D(int ei, int deg, 
    const vector<T>& U, const vector<vec3<T>> points, 
    const vector<T>& weights, Function fun) 
{
    
    auto guass = Guass1D<T>(4); //获取高斯点及权值
    int nep = deg + 1;
    VectorX<T> res = Eigen::VectorX<T>::Zero(nep);

    for (int g = 0; g < guass.rows(); g++) {
        T u = ((U[ei + 1] - U[ei]) * guass(g, 0) + (U[ei + 1] + U[ei])) / 2;
        vector<int> Li = LocalIndex1(deg);
        vector<int> Gi = GlobalIndex1(ei, deg);
        array2<T> nders = tinynurbs::bsplineDerBasis<T>(deg, ei, U, u, 1);

        T sumOfN0W = 0, sumOfN1W = 0;
        for (int k = 0; k < nep; k++) {
            int lk = Li[k];
            int gk = Gi[k];
            sumOfN0W += nders(0, lk) * weights[gk];
            sumOfN1W += nders(1, lk) * weights[gk];
        }

        VectorX<T> R0(nep), R1(nep) ; //下一步，后续计算可以采用矩阵乘法
        for (int k = 0; k < nep; k++) {
            int lk = Li[k];
            int gk = Gi[k];
            R0[k] = nders(0, lk) * weights[gk] / sumOfN0W;
            R1[k] = (nders(1, lk) * weights[gk] * sumOfN0W - nders(0, lk) * weights[gk] * sumOfN1W) 
                     / sumOfN0W / sumOfN0W;
        }
        T dx_du = 0.0, dy_du = 0.0;
        for (int k = 0; k < nep; k++) {
            int gk = Gi[k];
            dx_du += R1[k] * points[gk].x; 
            dy_du += R1[k] * points[gk].y;
        }

        T x = 0.0, y = 0.0;
        for (int k = 0; k < nep; k++) {
            int gk = Gi[k];
            x += R0[k] * points[gk].x; 
            y += R0[k] * points[gk].y;
        }
        T Jx = (U[ei + 1] - U[ei]) / 2.0;
        T J1 = sqrt(dx_du * dx_du + dy_du * dy_du);
        res += guass(g, 1) * R0 * fun(x, y) * J1 * Jx;
    }
    return res;
}

template<typename T, typename Function>
VectorX<T> calcEquivalentVals2D(int deg, 
    const vector<T>& U, const vector<vec3<T>> points, 
    const vector<T>& weights, Function fun) 
{
    VectorX<T> res = Eigen::VectorX<T>::Zero(points.size());
    for (int ei = deg; ei < points.size(); ei++) {
        VectorX<T> eleVals = calcElementeEquivalentVals2D(ei, deg, U, points, weights, fun);
        vector<int> dofsIdx = GlobalIndex1(ei, deg);
        assembly(eleVals, dofsIdx, res);
    }
    return res;
}


template<typename T, typename Function>
VectorX<T> calcEquivalentVals2D(const tinynurbs::RationalCurve<T>& crv, Function fun)
{
    return calcEquivalentVals2D(crv.degree, crv.knots, crv.control_points, crv.weights, fun);
}

template<typename T, typename Function>
VectorX<T> applyNewmannBdryVals(const tinynurbs::RationalSurface<T>& srf, Boundary b, Function fun)
{
    return calcEquivalentVals2D(getBoundary(srf, b), fun);
}

template<typename T, typename Function>
pair<VectorX<T>, vector<int>> applyNewmannBC(const tinynurbs::RationalSurface<T>& srf, Boundary b, Function fun, int lb)
{
    VectorX<T> vals = applyNewmannBdryVals(srf, b, fun);
    vector<int> idxs = getSurfaceBdryIdxs(srf, b);
    if (lb == 0) Index1ToDofX(idxs);
    if (lb == 1) Index1ToDofY(idxs);
   return {vals, idxs};
}

template<typename T, typename Function>
VectorX<T> applyDrchltBdryVals(const tinynurbs::RationalSurface<T>& srf, Boundary b, Function fun)
{
    return calcEquivalentVals2D(getBoundary(srf, b), fun);
}

template<typename T, typename Function>
pair<VectorX<T>, vector<int>> applyDrchltBC(const tinynurbs::RationalSurface<T>& srf, Boundary b, Function fun, int lb) 
{
    auto vals = applyDrchltBdryVals(srf, b, fun);
    auto idxs = getSurfaceBdryIdxs(srf, b);
    if (lb == 0) Index1ToDofX(idxs);
    if (lb == 1) Index1ToDofY(idxs);
    return {vals, idxs};
}

template<typename T>
vector<int> getSurfaceBdryIdxsU0(const tinynurbs::RationalSurface<T>& srf)
{
    vector<int> res;
    for (int i = 0; i < srf.control_points.cols(); i++) {
        res.push_back(i);
    }
    return res;
}

template<typename T>
vector<int> getSurfaceBdryIdxsU1(const tinynurbs::RationalSurface<T>& srf)
{
    vector<int> res;
    int beg = srf.control_points.size() - srf.control_points.cols();
    for (int i = 0; i < srf.control_points.cols(); i++) {
        res.push_back(beg + i);
    }
    return res;
}

template<typename T>
vector<int> getSurfaceBdryIdxsV0(const tinynurbs::RationalSurface<T>& srf)
{
    vector<int> res;
    int cols = srf.control_points.cols();
    for (int i = 0; i < srf.control_points.rows(); i++) {
        res.push_back(i * cols);
    }
    return res;
}

template<typename T>
vector<int> getSurfaceBdryIdxsV1(const tinynurbs::RationalSurface<T>& srf)
{
    vector<int> res;
    int cols = srf.control_points.cols();
    for (int i = 0; i < srf.control_points.rows(); i++) {
        res.push_back(i * cols + cols - 1);
    }
    return res;
}


template <typename T>
vector<int> getSurfaceBdryIdxs(const tinynurbs::RationalSurface<T>& srf, Boundary b)
{
    switch (b)
    {
    case U0: return getSurfaceBdryIdxsU0(srf);
    case U1: return getSurfaceBdryIdxsU1(srf);
    case V0: return getSurfaceBdryIdxsV0(srf);
    case V1: return getSurfaceBdryIdxsV1(srf);
    default: break;
    }
    return {};
}

template<typename T>
vector<int> getBoundaryDofs(const tinynurbs::RationalSurface<T>& srf, Boundary b, int d)
{
    vector<int> res = getSurfaceBdryIdxs(srf, b);
    if (d == 1) {
        for (int& a : res) {
            a = 2 * a;
        }
    }
    if (d == 2) {
        for (int& a : res) {
            a = 2 * a + 1;
        }
    }
    return res;
}


#endif //ABASYS_IMPOSEBC_H
