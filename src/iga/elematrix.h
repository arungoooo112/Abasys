#ifndef ELASTICCOEFFICIENTMATRIX_H
#define ELASTICCOEFFICIENTMATRIX_H

#include <Eigen/Dense>

#include "core/surface.h"
#include "core/basis.h"
#include "util/util.h"
#include "util/coord.h" 
#include "util/array2.h"

#include "guass.h"
#include "index.h"

using namespace Eigen;
using namespace std;

namespace abab {
enum ID {
    PlaneStress, 
    PlaneStrain
};


//获取平面问题弹性系数矩阵D
//@param[in] 平面问题性质指示参数，1为平面应力，2为平面应变
//@param[in] 弹性系数 E
//@param[in] 泊松比
//@return 平面问题弹性系数矩阵D
template<typename T>
Eigen::Matrix<T, 3, 3> getElastMatrix(T E, T nu, ID id) {
    if (id == ID::PlaneStress) 
    {
        Eigen::Matrix<T, 3, 3> D;
        D << 1, nu, 0,
            nu, 1, 0,
            0, 0, (1 - nu) / 2;
        D = D * (E / (1 - nu * nu));
        return D;
    }
    else //if (id == ID::PlaneStrain) 
    {
        Eigen::Matrix<T, 3, 3> D;
        D << 1 - nu, nu, 0,
            nu, 1 - nu, 0,
            0, 0, (1 - 2 * nu) / 2;
        D = D * (E / (1 + nu) / (1 - 2 * nu));
        return D;
    }
}

template<typename T>
MatrixX<T> getShapeMatrix(    
    int deg_u, int deg_v,
    const vector<T>& knots_u, const vector<T>& knots_v,
    const array2<T>& weights,
     T u, T v) 
{
    int nep = (deg_u + 1) * (deg_v + 1); //获取每个单元的控制点数 ,numbers of element's points 
    int ei = tinynurbs::findSpan(deg_u, knots_u, u);
    int ej = tinynurbs::findSpan(deg_v, knots_v, v);

    array2<int> Lij = LocalIndex2(deg_u, deg_v);          //获取局部单元的二维索引
    array2<int> Gij = GlobalIndex2(ei, ej, deg_u, deg_v); //获取所在单元的二维索引

    vector<T> nbasic = tinynurbs::bsplineBasis<T>(deg_u, ei, knots_u, u);
    vector<T> mbasic = tinynurbs::bsplineBasis<T>(deg_v, ej, knots_v, v);

    T sumOfNMW = 0;
    for (int k = 0; k < nep; k++) 
    {
        int li = Lij(k, 0), lj = Lij(k, 1);
        int gi = Gij(k, 0), gj = Gij(k, 1);
        sumOfNMW += nbasic[li] * mbasic[lj] * weights(gi, gj);
    }

    vector<T> R(nep);
    for (int k = 0; k < nep; k++)
    {
        int li = Lij(k, 0), lj = Lij(k, 1);
        int gi = Gij(k, 0), gj = Gij(k, 1);
        R[k] = nbasic[li] * mbasic[lj] * weights(gi, gj) / sumOfNMW;
    }
    
    MatrixX<T> shape;
    shape.setZero(2, 2 * nep);
    for (int k = 0; k < nep; k++) {
         Eigen::MatrixX<T> si(2, 2);
         si << R[k], 0, 0, R[k];
         shape.block(0, 2 * k, 2, 2) = si;
    }
    return shape;
}

template<typename T>
MatrixX<T> getShapeMatrix(    
    const tinynurbs::RationalSurface<T>& surf,
    T u, T v) 
{
    return getShapeMatrix(surf.degree_u, surf.degree_v, surf.knots_u, surf.knots_v, surf.weights, u, v);
}

template<typename T>
MatrixX<T> getStrainMatrix(    
    int deg_u, int deg_v,
    const vector<T>& knots_u, const vector<T>& knots_v,
    const array2<vec3<T>>& points, 
    const array2<T>& weights, 
    T u, T v) 
{   
    int nep = (deg_u + 1) * (deg_v + 1); //获取每个单元的控制点数 ,numbers of element's points 
    int ei = tinynurbs::findSpan(deg_u, knots_u, u);
    int ej = tinynurbs::findSpan(deg_v, knots_v, v);

    array2<int> Lij = LocalIndex2(deg_u, deg_v); //获取局部单元的二维索引
    array2<int> Gij = GlobalIndex2(ei, ej, deg_u, deg_v); //获取所在单元的二维索引

    //计算高斯点处的基函数及其导数nders,mders (2, deg_u + 1)，(2, deg_v + 1)
    array2<T> nders = tinynurbs::bsplineDerBasis<T>(deg_u, ei, knots_u, u, 1);
    array2<T> mders = tinynurbs::bsplineDerBasis<T>(deg_v, ej, knots_v, v, 1);

    //计算 dx/du, dx/dv, dy/du, dy/dv;
    T sumOfNMW = 0, sumOfN1MW = 0, sumOfNM1W = 0;
    for (int k = 0; k < nep; k++) {
        int li = Lij(k, 0), lj = Lij(k, 1);
        int gi = Gij(k, 0), gj = Gij(k, 1);
        sumOfNMW  += nders(0, li) * mders(0, lj) * weights(gi, gj);
        sumOfN1MW += nders(1, li) * mders(0, lj) * weights(gi, gj);
        sumOfNM1W += nders(0, li) * mders(1, lj) * weights(gi, gj);
    }
    std::vector<T> dR_du(nep), dR_dv(nep), R0(nep);
    for (int k = 0; k < nep; k++) {
        int li = Lij(k, 0), lj = Lij(k, 1);
        int gi = Gij(k, 0), gj = Gij(k, 1);
        R0[k] = nders(0, li) * mders(0, lj) * weights(gi, gj) / sumOfNMW;
        dR_du[k] = (nders(1, li) * mders(0, lj) * weights(gi, gj) * sumOfNMW
            - nders(0, li) * mders(0, lj) * weights(gi, gj) * sumOfN1MW) / sumOfNMW / sumOfNMW;
        dR_dv[k] = (nders(0, li) * mders(1, lj) * weights(gi, gj) * sumOfNMW
            - nders(0, li) * mders(0, lj) * weights(gi, gj) * sumOfNM1W) / sumOfNMW / sumOfNMW;
    }

    T dx_du = 0.0, dx_dv = 0.0, dy_du = 0.0, dy_dv = 0.0;
    for (int k = 0; k < nep; k++) {
        int gi = Gij(k, 0), gj = Gij(k, 1);
        dx_du += dR_du[k] * points(gi, gj).x;
        dx_dv += dR_dv[k] * points(gi, gj).x;
        dy_du += dR_du[k] * points(gi, gj).y;
        dy_dv += dR_dv[k] * points(gi, gj).y;
    }

    //计算雅可比行列式 |J|
    T J1 = (knots_u[ei + 1] - knots_u[ei]) / 2.0 * (knots_v[ej + 1] - knots_v[ej]) / 2.0;
    T J2 = dx_du * dy_dv - dx_dv * dy_du;

    //计算 dR/dx, dR/dy;
    std::vector<T> dR_dx(nep), dR_dy(nep);
    for (int k = 0; k < nep; k++) {
        dR_dx[k] = (dy_dv * dR_du[k] - dy_du * dR_dv[k]) / J2;
        dR_dy[k] = (dx_du * dR_dv[k] - dx_dv * dR_du[k]) / J2;
    }

    //计算几何矩阵B, 下一步，修改B矩阵的表示，与后续中自由度索引相联系
    Eigen::MatrixX<T> B(3, 2 * nep);
    for (int k = 0; k < nep; k++) {
        Eigen::MatrixX<T> Bi(3, 2);
        Bi <<   dR_dx[k],        0,
                        0, dR_dy[k], 
                dR_dy[k], dR_dx[k];
        B.block(0, 2 * k, 3, 2) = Bi;
    }
    return B;
}

template<typename T>
MatrixX<T> getStrainMatrix(    
    const tinynurbs::RationalSurface<T>& surf,
    T u, T v) 
{
    return getStrainMatrix(surf.degree_u, surf.degree_v, surf.knots_u, surf.knots_v, surf.control_points, surf.weights, u, v);
}
//} //namespace iga


// template<typename T>
// Eigen::Matrix<T, 3, 3> PlaneStressMatrix(T E, T nu) {
//     Eigen::Matrix<T, 3, 3> D;
//     D << 1, nu, 0,
//         nu, 1, 0,
//         0, 0, (1 - nu) / 2;
//     D = D * (E / (1 - nu * nu));
//     return D;
// }

// template<typename T>
// Eigen::Matrix<T, 3, 3> PlaneStrainMatrix(T E, T nu) {
//     Eigen::Matrix<T, 3, 3> D;
//     D << 1 - nu, nu, 0,
//         nu, 1 - nu, 0,
//         0, 0, (1 - 2 * nu) / 2;
//     D = D * (E / (1 + nu) / (1 - 2 * nu));
//     return D;
// }

} //namespace abab
#endif //ELASTICCOEFFICIENTMATRIX_H