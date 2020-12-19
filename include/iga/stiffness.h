#ifndef ABASYS_STIFFNESS_H
#define ABASYS_STIFFNESS_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "util/coord.h" 
#include "util/array2.h"
#include "nurbs/basis.h"
#include "nurbs/surface.h"

#include "guass.h"
#include "index.h"
#include "assembly.h"
#include "elematrix.h"

using std::vector;

namespace abab {

template<typename T>
Eigen::MatrixX<T> ElementStiffnessMatrix(
    int ei, int ej, //单元编号
    int udeg, int vdeg,
    const vector<T>& U, const vector<T>& V,
    const array2<vec3<T>>& points, //形状参数
    const array2<T>& weights, 
    T E, T NU, ID id)
{
    int nep = (udeg + 1) * (vdeg + 1); //获取每个单元的控制点数 ,numbers of element's points 
    array2<T> guass = Guass2D<T>(9); //获取高斯点及权值

    Eigen::MatrixX<T> Ke;
    Ke.setZero(2 * nep, 2 * nep);                //初始化单元刚度矩阵

    //遍历高斯点
    for (int g = 0; g < guass.rows(); g++) {
        T u = ((U[ei + 1] - U[ei]) * guass(g, 0) + (U[ei + 1] + U[ei])) / 2;//高斯点对应的u值
        T v = ((V[ej + 1] - V[ej]) * guass(g, 1) + (V[ej + 1] + V[ej])) / 2;//高斯点对应的v值

        array2<int> Lij = LocalIndex2(udeg, vdeg);//获取局部单元的二维索引
        array2<int> Gij = GlobalIndex2(ei, ej, udeg, vdeg); //获取所在单元的二维索引

        //计算高斯点处的基函数及其导数nders,mders (2, udeg + 1)，(2, vdeg + 1)
        array2<T> nders = tinynurbs::bsplineDerBasis<T>(udeg, ei, U, u, 1);
        array2<T> mders = tinynurbs::bsplineDerBasis<T>(vdeg, ej, V, v, 1);

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

        //T x = 0.0, y = 0.0;
        //for (int k = 0; k < nep; k++)
        //{
        //    int gi = Gij(k, 0), gj = Gij(k, 1);
        //    x += R0[k] * points(gi, gj).x;
        //    y += R0[k] * points(gi, gj).y;
        //}

        T dx_du = 0.0, dx_dv = 0.0, dy_du = 0.0, dy_dv = 0.0;
        for (int k = 0; k < nep; k++) {
            int gi = Gij(k, 0), gj = Gij(k, 1);
            dx_du += dR_du[k] * points(gi, gj).x;
            dx_dv += dR_dv[k] * points(gi, gj).x;
            dy_du += dR_du[k] * points(gi, gj).y;
            dy_dv += dR_dv[k] * points(gi, gj).y;
        }

        //计算雅可比行列式 |J|
        T J1 = (U[ei + 1] - U[ei]) / 2.0 * (V[ej + 1] - V[ej]) / 2.0;
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

        //弹性系数矩阵D
        Eigen::Matrix3<T> D = getElastMatrix(E, NU, id);


        //计算高斯点处的单元刚度矩阵，并求和
        Ke += B.transpose() * D * B * J2 * J1 * guass(g, 2);
    }
    return Ke;
}

template<typename T>
Eigen::MatrixX<T> ElementStiffnessMatrix(const tinynurbs::RationalSurface<T>& surf,
    int ei, int ej, T E, T NU, ID id)
{
    return ElementStiffnessMatrix(
        ei, ej,
        surf.degree_u, surf.degree_v,
        surf.knots_u, surf.knots_v, surf.control_points, surf.weights,
        E, NU, id);
}

template<typename T>
Eigen::MatrixX<T> assemblyMatrix(const tinynurbs::RationalSurface<T>& surf, T E, T nu, ID id) 
{
    int deg_u = surf.degree_u;
    int deg_v = surf.degree_v;

    int dofs = surf.control_points.size() * 2;
    int nep = (deg_u + 1) * (deg_v + 1);

    Eigen::MatrixX<T> KK;
    KK.setZero(dofs, dofs);
    Eigen::MatrixX<T> Ke = Eigen::MatrixX<T>::Zero(nep * 2, nep * 2);//初始化单元刚度矩阵Ke

    //求解所有的单元刚度矩阵，并装配
    for (int ei = deg_u; ei < surf.control_points.rows(); ++ei) {
        for (int ej = deg_v; ej < surf.control_points.cols(); ++ej) {
            Ke = ElementStiffnessMatrix(surf, ei, ej, E, nu, id);//求（ei，ej）单元的刚度矩阵
            using std::cout; using std::endl;
            cout << ei << " " << ej << ": \n" << Ke << endl << endl;
            array2<int> idx2s = GlobalIndex2(ei, ej, deg_u, deg_v);//单元（ei，ej）控制点的二维索引
            std::vector<int> idx1s = Index2To1ByV(idx2s, surf.control_points.cols());//控制点按行编号的一维索引
            std::vector<int> assIdxs = getDofsIndex(idx1s);//装配索引
            assembly(Ke, assIdxs, KK);//装配单元（ei，ej）的刚度矩阵
        }
    }
    return KK;
}

template<typename T>
Eigen::SparseMatrix<T> assemblySparse(const tinynurbs::RationalSurface<T>& surf, T E, T nu, ID id) 
{
    int deg_u = surf.degree_u;
    int deg_v = surf.degree_v;

    int dofs = surf.control_points.size() * 2;
    int nep = (deg_u + 1) * (deg_v + 1);

    Eigen::SparseMatrix<T> KK(dofs, dofs);
    vector<Eigen::Triplet<T>> coefficients;
    Eigen::MatrixX<T> Ke = Eigen::MatrixX<T>::Zero(nep * 2, nep * 2);//初始化单元刚度矩阵Ke

    //求解所有的单元刚度矩阵，并装配
    int ki = 1; // 单元编号
    for (int ei = deg_u; ei < surf.control_points.rows(); ++ei) {
        for (int ej = deg_v; ej < surf.control_points.cols(); ++ej) {
            Ke = ElementStiffnessMatrix(surf, ei, ej, E, nu, id);//求（ei，ej）单元的刚度矩阵

            std::cout << "ke " << ki++ << std::endl;
            std::cout << Ke << std::endl << std::endl;

            array2<int> idx2s = GlobalIndex2(ei, ej, deg_u, deg_v);//单元（ei，ej）控制点的二维索引
            std::vector<int> idx1s = Index2To1ByV(idx2s, surf.control_points.cols());//控制点按行编号的一维索引
            std::vector<int> assIdxs = getDofsIndex(idx1s);//装配索引
            for (int i = 0; i < assIdxs.size(); i++) {
                for (int j = 0; j < assIdxs.size(); j++) {
                    coefficients.emplace_back( assIdxs[i], assIdxs[j], Ke(i, j) );
                }
            }
        }
    }
    KK.setFromTriplets(coefficients.begin, coefficients.end());
    return KK;
}


} //namespace abab
#endif //ABASYS_STIFFNESS_H