#ifndef ABASYS_ASSEMBLY_H
#define ABASYS_ASSEMBLY_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

namespace abab {

//按系数表简化向量
//@param[in] 简化的系数表
//@param[in] 稀疏的向量
//@return 简化后的向量
template <typename T> Eigen::VectorX<T> simplyfy(const vector<int> &vec, const Eigen::VectorX<T> &FF)
{
    Eigen::VectorX<T> res = Eigen::VectorX<T>::Zero(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        res[i] = FF(vec[i]);
    }
    return res;
}

//简化矩阵，用于简化整体刚度矩阵
//@param[in] 简化的系数表
//@param[in] 稀疏矩阵
//@return 简化后的矩阵
template <typename T> Eigen::SparseMatrix<T> simplyfy(const vector<int> &vec, const Eigen::SparseMatrix<T> &KK)
{
    Eigen::SparseMatrix<T> res(vec.size(), vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            res.coeffRef(i, j) = res.coeffRef(j, i) = KK.coeff(vec[i], vec[j]);
        }
        res.coeffRef(i, i) = KK.coeff(vec[i], vec[i]);
    }
    res.makeCompressed();
    return res;
}

template <typename T> Eigen::MatrixX<T> simplyfy(const vector<int> &vec, const Eigen::MatrixX<T> &KK)
{
    Eigen::MatrixX<T> res = Eigen::MatrixX<T>::Zero(vec.size(), vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        for (int j = 0; j < i; j++)
        {
            res(i, j) = res(j, i) = KK(vec[i], vec[j]);
        }
        res(i, i) = KK(vec[i], vec[i]);
    }
    return res;
}
    //装配矩阵，用于装配单元刚度矩阵到整体刚度矩阵
    //@param[in] 单元矩阵Ke
    //@param[in] 装配的系数表
    //@param[in] 整体矩阵K
    //@param[out] 整体矩阵K
    template<typename T>
    void assembly(const Eigen::MatrixX<T>& Ke, const std::vector<int>& vec, Eigen::MatrixX<T>& KK) {
        assert(Ke.rows() == Ke.cols() && KK.rows() == KK.cols() && Ke.rows() == vec.size());
        for (int i = 0; i < vec.size(); i++) {
            for (int j = 0; j < i; j++) { 
                KK(vec[j], vec[i]) = ( KK(vec[i], vec[j]) += Ke(i, j) );
            }
            KK(vec[i], vec[i]) += Ke(i, i);
        }
    }

    // template<typename T>
    // void assembly(const Eigen::MatrixX<T>& Ke, const std::vector<int>& vec, Eigen::SparseMatrix<T>& KK) {
    //     assert(Ke.rows() == Ke.cols() && KK.rows() == KK.cols() && Ke.rows() == vec.size());
    //     for (int i = 0; i < vec.size(); i++) {
    //         for (int j = 0; j < i; j++) { 
    //             KK.coeffRef(vec[j], vec[i]) = ( KK.coeffRef(vec[i], vec[j]) += Ke(i, j) );
    //         }
    //         KK.coeffRef(vec[i], vec[i]) += Ke(i, i);
    //     }
    // }

    //装配向量，用于装配单元力或位移向量到整体力或位移向量
    //@param[in] 单元向量
    //@param[in] 装配的系数表
    //@param[in] 整体向量
    //@param[out] 整体向量
    template<typename T>
    void assembly(const Eigen::VectorX<T>& Fe, const std::vector<int>& vec, Eigen::VectorX<T>& FF) {
        assert(Fe.size() == vec.size());
        for (int i = 0; i < vec.size(); i++) {
            FF(vec[i]) += Fe(i);
        }
    }

    template <typename T>
    void assembly(const vector<std::pair<Eigen::VectorX<T>, vector<int>>>& containor, Eigen::VectorX<T> &FF)
    {
        for (const auto & c : containor)
        {
            assembly(c.first, c.second, FF);
        }
    }





} //namespace abab


#endif //ABASYS_ASSEMBLY_H