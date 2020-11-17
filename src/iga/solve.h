#ifndef TINYIGA_SOLVE_H
#define TINYIGA_SOLVE_H

#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "iga/index.h"
#include "iga/assembly.h"

using namespace Eigen;
using namespace std;

namespace abab
{



template<typename T>
void solve(const MatrixX<T>& KK, VectorX<T> &U, 
            VectorX<T>& F, vector<int> & knowns) 
{
    sort(knowns.begin(), knowns.end());
    knowns.erase(unique(knowns.begin(), knowns.end()), knowns.end());
    
    int n = KK.rows();
    vector<int> other;
    int a = 0;
    for (int i = 0; i < n; ++i) {
        if (i != knowns[a]) other.push_back(i);
        else a++;
    }

    int nsimp = other.size();
    MatrixX<T> Ktemp(nsimp, nsimp);//简化后的整体刚度矩阵
    VectorX<T> Ftemp(nsimp);//简化后的载荷向量
    VectorX<T> Utemp(nsimp);//简化后的位移向量

    Ktemp = simplyfy(other, KK);
    Ftemp = simplyfy(other, F);

    //求解
    //Eigen::SparseLU<MatrixX<T>> solver;
    //solver.compute(Ktemp);
    Utemp = Ktemp.fullPivLu().solve(Ftemp);
    Ftemp = Ktemp * Utemp;

    assembly(Utemp, other, U);    //装配位移向量
    assembly(Ftemp, other, F);    //装配载荷向量

}

template <typename T>
void solve(const SparseMatrix<T> &KK, VectorX<T> &U, VectorX<T> &F, vector<int> &knowns)
{
    sort(knowns.begin(), knowns.end());
    knowns.erase(unique(knowns.begin(), knowns.end()), knowns.end());

    int n = KK.rows();
    vector<int> other;
    int a = 0;
    for (int i = 0; i < n; ++i)
    {
        if (i != knowns[a])
            other.push_back(i);
        else
            a++;
    }

    int nsimp = other.size();
    SparseMatrix<T> Ktemp(nsimp, nsimp); //简化后的整体刚度矩阵
    VectorX<T> Ftemp(nsimp);        //简化后的载荷向量
    VectorX<T> Utemp(nsimp);        //简化后的位移向量

    Ktemp = simplyfy(other, KK);
    Ftemp = simplyfy(other, F);

    //求解
    Eigen::SparseLU<SparseMatrix<T>> solver;
    solver.compute(Ktemp);
    Utemp = solver.solve(Ftemp);
    Ftemp = Ktemp * Utemp;

    assembly(Utemp, other, U); //装配位移向量
    assembly(Ftemp, other, F); //装配载荷向量
}

} // namespace ababa
#endif //TINYIGA_SOLVE_H