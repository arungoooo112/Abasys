
#ifndef TINYNURBS_REFINE_H
#define TINYNURBS_REFINE_H

#include <vector>
#include "basis.h"
#include "util/array2.h"
#include "util/coord.h"
#include "util/util.h"

using std::vector;

namespace tinynurbs
{

template <typename T, template <typename> class vecType>
void CurveKnotIns(int p, const vector<T> &UP, const vector<vecType<T>> &Pw,
    T u, int k, int s, int r, vector<T> &UQ, vector<vecType<T>> &Qw)
{
    int mp = UP.size() - 1;
    int np = mp - p - 1;
    int nq = np + r;
    T Rw[p + 1];
    //update knotvector
    for (int i = 0; i <= k; ++i) UQ[i] = UP[i];
    for (int i = 1; i <= r; ++i) UQ[k + i] = u;
    for (int i = k + 1; i <= mp; i++) UQ[i + r] = UP[i];
    // save unchanged control points
    for (int i = 0; i <= k - p; ++i) Qw[i] = Pw[i];
    for (int i = k - s; i <= np; ++i) Qw[i + r] = Pw[i];
    for (int i = 0; i <= p - s; ++i) Rw[i] = Pw[k - p + i];
    // insert r-th
    int L = 0;
    for (int j = 1; j <= r; ++j) {
        L = k - p + j;
        for (int i = 0; i <= p -j -s; ++i) {
            T alpha = (u - UP[L + i]) / (UP[i + k + 1] - UP[L + i]);
            Rw[i] = alpha * Rw[i + 1] + (1.0 - alpha) * Rw[i];
        }
        Qw[L] = Rw[0];
        Qw[k + r - j - s] = Rw[p -j -s];
    }
    for (int i = L + 1; i< k -s; ++i) Qw[i] = Rw[i - L];
}


template <typename T, template <typename> class vecType>
void SurfaceKnotIns(int p, const vector<T> &UP, int q, const vector<T> &VP,
    const array2<vecType<T>> &Pw, bool dir, T uv,  int k, int s, int r, 
    vector<T> &UQ, vector<T> &VQ, array2<vecType<T>> &Qw) 
{
    int mu = UP.size() - 1;
    int nu = mu - p - 1;
    int mv = VP.size() - 1;
    int nv = mv - q - 1;
    T Rw[p + 1];
    if (0 == dir) {
        //update knotvector
        for (int i = 0; i <= k; ++i) UQ[i] = UP[i];
        for (int i = 1; i <= r; ++i) UQ[k + i] = uv;
        for (int i = k + 1; i <= mu; i++) UQ[i + r] = UP[i];
        VQ = VP;
        int L = 0;
        T alpha[][r];
        for (int j = 1; j <= r; ++j) {
            L = k - p +j;
            for (int i = 0; i <= p -j -s; ++i) 
                alpha[i][j] = (uv- UP[L+i]) / (UP[i + k + 1] - UP[L + i]);
        }

        for (int row = 0; row <= nv; ++row) {
            for (int i = 0; i <= k -p; ++i) Qw(i, row) = Pw(i, row);
            for (int i = k -s; i <= nu; ++i) Qw( i + r, row) = Pw(i, row);

            for (int i = 0; i <= p -s; ++i) Rw[i] = Pw(k -p + i, row);

            for (int j = 1; j <= r; ++j) {
                L = k - p +j;
                for (int i = 0; i <= p - j - s; ++i) 
                    Rw[i] = alpha[i][j] * Rw[i + 1] + (1.0 - alpha[i][j]) * Rw[i];
                Qw(L, row) = Rw[0];
                Qw[k + r - j - s, row] = Rw[p - j -s];
            }
            for (int i = L + 1; i <= k -s; ++i) Qw(i ,row) = Rw[i - L];

        }

    }
    if (1 == dir) {
        //update knotvector
        for (int i = 0; i <= k; ++i) VQ[i] = VP[i];
        for (int i = 1; i <= r; ++i) VQ[k + i] = uv;
        for (int i = k + 1; i <= mv; i++) VQ[i + r] = VP[i];
        UQ = UP;
        int L = 0;
        T alpha[][r];
        for (int j = 1; j <= r; ++j) {
            L = k - q +j;
            for (int i = 0; i <= q -j -s; ++i) 
                alpha[i][j] = (uv- VP[L+i]) / (VP[i + k + 1] - VP[L + i]);
        }

        for (int row = 0; row <= nu; ++row) {
            for (int i = 0; i <= k -q; ++i) Qw(i, row) = Pw(i, row);
            for (int i = k -s; i <= nv; ++i) Qw(i + r, row) = Pw(i, row);

            for (int i = 0; i <= q -s; ++i) Rw[i] = Pw(k -q + i, row);

            for (int j = 1; j <= r; ++j) {
                L = k - q +j;
                for (int i = 0; i <= q - j - s; ++i) 
                    Rw[i] = alpha[i][j] * Rw[i + 1] + (1.0 - alpha[i][j]) * Rw[i];
                Qw(L, row) = Rw[0];
                Qw[k + r - j - s, row] = Rw[q - j -s];
            }
            for (int i = L + 1; i <= k -s; ++i) Qw(i ,row) = Rw[i - L];
        }
    }
    
}

// 通过向原曲线的节点向量插入一组新的节点集合来细化曲线，
// 输入：
// 曲线参数：p，U(m+1)，Pw(n+1)
// 待插入的一组节点集合：X(r+1)
// 输出：
// 新的节点向量和新的控制点坐标：Ubar(m+r+2)，Qw(n+r+2)

template <typename T, template<typename> class VecType>
void RefineKnotVectCurve(int p, const vector<T> &U, 
    const vector<VecType<T>> &Pw, const vector<T> &X,
    vector<T> &Ubar, vector<VecType<T>> &Qw)
{
    if (X.size() == 0)
        return;
    int n = Pw.size() - 1; // 基函数，或者控制点的最大下标
    int r = X.size() - 1; //待插入的节点集合的最大下标
    int m = n + p + 1;

    Ubar.resize(m+r+2);
    Qw.resize(n+r+2);

    int a = findSpan(p, U, X[0]);
    int b = findSpan(p, U, X[r]);
    b = b + 1;
    for (int j = 0; j <= a - p; ++j) Qw[j] = Pw[j];
    for (int j = b - 1; j <= n; ++j) Qw[j + r + 1] = Pw[j];
    for (int j = 0; j <= a; ++j) Ubar[j] = U[j];
    for (int j = b + p; j <= m; ++j) Ubar[j + r + 1] = U[j];

    int i = b + p - 1;
    int k = b + p + r;
    for (int j = r; j >= 0; j--)
    {
        while(X[j] <= U[i] && i > a)
        {
            Qw[k - p -1] = Pw[i - p - 1];
            Ubar[k] = U[i];
            k--;
            i--;
        }
        Qw[k - p - 1] = Qw[k - p];
        for (int l = 1; l <= p; ++l) 
        {
            int ind = k - p + l;
            T alfa = Ubar[k + l] - X[j];
            if (abs(alfa) <= 1e-6)
                Qw[ind - 1] = Qw[ind];
            else {
                alfa /= (Ubar[k + l] - U[i - p + l]);
                Qw[ind - 1] = alfa * Qw[ind - 1] + (1.0 - alfa) * Qw[ind];        
            }
        }
        Ubar[k] = X[j];
        k--;
    }
}




// 通过插入一组新的节点，沿u方向来细化曲面

template <typename T, template<typename> class VecType>
void RefineKnotVectU(int p, const vector<T> &U,
    const array2<VecType<T>> &Pw, const vector<T> &X,
    vector<T> &Ubar, array2<VecType<T>> &Qw)
{   
    if(X.size() == 0) return;

    int n = Pw.rows()-1; 
    int m = n+p+1;
    int r = X.size()-1;

    Ubar.resize(m+r+2);
    Qw.resize(n + r + 2, Pw.cols());

    int a = findSpan(p, U, X[0]);
    int b = findSpan(p, U, X[r]);
    b = b + 1;
    for (int j = 0; j <= a; ++j) 
        Ubar[j] = U[j];
    for (int j = b + p; j <= m; ++j) 
        Ubar[j + r + 1] = U[j];
    for (int col = 0; col < Pw.cols(); ++col) 
    {
        for (int k = 0; k <= a - p; ++k) 
            Qw(k, col) = Pw(k, col);
        for (int k = b - 1; k <= n; ++k) 
            Qw(k + r + 1, col) = Pw(k, col);
    }

    int i = b+p-1 ; 
    int k = b+p+r ;
    
    for (int j = r; j >= 0; --j) 
    {
        while(X[j] <= U[i] && i > a) 
        {
            for (int col = 0; col < Pw.cols(); ++col)  
                Qw(k - p -1, col) = Pw(i - p - 1, col);  
            Ubar[k] = U[i];
            k = k - 1;
            i = i - 1;    
        }
        for (int col = 0; col < Pw.cols(); ++col)   
            Qw(k - p - 1, col) = Qw(k - p, col);
        for (int l = 1; l <= p; ++l)
        {
            int ind = k - p + l;
            T alfa = Ubar[k + l] - X[j];
            if (abs(alfa) <= 1e-6)
                for (int col = 0; col < Pw.cols(); ++col)  
                    Qw(ind - 1, col) = Qw(ind, col);
            else {
                alfa /= (Ubar[k + l] - U[i - p + l]);
                for (int col = 0; col < Pw.cols(); ++col)
                    Qw(ind - 1, col) = alfa * Qw(ind - 1, col) + (1.0 - alfa) * Qw(ind, col);    
            }
        }
        Ubar[k] = X[j];
        k = k - 1;
    }
}


// 通过插入一组新的节点，沿v方向来细化曲面

template <typename T, template<typename> class VecType>
void RefineKnotVectV(int q, const vector<T> &V,
    const array2<VecType<T>> &Pw, const vector<T> &X,
    vector<T> &Vbar, array2<VecType<T>> &Qw)
{   
    if(X.size() == 0) return;

    int n = Pw.cols() - 1; 
    int m = n+q+1;
    int r = X.size()-1;

    Vbar.resize(m + r + 2);
    Qw.resize(Pw.rows(), n + r + 2);

    int a = findSpan(q, V, X[0]);
    int b = findSpan(q, V, X[r]);
    b = b + 1;
    for (int j = 0; j <= a; ++j) 
        Vbar[j] = V[j];
    for (int j = b + q; j <= m; ++j) 
        Vbar[j + r + 1] = V[j];
    for (int row = 0; row < Pw.rows(); ++row) {
        for (int k = 0; k <= a - q; ++k) 
            Qw(row, k) = Pw(row, k);
        for (int k = b - 1; k <= n; ++k) 
            Qw(row, k + r + 1) = Pw(row, k);
    }

    int i = b+q-1 ; 
    int k = b+q+r ;
    
    for (int j = r; j >= 0; --j) 
    {
        while(X[j] <= V[i] && i > a) 
        {
            Vbar[k] = V[i];
            for (int row = 0; row < Pw.rows(); ++row)  
                Qw(row, k - q -1) = Pw(row, i - q - 1);  
            k = k - 1;
            i = i - 1;
            
        }
        for (int row = 0; row < Pw.rows(); ++row)   
            Qw(row, k - q - 1) = Qw(row, k - q);
        for (int l = 1; l <= q; ++l)
        {
            int ind = k - q + l;
            T alfa = Vbar[k + l] - X[j];
            if (abs(alfa) <= 1e-6)
                for (int row = 0; row < Pw.rows(); ++row)  
                    Qw(row, ind - 1) = Qw(row, ind);
            else {                
                alfa = alfa / (Vbar[k + l] - V[i - q + l]);
                for (int row = 0; row < Pw.rows(); ++row)
                    Qw(row, ind - 1) = alfa * Qw(row, ind - 1) + (1.0 - alfa) * Qw(row, ind);
            }
        }
        Vbar[k] = X[j];
        k = k - 1;
    }
}

template <typename T, template<typename> class VecType>
void RefineKnotVectSurface(int p, int q, const vector<T> &U, const vector<T> &V,
    const array2<VecType<T>> &Pw, const vector<T> &X, const vector<T> &Y,
                           vector<T> &Ubar, vector<T> &Vbar, array2<VecType<T>> &Qw)
{
    array2<VecType<T>> Mw;
    RefineKnotVectU(p, U, Pw, X, Ubar, Mw);
    RefineKnotVectV(q, V, Mw, Y, Vbar, Qw);
}

template<typename T>
void insertKnotVect(vector<T>& U, const vector<T>& X) {
    U.insert(U.end(), X.begin(), X.end());
    sort(U.begin(), U.end());
}

template<typename T>
vector<T> getinsertKnotVect(vector<T>& U, int n) {
    vector<T> tmp;
    for (int i = 0; i < U.size() - 1; i++) {
        if (!close(U[i], U[i+1], 1e-6)) {
            T h = (U[i+1] - U[i]) / (n+1);
            for (int j = 1; j <= n; j++) {
                tmp.push_back(U[i] + h * j); 
            }
        }
    }
    return tmp;
}


} //namespace tinynurbs

#endif