#ifndef TINYIGA_INDEX_H
#define TINYIGA_INDEX_H

#include <vector>

#include "./array2.h" //array2
#include "./coord.h"

using std::vector;
namespace abab {

//控制点二维索引转一维索引，按行线性化
//@param[in] 控制点二维索引
//@param[in] 每行元素个数，即总列数
//@return 一维索引
vector<int> Index2To1ByV(const array2<int>& idxs, int cols) {
    std::vector<int> res;
    for (int i = 0; i < idxs.rows(); i++) 
    {
        res.push_back(idxs(i, 0) * cols + idxs(i, 1));
    }
    return res;
}

//控制点二维索引转一维索引，按列线性化
//@param[in] 控制点二维索引
//@param[in] 每行元素个数，即总列数
//@return 一维索引
std::vector<int> Index2To1ByU(const array2<int>& idxs, int rows) {
    std::vector<int> res;
    for (int i = 0; i < idxs.rows(); i++)
    {
        res.push_back(idxs(i, 0) + idxs(i, 1) * rows);
    }
    return res;
}

//由一维索引生成对应的装配索引
//@param[in] 控制点一维索引
//@return 控制点装配索引
std::vector<int> getDofsIndex(const std::vector<int>& Index1) {
    std::vector<int> res;
    for (const auto& a : Index1) {
        res.push_back(2 * a);
        res.push_back(2 * a + 1);
    }
    return res;
}

void Index1ToDofX(std::vector<int>& Index1) {
    for (auto& a : Index1) 
        a = 2 * a;
}

void Index1ToDofY(std::vector<int>& Index1) {
    for (auto& a : Index1) 
        a = 2 * a + 1;
}


//获取局部的控制点二维索引
//@param[in] u，v方向的次数
//@return 二维数组（nep，2）nep：numbers of element's points 

array2<int> LocalIndex2(int deg_u, int deg_v)
{
    int nep = (deg_u + 1) * (deg_v + 1); //numbers of element's points 
    array2<int> res(nep, 2);
    int k = 0;
    for (int j = 0; j <= deg_v; ++j) {
        for (int i = 0; i <= deg_u; ++i)
        { //沿着u向顺序添加控制点的索引
            res(k, 0) = i;
            res(k, 1) = j;
            k++;
        }
    }
    return res;
}

//获取全局的控制点二维索引
//@param[in] 单元编号（ei，ej）
//@param[in] u，v方向的次数
//@return 二维数组（nep，2）nep：numbers of element's points 
array2<int> GlobalIndex2(int ei, int ej, int deg_u, int deg_v)
{
    array2<int> res = LocalIndex2(deg_u, deg_v);
    int iBegin = ei - deg_u;
    int jBegin = ej - deg_v;
    int nep = (deg_u + 1) * (deg_v + 1);
    for (int n = 0; n < nep; ++n) {
        res(n, 0) += iBegin;
        res(n, 1) += jBegin;
    }
    return res;
}

//获取局部的控制点二维索引
//@param[in] u，v方向的次数
//@return 二维数组（nep，2）nep：numbers of element's points 

vector<int> LocalIndex1(int degree) {
    vector<int> res;
    for (int i = 0; i <= degree; ++i) {
        res.push_back(i);     
    }
    return res;
}

//获取全局的控制点二维索引
//@param[in] 单元编号（ei，ej）
//@param[in] u，v方向的次数
//@return 二维数组（nep，2）nep：numbers of element's points 
vector<int> GlobalIndex1(int e, int degree) {
    vector<int> res = LocalIndex1(degree);
    int begin = e - degree;  
    for (int n = 0; n <= degree; ++n) {
        res[n] += begin;   
    }
    return res;
}
} //namespace abab
#endif //TINYIGA_INDEX_H