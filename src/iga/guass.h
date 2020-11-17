#ifndef TINYIGA_GUASS_H
#define TINYIGA_GUASS_H

#include "util/array2.h"
namespace abab {
//获取二维高斯积分点
//@param[in) 高斯点个数
//@return 数组，包含高斯点及权值
template<typename T>
array2<T> Guass2D(int n) {
    array2<T> gauss(n, 3);
    if (n == 9) {
        T gs = (T)0.774596669241483;
        gauss(0, 0) = -gs; gauss(0, 1) = -gs; gauss(0, 2) = (T)0.308641975308642;
        gauss(1, 0) = -gs; gauss(1, 1) =   0; gauss(1, 2) = (T)0.493827160493828;
        gauss(2, 0) = -gs; gauss(2, 1) =  gs; gauss(2, 2) = (T)0.308641975308642;
        gauss(3, 0) =   0; gauss(3, 1) =  gs; gauss(3, 2) = (T)0.493827160493828;
        gauss(4, 0) =  gs; gauss(4, 1) =  gs; gauss(4, 2) = (T)0.308641975308642;
        gauss(5, 0) =  gs; gauss(5, 1) =   0; gauss(5, 2) = (T)0.493827160493828;
        gauss(6, 0) =  gs; gauss(6, 1) = -gs; gauss(6, 2) = (T)0.308641975308642;
        gauss(7, 0) =   0; gauss(7, 1) = -gs; gauss(7, 2) = (T)0.493827160493828;
        gauss(8, 0) =   0; gauss(8, 1) =   0; gauss(8, 2) = (T)0.790123456790124;
    }
    return gauss;
}

//获取一维高斯积分点
//@param[in) 高斯点个数
//@return 数组，包含高斯点及权值
template<typename T>
array2<T> Guass1D(int n) {
    array2<T> gauss(n, 2);
    if (n == 4) {
        gauss(0, 0) = (T) 0.861136311594053; gauss(0, 1) = (T)0.347854845137454;
        gauss(1, 0) = (T) 0.339981043584856; gauss(1, 1) = (T)0.652145154862546;
        gauss(2, 0) = (T)-0.339981043584856; gauss(2, 1) = (T)0.652145154862546;
        gauss(3, 0) = (T)-0.861136311594053; gauss(3, 1) = (T)0.347854845137454;
    }
    return gauss;
}

} //namespace abab
#endif //TINYIGA_GUASS_H