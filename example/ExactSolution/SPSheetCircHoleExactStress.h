#ifndef SPSHEETCIRCHOLEEXACTSTRESS_H
#define SPSHEETCIRCHOLEEXACTSTRESS_H

#include <cmath>

template <typename T>
struct Stressxx 
{
    Stressxx(T a, T tx) : a(a), Tx(tx) {}

    T operator()(T x, T y) {
        T r = sqrt(x * x + y * y);
        T theta = atan(y / x);
        T c2t = cos(2 * theta);
        T c4t = cos(4 * theta);
        T term1 = pow(a / r, 2);
        T term2 = pow(a / r, 4);
        return Tx * (1 - term1 * (3.0 / 2.0 * c2t + c4t) + 3.0 / 2.0 * term2 * c4t);
    }
    T a;
    T Tx;
};

template <typename T>
struct Stressyy 
{
    Stressyy(T a, T tx) : a(a), Tx(tx) {}

    T operator()(T x, T y) {
        T r = sqrt(x * x + y * y);
        T theta = atan(y / x);
        T c2t = cos(2 * theta);
        T c4t = cos(4 * theta);
        T term1 = pow(a / r, 2);
        T term2 = pow(a / r, 4);
        return Tx * (- term1 * (1.0 / 2.0 * c2t - c4t) - 3.0 / 2.0 * term2 * c4t);
    }
    T a;
    T Tx;
};


template <typename T>
struct Stressxy 
{
    Stressxy(T a, T tx) : a(a), Tx(tx) {}

    T operator()(T x, T y) {
        T r = sqrt(x * x + y * y);
        T theta = atan(y / x);
        T s2t = sin(2 * theta);
        T s4t = sin(4 * theta);
        T term1 = pow(a / r, 2);
        T term2 = pow(a / r, 4);
        return Tx * (- term1 * (1.0 / 2.0 * s2t + s4t) - 3.0 / 2.0 * term2 * s4t);
    }
    T a;
    T Tx;
};

template <typename T>
struct Spsheetcircholeexactstress
{   
    Spsheetcircholeexactstress(T a, T tx) : xx(a, tx), yy(a, tx), xy(a, tx) {}
    Stressxx<T> xx;
    Stressyy<T> yy;
    Stressxy<T> xy;
};


#endif //SPSHEETCIRCHOLEEXACTSTRESS_H