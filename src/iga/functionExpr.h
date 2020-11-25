#ifndef ABASYS_FUNCTIONEXPR_H
#define ABASYS_FUNCTIONEXPR_H

#include <iostream>
#include <string>
#include <vector>

#include "exprtk/exprtk.hpp";

using std::cout;
using std::string;

template <typename T>
class FunctionExpr2
{
    typedef exprtk::symbol_table<T> symbol_table_t;
    typedef exprtk::expression<T>     expression_t;
    typedef exprtk::parser<T>             parser_t;

public:
    FunctionExpr2(string expr);
    FunctionExpr2 operator=(string expr);

    //FunctionExpr2(string v1, string v2, string expr);

    //void addVariety(string var);
    //void addConstant(string con);

public:
    //Array1<T> operator()(Array2<T> inMat);
    T operator()(T u, T v);

private:
    string expression_string;
    symbol_table_t symbol_table;
    expression_t expression;
    parser_t parser;
private:
    
    T _x, _y;
};

template <typename T>
FunctionExpr2<T>::FunctionExpr2(string expr) : expression_string(expr) {
    symbol_table.add_variable("x", _x);
    symbol_table.add_variable("y", _y);
    symbol_table.add_constants();
    expression.register_symbol_table(symbol_table);
    if (!parser.compile(expression_string, expression))
        printf("Compilation error...\n");
}

template <typename T>
FunctionExpr2<T> FunctionExpr2<T>::operator=(string expr){

}
//template <typename T>
//Array1<T> FunctionExpr2<T>::operator()(Array2<T> inMat)
//{
//    typedef exprtk::symbol_table<T> symbol_table_t;
//    typedef exprtk::expression<T>     expression_t;
//    typedef exprtk::parser<T>             parser_t;
//
//    T x = T(0);
//    T y = T(0);
//
//    symbol_table_t symbol_table;
//    symbol_table.add_variable("v1", x);
//    symbol_table.add_variable("v2", y);
//
//    expression_t expression;
//    expression.register_symbol_table(symbol_table);
//
//    parser_t parser;
//    parser.compile(expr_string, expression);
//
//    Array1 outVec(inMat.rows());
//    for (int i = 0; i < inMat.rows(); ++i)
//    {
//        x = inMat(i, 0);
//        y = inMat(i, 1);
//        outVec[i] = expression.value();
//    }
//    return outVec;
//}

template <typename T>
T FunctionExpr2<T>::operator()(T x, T y) {
    _x = x;
    _y = y;
    static const T pi = T(3.141592653589793238462643383279502);
    return expression.value();
}

using FunctionExpr2d = FunctionExpr2<double>;


#endif //ABASYS_FUNCTIONEXPR_H