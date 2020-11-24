#ifndef TINYNURBS_MODEL_H
#define TINYNURBS_MODEL_H

#include <cmath>
#include <vector>
#include <initializer_list>
#include "nurbs/curve.h"
#include "nurbs/surface.h"
#include "./coord.h"
#include "./array2.h"

static const double PI = 3.14159265358979323846;

namespace tinynurbs{


std::vector<double> createKnotVector(int p, int np) {
    int n = np - 1;
    int m = n + p + 1;
    int ne = np - p;
    double h = 1.0 / ne;
    std::vector<double> knots(p, 0.0);
    for (int i = 0; i < ne + 1; i++) {
        knots.push_back(i * h);
    }
    knots.resize(m + 1, 1.0);
    return knots;
}

//template<typename T, typename... Args>
//RationalSurface<T> createSurface(const RationalCurve<T>& c, const Args&... rest) 
//{
//    array2<T> weights;
//    array2<vec3<T>> points;
//    points.push_back(c.control_points);
//    weights.push_back(c.weights);
//    points.push_back((rest...).control_points);
//    weights.push_back((rest...).weights);
//    std::vector<T> knot_u = createKnotVector(2, weights.rows());
//    std::vector<T> knot_v = createKnotVector(2, weights.cols());
//    return RationalSurface<T>(2, 2, knot_u, knot_v, points, weights);
//}

template <typename T>
RationalSurface<T> createSurface(std::initializer_list<RationalCurve<T>> il) 
{
    array2<T> weights;
    array2<vec3<T>> points;
    for (const RationalCurve<T>& c : il) {
        points.push_back(c.control_points);
        weights.push_back(c.weights);
    }
    std::vector<T> knot_u = createKnotVector(2, weights.rows());
    std::vector<T> knot_v = createKnotVector(2, weights.cols());
    return RationalSurface<T>(2, 2, knot_u, knot_v, points, weights);
}

template<typename T>
RationalSurface<T> createSquare(T a = 1.0, T b = 1.0, T ox = 0.0, T oy = 0.0) {
    std::vector<T> knots_u = {0,0,1,1}; // 定义节点向量
    std::vector<T> knots_v = {0,0,1,1}; // 定义节点向量
    array2<T> weights(2, 2, 1); //定义权值
    array2<vec3<T>> points(2, 2);
    points(0,0) = vec3<T>(0, 0, 0);
    points(1,0) = vec3<T>(a, 0, 0);
    points(0,1) = vec3<T>(0, b, 0);
    points(1,1) = vec3<T>(a, b, 0);

    return RationalSurface<T>(1,1,knots_u, knots_v, points, weights);
}

template<typename T>
RationalCurve<T> createLine(std::initializer_list<vec3<T>> il) 
{
    int np = il.size();
    std::vector<T> knots = createKnotVector(1, np);
    return RationalCurve<T>(1, knots, il, std::vector<T>(np, 1.0));
}

template<typename T>
RationalCurve<T> createCircle(T radius, T sang, T eang, T ox, T oy)
{
    T sweep = eang - sang; 
    if (sweep < 0)
        sweep = 2 * PI + sweep;  //计算圆弧夹角

    int narcs = 0; //计算小圆弧段数
    std::vector<T> knots; // 定义节点向量
    std::vector<T> weights; //定义权值

    T dsweep = 0.0;
    T wm = 0.0;

    if (sweep <= PI / 2) {
        narcs = 1;
        knots = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
        dsweep = sweep / narcs / 2;
        wm = cos(dsweep);
        weights = {1.0, wm, 1.0};
    }
    else if (sweep <= PI)
    {
        narcs = 2;
        knots = {0, 0, 0, 0.5, 0.5, 1, 1, 1};
        dsweep = sweep / narcs / 2;
        wm = cos(dsweep);
        weights = {1.0, wm, 1.0, wm, 1.0};
    }
    else if (sweep <= 3 * PI / 2)
    {
        narcs = 3;
        knots = {0, 0, 0, 0.5, 0.5, 1, 1, 1};
        dsweep = sweep / narcs / 2;
        wm = cos(dsweep);
        weights = {1.0, wm, 1.0, wm, 1.0, wm, 1.0};
    }
    else {
        narcs = 4;
        knots = {0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0};
        dsweep = sweep / narcs / 2;
        wm = cos(dsweep);
        weights = {1.0, wm, 1.0, wm, 1.0, wm, 1.0, wm, 1.0};
    }


    T x  = radius*wm;
    T y  = radius*sin(dsweep);
    T xm = x+y*tan(dsweep);

    std::vector<vec3<T>> points = {vec3<T>(x, -y, 0.0), vec3<T>(xm, 0.0, 0.0), vec3<T>(x, y, 0.0)};
 
    for (auto &a : points) //对前三个点进行旋转,到开始角的位置
    {               
        //a = glm::rotateZ(a, glm::radians(sang + dsweep)); // rotate to start angle
        a = rotateZ(a, sang+dsweep);
        // T alpha = sang + dsweep;
        // T tx = a.x;
        // T ty = a.y;
        // a.x = tx * cos(alpha) - ty * sin(alpha);
        // a.y = tx * sin(alpha) + ty * cos(alpha);
    }

    for (int i = 3; i < 2 * narcs + 1; i++) { //通过旋转映射添加剩余节点
        //points.push_back(glm::rotateZ(points[i-2], glm::radians(2 * dsweep)));     
        //T alpha = dsweep + dsweep;
        //T tx = points[i - 2].x;
        //T ty = points[i - 2].y;
        //T x = tx * cos(alpha) - ty * sin(alpha);
        //T y = tx * sin(alpha) + ty * cos(alpha);
        //points.emplace_back(x, y, (T)0.0);
        points.push_back(rotateZ(points[i-2], 2 * dsweep));
    }

    //对所有坐标点进行偏移
    for (auto& a : points) {
        a += vec3<T>(ox, oy, 0.0);
    }

    return RationalCurve<T>(2, knots, points, weights);
}




} // namespace tinynurbs

#endif //TINYNURBS_MODEL_H