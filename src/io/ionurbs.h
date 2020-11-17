#ifndef TINYNURBS_IONURBS_H
#define TINYNURBS_IONURBS_H

#include <iostream>
#include <vector>
#include "../util/array2.h"
#include "../util/coord.h"
#include "../core/curve.h"
#include "../core/surface.h"


namespace tinynurbs {

template <typename T>
void outputCurve(std::ostream& os, int degree, const std::vector<T>& knots, 
                const std::vector<vec3<T>>& control_points, 
                const std::vector<T>& weights) 
{
    using std::endl;

    os << "degrees\t" << degree << endl;
    os << endl;

    os << "knots\t" << knots.size() << endl;
    for (const auto& u : knots)
        os << "    " << u << endl;
    os << endl; 

    os << "control_points\t" << control_points.size() << endl;
    for (int i = 0; i < control_points.size(); i++) {
         os << "    " << control_points[i].x << " \t" 
            << "    " << control_points[i].y << " \t" 
            << "    " << control_points[i].z << " \t" << endl;
    }
    os << endl; 

    os << "weights\t" << weights.size() << endl;
    for (int i = 0; i < weights.size(); i++) {
         os << "    " << weights[i] << endl;
    }
}

template <typename T> void outputCurve(std::ostream &os, const RationalCurve<T> &crv)
{   
    outputCurve(os, crv.degree, crv.knots, crv.control_points, crv.weights);  
}

template <typename T> std::ostream &operator<<(std::ostream & os, const RationalCurve<T> &crv) 
{
    outputCurve(os, crv);
    return os;
}

template <typename T>
void outputSurface(std::ostream &os, int deg_u, int deg_v,
                    const std::vector<T> &knots_u, const std::vector<T> &knots_v,
                    const array2<vec3<T>> &points, const array2<T>& weights) 
{
    using std::endl;

    os << "degrees\t" << deg_u << "\t" << deg_v << endl;
    os << endl;
    
    os << "knots_u\t" << deg_u << "\t" << knots_u.size() << endl;
    for (const auto& u : knots_u)
        os << "    " << u << endl;
    os << endl;

    os << "knots_v\t" << deg_v << "\t" << knots_v.size() << endl;
    for (const auto& v : knots_v)
        os << "    " << v << endl;
    os << endl;

    os << "points\t" << points.rows() << "\t" << points.cols() << endl;
    for (int j = 0; j < points.cols(); j++)
    {
        for (int i = 0; i < points.rows(); i++) {
             os << "    " 
                << points(i, j).x << " \t" 
                << points(i, j).y << " \t" 
                << points(i, j).z << endl;
        }
    }

    os << "weights\t" << weights.rows() << "\t" << weights.cols() << endl;
    for (int j = 0; j < weights.cols(); j++)
    {
        for (int i = 0; i < weights.rows(); i++) {
             os << "    " << weights(i, j) << endl;
        }
    }
    os << endl;
}

template <typename T>
void outputSurface(const RationalSurface<T>& srf) {
    outputSurface(std::cout, srf.degree_u, srf.degree_v, srf.knots_u, srf.knots_v, srf.control_points, srf.weights);
}

template <typename T> std::ostream &operator<<(std::ostream & os, const RationalSurface<T> &srf) 
{
    outputSurface(srf);
    return os;
}



}

#endif //TINYNURBS_IONURBS_H