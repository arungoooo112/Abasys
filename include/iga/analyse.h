#ifndef IGA_ANALYSE_H
#define IGA_ANALYSE_H

#include <Eigen/Dense>

#include "nurbs/surface.h"
#include "elematrix.h"
#include "assembly.h"

using namespace Eigen;
using namespace std;
using namespace tinynurbs;

namespace abab {

template<typename T>
Vector2<T> getDispltVals(const tinynurbs::RationalSurface<T>& surf, T u, T v,
    const VectorX<T>& U) 
{   
    int ei = tinynurbs::findSpan(surf.degree_u, surf.knots_u, u);
    int ej = tinynurbs::findSpan(surf.degree_v, surf.knots_v, v);
    
    array2<int> idx2s = GlobalIndex2(ei, ej, surf.degree_u, surf.degree_v);//单元（ei，ej）控制点的二维索引
    std::vector<int> idx1s = Index2To1ByV(idx2s, surf.control_points.cols());//控制点按行编号的一维索引
    std::vector<int> assIdxs = getDofsIndex(idx1s);//装配索引

    VectorX<T> Ue = simplyfy(assIdxs, U);
    return getShapeMatrix(surf, u, v) * Ue;
}

template<typename T>
Vector3<T> getStrainVals(const tinynurbs::RationalSurface<T>& surf, T u, T v, 
    const VectorX<T>& U) 
{   
    int ei = tinynurbs::findSpan(surf.degree_u, surf.knots_u, u);
    int ej = tinynurbs::findSpan(surf.degree_v, surf.knots_v, v);
    
    array2<int> idx2s = GlobalIndex2(ei, ej, surf.degree_u, surf.degree_v);//单元（ei，ej）控制点的二维索引
    std::vector<int> idx1s = Index2To1ByV(idx2s, surf.control_points.cols());//控制点按行编号的一维索引
    std::vector<int> assIdxs = getDofsIndex(idx1s);//装配索引

    VectorX<T> Ue = simplyfy(assIdxs, U);

    return getStrainMatrix(surf, u, v) * Ue;
}


} //namespace abab



#endif //IGA_ANALYSE_H