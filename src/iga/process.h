#ifndef TINYIGA_PROCESS_H
#define TINYIGA_PROCESS_H

#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "core/surface.h"
#include "core/boundary.h"
#include "imposebc.h"
#include "assembly.h"
#include "elematrix.h"
#include "solve.h"

using namespace std;
using namespace Eigen;
using namespace tinynurbs;
namespace abab {
enum Label {X, Y};

template<typename T>
class Process {
public: 
    Process(const RationalSurface<T>& s);

    void calStiffnessMatrix(T E, T nu, ID id);

    template<typename Function>
    void applyNewmannBC(Boundary b, Function fun, Label lb);
    template<typename Function>
    void applyDrchltBC(Boundary b, Function fun, Label lb);

    void solve();

private:
    VectorX<T> F;
    VectorX<T> U;
    MatrixX<T> KK;
    
    vector<int> knowns;
    vector<pair<VectorX<T>, vector<int>>> fContainors;
    vector<pair<VectorX<T>, vector<int>>> uContainors;
    RationalSurface<T> surf;
};

template<typename T>
Process<T>::Process(const RationalSurface<T>& s) : surf(s) {
        int np = surf.control_points.size();
        F.setZero(2 * np);
        U.setZero(2 * np);
        KK.resize(2 * np, 2 * np);
}

template<typename T>
void Process<T>::calStiffnessMatrix(T E, T nu, ID id)
{   
    KK = assemblySparse(surf, E, nu, id);
}

template<typename T>
template<typename Function>
void Process<T>::applyNewmannBC(Boundary b, Function fun, Label lb) {
    VectorX<T>  vals = applyNewmannBdryVals(surf, b, fun);
    vector<int> idxs = getSurfaceBdryIdxs(surf, b);
    if (lb == X) Index1ToDofX(idxs);
    if (lb == Y) Index1ToDofY(idxs);
    fContainors.emplace_back(vals, idxs);
}

template<typename T>
template<typename Function>
void Process<T>::applyDrchltBC(Boundary b, Function fun, Label lb) {
    auto vals = applyDrchltBdryVals(surf, b, fun);
    auto idxs = getSurfaceBdryIdxs(surf, b);
    if (lb == X) Index1ToDofX(idxs);
    if (lb == Y) Index1ToDofY(idxs);
    uContainors.emplace_back(vals, idxs);
    knowns.insert(knowns.end(), idxs.begin(), idxs.end());
}

template<typename T>
void Process<T>::solve() 
{
    sort(knowns.begin(), knowns.end());
    knowns.erase(unique(knowns.begin(), knowns.end()), knowns.end());
    
    for (const auto& f : fContainors) {
        assembly(f.first, f.second, F);
    }

    for (const auto& u : uContainors) {
        assembly(u.first, u.second, U);
    }

    cout << "ready to solve\n";
    //cout << "KK: \n" << KK.rows() << " " << KK.cols() << endl << KK << endl << endl;
    //cout << "U: \n" << U.rows() << " " << U.cols() << endl << U << endl << endl;
    //cout << "F: \n" << F.rows() << " " << F.cols() << endl << F << endl << endl;
    //std::for_each(knowns.begin(), knowns.end(), [](int x) { cout << x << "\n"; });

    cout << endl << endl;
    abab::solve(KK, U, F, knowns);

    cout << "solve U: \n" << U.rows() << " " << U.cols() << endl << U << endl << endl;
}

} //namespace abab
#endif //TINYIGA_PROCESS_H