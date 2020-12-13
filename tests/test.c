#include <iostream>
#include <fstream>
#include <cmath>
#include <initializer_list>
#include <Eigen/Dense>

#include "core/bases/check.h"
#include "core/model.h"
#include "core/check.h"
#include "core/surface.h"
#include "core/basis.h"
#include "core/evaluate.h"

#include "io/ionurbs.h"
#include "iga/stiffness.h"
#include "iga/assembly.h"
#include "iga/functionExpr.h"
#include "iga/process.h"
#include "iga/solve.h"
#include "iga/analyse.h"

#include "E:\tinynurbs\modules\ExactSolution\SPSheetCircHoleExactStress.h"


using namespace std;
using namespace Eigen;
using namespace tinynurbs;
using namespace abab;

bool testCurve() {
	auto circle = createCircle<double>(1, 0.0, 2 * PI / 3, 0, 0);
	//outputCurve(std::cout, circle);
	cout << circle;
	//cout << isValidCurve(circle.degree, circle.knots, circle.control_points);
	return true;
}

bool testSurface() {
	auto square = createSquare<double>();
	//outputCurve(std::cout, circle);
	cout << square << endl;
	cout << isValidSurface(square);
	return true;
}

bool teststiffness() {
	auto square = createSquare<double>();
	//outputCurve(std::cout, circle);
    cout << square << endl << endl;

	cout << ElementStiffnessMatrix(square, 1, 1, 1e6, 0.25, ID::PlaneStress);

	//cout << assembly(square, 1e6, 0.25, ID::PlaneStress);

	return true;
}

template<typename T, typename Fun>
T testFunctionExpr(T x, T y, Fun& g) {
	return g(x, y);
}

bool testFunctionExpr() {
	FunctionExpr2d fun("x + y");
	cout << fun(1, 2);
	cout << testFunctionExpr(1, 1, fun);
	return true;
}

bool testHRefine() {
    auto square = createSquare<double>();
    //outputCurve(std::cout, circle);
    cout << square << endl << endl;

    array2<hvec3d> hpoints = cartesianToHomogenous(square.control_points, square.weights);
    vector<double> Ubar;
    array2<hvec3d> Qw;
    vector<double> X = { 0.5 };
    RefineKnotVectU(1, square.knots_u, hpoints,
        X, Ubar, Qw);

    return 1;
}

bool SheetCircHole1()
{
    //---------CAD model----------------------------------------
    double a = 1; //圆弧半径
    double L = 5; //四边形边长

    vector<vec3d> ctrlpts1 = { vec3d(0, L, 0), vec3d(-L, L, 0), vec3d(-L, L, 0),  vec3d(-L, 0, 0)};
    vector<double> weights1 = { 1.0,1.0,1.0,1.0 };

    RCurve3d curv2 = createCircle((a + L) / 2, PI / 2, PI, 0.0, 0.0);
    curv2.HRefine({ 0.5 });
    vector<vec3d> ctrlpts2 = curv2.control_points;
    vector<double> weights2 = curv2.weights;

    RCurve3d curv3 = createCircle(a, PI / 2, PI, 0.0, 0.0);
    curv3.HRefine({ 0.5 });
    vector<vec3d> ctrlpts3 = curv3.control_points;
    vector<double> weights3 = curv3.weights;

    array2<vec3d> ctrlpts = { ctrlpts1, ctrlpts2, ctrlpts3 };
    array2<double> weights = { weights1, weights2, weights3 };

    vector<double> kntvect1 = { 0, 0, 0, 1, 1, 1 };
    vector<double> kntvect2 = { 0, 0, 0, 0.5, 1, 1, 1 };

    RSurface3d Surf(2, 2, kntvect1, kntvect2, ctrlpts, weights);

    cout << Surf << endl << endl << "-----------------------------" << endl;

    cout << isValidSurface(Surf) << "  valid?" << endl;

    Surf.HRefine(2, 2);
    cout << Surf << endl << endl << "-----------------------------" << endl;

    //----------------------------------------------------------
    //----------------------------------------------------------
    //----------------process-----------------------------------
    double E = 1e5;
    double nu = 0.3;
    //auto ke = ElementStiffnessMatrix(Surf, 2, 2, E, nu, PlaneStress);//求（ei，ej）单元的刚度矩阵
    //
    //cout << ke;
    Process<double> process(Surf);

    //----------Assembling the system--------------------------


    //MatrixX<double> K = assembly(Surf, E, nu, PlaneStress);
    process.calStiffnessMatrix(E, nu, PlaneStress);

    //----------boundary condition-----------------------------

    Spsheetcircholeexactstress<double> S(a, 1.0);
    
    auto Sxx = [&](double x, double y) {return -S.xx(x, y) * (abs(x + L) <= 1e-6); };
    auto Sxy = [&](double x, double y) {return -S.xy(x, y) * (abs(x + L) <= 1e-6); };
    auto Syy = [&](double x, double y) {return  S.yy(x, y) * (abs(y - L) <= 1e-6); };
    auto Syx = [&](double x, double y) {return  S.xy(x, y) * (abs(y - L) <= 1e-6); };

    process.applyNewmannBC(U0, Sxx, X);
    process.applyNewmannBC(U0, Sxy, Y);
    process.applyNewmannBC(U0, Syx, X);
    process.applyNewmannBC(U0, Syy, Y);

    auto h = [](double x, double y) {return 0.0; };
    process.applyDrchltBC(V1, h, Y);
    process.applyDrchltBC(V0, h, X);
    //----------solve-------------------------------------------
    process.solve();
    //-------------post processing------------------------------

    return true;
}


bool SheetCircHole()
{
    cout << "构造曲线模型 Surf: " << endl;
    //---------CAD model----------------------------------------
    double a = 1; //圆弧半径
    double L = 5; //四边形边长

    RCurve3d curv1 = createLine({ vec3d(0, L, 0), vec3d(-L, L, 0), vec3d(-L, L, 0),  vec3d(-L, 0, 0) });

    RCurve3d curv2 = createCircle((a + L) / 2, PI / 2, PI, 0.0, 0.0);
    curv2.HRefine({ 0.5 });

    RCurve3d curv3 = createCircle(a, PI / 2, PI, 0.0, 0.0);
    curv3.HRefine({ 0.5 });

    auto Surf = createSurface({ curv1, curv2, curv3 });

    cout << Surf << endl << endl;

    cout << "Is Surf valid?  " << isValidSurface(Surf) << endl;

    cout << "对Surf进行加密，（20，20）\n";
    Surf.HRefine(50, 50);
    cout << "加密后的Surf：\n";
    cout << Surf << endl << endl;

    double E = 1e5;
    double nu = 0.3;

    cout << "calculate K:\n";
    MatrixX<double> K = assemblyMatrix(Surf, E, nu, PlaneStress);
    
    cout << K;
    Spsheetcircholeexactstress<double> S(a, 1.0);

    auto Sxx = [&](double x, double y) {return -S.xx(x, y) * (abs(x + L) <= 1e-6); };
    auto Sxy = [&](double x, double y) {return -S.xy(x, y) * (abs(x + L) <= 1e-6); };
    auto Syy = [&](double x, double y) {return  S.yy(x, y) * (abs(y - L) <= 1e-6); };
    auto Syx = [&](double x, double y) {return  S.xy(x, y) * (abs(y - L) <= 1e-6); };

    vector<pair<VectorX<double>, vector<int>>> fContainor;
    fContainor.push_back(applyNewmannBC(Surf, U0, Sxx, X));
    fContainor.push_back(applyNewmannBC(Surf, U0, Sxy, Y));
    fContainor.push_back(applyNewmannBC(Surf, U0, Syx, X));
    fContainor.push_back(applyNewmannBC(Surf, U0, Syy, Y));

    auto h = [](double x, double y) {return 0.0; };

    vector<pair<VectorX<double>, vector<int>>> uContainor;
    uContainor.push_back(applyDrchltBC(Surf, V1, h, Y));
    uContainor.push_back(applyDrchltBC(Surf, V0, h, X));

    VectorX<double> F = Eigen::VectorXd::Zero(K.rows());
    VectorX<double> U = Eigen::VectorXd::Zero(K.rows());

    assembly(fContainor, F);
    assembly(uContainor, U);

    vector<int> knowns;
    for (const auto& ui : uContainor)
        knowns.insert(knowns.end(), ui.second.begin(), ui.second.end());
    //abab::solveByLarge(K, U, F, knowns);
    abab::solve(K, U, F, knowns);

    cout << "\npost process:\n"<< endl << endl;
    fstream out("post.csv");
    double u = 1.0, v = 0.0;
    int n = 180;
    for (int i = 0; i <= n; i++)
    {
        v = 1.0 / n * i;
        auto crd = surfacePoint(Surf, u, v);
        double x = crd[0], y = crd[1];

        double t = atan(-x / y) / PI * 180 + 90.0;

        Vector3d strain = getStrainVals(Surf, u, v, U);
        Matrix3d D = getElastMatrix(E, nu, PlaneStress);
        Vector3d stress = D * strain;
        out << t << ", " << S.xx(x, y) << ", " << -stress[0] << endl;
        printf("%3.2f : %f  %f\n", t, S.xx(x, y), -stress[0]);
    }

    //double u = 1.0, v = 0.0;
    //auto crd = surfacePoint(Surf, u, v);
    //double x = crd[0], y = crd[1];
    //cout << "\ncalcu stress on (u, v) = (1.0, 0.0):\n";
    //cout << "\ncalcu stress on (x, y) = (" << x << ", " << y << "):\n";

    //cout << S.xx(x, y) << endl;
    //cout << S.yy(x, y) << endl;
    //cout << S.xy(x, y) << endl;

    //cout << "\ncalcu stress:\n";
    //Vector3d strain = getStrainVals(Surf, u, v, U);
    //Matrix3d D = getElastMatrix(E, nu, PlaneStress);
    //Vector3d stress = D * strain;
    //cout << stress;

    return true;
}

bool SheetCircHole2()
{
    //---------CAD model----------------------------------------
    double a = 1; //圆弧半径
    double L = 5; //四边形边长

    RCurve3d curv1 = createLine({ vec3d(0, L, 0), vec3d(-L, L, 0), vec3d(-L, L, 0),  vec3d(-L, 0, 0) });

    RCurve3d curv2 = createCircle((a + L) / 2, PI / 2, PI, 0.0, 0.0);
    curv2.HRefine({ 0.5 });

    RCurve3d curv3 = createCircle(a, PI / 2, PI, 0.0, 0.0);
    curv3.HRefine({ 0.5 });

    auto Surf = createSurface({ curv1, curv2, curv3 });

    cout << Surf << endl << endl << "-----------------------------" << endl;

    cout << isValidSurface(Surf) << "  valid?" << endl;

    Surf.HRefine(5, 5);
    cout << Surf << endl << endl << "-----------------------------" << endl;

    double E = 1e5;
    double nu = 0.3;

    cout << "calculate K:\n";
    SparseMatrix<double> K = assemblySparse(Surf, E, nu, PlaneStress);

    cout << K;
    Spsheetcircholeexactstress<double> S(a, 1.0);

    auto Sxx = [&](double x, double y) {return -S.xx(x, y) * (abs(x + L) <= 1e-6); };
    auto Sxy = [&](double x, double y) {return -S.xy(x, y) * (abs(x + L) <= 1e-6); };
    auto Syy = [&](double x, double y) {return  S.yy(x, y) * (abs(y - L) <= 1e-6); };
    auto Syx = [&](double x, double y) {return  S.xy(x, y) * (abs(y - L) <= 1e-6); };

    vector<pair<VectorX<double>, vector<int>>> fContainor;
    fContainor.push_back(applyNewmannBC(Surf, U0, Sxx, X));
    fContainor.push_back(applyNewmannBC(Surf, U0, Sxy, Y));
    fContainor.push_back(applyNewmannBC(Surf, U0, Syx, X));
    fContainor.push_back(applyNewmannBC(Surf, U0, Syy, Y));

    auto h = [](double x, double y) {return 0.0; };

    vector<pair<VectorX<double>, vector<int>>> uContainor;
    uContainor.push_back(applyDrchltBC(Surf, V1, h, Y));
    uContainor.push_back(applyDrchltBC(Surf, V0, h, X));

    VectorX<double> F = Eigen::VectorXd::Zero(K.rows());
    VectorX<double> U = Eigen::VectorXd::Zero(K.rows());

    assembly(fContainor, F);
    assembly(uContainor, U);

    vector<int> knowns;
    for (const auto& ui : uContainor)
        knowns.insert(knowns.end(), ui.second.begin(), ui.second.end());
    abab::solve(K, U, F, knowns);

    
    double u = 1.0, v = 0.0;
    for (int i = 0;i <= 50; i++) 
    {
        v += 0.02 * i;
        auto crd = surfacePoint(Surf, u, v);
        double x = crd[0], y = crd[1];

        double t = atan(-x / y) / PI * 180;
        
 /*       cout << "\ncalcu stress on (u, v) = (1.0, 0.0):\n";
        cout << "\ncalcu stress on (x, y) = (" << x << ", " << y << "):\n";*/

        //cout << "\nexact stress:\n";
        //cout << S.xx(x, y) << endl;
 /*       cout << S.yy(x, y) << endl;
        cout << S.xy(x, y) << endl;*/

        //cout << "\ncalcu stress:\n";
        Vector3d strain = getStrainVals(Surf, u, v, U);
        Matrix3d D = getElastMatrix(E, nu, PlaneStress);
        Vector3d stress = D * strain;
        printf("%2d : %2d  &2d\n", t, S.xx(x, y), stress[0]);
        //cout << stress;
    }
    return true;
}

bool testCreateSurface() {

    double a = 1; //圆弧半径
    double L = 5; //四边形边长

    RCurve3d curv1 = createLine({ vec3d(0, L, 0), vec3d(-L, L, 0), vec3d(-L, L, 0),  vec3d(-L, 0, 0) });

    RCurve3d curv2 = createCircle((a + L) / 2, PI / 2, PI, 0.0, 0.0);
    curv2.HRefine({ 0.5 });

    RCurve3d curv3 = createCircle(a, PI / 2, PI, 0.0, 0.0);
    curv3.HRefine({ 0.5 });

    auto surf = createSurface({ curv1, curv2, curv3 });

    cout << surf;

    return 1;
}

int main() {
    //后续可以采用宏替换启动程序

	//testCurve();
	//testSurface();
	//teststiffness();
	//testFunctionExpr();
    SheetCircHole();
    //testHRefine();
    //testCreateSurface();


	return -1;
}