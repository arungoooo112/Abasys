
#include "abasys.h"
#include "ExactSolution/SPSheetCircHoleExactStress.h"
// 一个分析样例
using namespace abab;

enum {X,Y};
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
    Surf.HRefine(20, 20);
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

    return true;
}