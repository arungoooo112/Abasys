
#include <iostream>
#include "abasys.h"
#include "ExactSolution/SPSheetCircHoleExactStress.h"
// 一个分析样例
using namespace abab;
using namespace std;

enum {X,Y};
void SheetCircHole()
{

    //---------CAD model----------------------------------------
    //----------------------------------------------------------

    cout << "-------------CAD建模，构造参数曲面----------------\n\n";

    double a = 1; //圆弧半径
    double L = 5; //四边形边长

    RCurve3d curv1 = createLine({ vec3d(0, L, 0), vec3d(-L, L, 0), vec3d(-L, L, 0),  vec3d(-L, 0, 0) });

    RCurve3d curv2 = createCircle((a + L) / 2, PI / 2, PI, 0.0, 0.0);
    curv2.HRefine({ 0.5 });

    RCurve3d curv3 = createCircle(a, PI / 2, PI, 0.0, 0.0);
    curv3.HRefine({ 0.5 });

    auto Surf = createSurface({ curv1, curv2, curv3 });

    cout << "初始参数曲面：\n";
    cout << Surf << endl << endl;

    if (!isValidSurface(Surf)) {
        cerr << "生成参数曲面无效" << endl;
        return ;
    }

    cout << "对参数曲面加密\n";
    Surf.HRefine(20, 20);
    cout << "加密后的参数曲面：\n";
    cout << Surf << endl << endl;

    //---------assembly K----------------------------------------
    //-----------------------------------------------------------

    cout << "------------------calculate K---------------------\n\n";

    double E = 1e5;
    double nu = 0.3;

    auto K = assemblyMatrix(Surf, E, nu, PlaneStress); // K 为密集矩阵
    //auto K = assemblySparse(Surf, E, nu, PlaneStress); // K 为稀疏矩阵 
    cout << K;

    //---------apply boundary conditions-------------------------
    //-----------------------------------------------------------

    cout << "------------apply boundary conditions----------------\n\n";

    Spsheetcircholeexactstress<double> S(a, 1.0); // 精确解

    auto Sxx = [&](double x, double y) {return -S.xx(x, y) * (abs(x + L) <= 1e-6); };
    auto Sxy = [&](double x, double y) {return -S.xy(x, y) * (abs(x + L) <= 1e-6); };
    auto Syy = [&](double x, double y) {return  S.yy(x, y) * (abs(y - L) <= 1e-6); };
    auto Syx = [&](double x, double y) {return  S.xy(x, y) * (abs(y - L) <= 1e-6); };

    vector<pair<VectorX<double>, vector<int>>> fContainor;
    fContainor.push_back(applyNewmannBC(Surf, U0, Sxx, X));
    fContainor.push_back(applyNewmannBC(Surf, U0, Sxy, Y));
    fContainor.push_back(applyNewmannBC(Surf, U0, Syx, X));
    fContainor.push_back(applyNewmannBC(Surf, U0, Syy, Y));

    auto h = [](double , double ) {return 0.0; };

    vector<pair<VectorX<double>, vector<int>>> uContainor;
    uContainor.push_back(applyDrchltBC(Surf, V1, h, Y));
    uContainor.push_back(applyDrchltBC(Surf, V0, h, X));

    VectorX<double> F = Eigen::VectorXd::Zero(K.rows());
    VectorX<double> U = Eigen::VectorXd::Zero(K.rows());

    assembly(fContainor, F);
    assembly(uContainor, U);


    //---------solve---------------------------------------------
    //-----------------------------------------------------------

    cout << "---------------------solve---------------------\n\n";

    vector<int> knowns; //已知的位移点
    for (const auto& ui : uContainor)
        knowns.insert(knowns.end(), ui.second.begin(), ui.second.end());

    //abab::solveByLarge(K, U, F, knowns);
    abab::solve(K, U, F, knowns);

    //---------post process--------------------------------------
    //-----------------------------------------------------------

    cout << "-----------------post process-------------------\n\n";

    fstream out;
    out.open("post.csv"); //创建文件

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
    out.close();

}
