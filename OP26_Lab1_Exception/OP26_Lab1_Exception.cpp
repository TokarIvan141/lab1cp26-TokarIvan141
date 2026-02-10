#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <stdexcept>

using namespace std;

struct TablePoint {
    double x, t, u;
};

class LabSolver {
private:
    void saveToFile(string name, vector<TablePoint> data) {
        ofstream f(name);
        if (f.is_open()) {
            for (size_t i = 0; i < data.size(); ++i) {
                f << data[i].x << " " << data[i].t << " " << data[i].u << "\n";
            }
            f.close();
        }
    }

    void prepareEnvironment() {
        vector<TablePoint> t1 = {{-1.0, -4.935, 1.935}, {0.0, 0.0, 4.571}, {1.0, 0.0, 3.0}};
        saveToFile("dat_X_1_1.dat", t1);

        vector<TablePoint> t2 = {{0.0, -4.935, 1.935}, {0.5, 3.356, 1.388}, {1.0, 5.89, 0.377}};
        saveToFile("dat_X1_00.dat", t2);

        vector<TablePoint> t3 = {{0.0, -4.935, 1.935}, {-0.5, -0.141, 0.976}, {-1.0, 3.48, 0.252}};
        saveToFile("dat_X00_1.dat", t3);
    }

    struct DataResult { double t, u; };

    DataResult get_table_data(double x) {
        string fname;
        double tx = x;
        if (abs(x) <= 1.0) fname = "dat_X_1_1.dat";
        else { tx = 1.0 / x; fname = (x < -1.0) ? "dat_X00_1.dat" : "dat_X1_00.dat"; }

        ifstream file(fname);
        if (!file.is_open()) throw runtime_error("file_error");

        vector<TablePoint> pts;
        double rx, rt, ru;
        while (file >> rx >> rt >> ru) pts.push_back({rx, rt, ru});
        file.close();

        if (pts.empty()) throw runtime_error("data_empty");

        for (size_t i = 0; i < pts.size(); ++i)
            if (abs(tx - pts[i].x) < 1e-7) return {pts[i].t, pts[i].u};

        for (size_t i = 0; i < pts.size() - 1; ++i) {
            if ((tx >= pts[i].x && tx <= pts[i+1].x) || (tx <= pts[i].x && tx >= pts[i+1].x)) {
                double t = pts[i].t + (pts[i+1].t - pts[i].t) * (tx - pts[i].x) / (pts[i+1].x - pts[i].x);
                double u = pts[i].u + (pts[i+1].u - pts[i].u) * (tx - pts[i].x) / (pts[i+1].x - pts[i].x);
                return {t, u};
            }
        }
        throw runtime_error("range_error");
    }

    double srz(double x, double y, double z) {
        DataResult dx = get_table_data(x), dy = get_table_data(y), dz = get_table_data(z);
        return (x > y) ? (dx.t + dz.u - dy.t) : (dy.t + dy.u - dz.u);
    }

    double glr(double x, double y) {
        if (abs(x) < 1.0) return x;
        if (abs(y) < 1.0) return y;
        double v = x * x + y * y - 4.0;
        if (v < 0.1) throw 2;
        return y / sqrt(v);
    }

    double gold(double x, double y) {
        if (x > y && abs(y) > 1e-9) return x / y;
        if (x < y && abs(x) > 1e-9) return y / x;
        throw 2;
    }

    double grs(double x, double y) {
        return 0.1389 * srz(x + y, gold(x, y), glr(x, x * y)) +
               1.8389 * srz(x - y, gold(y, x / 5.0), glr(5.0 * x, x * y)) +
               0.83 * srz(x - 0.9, glr(y, x / 5.0), gold(5.0 * y, y));
    }

    double alg1(double x, double y, double z) {
        return x * x * grs(y, z) + y * y * grs(x, z) + 0.33 * x * y * grs(x, z);
    }

    double gold1(double x, double y) {
        if (x > y && abs(y) > 0.1) return x / y;
        if (x <= y && abs(x) > 0.1) return y / x;
        return (x < y && abs(x) > 0.1) ? 0.15 : 0.1;
    }

    double grsl(double x, double y) {
        return 0.14 * srz(x + y, gold1(x, y), (abs(x) < 1.0 ? x : x * y)) +
               1.83 * srz(x - y, gold1(y, x / 5.0), (abs(4 * x) < 1.0 ? 4 * x : x * y)) +
               0.83 * srz(x, (abs(y) < 1.0 ? y : x / 4.0), gold1(4 * y, y));
    }

    double alg2(double x, double y, double z) {
        return x * grsl(x, y) + y * grsl(y, z) + z * grsl(z, x);
    }

    double alg3(double x, double y, double z) {
        return 1.3498 * z + 2.2362 * y - 2.348 * x * y;
    }

public:
    double solve(double x, double y, double z) {
        prepareEnvironment();
        try {
            try {
                return alg1(x, y, z);
            } catch (int e) {
                return alg2(x, y, z);
            }
        } catch (...) {
            return alg3(x, y, z);
        }
    }
};

int main() {
    double x, y, z;
    if (!(cin >> x >> y >> z)) return 0;
    LabSolver solver;
    cout << fixed << setprecision(6) << solver.solve(x, y, z) << endl;
    return 0;
}