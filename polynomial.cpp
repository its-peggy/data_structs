#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>
#include <utility>

#define EPS 1E-7
#define MAX_ITER 100

using namespace std;

typedef vector<double> VD;

class Polynomial {
public:
    Polynomial(const VD& coeffs) : c(coeffs) {} // constructor
    Polynomial(const Polynomial& p2) : c(p2.c) {} // copy constructor
    double operator=(const Polynomial& p2) { // operator= overloading
        c = p2.c;
        return 0;
    };
    Polynomial find_derivative() {
        VD d(c.size()-1);
        for (int i = 0; i < d.size(); ++i) {
            d[i] = c[i+1]*(i+1);
        }
        return Polynomial(d);
    }
    double operator()(double x) const {
        double res = 0;
        for (int i = c.size()-1; i >= 0; --i) {
            res += c[i];
            res *= x;
        }
        return res;
    }
private:
    VD c;
};

template<typename F, typename D>
double newton_solve(double x0, const F& f, const D& df) {
    for (int i = 0; i < MAX_ITER; ++i) {
        double dx = f(x0)/df(x0);
        if (abs(df(x0)) < EPS * abs(f(x0))) {
            cerr << "Slope too close to 0! \n";
            exit(-1);
        }
        if (abs(dx) < EPS) {
            return x0;
        }
        x0 -= dx;
    }
    cerr << "Max iterations exceeded! \n";
    return 0;
}

int main() {
    int deg;
    cout << "Enter polynomial degree:\n";
    cin >> deg;
    VD v(deg+1);
    cout << "Enter coefficients (highest degree term first):\n";
    for (int i = deg; i >= 0; --i) {
        cin >> v[i];
    }
    double x0;
    cout << "Enter x_0:\n";
    cin >> x0;
    Polynomial f(v);
    Polynomial df = f.find_derivative();
    double root = newton_solve(x0, f, df);
    cout << "Root = " << root << "\n";
    return 0;
}
