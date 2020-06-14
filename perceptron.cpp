#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <algorithm>
#include <cassert>
#include <chrono>
using namespace std;
using namespace std::chrono;
using VD = vector<double>;

double calc(VD& set, double y_cur, VD& w) {
    double res = 0;
    // dot prod
    for (int i = 0; i < w.size()-1; ++i) {
        res += w[i]*set[i];
    }
    res += w[w.size()-1]; // bias, pretend that there is an extra 1 at the end of "set"
    // apply step function
    if (res < 1E-4 && res > -1E-4) { // very close to 0
        res = (-1)*y_cur; // give it an incentive to change
    }
    else if (res > 0) {
        res = 1;
    }
    else { // res < 0
        res = -1;
    }
    return res;
}

VD perceptron(vector<VD>& a, const VD& y, const int IT_MAX) {
    int m = a.size(); // y.size() == m
    int n = a[0].size();
    VD w(n+1, 0);
    for (int i = 0; i < IT_MAX; ++i) { // specified # of iterations
        for (int j = 0; j < m; ++j) { // for each training set
            double res = calc(a[j], y[j], w); // calculate output
            for (int k = 0; k < w.size()-1; ++k) { // update all weights
                w[k] += 0.1*(y[j]-res)*a[j][k]; // learning rate = 0.1
            }
            w[w.size()-1] += 0.1*(y[j]-res)*1; // update last weight (bias)
        }
    }
    // check if not learned
    for (int i = 0; i < m; ++i) {
        if (y[i] != calc(a[i], y[i], w)) {
            return VD();
        }
    }
    return w;
}
