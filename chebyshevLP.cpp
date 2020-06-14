#include <vector>
#include <iostream>
#include <utility>
#include <algorithm>
#include <cassert>
#include <cmath>
using namespace std;
using VD = vector<double>;

void f(vector<VD>& A, VD& b, VD& c, const VD& x, const VD& y) {
    int n = x.size(); // n = number of equations
    vector<VD> segments;
    VD row(5, 0); // A is n x 5
    for (int i = 0; i < n; ++i) {
        A.push_back(row);
    }
    for (int i = 0; i < n; ++i) {
        // find segment equation
        VD res;
        double x1 = x.at(i);
        double x2 = x.at((i+1)%n);
        double y1 = y.at(i);
        double y2 = y.at((i+1)%n);
        double aa = (-1) * (y1-y2);
        double bb = (-1) * (x2-x1);
        double cc = (-1) * ((y1-y2)*x1 - (x1-x2)*y1);
        res.push_back(aa);
        res.push_back(bb);
        res.push_back(cc);
        b.push_back(cc);
        segments.push_back(res);
    }
    // c
    for (int i = 0; i < 4; ++i) {
        c.push_back(0);
    }
    c.push_back(1);
    // A
    for (int i = 0; i < n; ++i) {
        A[i][0] = segments[i][0];
        A[i][1] = (-1)*segments[i][0];
        A[i][2] = segments[i][1];
        A[i][3] = (-1)*segments[i][1];
        A[i][4] = sqrt(segments[i][0]*segments[i][0] + segments[i][1]*segments[i][1]);
    }
}
