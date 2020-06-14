#include <vector>
#include <iostream>
#include <climits>
#include <stack>
using namespace std;

typedef vector<int> VI;
typedef vector<long> VL;

vector<VL> n_hanoi;
vector<VI> k_hanoi;

// Initialize n_hanoi and k_hanoi
void hanoi_init() {
    // initialize n_hanoi
    VL init_row(10001, -1); // vector of 10001 -1's
    init_row[1] = 1;
    for (int i = 0; i <= 10; ++i) n_hanoi.push_back(init_row);
    n_hanoi[0][1] = -1; // twr = 0 is all impossible
    n_hanoi[1][1] = -1; // twr = 1 is all impossible
    // initialize k_hanoi
    VI init_row_k(10001, -1); // vector of 10001 -1's
    for (int i = 0; i <= 10; ++i) k_hanoi.push_back(init_row_k);
    k_hanoi[0][1] = -1; // twr = 0 is all impossible
    k_hanoi[1][1] = -1;
    // dp to fill n_hanoi
    for (int twrs = 3; twrs <= 10; ++twrs) {
        for (int dsks = 2; dsks <= 10000; ++dsks) {
            long cur_min = LONG_MAX;
            int best_k = -1;
            bool changed = false;
            for (int k = 1; k <= dsks-1; ++k) {
                long x = n_hanoi[twrs-1][dsks-k];
                long y = n_hanoi[twrs][k];
                if ((x != -1) && (y != -1)) { // actually valid
                    if (y <= (LONG_MAX-x)/2) { // no overflow
                        if (2*y + x <= cur_min) { // better solution found
                            changed = true;
                            if (2*y + x < cur_min) best_k = k;
                            cur_min = 2*y + x;
                        }
                    }
                }
            }
            n_hanoi[twrs][dsks] = (changed) ? cur_min : -1;
            k_hanoi[twrs][dsks] = best_k;
        }
    }
}

// Fill moves so that each element is a two-integer VI describing the move
// You may assume that, initially, aux[] = {0, 1, 2, ..., n_twrs - 1}
void hanoi(vector<VI>& moves, int n_twrs, int n_dsks, VI& aux) {
    if (n_dsks == 1) {
        vector<int> v;
        v.push_back(aux[0]);
        v.push_back(aux[1]);
        moves.push_back(v);
    }
    if (n_dsks > 1) {
        int temp_total = n_dsks;
        n_dsks = k_hanoi[n_twrs][n_dsks];
        swap(aux[1], aux[n_twrs-1]);
        hanoi(moves, n_twrs, n_dsks, aux);
        swap(aux[1], aux[n_twrs-1]);
        hanoi(moves, n_twrs-1, temp_total-n_dsks, aux);
        swap(aux[0], aux[n_twrs-1]);
        hanoi(moves, n_twrs, n_dsks, aux);
        swap(aux[0], aux[n_twrs-1]);
    }
}
