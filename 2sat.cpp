/*
Jan 2020
Solves 2-SAT using CNF and SCCs
*/

#include <vector>
#include <iostream>
#include <utility>
#include <stack>

using namespace std;

typedef vector<vector<int> > VVI;
typedef vector<int> VI;
typedef vector<int> VB;
typedef stack<int> SI;

int V, C, a, b, ai, bi, aneg, bneg;
VVI adj, adj_T;
VB vst, vst_T, res;
SI s;
VI scc;
int cur_scc = 0;

void dfs1(int v) {
    vst[v] = true;
    for (int node : adj[v]) {
        if (!vst[node]) {
            dfs1(node);
        }
    }
    s.push(v);
}

void dfs2(int v, int scc_index) {
    scc[v] = scc_index;
    for (int node : adj_T[v]) {
        if (scc[node] == -1) {
            dfs2(node, scc_index);
        }
    }
}

int main() {
    cin >> V >> C;
    adj = VVI(2*V);
    adj_T = VVI(2*V);
    vst = VB(2*V);
    vst_T = VB(2*V);
    scc = VI(2*V, -1);
    res = VB(V);

    // 1  -1   2  -2   3  -3
    // 0   1   2   3   4   5
    // positive x --> 2x-2
    // negative x --> 2|x|-1

    // make graphs
    for (int i = 0; i < C; ++i) {
        cin >> a >> b; // (a or b) == (~a -> b) and (~b -> a)
        ai = (a < 0) ? -2*a-1 : 2*a-2;
        bi = (b < 0) ? -2*b-1 : 2*b-2;
        // xneg: (4|x|-3) - xi
        aneg = (4*abs(a)-3) - ai;
        bneg = (4*abs(b)-3) - bi;
        // add to G
        adj[aneg].push_back(bi);
        adj[bneg].push_back(ai);
        // add to G_T
        adj_T[bi].push_back(aneg);
        adj_T[ai].push_back(bneg);
    }

    // kosaraju - step 1
    for (int i = 0; i < vst.size(); ++i) {
        if (!vst[i]) {
            dfs1(i);
        }
    }

    // kosaraju - step 2
    for (int i = 0; i < 2*V; ++i) {
        int v = s.top();
        s.pop();
        if (scc[v] == -1) {
            dfs2(v, cur_scc++);
        }
    }

    // check if unsatisfiable, or assign truth values
    for (int i = 0; i < scc.size(); i += 2) {
        if (scc[i] == scc[i+1]) {
            cout << "Not satisfiable\n";
            exit(0);
        }
        res[i/2] = (scc[i] > scc[i+1]);
    }

    // print
    for (int i = 0; i < res.size(); ++i) {
        if (res[i]) {
            cout << "T\n";
        }
        else {
            cout << "F\n";
        }
    }
}
