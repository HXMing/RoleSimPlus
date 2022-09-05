#include "rolesim.h"
#include "dyhungarian.h"
#include <vector>
#include <iostream>
using namespace std;
#define MAX_RS 10000
int iterative_time = 0;
vector<vector<int>> OriginalRoleSim(vector<vector<int>>& G, int k, double beta)
{
    int n = G.size();
    vector<vector<int>> H(n, vector<int>(n, MAX_RS));
    vector<vector<int>> I(n);
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < n;++j) {
            if (i == j) continue;
            if (G[i][j] == 1)
                I[j].push_back(i);
        }
    }
    for (int i = 0;i < k;++i) {
        H = IterateRoleSim(H, I, beta);
        ++iterative_time;
    }
    return H;
}

vector<vector<int>> IterateRoleSim(vector<vector<int>>& H, vector<vector<int>>& I, double beta)
{
    int n = H.size();
    vector<vector<int>> H1(n, vector<int>(n));
    for (int v = 0;v < n;++v) {
        for (int u = 0;u < n;++u) {
            //cout << "iterative time: " << iterative_time << " || ";
            //cout << "nodes pair: " << v << " " << u << " || ";
            //cout << "size: " << I[v].size() << " " << I[u].size() << endl;
            int Nv = I[v].size(), Nu = I[u].size();
            int maxn = (Nv > Nu ? Nv : Nu);
            if (maxn == 0) {
                H1[v][u] = beta * MAX_RS;
                continue;
            }
            int** wei = new int* [maxn];
            for (int i = 0;i < maxn;++i)wei[i] = new int[maxn];
            for (int i = 0;i < maxn;++i)for (int j = 0;j < maxn;++j)wei[i][j] = 0;
            int* rd = new int[maxn];
            int* cd = new int[maxn];
            int* rs = new int[maxn];
            int* cs = new int[maxn];
            for (int i = 0;i < Nv;++i) {
                for (int j = 0;j < Nu;++j) {
                    wei[i][j] = H[I[v][i]][I[u][j]];
                }
            }
            int maxwei = Hungarian(maxn, wei, rd, cd, rs, cs);
            H1[v][u] = (1 - beta) * maxwei / maxn + beta * MAX_RS;
            for (int i = 0;i < maxn;++i) delete[] wei[i];
            delete[] wei;
            delete[] rd;
            delete[] cd;
            delete[] rs;
            delete[] cs;
        }
    }
    return H1;
}