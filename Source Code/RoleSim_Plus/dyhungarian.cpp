#include "dyhungarian.h"
#include <cstring>
#include <cmath>
#include <iostream>
#include <set>
using namespace std;

static const int MAXN = 10000;
static const int INF = 0x3f3f3f3f;
static const double eps = 1e-6;
static int visr[MAXN], visc[MAXN];
static int slack[MAXN];


double Hungarian(int dim, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs)
{
    for (int j = 0;j < dim;++j) cs[j] = -1;
    for (int j = 0;j < dim;++j) cd[j] = 0.0;
    for (int i = 0;i < dim;++i) {
        rd[i] = wei[i][0];
        for (int j = 1;j < dim;++j) {
            if (rd[i] < wei[i][j]) rd[i] = wei[i][j];
        }
    }
    for (int i = 0;i < dim;++i) {
        for (int j = 0;j < dim;++j) slack[j] = INF;
        while (1) {
            //memset(visr, 0, sizeof(visr));
           // memset(visc, 0, sizeof(visc));
            for (int i = 0;i < dim;++i) visr[i] = visc[i] = 0;
            if (dfs(i,dim,wei,rd,cd,rs,cs)) break;
            double d = INF;
            for (int j = 0;j < dim;++j) {
                if (!visc[j] && d > slack[j])
                    d = slack[j];
            }
            for (int j = 0;j < dim;++j) {
                if (visr[j]) rd[j] -= d;
                if (visc[j]) cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    double max_weight = 0.0;
    for (int i = 0;i < dim;++i) {
        max_weight += wei[cs[i]][i];
    }
    return max_weight;
}

double dyHungarian_r(int dim, int r, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs)
{
    rd[r] = wei[r][0] - cd[0];
    for (int j = 1;j < dim;++j) {
        if (wei[r][j] - cd[j] > rd[r])
            rd[r] = wei[r][j] - cd[j];
    }
    cs[rs[r]] = -1;
    for (int i = r;i <=r;++i) {
        for (int j = 0;j < dim;++j) slack[j] = INF;
        while (1) {
            //memset(visr, 0, sizeof(visr));
           // memset(visc, 0, sizeof(visc));
            for (int i = 0;i < dim;++i) visr[i] = visc[i] = 0;
            if (dfs(i, dim, wei, rd, cd, rs, cs)) break;
            double d = INF;
            for (int j = 0;j < dim;++j) {
                if (!visc[j] && d > slack[j])
                    d = slack[j];
            }
            for (int j = 0;j < dim;++j) {
                if (visr[j]) rd[j] -= d;
                if (visc[j]) cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    double max_weight = 0.0;
    for (int i = 0;i < dim;++i) {
        max_weight += wei[cs[i]][i];
    }
    return max_weight;
}

double dyHungarian_c(int dim, int c, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs) {
    cd[c] = wei[0][c] - rd[0];
    for (int i = 1;i < dim;++i) {
        if (wei[i][c] - rd[i] > cd[c])
            cd[c] = wei[i][c] - rd[i];
    }
    int r = cs[c];
    cs[c] = -1;
    for (int i = r;i <= r;++i) {
        for (int j = 0;j < dim;++j) slack[j] = INF;
        while (1) {
            //memset(visr, 0, sizeof(visr));
           // memset(visc, 0, sizeof(visc));
            for (int i = 0;i < dim;++i) visr[i] = visc[i] = 0;
            if (dfs(i, dim, wei, rd, cd, rs, cs)) break;
            double d = INF;
            for (int j = 0;j < dim;++j) {
                if (!visc[j] && d > slack[j])
                    d = slack[j];
            }
            for (int j = 0;j < dim;++j) {
                if (visr[j]) rd[j] -= d;
                if (visc[j]) cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    double max_weight = 0.0;
    for (int i = 0;i < dim;++i) {
        max_weight += wei[cs[i]][i];
    }
    return max_weight;
}

int dfs(int r, int dim, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs) {
    visr[r] = 1;
    for (int j = 0;j < dim;++j) {
        if (visc[j]) continue;
        double theta = rd[r] + cd[j] - wei[r][j];
        if (fabs(theta) <= eps) {
            visc[j] = 1;
            if (cs[j] == -1 || dfs(cs[j],dim, wei,rd,cd,rs,cs)) {
                cs[j] = r, rs[r] = j;
                return 1;
            }
        }
        else {
            if (theta < slack[j])
                slack[j] = theta;
        }
    }
    return 0;
}

double Hungarian(Mi& mi)
{
    for (int j = 0;j < mi.dim;++j) mi.cs[j] = -1;
    for (int j = 0;j < mi.dim;++j) mi.cd[j] = 0.0;
    for (int i = 0;i < mi.dim;++i) {
        mi.rd[i] = mi.wei[i][0];
        for (int j = 1;j < mi.dim;++j) {
            if (mi.rd[i] < mi.wei[i][j]) mi.rd[i] = mi.wei[i][j];
        }
    }
    for (int i = 0;i < mi.dim;++i) {
        for (int j = 0;j < mi.dim;++j) slack[j] = INF;
        while (1) {
            //memset(visr, 0, sizeof(visr));
           // memset(visc, 0, sizeof(visc));
            for (int i = 0;i < mi.dim;++i) visr[i] = visc[i] = 0;
            if (dfs(i, mi)) break;
            double d = INF;
            for (int j = 0;j < mi.dim;++j) {
                if (!visc[j] && d > slack[j])
                    d = slack[j];
            }
            for (int j = 0;j < mi.dim;++j) {
                if (visr[j]) mi.rd[j] -= d;
                if (visc[j]) mi.cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    double max_weight = 0.0;
    for (int i = 0;i < mi.dim;++i) {
        max_weight += mi.wei[mi.cs[i]][i];
    }
    return max_weight;
}


int dfs(int r, Mi& mi) {
    visr[r] = 1;
    for (int j = 0;j < mi.dim;++j) {
        if (visc[j]) continue;
        double theta = mi.rd[r] + mi.cd[j] - mi.wei[r][j];
        if (fabs(theta) <= eps) {
            visc[j] = 1;
            if (mi.cs[j] == -1 || dfs(mi.cs[j], mi)) {
                mi.cs[j] = r, mi.rs[r] = j;
                return 1;
            }
        }
        else {
            if (theta < slack[j])
                slack[j] = theta;
        }
    }
    return 0;
}

//double dymatch_col(int idx, Mi& mf, Mi& mc)
//{
//    for (int j = idx;j < mf.dim;++j) {
//        double maxval = -INF;
//        for (int i = 0;i < mf.dim;++i) {
//            maxval = max(maxval, mc.wei[i][j] - mc.rd[i]);
//        }
//        mc.cd[j] = maxval;
//    }
//    set<int> matched;
//    for (int j = 0;j < idx;++j) {
//        matched.insert(mc.cs[j]);
//    }
//    for (int i = 0;i < mf.dim;++i) {
//        if (matched.find(i) != matched.end()) continue;
//        for (int j = 0;j < mf.dim;++j) slack[j] = INF;
//        while (1) {
//            for (int j = 0;j < mf.dim;++j) visr[j] = 0, visc[j] = 0;
//            if (dfs(i, mf.dim, mc.wei, mc.rs, mc.cs, mc.rd, mc.cd))break;
//            double d = INF;
//            for (int j = 0;j < mf.dim;++j) {
//                if (!visc[j] && d > slack[j])
//                    d = slack[j];
//            }
//            for (int j = 0;j < mf.dim;++j) {
//                if (visr[j]) mc.rd[j] -= d;
//                if (visc[j]) mc.cd[j] += d;
//                else slack[j] -= d;
//            }
//        }
//    }
//    for (int k = mf.dim;k < mc.dim;++k) {
//        double maxval = -INF;
//        for (int i = 0;i < k;++i) {
//            maxval = max(maxval, mc.wei[i][k] - mc.rd[i]);
//        }
//        mc.cd[k] = max(maxval, mc.wei[k][k]);
//        maxval = -INF;
//        for (int i = 0; i <= k; ++i) {
//            maxval = max(maxval, mc.wei[k][i] - mc.cd[i]);
//        }
//        mc.rd[k] = maxval;
//        for (int i = 0; i <= mc.dim; ++i) slack[i] = INF;
//        while (1) {
//            for (int j = 0; j <= mc.dim; ++j) visr[j] = 0, visc[j] = 0;
//            if (dfs(k, k + 1, mc.wei, mc.rs, mc.cs, mc.rd, mc.cd)) break;
//            double d = INF;
//            for (int j = 0; j <= k; ++j) {
//                if (!visc[j]&&d>slack[j]) {
//                    d = slack[j];
//                }
//            }
//            for (int j = 0; j <= k; ++j) {
//                if (visr[j]) mc.rd[j] -= d;
//                if (visc[j]) mc.cd[j] += d;
//                else slack[j] -= d;
//            }
//        }
//    }
//    double maxwei = 0.0;
//    for (int i = 0;i < mc.dim;++i) {
//        maxwei += mc.wei[mc.cs[i]][i];
//    }
//    return maxwei;
//}

int dfs(int r, int dimx, double** wei, int* rs, int* cs, double* rd, double* cd) {
    visr[r] = 1;
    for (int j = 0; j < dimx; ++j) {
        if (visc[j]) continue;
        double theta = rd[r] + cd[j] - wei[r][j];
        if (fabs(theta) <= eps) {
            visc[j] = 1;
            if (cs[j] == -1 || dfs(cs[j], dimx, wei, rs, cs, rd, cd)) {
                cs[j] = r, rs[r] = j;
                return 1;
            }
        }
        else {
            if (theta < slack[j])
                slack[j] = theta;
        }
    }
    return 0;
}

int Hungarian(int dim, int** wei, int* rd, int* cd, int* rs, int* cs)
{
    for (int j = 0;j < dim;++j) cs[j] = -1;
    for (int j = 0;j < dim;++j) cd[j] = 0;
    for (int i = 0;i < dim;++i) {
        rd[i] = wei[i][0];
        for (int j = 1;j < dim;++j) {
            if (rd[i] < wei[i][j]) rd[i] = wei[i][j];
        }
    }
    for (int i = 0;i < dim;++i) {
        for (int j = 0;j < dim;++j) slack[j] = INF;
        while (1) {
            //memset(visr, 0, sizeof(visr));
           // memset(visc, 0, sizeof(visc));
            for (int i = 0;i < dim;++i) visr[i] = visc[i] = 0;
            if (dfs(i, dim, wei, rd, cd, rs, cs)) break;
            int d = INF;
            for (int j = 0;j < dim;++j) {
                if (!visc[j] && d > slack[j])
                    d = slack[j];
            }
            for (int j = 0;j < dim;++j) {
                if (visr[j]) rd[j] -= d;
                if (visc[j]) cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    int max_weight = 0;
    for (int i = 0;i < dim;++i) {
        max_weight += wei[cs[i]][i];
    }
    return max_weight;
}

int dymatch_col(int f_dim, int dim, int idx, int** wei, int* rd, int* cd, int* rs, int* cs)
{
    for (int j = idx;j < f_dim;++j) {
        int maxval = -INF;
        for (int i = 0;i < f_dim;++i) {
            maxval = max(maxval, wei[i][j] - rd[i]);
        }
        cd[j] = maxval;
    }
    set<int> matched;
    for (int j = 0;j < idx;++j) {
        matched.insert(cs[j]);
    }
    for (int i = 0;i < f_dim;++i) {
        if (matched.find(i) != matched.end()) continue;
        for (int j = 0;j < f_dim;++j) slack[j] = INF;
        while (1) {
            for (int j = 0;j < f_dim;++j) visr[j] = 0, visc[j] = 0;
            if (dfs(i, f_dim, wei, rd, cd, rs, cs))break;
            double d = INF;
            for (int j = 0;j < f_dim;++j) {
                if (!visc[j] && d > slack[j])
                    d = slack[j];
            }
            for (int j = 0;j < f_dim;++j) {
                if (visr[j]) rd[j] -= d;
                if (visc[j]) cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    for (int k = f_dim;k < dim;++k) {
        int maxval = -INF;
        for (int i = 0;i < k;++i) {
            maxval = max(maxval, wei[i][k] - rd[i]);
        }
        cd[k] = max(maxval, wei[k][k]);
        maxval = -INF;
        for (int i = 0; i <= k; ++i) {
            maxval = max(maxval, wei[k][i] - cd[i]);
        }
        rd[k] = maxval;
        for (int i = 0; i <= dim; ++i) slack[i] = INF;
        while (1) {
            for (int j = 0; j <= dim; ++j) visr[j] = 0, visc[j] = 0;
            if (dfs(k, k + 1, wei, rd, cd, rs, cs)) break;
            double d = INF;
            for (int j = 0; j <= k; ++j) {
                if (!visc[j] && d > slack[j]) {
                    d = slack[j];
                }
            }
            for (int j = 0; j <= k; ++j) {
                if (visr[j]) rd[j] -= d;
                if (visc[j]) cd[j] += d;
                else slack[j] -= d;
            }
        }
    }
    int maxwei = 0;
    for (int i = 0;i < dim;++i) {
        maxwei += wei[cs[i]][i];
    }
    return maxwei;
}

int dfs(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs){
    visr[r] = 1;
    for (int j = 0; j < dimx; ++j) {
        if (visc[j]) continue;
        int theta = rd[r] + cd[j] - wei[r][j];
        if (theta == 0) {
            visc[j] = 1;
            if (cs[j] == -1 || dfs(cs[j], dimx, wei, rd, cd, rs, cs)) {
                cs[j] = r, rs[r] = j;
                return 1;
            }
        }
        else {
            if (theta < slack[j])
                slack[j] = theta;
        }
    }
    return 0;
}