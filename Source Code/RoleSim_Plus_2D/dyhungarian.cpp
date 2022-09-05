#include "dyhungarian.h"
#include <algorithm>
using namespace std;

static const int MAXN = 10000;
static const int INF = 0x3f3f3f3f;
static int visr[MAXN], visc[MAXN];
static int slack[MAXN];

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
            if (dfs_r(i, dim, wei, rd, cd, rs, cs)) break;
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

int dfs_r(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs) {
    visr[r] = 1;
    for (int j = 0; j < dimx; ++j) {
        if (visc[j]) continue;
        int theta = rd[r] + cd[j] - wei[r][j];
        if (theta == 0) {
            visc[j] = 1;
            if (cs[j] == -1 || dfs_r(cs[j], dimx, wei, rd, cd, rs, cs)) {
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

int dfs_c(int c, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs)
{
    visc[c] = 1;
    for (int i = 0;i < dimx;++i) {
        if (visr[i]) continue;
        int theta = rd[i] + cd[c] - wei[i][c];
        if (theta == 0) {
            visr[i] = 1;
            if (rs[i] == -1 || dfs_c(rs[i], dimx, wei, rd, cd, rs, cs)) {
                rs[i] = c, cs[c] = i;
                return 1;
            }
        }
        else {
            if (theta < slack[i])
                slack[i] = theta;
        }
    }
    return 0;
}

int dymatch_row(int start, int end, int f_dim, int dim, int** wei, int* rd, int* cd, int* rs, int* cs)
{
    for (int i = start;i < end;++i) {
        int maxval = -INF;
        for (int j = 0;j < f_dim;++j) {
            maxval = max(maxval, wei[i][j] - cd[j]);
        }
        rd[i] = maxval;
    }
    for (int i = start;i < end;++i) {
        for (int j = 0;j < f_dim;++j) slack[j] = INF;
        while (1)
        {
            for (int j = 0;j < f_dim;++j) visr[j] = 0, visc[j] = 0;
            if (dfs_r(i, f_dim, wei, rd, cd, rs, cs))break;
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
        cs[k] = -1;
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
            if (dfs_r(k, k + 1, wei, rd, cd, rs, cs)) break;
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

int dymatch_col(int start, int end, int f_dim, int dim, int** wei, int* rd, int* cd, int* rs, int* cs)
{
    for (int j = start;j < end;++j) {
        int maxval = -INF;
        for (int i = 0;i < f_dim;++i) {
            maxval = max(maxval, wei[i][j] - rd[i]);
        }
        cd[j] = maxval;
    }
    for (int j = start;j < end;++j) {
        for (int idx = 0;idx < f_dim;++idx)slack[idx] = INF;
        while (1)
        {
            for (int idx = 0;idx < f_dim;++idx) visr[idx] = 0, visc[idx] = 0;
            if (dfs_c(j, f_dim, wei, rd, cd, rs, cs))break;
            double d = INF;
            for (int i = 0;i < f_dim;++i) {
                if (!visr[i] && d > slack[i])
                    d = slack[i];
            }
            for (int i = 0;i < f_dim;++i) {
                if (visc[i])cd[i] -= d;
                if (visr[i])rd[i] += d;
                else slack[i] -= d;
            }
        }
    }
    for (int k = f_dim;k < dim;++k) {
        cs[k] = -1;
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
            if (dfs_r(k, k + 1, wei, rd, cd, rs, cs)) break;
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