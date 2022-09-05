#include "mst.h"
#include <algorithm>
#include <iostream>
using namespace std;

#define oo 0x3f3f3f3f

static const int MAXN = 10000;
int father[MAXN];

bool cmp(vector<int>& a, vector<int>& b) {
    return a[2] < b[2];
}


vector<edge> prim(vector<vector<int>>& weight)
{
    int n = weight.size();
    int* lowcost = new int[n];
    int* nearvex = new int[n];
    vector<edge> res;
    lowcost[0] = 0;
    nearvex[0] = -1;
    for (int i = 1;i < n;++i) {
        lowcost[i] = weight[i][0];
        nearvex[i] = 0;
    }
    int k = n - 1;
    while (k--) {//循环n-1次找到n-1条边
        int minval = oo, v = -1;
        for (int i = 0;i < n;++i) {
            if (nearvex[i] != -1 && lowcost[i] < minval) {
                minval = lowcost[i];
                v = i;
            }
        }
        res.push_back(edge(nearvex[v], v, weight[nearvex[v]][v]));
        for (int i = 0;i < n;++i) {
            if (nearvex[i] != -1 && lowcost[i] > weight[i][v])
                lowcost[i] = weight[i][v], nearvex[i] = v;
        }
        nearvex[v] = -1;//加入到最小生成树
    }
    delete[] lowcost;
    delete[] nearvex;
    return res;
}


void init(int n) {
    for (int i = 0;i < n;++i) father[i] = i;
}

int find(int i) {
    if (father[i] != i)
        father[i] = find(father[i]);
    return father[i];
}

void join(int i, int j) {
    int fi = find(i), fj = find(j);
    father[fi] = fj;
}

vector<edge> kruskal(int n, vector<vector<int>>& edges)
{
    vector<edge> res;
    sort(edges.begin(), edges.end(), cmp);
    init(n);
    int k = n - 1;
    for (int i = 0;i < edges.size();++i) {
        int a = edges[i][0], b = edges[i][1], w = edges[i][2];
        int fa = find(a), fb = find(b);
        if (fa == fb) continue;
        join(a, b);
        res.push_back(edge(a, b, w));
        k -= 1;
        if (k == 0)break;
    }
    return res;
}


vector<vector<int>> kruskal2(int n, vector<vector<int>>& edges)
{
    vector<vector<int>> res(n);
    sort(edges.begin(), edges.end(), cmp);
    init(n);
    int k = n - 1;
    for (int i = 0;i < edges.size();++i) {
        int a = edges[i][0], b = edges[i][1];
        int fa = find(a), fb = find(b);
        if (fa == fb) continue;
        join(a, b);
        res[a].push_back(b);
        res[b].push_back(a);
        k -= 1;
        if (k == 0)break;
    }
    return res;
}
