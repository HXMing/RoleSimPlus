#include "rolesimplus.h"
#include "mst.h"
#include <unordered_map>
#include <set>
#include <queue>
#include <iostream>
#include "matchinfo.h"
#include "dyhungarian.h"
#include <ctime>
using namespace std;
#define MAX_RS 10000
static const int oo = 0x3f3f3f3f;

vector<vector<int>> RoleSim_Plus(vector<vector<int>>& G, int k, double beta)
{
    int n = G.size();
    vector<vector<int>> H(n, vector<int>(n,MAX_RS));


    vector<vector<int>> W = get_W(G);
    pathnode* root = get_mst(W, G);
    dfs_init_v_common_index(root);

    for (int i = 0;i < k;++i) {
        H = Iterative(root, G, H, beta);
    }
    freeupspace(root);
    return H;
}

vector<vector<int>> get_W(vector<vector<int>>& G)
{
    int n = G.size();
    vector<vector<int>> W(n + 1, vector<int>(n + 1));
    for (int i = 0;i < n;++i) {
        W[i + 1][i + 1] = oo;
        for (int j = i + 1;j < n;++j) {
            W[i + 1][j + 1] = W[j + 1][i + 1] = distance(G[i], G[j]);
        }
    }
    W[0][0] = oo;
    for (int i = 1;i <= n;++i) {
        W[0][i] = W[i][0] = G[i - 1].size();
    }
    return W;
}

int distance(vector<int>& set1, vector<int>& set2)
{
    int n = set1.size();
    int m = set2.size();

    unordered_map<int, int> nums;
    for (int i = 0; i < n; ++i) {
        ++nums[set1[i]];
    }

    for (int i = 0; i < m; ++i) {
        ++nums[set2[i]];
    }

    int sum = 0;
    unordered_map<int, int>::iterator it;
    for (it = nums.begin(); it != nums.end(); ++it) {
        if (it->second == 2) {
            ++sum;
        }
    }
    return (n - sum) + (m - sum);
}



/*
 * 1.使用kruskal算法求出边集合
 * 2.使用边集合构造路径
 **/
pathnode* get_mst(vector<vector<int>>& W, vector<vector<int>>& graph) //W是转移代价矩阵(n+1 X n+1)
{
    int n = W.size();
    vector<vector<int>> G;
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < n;++j) {
            if (W[i][j] != oo)
                G.push_back({ i, j, W[i][j] });
        }
    }
    vector<vector<int>> I;
    I.push_back({});
    for (int i = 0;i < graph.size();++i) I.push_back(graph[i]);
    vector<vector<int>> g = kruskal2(n, G);
    pathnode* root = new pathnode(0, nullptr);
    queue<int> q;
    set<int> visited;
    unordered_map<int, pathnode*> pos;
    q.push(0);
    visited.insert(0);
    pos[0] = root;
    while (!q.empty()) {
        int v = q.front();
        pathnode* p = pos[v];
        q.pop();
        for (int i = 0;i < g[v].size();++i) {
            int x = g[v][i];
            if (visited.find(x) != visited.end()) continue;
            pathnode* next = new pathnode(x, pos[v]);
            p->nexts.push_back(next);
            visited.insert(x);
            pos[x] = next;
            q.push(x);

            set<int> fat, chi;
            for (int i = 0; i < I[v].size(); ++i) {
                fat.insert(I[v][i]);
            }
            for (int i = 0; i < I[x].size(); ++i) {
                chi.insert(I[x][i]);
            }
            set_intersection(fat.begin(), fat.end(), chi.begin(), chi.end(), inserter(next->commom, next->commom.begin()));
            set_difference(chi.begin(), chi.end(), fat.begin(), fat.end(), inserter(next->diff, next->diff.begin()));
        }
    }
    return root;
}

vector<vector<int>> Iterative(pathnode* root, vector<vector<int>>& G, vector<vector<int>>& H, double beta)
{
    int n = H.size();
    vector<vector<int>> H1(n, vector<int>(n));
    for (int u = 0;u < n;++u) {
        unordered_map<int, Mi*> forlapjv;
        dfs(u, root, G, forlapjv, H, H1, beta);
    }
    return H1;
}

int create_first_level(Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& I, vector<vector<int>>& H) {
    int vid = curent_node->val - 1;
    int Nu = I[uid].size(), Nv = I[vid].size();
    int mcdim = max(Nu, Nv);
    mcptr->dim = mcdim;
    mcptr->wei = new int* [mcdim];
    for (int i = 0; i < mcdim; ++i) mcptr->wei[i] = new int[mcdim];
    mcptr->rs = new int[mcdim];
    mcptr->cs = new int[mcdim];
    mcptr->rd = new int[mcdim];
    mcptr->cd = new int[mcdim];

    for (int j = 0; j < Nv; ++j) {
       int jj = curent_node->nodes_seq[j];
         //int jj = I[vid][j];
        for (int i = 0; i < Nu; ++i) {
            mcptr->wei[i][j] = H[I[uid][i]][jj];
        }
    }
    if (Nv > Nu) {
        for (int i = Nu; i < mcdim; ++i) {
            for (int j = 0; j < mcdim; ++j) {
                mcptr->wei[i][j] = 0;
            }
        }
    }
    if (Nv < Nu) {
        for (int j = Nv; j < mcdim; ++j) {
            for (int i = 0; i < mcdim; ++i) {
                mcptr->wei[i][j] = 0;
            }
        }
    }
    int idx = curent_node->common_size;
    return idx;
}

int create_v2(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& I, vector<vector<int>>& H) {
    int vid = curent_node->val - 1;
    int Nu = I[uid].size(), Nv = I[vid].size();
    int mcdim = max(mfptr->dim, Nv); // when Nv > mf.dim > Nu, mc.dim > mf.dim

    //assigment space for current node(mc)
    mcptr->dim = mcdim;
    mcptr->wei = new int* [mcdim];
    for (int i = 0; i < mcdim; ++i) mcptr->wei[i] = new int[mcdim];
    mcptr->rs = new int[mcdim];
    mcptr->cs = new int[mcdim];
    mcptr->rd = new int[mcdim];
    mcptr->cd = new int[mcdim];

    //Initialize the current node information
    int mfdim = mfptr->dim;
    for (int j = 0; j < curent_node->common_size; ++j) {
        for (int i = 0; i < Nu; ++i) {
            mcptr->wei[i][j] = mfptr->wei[i][curent_node->index[j]];
        }
        mcptr->cd[j] = mfptr->cd[curent_node->index[j]];
        mcptr->cs[j] = mfptr->cs[curent_node->index[j]];
        mcptr->rs[mcptr->cs[j]] = j;
    }

    for (int j = curent_node->common_size; j < curent_node->common_size + curent_node->diff_size; ++j) {
        for (int i = 0; i < Nu; ++i) {
            mcptr->wei[i][j] =  H[I[uid][i]][curent_node->nodes_seq[j]];
        }
        mcptr->cs[j] = -1;
    }
    for (int j = Nv; j < mcdim; ++j) {
        mcptr->cs[j] = -1;
        for (int i = 0; i < Nu; ++i) {
            mcptr->wei[i][j] = 0;
        }
    }

    for (int i = Nu; i < mcdim; ++i) {
        for (int j = 0; j < mcdim; ++j) {
            mcptr->wei[i][j] = 0;
        }
    }

    for (int i = 0; i < mfptr->dim; ++i) mcptr->rd[i] = mfptr->rd[i];

    int idx = curent_node->common_size;
    return idx;
}

void dfs(int u, pathnode* root, vector<vector<int>>& I, unordered_map<int, Mi*>& forlapjv, vector<vector<int>>& H, vector<vector<int>>& H1, double beta) {
    if (root == nullptr)
        return;
    Mi mc, * mcptr;
    int N = 0;
    mcptr = &mc;
    if (root->val == 0) {   //空集，起点
        mc.dim = 0;
        forlapjv[0] = mcptr;
    }
    else if (root->fa->val == 0) {
        int vid = root->val - 1;
        int Nu = I[u].size(), Nv = I[vid].size();
        N = max(Nu, Nv);
        create_first_level(mcptr, root, u, I, H);
        //double weight = Hungarian(mc);
        int weight = Hungarian(mc.dim, mc.wei, mc.rd, mc.cd, mc.rs, mc.cs);
        forlapjv[vid + 1] = mcptr;
        if (N == 0)
            H1[u][vid] = beta * MAX_RS;
        else
            H1[u][vid] = (1 - beta) * weight / mc.dim + beta * MAX_RS;
    }
    else {
        Mi mf = *(forlapjv[root->fa->val]);
        //==========================
        Mi* mfptr = forlapjv[root->fa->val];
        int f_dim = mfptr->dim;
        int idx = create_v2(mfptr, mcptr, root, u, I, H);

        //double weight = dymatch_col(idx, mf, mc);
        int weight = dymatch_col(mfptr->dim, mcptr->dim, idx, mcptr->wei, mcptr->rd, mcptr->cd, mcptr->rs, mcptr->cs);
        forlapjv[root->val] = &mc;
        if (root->val <= H1.size() && root->val > 0) {
            int Nu = I[u].size(), Nv = I[root->val - 1].size();
            int maxn = max(Nu, Nv);
            if (maxn == 0)
                H1[u][root->val - 1] = beta * MAX_RS;
            else
                H1[u][root->val - 1] = (1 - beta) * weight / maxn + beta * MAX_RS;
        }

    }
    for (int i = 0; i < root->nexts.size(); ++i) {
        dfs(u, root->nexts[i], I, forlapjv, H, H1, beta);
    }
    if (root->val != 0) {
        delete[] mc.cs;
        delete[] mc.rs;
        delete[] mc.rd;
        delete[] mc.cd;
        for (int i = 0; i < mc.dim; ++i) {
            delete[] mc.wei[i];
        }
        delete[] mc.wei;
        //delete[] mc.index;
    }
}

void dfs_init_v_common_index(pathnode* root) {
    if (root == nullptr)
        return;
    if (root->val == 0) {
        // empty set,do nothing
        root->common_size = root->commom.size();
        root->diff_size = root->diff.size();
    }
    else if (root->fa->val == 0) {//first level the in tree,record the sequence of node
        root->common_size = root->commom.size();
        root->diff_size = root->diff.size();
        root->nodes_seq = new int[root->diff.size()];
        set<int>::iterator it;
        int x = 0;
        for (it = root->diff.begin(); it != root->diff.end(); ++it) {
            root->nodes_seq[x++] = *it;
        }
    }
    else {
        root->common_size = root->commom.size();
        root->diff_size = root->diff.size();
        unordered_map<int, int> node_to_idx;
        for (int i = 0; i < root->fa->commom.size() + root->fa->diff.size(); ++i) {
            node_to_idx[root->fa->nodes_seq[i]] = i;
        }
        root->nodes_seq = new int[root->diff.size() + root->commom.size()];
        root->index = new int[root->commom.size()];
        set<int>::iterator it;
        int x = 0;
        for (it = root->commom.begin(); it != root->commom.end(); ++it) {
            root->nodes_seq[x] = *it;
            root->index[x] = node_to_idx[*it]; //get common elem in father's index
            ++x;
        }
        for (it = root->diff.begin(); it != root->diff.end(); ++it) {
            root->nodes_seq[x++] = *it;
        }
    }
    for (int i = 0; i < root->nexts.size(); ++i) {
        dfs_init_v_common_index(root->nexts[i]);
    }
}

void freeupspace(pathnode* root) {
    if (root == nullptr)
        return;
    for (int i = 0; i < root->nexts.size(); ++i) {
        freeupspace(root->nexts[i]);
    }
    delete root;
}