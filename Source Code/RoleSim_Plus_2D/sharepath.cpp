#include "sharepath.h"
#include <iostream>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <set>

using namespace std;

#define oo 0x3f3f3f3f

vector<vector<int>> get_W(vector<vector<int>>& G) 
//input : G --- 图，邻接表形式存储
//outout: 代价转移矩阵，id偏移1，0为起点
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

pathnode* get_mst(vector<vector<int>>& g,vector<vector<int>>& W)
{
    
    int n = g.size() + 1;
    vector<vector<int>> id_to_pos(n * n);

    vector<edge> edges;
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < g.size();++j) {
            int u = g[j][0], v = g[j][1], w = W[u][v];
            edges.push_back(edge(i, u, i, v, w));
            edges.push_back(edge(u, i, v, i, w));
        }
    }
    vector<edge> g2 = kruskal2(n, edges);
    
    int id1, id2;
    for (int i = 0;i < g2.size();++i) {
        id1 = g2[i].r1 * n + g2[i].c1, id2 = g2[i].r2 * n + g2[i].c2;
        id_to_pos[id1].push_back(i);
        id_to_pos[id2].push_back(i);
    }
    queue<int> q;
    vector<pathnode*> pos(n * n);
    vector<bool> visited(n * n, false);
    pathnode* root = new pathnode(0, 0, nullptr);
    q.push(0);
    visited[0] = true;
    pos[0] = root;
    while (!q.empty()) {
        int id = q.front();
        pathnode* p = pos[id];
        q.pop();
        for (int i = 0;i < id_to_pos[id].size();++i) {
            id1 = g2[id_to_pos[id][i]].r1 * n + g2[id_to_pos[id][i]].c1;
            id2 = g2[id_to_pos[id][i]].r2 * n + g2[id_to_pos[id][i]].c2;
            if (id1 == id) id1 = id2;
            if (visited[id1]) continue;
            pathnode* next = new pathnode(id1 / n, id1 % n, p);
            p->nexts.push_back(next);
            visited[id1] = true;
            pos[id1] = next;
            q.push(id1);
        }
    }
    return root;
}

vector<vector<int>> get_edges(vector<vector<int>>& W)
{
    int n = W.size();
    vector<vector<int>> G;
    for (int i = 0;i < n;++i) {
        for (int j = 0;j < n;++j) {
            if (W[i][j] != oo)
                G.push_back({ i, j, W[i][j] });
        }
    }
    vector<vector<int>> g = kruskal(n, G);
    return g;
}

void common_elem_in_edge(unordered_map<int,int*>& common_elem,unordered_map<int,int*>& diff_elem,vector<vector<int>>& edges, vector<vector<int>>& I)
{   
    int n = edges.size() + 1;
    vector<set<int>> elems(n - 1);
    for (int i = 0;i < I.size();++i) {
        for (int j = 0;j < I[i].size();++j) {
            elems[i].insert(I[i][j]);
        }
    }
    for (int i = 0;i < edges.size();++i) {
        int u = edges[i][0], v = edges[i][1];
        if (u == 0 && v == 0) {
            int* common = new int[1];
            common[0] = 0;
            common_elem[0] = common;
            int* diff = new int[1];
            diff[0] = 0;
            diff_elem[0] = diff;
        }
        else if (u == 0 && v != 0) {
            int* common = new int[1];
            common[0] = 0;
            common_elem[u * n + v] = common;
            common_elem[v * n + u] = common;
            int* diff_u = new int[1];
            diff_u[0] = 0;
            diff_elem[u * n + v] = diff_u;
            int* diff_v = new int[elems[v - 1].size() + 1];
            diff_v[0] = elems[v - 1].size();
            set<int>::iterator it;
            int idx = 1;
            for (it = elems[v - 1].begin();it != elems[v - 1].end();++it) {
                diff_v[idx++] = *it;
            }
            diff_elem[v * n + u] = diff_v;
        }
        else if (u != 0 && v == 0) {
            int* common = new int[1];
            common[0] = 0;
            common_elem[u * n + v] = common;
            common_elem[v * n + u] = common;
            int* diff_u = new int[elems[u - 1].size() + 1];
            diff_u[0] = elems[u - 1].size();
            set<int>::iterator it;
            int idx = 1;
            for (it = elems[u - 1].begin();it != elems[u - 1].end();++it) {
                diff_u[idx++] = *it;
            }
            diff_elem[u * n + v] = diff_u;
            int* diff_v = new int[1];
            diff_v[0] = 1;
            diff_elem[v * n + u] = diff_v;
        }
        else {
            set<int> common, df_u_v, df_v_u;
            set_intersection(elems[u - 1].begin(), elems[u - 1].end(), elems[v - 1].begin(), elems[v - 1].end(), inserter(common, common.begin()));
            int* com = new int[common.size() + 1];
            com[0] = common.size();
            set<int>::iterator it;
            int idx;
            idx = 1;
            for (it = common.begin();it != common.end();++it) {
                com[idx++] = *it;
            }
            common_elem[u * n + v] = com;
            common_elem[v * n + u] = com;
            set_difference(elems[u - 1].begin(), elems[u - 1].end(), elems[v - 1].begin(), elems[v - 1].end(), inserter(df_u_v, df_u_v.begin()));
            int* diff_u = new int[df_u_v.size() + 1];
            diff_u[0] = df_u_v.size();
            idx = 1;
            for (it = df_u_v.begin();it != df_u_v.end();++it) {
                diff_u[idx++] = *it;
            }
            diff_elem[u * n + v] = diff_u;

            set_difference(elems[v - 1].begin(), elems[v - 1].end(), elems[u - 1].begin(), elems[u - 1].end(), inserter(df_v_u, df_v_u.begin()));
            int* diff_v = new int[df_v_u.size() + 1];
            diff_v[0] = df_v_u.size();
            idx = 1;
            for (it = df_v_u.begin();it != df_v_u.end();++it) {
                diff_v[idx++] = *it;
            }
            diff_elem[v * n + u] = diff_v;
        }
    }
}

void dfs_create_common_index(int n, pathnode* root, unordered_map<int, int*>& common_elem, unordered_map<int, int*>& diff_elem)
{
    if (root == nullptr)
        return;
    int uid = -1, vid = -1, fuid = -1, fvid = -1;
    uid = root->row, vid = root->col;
    if (root->fa != nullptr)
        fuid = root->fa->row, fvid = root->fa->col;
    if (uid == 0 && vid == 0) {
        root->common_size = 0;
        root->diff_size = 0;
        root->col_node_seq = nullptr;
        root->row_node_seq = nullptr;
        root->col_common_elem_index = nullptr;
        root->row_common_elem_index = nullptr;
    }
    else if (fuid == uid && fvid == 0) {
        int index = vid * n + fvid;
        root->common_size = common_elem[index][0];
        root->diff_size = diff_elem[index][0];
        root->row_node_seq = root->fa->row_node_seq;
        //root->row_common_elem_index = root->fa->row_common_elem_index;
        root->row_common_elem_index = nullptr;
        root->col_node_seq = new short[root->common_size + root->diff_size];
        root->col_common_elem_index = nullptr;
        for (int i = 1;i <= diff_elem[index][0];++i) {
            root->col_node_seq[i - 1] = (short)diff_elem[index][i];
        }
    }
    else if (fuid == uid && fvid != 0) {
        int index = vid * n + fvid;
        root->common_size = common_elem[index][0];
        root->diff_size = diff_elem[index][0];
        root->row_node_seq = root->fa->row_node_seq;
        //root->row_common_elem_index = root->fa->row_common_elem_index;
        root->row_common_elem_index = nullptr;
        root->col_node_seq = new short[root->common_size + root->diff_size];
        root->col_common_elem_index = new short[root->common_size];
        unordered_map<int, int> node_to_idx;
        for (int i = 0;i < diff_elem[fvid * n + vid][0] + common_elem[fvid * n + vid][0];++i) {
            node_to_idx[root->fa->col_node_seq[i]] = i;
        }
        int x = 0;
        for (int i = 0;i < root->common_size;++i) {
            root->col_node_seq[x] = common_elem[index][i + 1];
            root->col_common_elem_index[x] = node_to_idx[root->col_node_seq[x]];
            ++x;
        }
        for (int i = 0;i < root->diff_size;++i) {
            root->col_node_seq[x] = diff_elem[index][i + 1];
            ++x;
        }
    }
    else if (fvid == vid && fuid == 0) {
        int index = uid * n + fuid;
        root->common_size = common_elem[index][0];
        root->diff_size = diff_elem[index][0];
        root->col_node_seq = root->fa->col_node_seq;
        //root->col_common_elem_index = root->fa->col_common_elem_index;
        root->col_common_elem_index = nullptr;
        root->row_node_seq=new short[root->common_size + root->diff_size];
        root->row_common_elem_index = nullptr;
        for (int i = 1;i <= diff_elem[index][0];++i) {
            root->row_node_seq[i - 1] = (short)diff_elem[index][i];
        }
    }
    else if (fvid == vid && fuid != 0) {
        int index = uid * n + fuid;
        root->common_size = common_elem[index][0];
        root->diff_size = diff_elem[index][0];
        root->col_node_seq = root->fa->col_node_seq;
        //root->col_common_elem_index = root->fa->col_common_elem_index;
        root->col_common_elem_index = nullptr;
        root->row_node_seq = new short[root->common_size + root->diff_size];
        root->row_common_elem_index = new short[root->common_size];
        unordered_map<int, int> node_to_idx;
        for (int i = 0;i < diff_elem[fuid * n + uid][0] + common_elem[fuid * n + uid][0];++i) {
            node_to_idx[root->fa->row_node_seq[i]] = i;
        }
        int x = 0;
        for (int i = 0;i < root->common_size;++i) {
            root->row_node_seq[x] = common_elem[index][i + 1];
            root->row_common_elem_index[x] = node_to_idx[root->row_node_seq[x]];
            ++x;
        }
        for (int i = 0;i < root->diff_size;++i) {
            root->row_node_seq[x] = diff_elem[index][i + 1];
            ++x;
        }
    }
    for (int i = 0;i < root->nexts.size();++i) {
        dfs_create_common_index(n, root->nexts[i], common_elem, diff_elem);
    }
}


pathnode* get_shared_path(vector<vector<int>>& I)
{
    vector<vector<int>> W = get_W(I);
    vector<vector<int>> g = get_edges(W);
    pathnode* root = get_mst(g, W);
    unordered_map<int, int*> common_elems, diff_elems;
    common_elem_in_edge(common_elems, diff_elems, g, I);
    dfs_create_common_index(W.size(), root, common_elems, diff_elems);
    unordered_map<int, int*>::iterator it;
    unordered_map<int*, int> remap;
    for (it = common_elems.begin();it != common_elems.end();++it) {
        remap[it->second] = it->first;
    }
    for (unordered_map<int*, int>::iterator it = remap.begin();it != remap.end();++it) {
        delete[] it->first;
    }
    for (it = diff_elems.begin();it != diff_elems.end();++it) {
        delete[] it->second;
    }
    return root;
}

void freeupspace(pathnode* root) {
    if (root == nullptr)
        return;
    int uid = -1, vid = -1, fuid = -1, fvid = -1;
    uid = root->row, vid = root->col;
    if (root->fa != nullptr)
        fuid = root->fa->row, fvid = root->fa->col;
    if(uid==0&&vid==0){}
    else if (fuid == uid && fvid == 0) {
        delete[] root->col_node_seq;
    }
    else if (fuid == uid && fvid != 0) {
        delete[] root->col_node_seq;
        delete[] root->col_common_elem_index;
    }
    else if (fvid == vid && fuid == 0) {
        delete[] root->row_node_seq;
    }
    else if (fvid == 0 && fuid != 0) {
        delete[] root->row_node_seq;
        delete[] root->row_common_elem_index;
    }
    for (int i = 0;i < root->nexts.size();++i) {
        freeupspace(root->nexts[i]);
    }
    delete root;
}



//=============== mst ===================
void init(vector<int>& fa) {
    for (int i = 0;i < fa.size();++i)
        fa[i] = i;
}

int find(int i, vector<int>& father) {
    if (father[i] != i)
        father[i] = find(father[i], father);
    return father[i];
}

void join(int i, int j, vector<int>& father) {
    int fi = find(i, father), fj = find(j, father);
    father[fi] = fj;
}

vector<vector<int>> kruskal(int n, vector<vector<int>>& edges)
{
    int es = edges.size();
    vector<int> father(n);
    init(father);
    vector<vector<int>> res;
    sort(edges.begin(), edges.end(), [](vector<int>& a, vector<int>& b)->bool {return a[2] < b[2];});
    int k = n - 1;
    for (int i = 0;i < edges.size();++i) {
        int a = edges[i][0], b = edges[i][1];
        int fa = find(a, father), fb = find(b, father);
        if (fa == fb) continue;
        join(a, b, father);
        res.push_back({ a,b });
        k -= 1;
        if (k == 0)break;
    }
    return res;
}

vector<edge> kruskal2(int n, vector<edge>& edges)
{
    sort(edges.begin(), edges.end(), [](edge& a, edge& b)->bool {return a.w < b.w;});
    int es = edges.size();
    vector<int> father(n * n);
    init(father);
    vector<edge> res;
    int k = n * n - 1;
    for (int i = 0;i < edges.size();++i) {
        int a = edges[i].r1 * n + edges[i].c1, b = edges[i].r2 * n + edges[i].c2;
        int fa = find(a, father), fb = find(b, father);
        if (fa == fb) continue;
        join(a, b, father);
        res.push_back(edges[i]);
        k -= 1;
        if (k == 0)break;
    }
    return res;
}
