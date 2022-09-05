#pragma once
#include <vector>
#include "pathnode.h"
#include <unordered_map>
using namespace std;

typedef struct edge {
	int r1, c1, r2, c2, w;
	edge(int a, int b, int c, int d, int e):r1(a), c1(b), r2(c), c2(d), w(e) {}
}edge;

vector<vector<int>> get_W(vector<vector<int>>& G);
int distance(vector<int>& set1, vector<int>& set2);
pathnode* get_mst(vector<vector<int>>& edges,vector<vector<int>>& W);
vector<vector<int>> get_edges(vector<vector<int>>& W);
void common_elem_in_edge(unordered_map<int, int*>& common_elem, unordered_map<int, int*>& diff_elem, vector<vector<int>>& edges, vector<vector<int>>& I);
void dfs_create_common_index(int n, pathnode* root, unordered_map<int, int*>& common_elem, unordered_map<int, int*>& diff_elem);
pathnode* get_shared_path(vector<vector<int>>& I);
void freeupspace(pathnode* root);

//========== mst ===========
void init(vector<int>& fa);
int find(int i, vector<int>& fa);
void join(int i, int j, vector<int>& fa);
vector<vector<int>> kruskal(int n, vector<vector<int>>& edges);
vector<edge> kruskal2(int n, vector<edge>& edges);