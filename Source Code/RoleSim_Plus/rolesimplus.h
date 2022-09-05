#pragma once
#include <vector>
#include "pathnode.h"
#include <unordered_map>
#include "matchinfo.h"
using namespace std;

vector<vector<int>> RoleSim_Plus(vector<vector<int>>& G, int k, double beta);

vector<vector<int>> get_W(vector<vector<int>>& G);

int distance(vector<int>& set1, vector<int>& set2);

pathnode* get_mst(vector<vector<int>>& W, vector<vector<int>>& graph);

vector<vector<int>> Iterative(pathnode* root, vector<vector<int>>& G, vector<vector<int>>& H, double beta);

void dfs(int u, pathnode* root, vector<vector<int>>& I, unordered_map<int, Mi*>& forlapjv, vector<vector<int>>& H, vector<vector<int>>& H1, double beta);

void dfs_init_v_common_index(pathnode* root);

int create_first_level(Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& I, vector<vector<int>>& H);

int create_v2(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& I, vector<vector<int>>& H);

void freeupspace(pathnode* root);