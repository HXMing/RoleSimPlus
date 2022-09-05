#pragma once
#include <vector>
#include "matchinfo.h"
#include "pathnode.h"
#include <unordered_map>
using namespace std;

vector<vector<int>> RoleSim_Plus2(vector<vector<int>>& I, int k, double beta);
vector<vector<int>> Iterative(pathnode* root, vector<vector<int>>& I, vector<vector<int>>& H, double beta);
void dfs(pathnode* root, vector<vector<int>>& I, unordered_map<int, Mi*>& forlapjv, vector<vector<int>>& H, vector<vector<int>>& H1, double beta);
void create_hungarian_space(Mi* mcptr, pathnode* curent_node, int Nu, int Nv, vector<vector<int>>& H);
pair<int, int> create_dyhungarian_col_space(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int Nu, int Nv, int Nfv, vector<vector<int>>& H);
pair<int, int> create_dyhungarian_row_space(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int Nu, int Nv, int Nfu, vector<vector<int>>& H);