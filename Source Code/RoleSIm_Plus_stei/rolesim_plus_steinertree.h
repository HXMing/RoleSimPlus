#pragma once
#include "matchinfo.h"
#include <vector>
#include "pathtree.h"
#include <unordered_map>
using namespace std;

vector<vector<int>> Rolesim_plus_stei(vector<vector<int>>& G, double beta, int k);

void Iterative(pathnode* root, vector<vector<int>>& nei, vector<vector<int>>& Ha, vector<vector<int>>& Hb, double beta);
//void dfs(int u, pathnode* root, vector<vector<int>>& nei, unordered_map<int, Mi*>& forlapjv, vector<vector<double>>& Ha, vector<vector<double>>& Hb, double beta);
void freeupspace(pathnode* root);

void dfs2(pathnode* root);

int create(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<int>>& Ha);
int create_v2(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<int>>& Ha);
pair<int, int> create_dyhungarian_col_space(Mi* mfptr, Mi* mcptr, pathnode* curent_node, vector<vector<int>>& nei, int uid, vector<vector<int>>& H);
int create_first_level(Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<int>>& Ha);
void dfs_init_v_common_index(pathnode* root);