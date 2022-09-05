#pragma once
#include <vector>
using namespace std;

vector<vector<int>> OriginalRoleSim(vector<vector<int>>& G, int k, double beta);
vector<vector<int>> IterateRoleSim(vector<vector<int>>& H, vector<vector<int>>& I, double beta);