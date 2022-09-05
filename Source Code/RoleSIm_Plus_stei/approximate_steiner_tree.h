#pragma once
#include <vector>
#include "Solver.h"
#include "Structure.h"
#include <iostream>
#include <unordered_set>
#include "pathtree.h"
#include <unordered_map>
using namespace std;


pathnode* approximate_steiner_tree(vector<vector<int>>& graph, vector<int>& ids, vector<vector<int>>& sets);
pathnode* getsharepath(unordered_map<int, vector<int>>& nei, vector<vector<int>>& sets);