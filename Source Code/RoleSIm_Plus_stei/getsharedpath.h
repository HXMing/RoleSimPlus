#pragma once
#include <vector>
#include "approximate_steiner_tree.h"
using namespace std;


vector<vector<int>> generategraph(vector<vector<int>>& sets, vector<vector<int>>& aftersets);

void count_pruning_edges(vector<vector<int>>& graph);
vector<vector<int>> mapreduce(vector<vector<int>>& sets, vector<vector<int>>& aftersets);
vector<int> distance(vector<int>& set1, vector<int>& set2);
void purning(vector<vector<int>>& graph, vector<vector<int>>& W, int n);

//pathnode* getsharepath(vector<vector<int>>& G, vector<vector<int>>& nei);
pathnode* getsharepath(vector<vector<int>>& G);