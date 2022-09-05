#pragma once
#include <vector>
using namespace std;

typedef struct edge {
	int start, end;
	int weight;
	edge(int s,int e,int w):start(s),end(e),weight(w){}
}edge;

vector<edge> prim(vector<vector<int>>& weight);
vector<edge> kruskal(int n, vector<vector<int>>& edges);
bool cmp(vector<int>& a, vector<int>& b);
vector<vector<int>> kruskal2(int n, vector<vector<int>>& edges);

int find(int x);
void join(int x, int y);
void init(int n);