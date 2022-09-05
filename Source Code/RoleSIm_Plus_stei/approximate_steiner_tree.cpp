#include "approximate_steiner_tree.h"
#include <unordered_map>
#include <queue>
#include <set>
using namespace std;
//https://github.com/jacky860226/steiner-tree-2-approximation
#define oo 10000
pathnode* approximate_steiner_tree(vector<vector<int>>& graph, vector<int>& ids, vector<vector<int>>& sets)
{
    size_t n, m;
    n = graph.size();
    steiner_tree::UndirectedGraph<double> G;
    G.setVertixNum(n);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (graph[i][j] != oo)
                G.addEdge(i, j, graph[i][j]);
        }
    }

    size_t terminalNum;

    terminalNum = ids.size();
    unordered_set<size_t> terminals;
    while (terminalNum--)
    {
        size_t terminal;
        terminal = ids[terminalNum];

        terminals.emplace(terminal);
    }
    steiner_tree::Solver<double> solver(G);

    auto res = solver.solve(terminals);

    unordered_map<int, vector<int>> nei;
    for (auto eid : *res)
    {
        nei[G.getEdge(eid).v1].push_back(G.getEdge(eid).v2);
        nei[G.getEdge(eid).v2].push_back(G.getEdge(eid).v1);
    }

    return getsharepath(nei, sets);
}

pathnode* getsharepath(unordered_map<int, vector<int>>& nei, vector<vector<int>>& sets)
{
    queue<int> q;
    set<int> vised;
    q.push(0);
    vised.insert(0);
    pathnode* root = new pathnode(0, NULL);
    unordered_map<int, pathnode*> pos;
    pos[0] = root;
    while (!q.empty()) {
        int v = q.front();
        pathnode* p = pos[v];
        q.pop();
        for (int i = 0; i < nei[v].size(); ++i) {
            int x = nei[v][i];
            if (vised.find(x) == vised.end()) {
                pathnode* next = new pathnode(x, pos[v]);
                p->nexts.push_back(next);
                vised.insert(x);
                pos[x] = next;
                q.push(x);

                set<int> fat, chi;
                for (int i = 0; i < sets[v].size(); ++i) {
                    fat.insert(sets[v][i]);
                }
                for (int i = 0; i < sets[x].size(); ++i) {
                    chi.insert(sets[x][i]);
                }

                set_intersection(fat.begin(), fat.end(), chi.begin(), chi.end(), inserter(next->commom, next->commom.begin()));
                set_difference(chi.begin(), chi.end(), fat.begin(), fat.end(), inserter(next->diff, next->diff.begin()));
            }
        }
    }
    return root;
}