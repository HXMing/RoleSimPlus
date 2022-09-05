#include "getsharedpath.h"
#include <string>
#include <ctime>
#include <iostream>
using namespace std;
#define oo 10000
pathnode* getsharepath(vector<vector<int>>& G) {
    //count_pruning_edges(G);
    int n = G.size();
    vector<vector<int>> sets;
    sets.push_back({});
    for (int i = 0; i < n; ++i) {
        sets.push_back(G[i]);
    }
    vector<vector<int>> aftermapsets;

    vector<vector<int>> W = generategraph(sets, aftermapsets);
    purning(aftermapsets, W, n + 1);

    vector<int> ids;
    for (int i = 0; i <= n; ++i)
        ids.push_back(i);
    pathnode* root = approximate_steiner_tree(W, ids, aftermapsets);
    return root;
}

vector<vector<int>> generategraph(vector<vector<int>>& sets, vector<vector<int>>& aftersets)
{
    vector<vector<int>> newsets = mapreduce(sets, aftersets);
    int n = newsets.size();
    vector<vector<int>> G(n, vector<int>(n, 0));
    for (int i = 0; i < n; ++i) {
        G[i][i] = oo;
        for (int j = i + 1; j < n; ++j) {
            vector<int> d = distance(newsets[i], newsets[j]);
            G[i][j] = d[0], G[j][i] = d[1];
        }
    }
    return G;
}

vector<vector<int>> mapreduce(vector<vector<int>>& sets, vector<vector<int>>& aftersets)
{
    vector<vector<int>> newsets;
    unordered_map<int, string> mapped;
    for (int i = 0; i < sets.size(); ++i) {
        string tag = std::to_string(i);
        tag += '-';
        for (int j = 0; j < sets[i].size(); ++j) {
            mapped[sets[i][j]] += tag;
        }
    }
    unordered_map<string, vector<int>> remapped; //reverse
    unordered_map<int, string>::iterator it;
    for (it = mapped.begin(); it != mapped.end(); ++it) {
        if (it->second.size() != 0) {
            remapped[it->second].push_back(it->first);
        }
    }
    unordered_map<string, vector<int>>::iterator it2;
    for (it2 = remapped.begin(); it2 != remapped.end(); ++it2) {
        string s = it2->first;
        int cnt = 0;
        for (int i = 0; i < s.size(); ++i) {
            if (s[i] == '-') {
                ++cnt;
            }
        }
        if (cnt > 1 && cnt <= 2)
            newsets.push_back(it2->second);
    }
    vector<vector<int>> resset = sets;
    for (int i = 0; i < newsets.size(); ++i) {
        resset.push_back(newsets[i]);
    }
    aftersets = resset;
    return resset;
}

vector<int> distance(vector<int>& set1, vector<int>& set2)
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
    vector<int> ans(2, 0);
    if (sum < n && n <= min(2 * sum, m)) {
        ans[0] = n + m - 2 * sum;
    }
    else if (sum == n) {
        ans[0] = m - n;
    }
    else {
        ans[0] = oo;
    }

    if (sum < m && m <= min(2 * sum, n)) {
        ans[1] = n + m - 2 * sum;
    }
    else if (sum == m) {
        ans[1] = n - m;
    }
    else {
        ans[1] = oo;
    }
    return ans;
}

void purning(vector<vector<int>>& graph, vector<vector<int>>& W, int L)
{
    int n = graph.size();
    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (W[i][j] == oo) continue;
            if (W[i][j] >= graph[j].size() && graph[i].size() > 0) {
                W[i][j] = oo;
            }
        }
    }
    vector<int> minval(n, oo);
    for (int j = 0; j < n; ++j) {
        for (int i = L; i < n; ++i) {
            if (W[i][j] != oo && W[i][j] < minval[j]) {
                minval[j] = W[i][j];
            }
        }
    }
    for (int j = 0; j < n; ++j) {
        for (int i = L; i < n; ++i) {
            if (W[i][j] == oo) continue;
            if (minval[j] < W[i][j]) W[i][j] = oo;
        }
    }
}

void count_pruning_edges(vector<vector<int>>& graph) {
    int n = graph.size();
    vector<vector<int>> sets;
    sets.push_back({});
    for (int i = 0; i < n; ++i) {
        sets.push_back(graph[i]);
    }
    vector<vector<int>> aftermapsets;

    vector<vector<int>> W = generategraph(sets, aftermapsets);
    int L = n + 1;
    n = aftermapsets.size();

    int R1 = 0, R2 = 0, R3 = 0;
    for (int j = 0; j < L; ++j) {
        for (int i = 1; i < n; ++i) {
            if (W[i][j] == oo)continue;
            if (W[i][j] <= W[0][j]) {
                W[0][j] = oo;
                ++R2;
                break;
            }
        }
    }

    for (int i = 1; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (W[i][j] == oo) continue;
            if (W[i][j] >= aftermapsets[j].size() && aftermapsets[i].size() > 0) {
                W[i][j] = oo;
                ++R1;
            }
        }
    }

    vector<int> minval(n, oo);
    for (int j = 0; j < n; ++j) {
        for (int i = L; i < n; ++i) {
            if (W[i][j] != oo && W[i][j] < minval[j]) {
                minval[j] = W[i][j];
            }
        }
    }
    for (int j = 0; j < n; ++j) {
        for (int i = L; i < n; ++i) {
            if (W[i][j] == oo) continue;
            if (minval[j] < W[i][j]) W[i][j] = oo, ++R3;
        }
    }
    cout << "the number of pruned edges\n";
    cout << "R1 = " << R1 << endl;
    cout << "R2 = " << R2 << endl;
    cout << "R3 = " << R3 << endl;
}