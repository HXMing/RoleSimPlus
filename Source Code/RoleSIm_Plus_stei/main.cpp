#include <iostream>
#include <vector>
#include <fstream>
#include <unordered_map>
#include <set>
#include "rolesim_plus_steinertree.h"
#include <ctime>
#include "getsharedpath.h"
using namespace std;

vector<vector<int>> Load_dataset() {
	int nodes, edges;
	ofstream ofs;
	FILE* fn1;
	//fopen_s(&fn1, "D:\\datasets\\CA-GrQc.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\Wiki-Vote.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\p2p-Gnutella08.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\CA-HepTh.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\p2p-Gnutella09.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\soc-sign-bitcoinotc.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\facebook_combined.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\p2p-Gnutella06.txt", "r");
	fopen_s(&fn1, "D:\\datasets\\soc-sign-bitcoinalpha.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\email-Eu-core.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\CollegeMsg.txt", "r");
	//fopen_s(&fn1, "D:\\datasets\\as20000102.txt", "r");
	unordered_map<int, set<int>> mp;
	unordered_map<int, int> id;
	int from, to,nouse;
	int n = 0;
	int k = 0;
	while (fscanf_s(fn1, "%d\t%d\t\%d\n", &from, &to,&nouse) != EOF) {
		if (id.find(from) == id.end()) {
			id[from] = k++;
		}
		if (id.find(to) == id.end()) {
			id[to] = k++;
		}
		if (from != to)
			mp[id[from]].insert(id[to]);
		n = max(n, id[from]);
		n = max(n, id[to]);
	}

	unordered_map<int, set<int>>::iterator it;
	set<int>::iterator it2;
	vector<vector<int>> G(n + 1, vector<int>(n + 1, 0));
	for (it = mp.begin();it != mp.end();++it) {
		from = it->first;
		set<int> in_nodes = it->second;
		for (it2 = in_nodes.begin();it2 != in_nodes.end();++it2) {
			to = *it2;
			if (from != to)
				G[from][to] = 1;
		}
	}
	fclose(fn1);
	return G;
}

vector<vector<int>> Load_graph(vector<vector<int>>& g) {
	int n = g.size();
	vector<vector<int>> I(n);
	for (int i = 0;i < g.size();++i) {
		for (int j = 0;j < g.size();++j) {
			if (i == j) continue;
			if (g[i][j] == 1)
				I[j].push_back(i);
		}
	}
	return I;
}

int main()
{
	vector<vector<int>> G1 = {
		{0,0,1,0,0,1,0},//1
		{0,0,0,1,0,0,1},//2
		{1,1,0,1,1,0,1},//3
		{1,0,1,0,0,1,1},//4
		{1,0,0,0,0,0,1},//5
		{1,0,1,1,0,0,1},//6
		{1,0,0,0,0,0,0} //7
	};
	vector<vector<int>> G = Load_dataset();
	vector<vector<int>> I = Load_graph(G);		//邻接表的形式储存图
	//count_pruning_edges(I);
	cout << "start\n";
	clock_t start, end;
	start = clock();
	vector<vector<int>> S = Rolesim_plus_stei(I, 0.8, 5);
	end = clock();
	cout << "time = " << (end - start) / 1000.0 << endl;
	return 0;
}