#include "rolesim_plus_steinertree.h"
#include "getsharedpath.h"
#include "dyhungarian.h"
#include <ctime>
#include <fstream>
#include <string>
using namespace std;
#define MAX_RS 10000

vector<vector<int>> Rolesim_plus_stei(vector<vector<int>>& G, double beta, int k)
{
	
	int n = G.size();
	vector<vector<int>> H1(n, vector<int>(n, MAX_RS));
	vector<vector<int>> H2(n, vector<int>(n, MAX_RS));

	
	clock_t start, end;
	//start = clock();
	pathnode* root = getsharepath(G);
	dfs_init_v_common_index(root);
	//end = clock();
	//cout << "Phrase 1 time = " << (end - start) / 1000.0 << endl;

	for (int i = 0; i < k; ++i) {
		if (i % 2 == 0) {
			Iterative(root, G, H1, H2, beta);
		}
		else {
			Iterative(root, G, H2, H1, beta);
		}
		//end = clock();
		//cout << "on iterate " << i + 1 << " time = " << (end - start) / 1000.0 << endl;
	}
	freeupspace(root);
	if (k % 2 == 1) {
		return H2;
	}
	return H1;
}

void dfs_init_v_common_index(pathnode* root) {
	if (root == nullptr)
		return;
	if (root->val == 0) {
		// empty set,do nothing
		root->common_size = root->commom.size();
		root->diff_size = root->diff.size();
	}
	else if (root->fa->val == 0) {//first level in the tree,record the sequence of node
		root->common_size = root->commom.size();
		root->diff_size = root->diff.size();
		root->nodes_seq = new int[root->diff.size()];
		set<int>::iterator it;
		int x = 0;
		for (it = root->diff.begin(); it != root->diff.end(); ++it) {
			root->nodes_seq[x++] = *it;
		}
	}
	else {
		root->common_size = root->commom.size();
		root->diff_size = root->diff.size();
		unordered_map<int, int> node_to_idx;
		for (int i = 0; i < root->fa->commom.size() + root->fa->diff.size(); ++i) {
			node_to_idx[root->fa->nodes_seq[i]] = i;
		}
		root->nodes_seq = new int[root->diff.size() + root->commom.size()];
		root->index = new int[root->commom.size()];
		set<int>::iterator it;
		int x = 0;
		for (it = root->commom.begin(); it != root->commom.end(); ++it) {
			root->nodes_seq[x] = *it;
			root->index[x] = node_to_idx[*it]; //获取父节点中的公共元素
			++x;
		}
		for (it = root->diff.begin(); it != root->diff.end(); ++it) {
			root->nodes_seq[x++] = *it;
		}
	}
	for (int i = 0; i < root->nexts.size(); ++i) {
		dfs_init_v_common_index(root->nexts[i]);
	}
}

int create_first_level(Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<int>>& Ha) {
	int vid = curent_node->val - 1;
	int Nu = nei[uid].size(), Nv = nei[vid].size();
	int mcdim = max(Nu, Nv);
	mcptr->dim = mcdim;
	mcptr->cost = new int* [mcdim];
	for (int i = 0; i < mcdim; ++i) mcptr->cost[i] = new int[mcdim];
	mcptr->rowsol = new int[mcdim];
	mcptr->colsol = new int[mcdim];
	mcptr->u = new int[mcdim];
	mcptr->v = new int[mcdim];

	for (int j = 0; j < Nv; ++j) {
		int jj = curent_node->nodes_seq[j];
		for (int i = 0; i < Nu; ++i) {
			mcptr->cost[i][j] = Ha[nei[uid][i]][jj];
		}
	}
	if (Nv > Nu) {
		for (int i = Nu; i < mcdim; ++i) {
			for (int j = 0; j < mcdim; ++j) {
				mcptr->cost[i][j] = 0;
			}
		}
	}
	if (Nv < Nu) {
		for (int j = Nv; j < mcdim; ++j) {
			for (int i = 0; i < mcdim; ++i) {
				mcptr->cost[i][j] = 0;
			}
		}
	}
	int idx = curent_node->common_size;
	return idx;
}

int create_v2(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<int>>& Ha) {
	int vid = curent_node->val - 1;
	int Nu = nei[uid].size(), Nv = nei[vid].size();
	int mcdim = max(mfptr->dim, Nv); // when Nv > mf.dim > Nu, mc.dim > mf.dim

	//assigment space for current node(mc)
	mcptr->dim = mcdim;
	mcptr->cost = new int* [mcdim];
	for (int i = 0; i < mcdim; ++i) mcptr->cost[i] = new int[mcdim];
	mcptr->rowsol = new int[mcdim];
	mcptr->colsol = new int[mcdim];
	mcptr->u = new int[mcdim];
	mcptr->v = new int[mcdim];

	//Initialize the current node information
	int mfdim = mfptr->dim;
	for (int j = 0; j < curent_node->common_size; ++j) {
		for (int i = 0; i < Nu; ++i) {
			mcptr->cost[i][j] = mfptr->cost[i][curent_node->index[j]];
		}
		mcptr->v[j] = mfptr->v[curent_node->index[j]];
		mcptr->colsol[j] = mfptr->colsol[curent_node->index[j]];
		mcptr->rowsol[mcptr->colsol[j]] = j;
	}

	for (int j = curent_node->common_size; j < curent_node->common_size + curent_node->diff_size; ++j) {
		for (int i = 0; i < Nu; ++i) {
			mcptr->cost[i][j] = Ha[nei[uid][i]][curent_node->nodes_seq[j]];
		}
		mcptr->colsol[j] = -1;
	}
	for (int j = Nv; j < mcdim; ++j) {
		mcptr->colsol[j] = -1;
		for (int i = 0; i < Nu; ++i) {
			mcptr->cost[i][j] = 0;
		}
	}

	for (int i = Nu; i < mcdim; ++i) {
		for (int j = 0; j < mcdim; ++j) {
			mcptr->cost[i][j] = 0;
		}
	}

	for (int i = 0; i < mfptr->dim; ++i) mcptr->u[i] = mfptr->u[i];

	int idx = curent_node->common_size;
	return idx;
}

pair<int, int> create_dyhungarian_col_space(Mi* mfptr, Mi* mcptr, pathnode* curent_node, vector<vector<int>>& nei, int uid, vector<vector<int>>& H)
{
	int vid = curent_node->val - 1;
	int Nu = nei[uid].size(), Nv = nei[vid].size(), Nfv = nei[curent_node->fa->val - 1].size();
	int mcdim = max(mfptr->dim, Nv); // when Nv > mf.dim > Nu, mc.dim > mf.dim
	//assigment space for current node(mc)
	mcptr->dim = mcdim;
	mcptr->cost = new int* [mcdim];
	for (int i = 0; i < mcdim; ++i) mcptr->cost[i] = new int[mcdim];
	mcptr->rowsol = new int[mcdim];
	mcptr->colsol = new int[mcdim];
	mcptr->u = new int[mcdim];
	mcptr->v = new int[mcdim];

	//Initialize the current node information
	int mfdim = mfptr->dim;
	for (int j = 0;j < curent_node->common_size;++j) {
		for (int i = 0;i < Nu;++i) {
			mcptr->cost[i][j] = mfptr->cost[i][curent_node->index[j]];
		}
	}
	for (int j = curent_node->common_size;j < Nv;++j) {
		for (int i = 0;i < Nu;++i) {
			mcptr->cost[i][j] = H[nei[uid][i]][curent_node->nodes_seq[j]];
		}
	}
	for (int j = Nv; j < mcdim; ++j) {
		for (int i = 0; i < Nu; ++i) {
			mcptr->cost[i][j] = 0;
		}
	}
	for (int i = Nu; i < mcdim; ++i) {
		for (int j = 0; j < mcdim; ++j) {
			mcptr->cost[i][j] = 0;
		}
	}
	int start, end;
	start = curent_node->common_size;
	end = max(min(Nv, mfdim), Nfv);
	for (int i = 0;i < mfdim;++i)mcptr->u[i] = mfptr->u[i], mcptr->rowsol[i] = -1;
	for (int j = 0;j < start;++j) {
		mcptr->v[j] = mfptr->v[curent_node->index[j]];
		mcptr->colsol[j] = mfptr->colsol[curent_node->index[j]];
		mcptr->rowsol[mcptr->colsol[j]] = j;
	}
	for (int j = end;j < mfdim;++j) {
		mcptr->v[j] = mfptr->v[j];
		mcptr->colsol[j] = mfptr->colsol[j];
		mcptr->rowsol[mcptr->colsol[j]] = j;
	}
	return make_pair(start, end);
}


void dfs_v2(int u, pathnode* root, vector<vector<int>>& nei, unordered_map<int, Mi*>& forlapjv, vector<vector<int>>& Ha, vector<vector<int>>& Hb, double beta)
{
	if (root == nullptr)
		return;
	Mi mc, * mcptr;
	mcptr = &mc;
	if (root->val == 0) {
		mcptr->dim = 0;
		forlapjv[0] = mcptr;
	}
	else if (root->fa->val == 0) {
		int vid = root->val - 1;
		int Nu = nei[u].size(), Nv = nei[vid].size();
		create_first_level(mcptr, root, u, nei, Ha);
		int weight = Hungarian(mc.dim, mc.cost, mc.u, mc.v, mc.rowsol, mc.colsol);
		//double weight = mcptr->dim - lap(mc.dim, mc.cost, mc.rowsol, mc.colsol, mc.u, mc.v);
		forlapjv[vid + 1] = mcptr; // fixed u, store info of N(v)

		if (vid < Hb.size() && vid > -1) {
			if (mcptr->dim == 0)
				Hb[u][vid] = beta * MAX_RS;
			else
				Hb[u][vid] = (1 - beta) * weight / max(Nu, Nv) + beta * MAX_RS;
		}

	}
	else {
		Mi* mfptr = forlapjv[root->fa->val];
		int f_dim = mfptr->dim;
		//int idx = create_v2(mfptr, mcptr, root, u, nei, Ha);
		//int weight = dymatch_col(f_dim, mcptr->dim, idx, mcptr->cost, mcptr->u, mcptr->v, mcptr->rowsol, mcptr->colsol);

		pair<int, int> diff = create_dyhungarian_col_space(mfptr, mcptr, root, nei, u, Ha);
		int weight = dymatch_col(diff.first, diff.second, mfptr->dim, mcptr->dim, mcptr->cost, mcptr->u, mcptr->v, mcptr->rowsol, mcptr->colsol);
		forlapjv[root->val] = &mc;
		if (root->val <= Hb.size() && root->val > 0) {
			int Nu = nei[u].size(), Nv = nei[root->val - 1].size();
			int maxn = max(Nu, Nv);
			if (maxn == 0)
				Hb[u][root->val - 1] = beta * MAX_RS;
			else
				Hb[u][root->val - 1] = (1 - beta) * weight / maxn + beta * MAX_RS;
		}
	}
	for (int i = 0; i < root->nexts.size(); ++i) {
		dfs_v2(u, root->nexts[i], nei, forlapjv, Ha, Hb, beta);
	}
	if (root->val != 0) {
		delete[] mc.colsol;
		delete[] mc.rowsol;
		delete[] mc.u;
		delete[] mc.v;
		for (int i = 0; i < mc.dim; ++i) {
			delete[] mc.cost[i];
		}
		delete[] mc.cost;
	}
}

void Iterative(pathnode* root, vector<vector<int>>& nei, vector<vector<int>>& Ha, vector<vector<int>>& Hb, double beta) {
	int n = Ha.size();
	for (int i = 0; i < n; ++i) {
		unordered_map<int, Mi*> forlapjv;
		dfs_v2(i, root, nei, forlapjv, Ha, Hb, beta);
	}
}

int create(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int uid, vector<vector<int>>& nei, vector<vector<int>>& Ha) {
	set<int>& comelem = curent_node->commom;
	set<int>& diffelem = curent_node->diff;

	int vid = curent_node->val - 1;
	int Nu = nei[uid].size(), Nv = nei[vid].size();
	int mcdim = max(mfptr->dim, Nv); // when Nv > mf.dim > Nu, mc.dim > mf.dim

	//assigment space for current node(mc)
	mcptr->dim = mcdim;
	mcptr->cost = new int* [mcdim];
	for (int i = 0; i < mcdim; ++i) mcptr->cost[i] = new int[mcdim];
	mcptr->rowsol = new int[mcdim];
	mcptr->colsol = new int[mcdim];
	mcptr->u = new int[mcdim];
	mcptr->v = new int[mcdim];
	mcptr->index = new int[mcdim];

	//Initialize the current node information
	int mfdim = mfptr->dim;
	int x = 0;//x is mc column index in costmatrix
	for (int j = 0; j < mfdim; ++j) {
		int node = mfptr->index[j];
		if (comelem.find(node) != comelem.end()) {
			mcptr->index[x] = node;
			for (int i = 0; i < Nu; ++i) {
				mcptr->cost[i][x] = mfptr->cost[i][j];
			}
			mcptr->v[x] = mfptr->v[j];
			mcptr->colsol[x] = mfptr->colsol[j];
			mcptr->rowsol[mcptr->colsol[x]] = x;
			++x;
		}
		else {
			mcptr->rowsol[mfptr->colsol[j]] = -1; //let the row of mc's costmatrix become unmatchecd
		}
	}
	set<int>::iterator it;
	for (it = diffelem.begin(); it != diffelem.end(); ++it) {
		int vn = *it; // v's neighbor
		mcptr->index[x] = vn;
		mcptr->colsol[x] = -1;
		for (int i = 0; i < Nu; ++i) {
			mcptr->cost[i][x] = 1 - Ha[nei[uid][i]][vn];
		}
		++x;
	}
	for (int i = x; i < mcdim; ++i) mcptr->colsol[i] = -1;

	for (int i = 0; i < mfdim; ++i) //save mc.u from mf.u
		mcptr->u[i] = mfptr->u[i];


	if (Nv < Nu) { //Fill the dummy node with 1.0
		for (int i = Nv; i < mcdim; ++i) {
			mcptr->index[i] = -1;
			for (int j = 0; j < Nu; ++j) {
				mcptr->cost[j][i] = 1.0;
			}
		}
	}
	for (int i = Nu; i < mcdim; ++i) {//patch the costmatrix into a square with 1.0
		for (int j = 0; j < mcdim; ++j) {
			mcptr->cost[i][j] = 1.0;
		}
	}

	int idx = comelem.size();// divied the matched column and unmatched column in mc.cost
	return idx;
}

//void dfs(int u, pathnode* root, vector<vector<int>>& nei, unordered_map<int, Mi*>& forlapjv, vector<vector<double>>& Ha, vector<vector<double>>& Hb, double beta) {
//	if (root == nullptr)
//		return;
//	Mi mc, * mcptr;
//	int N = 0;
//	mcptr = &mc;
//
//	if (root->val == 0) {
//		mc.dim = 0;
//		forlapjv[0] = mcptr;
//		//cout << "u = " << u << " root->val = " << root->val << endl;
//	}
//	else if (root->fa->val == 0) {
//		int vid = root->val - 1;
//
//		int Nu = nei[u].size(), Nv = nei[vid].size();
//		N = max(Nu, Nv);
//
//		mc.dim = N;
//
//
//		// initilse mc.cost matrix to all 1s
//
//		mc.cost = new double* [N];
//		for (int i = 0; i < N; ++i) {
//			mc.cost[i] = new double[N];
//		}
//
//		mc.rowsol = new int[N];
//		mc.colsol = new int[N];
//		mc.u = new double[N];
//		mc.v = new double[N];
//		mc.index = new int[N];
//		// initilse mc.cost matrix to s(x,y) for x in N(u) and y in N(v)
//		for (int j = 0; j < Nv; ++j) {
//			int jj = nei[vid][j];
//			for (int i = 0; i < Nu; ++i) {
//				mc.cost[i][j] = 1 - Ha[nei[u][i]][jj];
//			}
//			mc.index[j] = jj;
//		}
//		if (Nv > Nu) {
//			for (int i = Nu; i < N; ++i) {
//				for (int j = 0; j < N; ++j) {
//					mc.cost[i][j] = 1.0;
//				}
//			}
//		}
//		if (Nv < Nu) {
//			for (int j = Nv; j < N; ++j) {
//				for (int i = 0; i < N; ++i) {
//					mc.cost[i][j] = 1.0;
//				}
//			}
//		}
//
//		// initilse  mc.index
//		for (int j = Nv; j < N; ++j) {
//			mc.index[j] = -1;
//		}
//		// N==0 ?
//
//		// use lapjv to comptue max weighted matching
//		//double weight = N - lap(N, mc.cost, mc.rowsol, mc.colsol, mc.u, mc.v);
//		double weight = N - hungarian(mc.dim, mc.cost, mc.rowsol, mc.colsol, mc.u, mc.v);
//		//cout << "weight = " << weight << endl;
//
//		forlapjv[vid + 1] = mcptr; // fixed u, store info of N(v)
//
//		if (vid < Hb.size() && vid > -1) {
//			if (vid == u)
//				Hb[u][vid] = 1;
//			else if (N == 0)
//				Hb[u][vid] = 0;
//			else
//				Hb[u][vid] = (1 - beta) * weight / N + beta;
//		}
//		//delete[] rowsol;
//	}
//	else {
//		Mi* mfptr = forlapjv[root->fa->val];
//		int f_dim = mfptr->dim;
//		int idx = create(mfptr, mcptr, root, u, nei, Ha);
//
//		double weight = dymatch_col(f_dim, mcptr->dim, idx, mcptr->cost, mcptr->rowsol, mcptr->colsol, mcptr->u, mcptr->v);
//		weight = mc.dim - weight;
//
//		forlapjv[root->val] = &mc;
//		if (root->val <= Hb.size() && root->val > 0) {
//			int Nu = nei[u].size(), Nv = nei[root->val - 1].size();
//			int maxn = max(Nu, Nv);
//			if (root->val - 1 == u)
//				Hb[u][root->val - 1] = 1;
//			else if (maxn == 0)
//				Hb[u][root->val - 1] = 0;
//			else
//				Hb[u][root->val - 1] = (1 - beta) * weight / maxn + beta;
//		}
//	}
//
//	for (int i = 0; i < root->nexts.size(); ++i) {
//		dfs(u, root->nexts[i], nei, forlapjv, Ha, Hb, beta);
//	}
//	if (root->val != 0) {
//		delete[] mc.colsol;
//		delete[] mc.rowsol;
//		delete[] mc.u;
//		delete[] mc.v;
//		for (int i = 0; i < N; ++i) {
//			delete[] mc.cost[i];
//		}
//		delete[] mc.cost;
//		delete[] mc.index;
//	}
//}

void freeupspace(pathnode* root) {
	if (root == nullptr)
		return;

	for (int i = 0; i < root->nexts.size(); ++i) {
		freeupspace(root->nexts[i]);
	}
	delete root;
}

void dfs2(pathnode* root) {
	if (root == nullptr)
		return;
	cout << "root->val = " << root->val << endl;
	if (root->fa != nullptr)
		cout << "root->fa->val = " << root->fa->val << endl;
	set<int> com = root->commom, diff = root->diff;
	set<int>::iterator it;
	cout << "common elem: ";
	for (it = com.begin(); it != com.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
	cout << "diff elem: ";
	for (it = diff.begin(); it != diff.end(); ++it) {
		cout << *it << " ";
	}
	cout << endl;
	cout << endl << endl;
	for (int i = 0; i < root->nexts.size(); ++i) {
		dfs2(root->nexts[i]);
	}
}
