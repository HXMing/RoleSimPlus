#include "rolesim_plus2.h"
#include "sharepath.h"
#include "dyhungarian.h"
#include <iostream>
#include <ctime>

#define MAX_RS 10000
vector<vector<int>> RoleSim_Plus2(vector<vector<int>>& I, int k, double beta)
{
	clock_t start, end;
	start = clock();
    int n = I.size();
	vector<vector<int>> H(n, vector<int>(n, MAX_RS));
    pathnode* root = get_shared_path(I);
	end = clock();
	cout << "2-DT Phrase 1 time = " << (end - start) / 1000.0 << endl;
	for (int i = 0;i < k;++i) {
        H = Iterative(root, I, H, beta);
    }
	freeupspace(root);
	return H;
}

vector<vector<int>> Iterative(pathnode* root, vector<vector<int>>& I, vector<vector<int>>& H, double beta)
{
	int n = H.size();
	vector<vector<int>> H1(n, vector<int>(n));
	unordered_map<int, Mi*> forlapjv;
	dfs(root, I, forlapjv, H, H1, beta);
    return H1;
}





void dfs(pathnode* root, vector<vector<int>>& I, unordered_map<int, Mi*>& forlapjv, vector<vector<int>>& H, vector<vector<int>>& H1, double beta)
{
	if (root == nullptr)
		return;
	int n = I.size() + 1;
	int uid = root->row, vid = root->col;
	Mi mc, * mcptr;
	mcptr = &mc;
	if (uid == 0 || vid == 0) {
		mc.dim = 0;
		forlapjv[uid * n + vid] = mcptr;
	}
	else {
		int fuid = root->fa->row, fvid = root->fa->col;
		int Nu = I[uid - 1].size(), Nv = I[vid - 1].size();
		if (fuid == 0 || fvid == 0) {	//匈牙利算法计算
			create_hungarian_space(mcptr, root, Nu, Nv, H);
			int weight = Hungarian(mc.dim, mc.wei, mc.rd, mc.cd, mc.rs, mc.cs);
			forlapjv[uid * n + vid] = mcptr;
			int maxn = max(Nu, Nv);
			if (maxn == 0) 
				H1[uid - 1][vid - 1] = beta * MAX_RS;
			else
				H1[uid - 1][vid - 1] = (1 - beta) * weight / maxn + beta * MAX_RS;
		}
		else if (fuid == uid) {		//列转移
			Mi* mfptr = forlapjv[fuid * n + fvid];
			int Nfv = I[fvid - 1].size();
			pair<int, int> diff = create_dyhungarian_col_space(mfptr, mcptr, root, Nu, Nv, Nfv, H);
			int weight = dymatch_col(diff.first, diff.second, mfptr->dim, mcptr->dim, mcptr->wei, mcptr->rd, mcptr->cd, mcptr->rs, mcptr->cs);
			forlapjv[uid * n + vid] = mcptr;
			int maxn = max(Nu, Nv);
			if (maxn == 0)
				H1[uid - 1][vid - 1] = beta * MAX_RS;
			else
				H1[uid - 1][vid - 1] = (1 - beta) * weight / maxn + beta * MAX_RS;
			
		}
		else if (fvid == vid) {		//行转移
			Mi* mfptr = forlapjv[fuid * n + fvid];
			int Nfu = I[fuid - 1].size();
			pair<int, int> diff = create_dyhungarian_row_space(mfptr, mcptr, root, Nu, Nv, Nfu, H);
			int weight = dymatch_row(diff.first, diff.second, mfptr->dim, mcptr->dim, mcptr->wei, mcptr->rd, mcptr->cd, mcptr->rs, mcptr->cs);
			forlapjv[uid * n + vid] = mcptr;
			int maxn = max(Nu, Nv);
			if (mcptr->dim == 0)
				H1[uid - 1][vid - 1] = beta * MAX_RS;
			else
				H1[uid - 1][vid - 1] = (1 - beta) * weight / maxn + beta * MAX_RS;
		}
		
	}
	for (int i = 0;i < root->nexts.size();++i) {
		dfs(root->nexts[i], I, forlapjv, H, H1, beta);
	}
	if (uid != 0 && vid != 0) {
		delete[] mc.cs;
		delete[] mc.rs;
		delete[] mc.rd;
		delete[] mc.cd;
		for (int i = 0; i < mc.dim; ++i) {
			delete[] mc.wei[i];
		}
		delete[] mc.wei;
	}
}

void create_hungarian_space(Mi* mcptr, pathnode* curent_node, int Nu, int Nv, vector<vector<int>>& H)
{
	int mcdim = max(Nu, Nv);
	mcptr->dim = mcdim;
	mcptr->wei = new int* [mcdim];
	for (int i = 0; i < mcdim; ++i) mcptr->wei[i] = new int[mcdim];
	mcptr->rs = new int[mcdim];
	mcptr->cs = new int[mcdim];
	mcptr->rd = new int[mcdim];
	mcptr->cd = new int[mcdim];
	for (int i = 0;i < Nu;++i) {
		int ii = curent_node->row_node_seq[i];
		for (int j = 0;j < Nv;++j) {
			int jj = curent_node->col_node_seq[j];
			mcptr->wei[i][j] = H[ii][jj];
		}
	}
	if (Nv > Nu) {
		for (int i = Nu;i < mcdim;++i) {
			for (int j = 0;j < mcdim;++j) {
				mcptr->wei[i][j] = 0;
			}
		}
	}
	if (Nu > Nv) {
		for (int j = Nv;j < mcdim;++j) {
			for (int i = 0;i < mcdim;++i) {
				mcptr->wei[i][j] = 0;
			}
		}
	}
}

pair<int, int> create_dyhungarian_col_space(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int Nu, int Nv, int Nfv, vector<vector<int>>& H)
{
	int mcdim = max(mfptr->dim, Nv); // when Nv > mf.dim > Nu, mc.dim > mf.dim
	//assigment space for current node(mc)
	mcptr->dim = mcdim;
	mcptr->wei = new int* [mcdim]; 
	for (int i = 0; i < mcdim; ++i) mcptr->wei[i] = new int[mcdim];
	mcptr->rs = new int[mcdim];
	mcptr->cs = new int[mcdim];
	mcptr->rd = new int[mcdim];
	mcptr->cd = new int[mcdim];

	//Initialize the current node information
	int mfdim = mfptr->dim;
	for (int j = 0;j < curent_node->common_size;++j) {
		for (int i = 0;i < Nu;++i) {
			mcptr->wei[i][j] = mfptr->wei[i][curent_node->col_common_elem_index[j]];
		}
	}
	for (int j = curent_node->common_size;j < Nv;++j) {
		for (int i = 0;i < Nu;++i) {
			mcptr->wei[i][j] = H[curent_node->row_node_seq[i]][curent_node->col_node_seq[j]];
		}
	}
	for (int j = Nv; j < mcdim; ++j) {
		for (int i = 0; i < Nu; ++i) {
			mcptr->wei[i][j] = 0;
		}
	}
	for (int i = Nu; i < mcdim; ++i) {
		for (int j = 0; j < mcdim; ++j) {
			mcptr->wei[i][j] = 0;
		}
	}
	int start, end;
	start = curent_node->common_size;
	end = max(min(Nv, mfdim), Nfv);
	for (int i = 0;i < mfdim;++i)mcptr->rd[i] = mfptr->rd[i], mcptr->rs[i] = -1;
	for (int j = 0;j < start;++j) {
		mcptr->cd[j] = mfptr->cd[curent_node->col_common_elem_index[j]];
		mcptr->cs[j] = mfptr->cs[curent_node->col_common_elem_index[j]];
		mcptr->rs[mcptr->cs[j]] = j;
	}
	for (int j = end;j < mfdim;++j) {
		mcptr->cd[j] = mfptr->cd[j];
		mcptr->cs[j] = mfptr->cs[j];
		mcptr->rs[mcptr->cs[j]] = j;
	}
	return make_pair(start, end);
}

pair<int, int> create_dyhungarian_row_space(Mi* mfptr, Mi* mcptr, pathnode* curent_node, int Nu, int Nv, int Nfu, vector<vector<int>>& H)
{
	int mcdim = max(mfptr->dim, Nu);
	//assigment space for current node(mc)
	mcptr->dim = mcdim;
	mcptr->wei = new int* [mcdim];
	for (int i = 0; i < mcdim; ++i) mcptr->wei[i] = new int[mcdim];
	mcptr->rs = new int[mcdim];
	mcptr->cs = new int[mcdim];
	mcptr->rd = new int[mcdim];
	mcptr->cd = new int[mcdim];

	//Initialize the current node information
	int mfdim = mfptr->dim;
	for (int i = 0;i < curent_node->common_size;++i) {
		for (int j = 0;j < Nv;++j) {
			mcptr->wei[i][j] = mfptr->wei[curent_node->row_common_elem_index[i]][j];
		}
	}
	for (int i = curent_node->common_size;i < Nu;++i) {
		for (int j = 0;j < Nv;++j) {
			mcptr->wei[i][j] = H[curent_node->row_node_seq[i]][curent_node->col_node_seq[j]];
		}
	}
	for (int i = Nu; i < mcdim; ++i) {
		for (int j = 0; j < Nv; ++j) {
			mcptr->wei[i][j] = 0;
		}
	}
	for (int j = Nv; j < mcdim; ++j) {
		for (int i = 0; i < mcdim; ++i) {
			mcptr->wei[i][j] = 0;
		}
	}
	int start, end;
	start = curent_node->common_size;
	end = max(min(Nu, mfdim), Nfu);
	for (int j = 0;j < mfdim;++j)mcptr->cd[j] = mfptr->cd[j], mcptr->cs[j] = -1;
	for (int i = 0;i < start;++i) {
		mcptr->rd[i] = mfptr->rd[curent_node->row_common_elem_index[i]];
		mcptr->rs[i] = mfptr->rs[curent_node->row_common_elem_index[i]];
		mcptr->cs[mcptr->rs[i]] = i;
	}
	for (int i = end;i < mfdim;++i) {
		mcptr->rd[i] = mfptr->rd[i];
		mcptr->rs[i] = mfptr->rs[i];
		mcptr->cs[mcptr->rs[i]] = i;
	}
	return make_pair(start, end);
}
