#pragma once
#include <vector>

//=================hungarian and dyhungarian 处理int类型====================
int Hungarian(int dim, int** wei, int* rd, int* cd, int* rs, int* cs);
int dymatch_col(int f_dim, int dim, int idx, int** wei, int* rd, int* cd, int* rs, int* cs);
int dfs(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs);
int dymatch_col(int start, int end, int f_dim, int dim, int** wei, int* rd, int* cd, int* rs, int* cs);
int dfs_c(int c, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs);
int dfs_r(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs);