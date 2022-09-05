#pragma once
int Hungarian(int dim, int** wei, int* rd, int* cd, int* rs, int* cs);
int dfs_r(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs);
int dfs_c(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs);
int dymatch_row(int start, int end, int fdim, int dim, int** wei, int* rd, int* cd, int* rs, int* cs);
int dymatch_col(int start, int end, int fdim, int dim, int** wei, int* rd, int* cd, int* rs, int* cs);