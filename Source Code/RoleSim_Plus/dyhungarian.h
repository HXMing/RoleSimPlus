#pragma once
#include <vector>
#include "matchinfo.h"
using namespace std;
double Hungarian(int dim, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs);
double dyHungarian_r(int dim, int r, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs);
double dyHungarian_c(int dim, int c, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs);
int dfs(int r, int dim, vector<vector<double>>& wei, vector<double>& rd, vector<double>& cd, vector<int>& rs, vector<int>& cs);

//--------------------------------------
//--------------------------------------

double Hungarian(Mi& mi);
int dfs(int r, Mi& mi);
//double dymatch_col(int idx, Mi& mf, Mi& mc);
int dfs(int r, int dimx, double** wei, int* rs, int* cs, double* rd, double* cd);

//=================hungarian and dyhungarian 处理int类型====================
int Hungarian(int dim, int** wei, int* rd, int* cd, int* rs, int* cs);
int dymatch_col(int f_dim, int dim, int idx, int** wei, int* rd, int* cd, int* rs, int* cs);
int dfs(int r, int dimx, int** wei, int* rd, int* cd, int* rs, int* cs);
