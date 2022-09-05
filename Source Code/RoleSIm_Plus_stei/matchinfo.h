#pragma once
#include <unordered_map>
using namespace std;
typedef struct Matchinfo {
	int dim;
	int** cost;
	int* u;
	int* v;
	int* colsol;
	int* rowsol;
	int* index;
	//unordered_map<int, int> index;
} Mi;