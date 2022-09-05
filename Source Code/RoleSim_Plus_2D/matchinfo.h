#pragma once
typedef struct Matchinfo {
	int dim;
	int** wei;
	int* rd;
	int* cd;
	int* rs;
	int* cs;
	int* index;
	//unordered_map<int, int> index;
} Mi;