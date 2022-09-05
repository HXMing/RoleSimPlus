#pragma once
#include <vector>
#include <set>
using namespace std;
typedef struct pathnode {
	vector<pathnode*> nexts;
	int val; //指示在斯坦纳树中的节点编号，与图中的编号偏移为1 ，例如斯坦树中 1-->4 == 5 表示图中0-->3的转移代价是5
	pathnode* fa;

	set<int> commom;
	set<int> diff;

	int* nodes_seq;
	int* index; //common elem index,get fahter's common elem index in cost matrix
	int common_size;
	int diff_size;

	pathnode(int v, pathnode* f) :val(v), fa(f) {}
}pathnode;