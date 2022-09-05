#pragma once
#include <vector>
#include <set>
using namespace std;

typedef struct pathnode {
	int row,col;
	pathnode* fa;
	vector<pathnode*> nexts;
	short* row_node_seq;			//行节点排列次序
	short* col_node_seq;			//列节点排列次序
	short* row_common_elem_index;	//行共同节点索引，相对于父节点的第几行
	short* col_common_elem_index;	//列共同节点索引，相对于父节点的第几列
	int common_size, diff_size;		
	pathnode(int r, int c,pathnode* f): row(r),col(c), fa(f) {}
}pathnode;