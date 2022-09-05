#pragma once
#include <vector>
#include <set>
using namespace std;

typedef struct pathnode {
	int row,col;
	pathnode* fa;
	vector<pathnode*> nexts;
	short* row_node_seq;			//�нڵ����д���
	short* col_node_seq;			//�нڵ����д���
	short* row_common_elem_index;	//�й�ͬ�ڵ�����������ڸ��ڵ�ĵڼ���
	short* col_common_elem_index;	//�й�ͬ�ڵ�����������ڸ��ڵ�ĵڼ���
	int common_size, diff_size;		
	pathnode(int r, int c,pathnode* f): row(r),col(c), fa(f) {}
}pathnode;