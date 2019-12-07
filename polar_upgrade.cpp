#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <queue>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
//#include <mpi.h>
using namespace std;
int mil = 32;
typedef struct heapnode {
	double data;
	int index;
} heapnode;
bool compare_heap(heapnode a, heapnode b) {
	return a.data < b.data;
}

class my_heap {

private:
	heapnode* data;
	int* position;
	int len;

public:

	my_heap(int n, double** w) {
		len = n - 1;
		data = new heapnode[n];
		position = new int[n - 1];
		for (int i = 0; i < n - 1; i++) {
			data[i].data = capacity(w[0][i], w[1][i]) + capacity(w[0][i + 1], w[1][i + 1]) - capacity(w[0][i] + w[0][i + 1], w[1][i] + w[1][i + 1]);
			data[i].index = i;
		}
		sort(data, data + n - 1, compare_heap);
		for (int i = 0; i < n - 1; i++)
		{
			position[data[i].index] = i;
		}
	}
	~my_heap() {
		delete[] data;
		delete[] position;
	}

	double capacity(double x, double y) {
		if (x != 0 && y != 0)
			return -(x + y) * log2((x + y) / 2) + x * log2(x) + y * log2(y);
		else {
			if (x == 0 && y != 0)
				return y;
			else
				if (y == 0 && x != 0)
					return x;
				else
					return 0;
		}
	}

	void swap(int i, int j) {
		heapnode hp;
		position[data[i].index] = j;
		position[data[j].index] = i;
		hp = data[i];
		data[i] = data[j];
		data[j] = hp;
	}

	void update(int pos, double delta) {
		int ind = position[pos];

		if (data[ind].data > delta) {
			data[ind].data = delta;
			int father = (ind + 1) / 2 - 1;
			while (father >= 0 && data[father].data > delta) {
				swap(father, ind);
				ind = father;
				father = (father + 1) / 2 - 1;
			}
		}
		else {
			data[ind].data = delta;
			int ls = (ind + 1) * 2 - 1;
			int rs = (ind + 1) * 2;
			while (true)
			{
				ls = (ind + 1) * 2 - 1;
				rs = (ind + 1) * 2;
				if (ls < len && (rs >= len || (data[ls].data < delta && data[ls].data <= data[rs].data))) {
					swap(ls, ind);
					ind = ls;
				}
				else if (rs < len && data[rs].data < delta && data[rs].data < data[ls].data) {
					swap(rs, ind);
					ind = rs;
				}
				else {
					break;
				}
			}
		}

	}

	void rm(int pos) {
		int ind = position[pos];
		data[ind].index = data[len - 1].index;
		position[data[len - 1].index] = ind;
		len--;
		update(pos, data[len].data);

	}


	heapnode getTop() {
		return data[0];
	}

};
double capacity(double x, double y) {
	if (x != 0 && y != 0)
		return -(x + y) * log2((x + y) / 2) + x * log2(x) + y * log2(y);
	else {
		if (x == 0 && y != 0)
			return y;
		else
			if (y == 0 && x != 0)
				return x;
			else
				return 0;
	}
}



enum Enumcomp { ASC, DESC };
class compare
{
private:
	Enumcomp comp;
public:
	compare(Enumcomp c) :comp(c) {};
	bool operator () (int num1, int num2)
	{
		switch (comp)
		{
		case ASC:
			return num1 < num2;
		case DESC:
			return num1 > num2;
		}
	}
};

void get_w_up(double** w, double** w_up, int len) {
	for (int u = 0; u < 2; u++)
		for (int i = 0; i < len; i++)
		{
			for (int j = 0; j < len; j++)
			{
				//W_up(u1 + 1, N * (y1 - 1) + y2) = 0.5 * (W(u1 + 1, y1) * W(1, y2) + W(mod(u1 + 1, 2) + 1, y1) * W(2, y2));
				w_up[u][i * len + j] = 0.5 * (w[u][i] * w[0][j] + w[u ^ 1][i] * w[1][j]);
			}
			//for (int j = 0; j < len; j++)
			//	w_up[u][i * len+len*len + j] = 0.5 * (w[u^1][i] * w[1][j] + w[u][i] * w[0][j]);
		}
}
void get_w_down_o(double** w, double** w_down, int len) {
	for (int u1 = 0; u1 < 2; u1++)
		for (int u2 = 0; u2 < 2; u2++)
			for (int i = 0; i < len; i++)
			{
				for (int j = 0; j < len; j++)
				{
					//W_down(u2 + 1, 2 * N * (y1 - 1) + 2 * (y2 - 1) + u1 + 1) 
					//= 0.5 * W(mod(u1 + u2, 2) + 1, y1) * W(u2 + 1, y2);
					w_down[u2][2 * len * i + 2 * j + u1] =
						0.5 * w[u1 ^ u2][i] * w[u2][j];
				}
			}
}

void get_w_down(double** w, double** w_down, int len) {
	for (int u2 = 0; u2 < 2; u2++)
		for (int i = 0; i < len; i++)
		{
			for (int j = 0; j < len; j++)
			{
				//W_down(u2 + 1, 2 * N * (y1 - 1) + 2 * (y2 - 1) + u1 + 1) 
				//= 0.5 * W(mod(u1 + u2, 2) + 1, y1) * W(u2 + 1, y2);
				w_down[u2][2 * len * i + 2 * j] =
					0.5 * w[u2][i] * w[u2][j];
				if (i < j)
					w_down[u2][2 * len * i + 2 * j + 1] =
					0.5 * w[u2][i] * w[u2 ^ 1][j];
				else
					w_down[u2][2 * len * i + 2 * j + 1] =
					0.5 * w[u2 ^ 1][i] * w[u2][j];
			}
		}
	double tmp;
	for (int i = 0; i < len * len; i++)
	{
		if (w_down[0][i] < w_down[1][i]) {
			tmp = w_down[0][i];
			w_down[0][i] = w_down[1][i];
			w_down[1][i] = tmp;
		}
	}
}


void degrading_merge_n(double** w_up, int& len, int miu, double** w) {
	if (len <= miu) {
		for (int i = 0; i < len; i++)
		{
			w[0][i] = w_up[0][i];
			w[1][i] = w_up[1][i];
		}
		return;
	}
	int* pre_g = new int[len];
	int* aft_g = new int[len];
	for (int i = 0; i < len - 1; i++)
	{
		pre_g[i] = i - 1;
		aft_g[i] = i + 1;
	}
	my_heap mh = my_heap(len, w_up);
	heapnode minnode;
	int minI = 0, final_node = len - 2;
	int len1 = len;

	while (len > miu)
	{
		minnode = mh.getTop();
		minI = minnode.index;
		int next = aft_g[minI];

		//update the w
		w_up[1][minI] = w_up[1][minI] + w_up[1][aft_g[minI]];
		w_up[0][minI] = w_up[0][minI] + w_up[0][aft_g[minI]];

		//update the previous deltaI
		if (minI > 0)
			mh.update(pre_g[minI], capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][pre_g[minI]], w_up[1][pre_g[minI]]) - capacity(w_up[0][minI] + w_up[0][pre_g[minI]], w_up[1][pre_g[minI]] + w_up[1][minI]));
		//deltaI[pre_g[minI]] = capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][pre_g[minI]], w_up[1][pre_g[minI]]) - capacity(w_up[0][minI] + w_up[0][pre_g[minI]], w_up[1][pre_g[minI]] + w_up[1][minI]);

	//update the next node
		if (minI == final_node) {
			mh.rm(minI);
			final_node = pre_g[minI];
		}
		else if (next == final_node) {
			final_node = minI;
			mh.update(minI, capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][aft_g[next]], w_up[1][aft_g[next]]) - capacity(w_up[0][minI] + w_up[0][aft_g[next]], w_up[1][aft_g[next]] + w_up[1][minI]));
			//deltaI[minI] = capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][aft_g[next]], w_up[1][aft_g[next]]) - capacity(w_up[0][minI] + w_up[0][aft_g[next]], w_up[1][aft_g[next]] + w_up[1][minI]);
			aft_g[minI] = aft_g[next];
			pre_g[aft_g[next]] = minI;
			mh.rm(next);
		}
		else {
			//update the deltaI
			mh.update(minI, capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][aft_g[next]], w_up[1][aft_g[next]]) - capacity(w_up[0][minI] + w_up[0][aft_g[next]], w_up[1][aft_g[next]] + w_up[1][minI]));
			//			deltaI[minI] = capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][aft_g[next]], w_up[1][aft_g[next]]) - capacity(w_up[0][minI] + w_up[0][aft_g[next]], w_up[1][aft_g[next]] + w_up[1][minI]);
			aft_g[minI] = aft_g[next];
			pre_g[aft_g[next]] = minI;
			mh.rm(next);
		}
		len -= 1;
	}

	for (int j = 0, i = 0; i < miu; i++)
	{
		w[0][i] = w_up[0][j];
		w[1][i] = w_up[1][j];
		j = aft_g[j];
	}
	delete[] pre_g;
	delete[] aft_g;
	return;
}

void erasure_symbol_merge(double** w, int& len) {
	int n = len;
	int offset = 0;
	double sump = 0;
	for (int i = n - 1; i >= 0; i--)
	{
		if (w[0][i] == w[1][i]) {
			sump += w[0][i];
			//sump += w[0][len - 1 - i];
			offset++;
		}
		else {
			break;
		}
	}
	if (offset > 1) {
		w[0][n - offset] = sump;
		w[1][n - offset] = sump;
		len -= (offset - 1);
	}

}

//double capacity(double x, double y) {
//	if (x != 0 && y != 0)
//		return -(x + y) * log2((x + y) / 2) + x * log2(x) + y * log2(y);
//	else {
//		if (x == 0 && y != 0)
//			return y;
//		else
//			if (y == 0 && x != 0)
//				return x;
//			else
//				return 0;
//	}
//}

double deltaIf(double a1, double b1, double a2, double b2) {
	if (a1 == 0 || a2 == 0 || b1 == 0 || b2 == 0) {
		return capacity(a1, b1) + capacity(a2, b2) - capacity(a1 + a2, b1 + b2);
	}
	double s1 = a1 * log2(a1 * (a1 + a2 + b1 + b2)) + a2 * log2(a2 * (a1 + a2 + b1 + b2))
		+ b1 * log2(b1 * (a1 + a2 + b1 + b2)) + b2 * log2(b2 * (a1 + a2 + b1 + b2));
	double s2 = -((a1 + b1) * log2(a1 + b1) + (a2 + b2) * log2(b2 + a2) + (a1 + a2) * log2(a1 + a2) + (b1 + b2) * log2(b1 + b2));
	return s1 + s2;

}

typedef struct mergenode {
	int index;
	double deltaI;
}mergenode;

/*qsort\u8981\u81ea\u5df1\u5b9a\u4e49\u51fd\u6570\uff0c\u4e0b\u9762\u662f\u5b9a\u4e49\u7684\u51fd\u6570*/
int cmp_mergenode(const void* a, const void* b)
{

	return ((mergenode*)a)->deltaI > ((mergenode*)b)->deltaI ? 1 : 0;
	/*\u5982\u679ca\u5927\u4e8eb\uff0c\u8fd4\u56de\u771f\uff0c\u53cd\u4e4b\u8fd4\u56deb
	\u8fd9\u4e2a\u51fd\u6570\u662f\u5347\u5e8f\u6392\u5217\u7684\uff0c\u5982\u679c\u8981\u964d\u5e8f\u6392\u5217\u5c31\u53cd\u8fc7\u6765\u8f93\u51fa*/
}

void degrading_merge(double** w_up, int& len, int miu, double** w) {
	if (len <= miu) {
		for (int i = 0; i < len; i++)
		{
			w[0][i] = w_up[0][i];
			w[1][i] = w_up[1][i];
		}
		return;
	}
	int* pre_g = new int[len];
	int* aft_g = new int[len];
	double* deltaI = new double[len];
	for (int i = 0; i < len - 1; i++)
	{
		pre_g[i] = i - 1;
		aft_g[i] = i + 1;
		deltaI[i] = capacity(w_up[0][i], w_up[1][i]) + capacity(w_up[0][i + 1], w_up[1][i + 1]) - capacity(w_up[0][i] + w_up[0][i + 1], w_up[1][i] + w_up[1][i + 1]);
		//if(i==0)
		/* code */
	}
	int minI = 0, final_node = len - 2;
	int len1 = len;
	while (len > miu)
	{
		minI = 0;
		double mindelta = deltaI[0], sc = 0;
		for (int i = 0, j = 0; j < len - 1; j++)
		{
			if (deltaI[i] < mindelta)
			{
				mindelta = deltaI[i];
				minI = i;
			}
			i = aft_g[i];
		}

		int next = aft_g[minI];

		//update the w
		w_up[1][minI] = w_up[1][minI] + w_up[1][aft_g[minI]];
		w_up[0][minI] = w_up[0][minI] + w_up[0][aft_g[minI]];

		//update the previous deltaI
		if (minI > 0)
			deltaI[pre_g[minI]] = capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][pre_g[minI]], w_up[1][pre_g[minI]]) - capacity(w_up[0][minI] + w_up[0][pre_g[minI]], w_up[1][pre_g[minI]] + w_up[1][minI]);

		//update the next node
		if (minI == final_node) {
			final_node = pre_g[minI];

		}
		else if (next == final_node) {
			final_node = minI;
			deltaI[minI] = capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][aft_g[next]], w_up[1][aft_g[next]]) - capacity(w_up[0][minI] + w_up[0][aft_g[next]], w_up[1][aft_g[next]] + w_up[1][minI]);
			aft_g[minI] = aft_g[next];
			pre_g[aft_g[next]] = minI;
		}
		else {
			//update the deltaI
			deltaI[minI] = capacity(w_up[0][minI], w_up[1][minI]) + capacity(w_up[0][aft_g[next]], w_up[1][aft_g[next]]) - capacity(w_up[0][minI] + w_up[0][aft_g[next]], w_up[1][aft_g[next]] + w_up[1][minI]);
			aft_g[minI] = aft_g[next];
			pre_g[aft_g[next]] = minI;
		}
		len -= 1;
		for (int si = 0, sj = 0; si < len; si++) {
			sc += w_up[0][sj];
			sc += w_up[1][sj];
			sj = aft_g[sj];
		}
		sc += 0;
	}

	for (int j = 0, i = 0; i < miu; i++)
	{
		w[0][i] = w_up[0][j];
		w[1][i] = w_up[1][j];
		j = aft_g[j];
	}
	delete[] pre_g;
	delete[] aft_g;
	delete[] deltaI;
	return;
}



typedef struct sortnode {
	double lr;
	int index;
	double p, q;
} sortnode;


/*qsort \u964d\u5e8f\u51fd\u6570*/
int cmp_func_decend(const void* a, const void* b)
{
	sortnode* a1 = (sortnode*)a;
	sortnode* b1 = (sortnode*)b;
	if (a1->lr > b1->lr)return -1;
	else if (a1->lr < b1->lr)return 1;
	else if (a1->p == b1->p)return 0;
	else if (a1->p > b1->p)return -1;
	else return 1;
	//return ((sortnode*)a)->lr < ((sortnode*)b)->lr ? 1 : -1;
	/*a>b \u4e14\u8fd4\u56de\u4e3a\u6b63\uff0c\u5219\u5347\u5e8f\u6392\u5217\uff1b\u8fd4\u56de\u4e3a0\uff0c\u5219\u4e0d\u786e\u5b9a\u987a\u5e8f*/
}

bool compare(sortnode a, sortnode b)
{
	if (a.lr > b.lr)return true;
	else return false;
}


/*qsort \u5347\u5e8f\u51fd\u6570*/
int cmp_func_acend(const void* a, const void* b)
{
	return ((sortnode*)a)->lr > ((sortnode*)b)->lr ? 1 : -1;
}

void LR_sort(double** w, int len) {
	sortnode* sn = new sortnode[len];
	for (int i = 0; i < len; i++)
	{
		sn[i].lr = log(w[0][i]) - log(w[1][i]);
		sn[i].index = i;
		sn[i].p = w[0][i];
		sn[i].q = w[1][i];
	}
	sort(sn, sn + len, compare);
	//qsort(sn, len, sizeof(sortnode), cmp_func_decend);
	for (int i = 0; i < len; i++)
	{
		w[0][i] = sn[i].p;
		w[1][i] = sn[i].q;
	}
	delete[] sn;
}

void LR_sort_down(double** w, int len) {
	sortnode* sn = new sortnode[len];
	for (int i = 0; i < len; i++)
	{
		sn[i].lr = log(w[0][i]) - log(w[1][i]);
		sn[i].index = i;
		sn[i].p = w[0][i];
		sn[i].q = w[1][i];
	}
	sort(sn, sn + len, compare);
	//qsort(sn, len, sizeof(sortnode), cmp_func_decend);
	for (int i = 0; i < len; i++)
	{
		w[0][2 * i] = sn[i].p;
		w[1][2 * i] = sn[i].q;
		w[0][2 * i + 1] = sn[i].p;
		w[1][2 * i + 1] = sn[i].q;
	}
	delete[] sn;
}

void get_Pe(double** w, double** w_up, double** w_down, int len, int bgn, int n, double* pe) {
	if (n == 1) {
		double sumw = 0;
		for (int i = 0; i < len; i++)
		{
			sumw += (w[0][i] > w[1][i] ? w[1][i] : w[0][i]);
		}
		pe[bgn] = sumw;
		//cout << "index: " << bgn << " error rate: " << double(pe[bgn]) << endl;
		cout << double(pe[bgn]) << endl;
		return;
	}
	double** w_tmp = new double* [2];
	for (int i = 0; i < 2; i++)
		w_tmp[i] = new double[mil / 2];
	//w_tmp[1] = new double[mil/2];
	int len1 = len * len * 2;
	int len2 = len1 * 2;
	//get_w_up(w, w_up, len);
	//LR_sort(w_up, len1);
	//erasure_symbol_merge(w_up, len1);
	//degrading_merge(w_up, len1, mil/2, w_tmp);
	//get_Pe(w_tmp, w_up,w_down,len1, bgn, n / 2, pe);



	//for (int i = 0; i < len2; i++)
	//{
	//	if (w_down[0][i] != w_down[1][len2 - 1 - i]) {
	//		cout << 1 << endl;
	//	}
	//}
	double sc = 0;
	get_w_down(w, w_down, len);
	LR_sort_down(w_down, len2 / 2);
	erasure_symbol_merge(w_down, len2);
	degrading_merge_n(w_down, len2, mil / 2, w_tmp);
	for (int i = 0; i < len2; i++) {
		sc += w_tmp[0][i] + w_tmp[1][i];
	}
	cout << sc << endl;
	get_Pe(w_tmp, w_up, w_down, len2, bgn + n / 2, n / 2, pe);
	for (int i = 0; i < 2; i++)
	{
		delete[] w_tmp[i];
	}
	delete[] w_tmp;


}


void degrade_find_frozen_bits(double p, int N, int rank, int numpro) {
	int mil1 = mil / 2;
	double** w = new double* [2];
	w[0] = new double[mil1];
	w[1] = new double[mil1];
	double** w_up = new double* [2];
	w_up[0] = new double[mil * mil / 2];
	w_up[1] = new double[mil * mil / 2];
	double* pe = new double[N];
	double** w_down = new double* [2];
	w_down[0] = new double[mil * mil];
	w_down[1] = new double[mil * mil];


	w[0][0] = 1 - p;
	w[0][1] = p;
	w[1][0] = p;
	w[1][1] = 1 - p;
	get_Pe(w, w_up, w_down, 1, 0, N, pe);
	delete[] w_down[0];
	delete[] w_down[1];
	delete[] w_up[0];
	delete[] w[0];
	delete[] w[1];
	delete[] w_up[1];
	delete[] w_down;
	delete[] w_up;
	delete[] w;
	delete[] pe;
	//sortnode* sn = new sortnode[N];
	//for (int i = 0; i < N; i++)
	//{
	//	sn[i].lr = pe[i];
	//	sn[i].index = i;
	//}
	//qsort(sn, N, sizeof(sortnode), cmp_func_acend);
	//for (int i = 0; i < N; i++)
	//{
	//	if (i > k)frozenbits[sn[i].index] = 0;
	//	else frozenbits[sn[i].index] = -1;
	//}
}

void updateW(double** w, int layers, int bgn, double** w_tmp, int* lens) {
	int len;
	if (bgn == -1) {
		bgn = 0;
		len = lens[bgn] * lens[bgn] * 2;
		get_w_up(w, w_tmp, lens[bgn]);
		LR_sort_down(w_tmp, len / 2);
		erasure_symbol_merge(w_tmp, len);
		degrading_merge_n(w_tmp, len, mil / 2, w + (bgn + 1) * 2);
		lens[bgn + 1] = len;
	}
	else {
		len = lens[bgn] * lens[bgn] * 2 * 2;
		get_w_down(w + 2 * bgn, w_tmp, lens[bgn]);
		LR_sort_down(w_tmp, len / 2);
		erasure_symbol_merge(w_tmp, len);
		degrading_merge_n(w_tmp, len, mil / 2, w + (bgn + 1) * 2);
		lens[bgn + 1] = len;
	}

	for (int i = bgn + 1; i < layers - 1; i++)
	{
		len = lens[i] * lens[i] * 2;
		get_w_up(w + 2 * i, w_tmp, lens[i]);
		LR_sort_down(w_tmp, len / 2);
		erasure_symbol_merge(w_tmp, len);
		degrading_merge_n(w_tmp, len, mil / 2, w + (i + 1) * 2);
		lens[i + 1] = len;
	}
}



void updateW(double** w, int layers, int bgn, int num, double** w_tmp, int* lens) {
	int len;
	if (bgn == -1) {
		bgn = 0;
	}

	for (int i = bgn; i < layers - 1; i++)
	{
		if ((num >> (layers - 2 - i)) & 1 > 0) {
			len = lens[i] * lens[i] * 2 * 2;
			get_w_down(w + 2 * i, w_tmp, lens[i]);
			LR_sort_down(w_tmp, len / 2);
			erasure_symbol_merge(w_tmp, len);
			degrading_merge_n(w_tmp, len, mil / 2, w + (i + 1) * 2);
			lens[i + 1] = len;
		}
		else {
			len = lens[i] * lens[i] * 2;
			get_w_up(w + 2 * i, w_tmp, lens[i]);
			LR_sort_down(w_tmp, len / 2);
			erasure_symbol_merge(w_tmp, len);
			degrading_merge_n(w_tmp, len, mil / 2, w + (i + 1) * 2);
			lens[i + 1] = len;
		}


	}
}


void degrade_find_frozen_bits(double p, int N, int nums, int rank, string outfile) {
	ofstream outf;
	string outfile1 = outfile + "_back_" + to_string(nums) + "_" + to_string(rank) + ".txt";
	FILE* fp = fopen(outfile1.c_str(), "r");
	int bbb = 0;
	if (fp == NULL) {
		
	}
	else {
		outf.open(outfile + "_back_" + to_string(nums) + "_" + to_string(rank) + ".txt_tmp");
		double pre = -1, now;
		while (true) {
			if (fscanf(fp, " %lf ", &now) <= 0)break;
			if (now > 0.5 || now < 0)break;
			if (pre > 0) { outf << pre << endl; bbb++; }
			pre = now;
		}
		outf.flush();
		fclose(fp);
		outf.close();
		remove(outfile1.c_str());
		rename((outfile + "_back_" + to_string(nums) + "_" + to_string(rank) + ".txt_tmp").c_str(), outfile1.c_str());
	}


	outf.open(outfile + "_back_" + to_string(nums) + "_" + to_string(rank) + ".txt", ios::app);
	int tmp = 0;
	tmp = N / nums * rank+bbb;

	int layers = int(log2(N)) + 1;
	double** w = new double* [layers * 2];
	int* lens = new int[layers];
	lens[0] = 1;
	for (int i = 0; i < layers * 2; i++)
	{
		w[i] = new double[mil / 2];
	}
	w[0][0] = 1 - p;
	w[1][0] = p;
	double** w_tmp = new double* [2];
	w_tmp[0] = new double[mil * mil];
	w_tmp[1] = new double[mil * mil];
	if (tmp > 0)
		updateW(w, layers, -1, tmp - 1, w_tmp, lens);
	for (int i = tmp; i < tmp+N/nums; i++)
	{
		int tmpi = i, layer = 1;
		while (!(tmpi & 1) && layer < layers) { tmpi >>= 1; layer++; }

		updateW(w, layers, layers - layer - 1, w_tmp, lens);

		double sumw = 0;
		for (int t = 0; t < lens[layers - 1]; t++)
		{
			//sumw += (w[2*(layers-1)][t] > w[2 * (layers - 1)+1][t] ? w[2 * (layers - 1)+1][t] : w[2 * (layers - 1)][t]);
			sumw += w[2 * (layers - 1) + 1][t];
		}
		//pe[bgn] = sumw;
		//cout << "index: " << bgn << " error rate: " << double(pe[bgn]) << endl;
		//cout << double(sumw) << endl;
		outf << i << "\t"<<sumw << endl;
	}

	for (int i = 0; i < layers * 2; i++)
	{
		delete[] w[i];
	}
	delete[] w_tmp[0];
	delete[] w_tmp[1];
	delete[] w_tmp;
	delete[] w;
	delete[] lens;
	outf.close();
}


void degrade_find_frozen_bits(double p, int N, string outfile) {
	ofstream outf;
	outf.open(outfile);
	int layers = int(log2(N)) + 1;
	double** w = new double* [layers * 2];
	int* lens = new int[layers];
	lens[0] = 1;
	for (int i = 0; i < layers * 2; i++)
	{
		w[i] = new double[mil / 2];
	}
	w[0][0] = 1 - p;
	w[1][0] = p;
	double** w_tmp = new double* [2];
	w_tmp[0] = new double[mil * mil];
	w_tmp[1] = new double[mil * mil];

	for (int i = 0; i < N; i++)
	{
		int tmpi = i, layer = 1;
		while (!(tmpi & 1) && layer < layers) { tmpi >>= 1; layer++; }

		updateW(w, layers, layers - layer - 1, w_tmp, lens);

		double sumw = 0;
		for (int t = 0; t < lens[layers - 1]; t++)
		{
			//sumw += (w[2*(layers-1)][t] > w[2 * (layers - 1)+1][t] ? w[2 * (layers - 1)+1][t] : w[2 * (layers - 1)][t]);
			sumw += w[2 * (layers - 1) + 1][t];
		}
		//pe[bgn] = sumw;
		//cout << "index: " << bgn << " error rate: " << double(pe[bgn]) << endl;
		cout << double(sumw) << endl;
		//outf << sumw << endl;
	}

	for (int i = 0; i < layers * 2; i++)
	{
		delete[] w[i];
	}
	delete[] w_tmp[0];
	delete[] w_tmp[1];
	delete[] w_tmp;
	delete[] w;
	delete[] lens;
	outf.close();
}


int main(int argc, char** argv) {
	int length, logn = 5;
	double p = 0.02;
	int pi = 2;

	int nums = 4;
	int rank = 1;

	for (int i = 1; i < argc; i += 2)
	{
		if (string(argv[i]) == string("-N") && i + 1 < argc) {
			logn = atoi(argv[i + 1]);
		}
		else if (string(argv[i]) == string("-mil") && i + 1 < argc) {
			mil = atoi(argv[i + 1]);
		}
		else if (string(argv[i]) == string("-p") && i + 1 < argc) {
			pi = atoi(argv[i + 1]);
			p = pi * 1.0 / 100;
		}
		else if (string(argv[i]) == string("-ps") && i + 1 < argc) {
			pi = atoi(argv[i + 1]);
			p = pi * 1.0 / 1000;
		}
		else if (string(argv[i]) == string("-nums") && i + 1 < argc) {
			nums = atoi(argv[i + 1]);
		}
		else if (string(argv[i]) == string("-rank") && i + 1 < argc) {
			rank = atoi(argv[i + 1]);
		}

	}

	string filename = string("frozen_bits_mil=") + to_string(mil) + string("_N=") + to_string(logn)
		+ string("_p=") + to_string(pi) + string(".txt");
	length = (1 << logn);

	//MPI_Init(&argc, &argv);
	//MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	//MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	//	degrade_find_frozen_bits(0.02, length, 0, 0);

	degrade_find_frozen_bits(p, length, nums, rank, filename);
	/*
		MPI_Finalize();*/

		//double** w = new double* [2];
		//w[0] = new double[20];
		//w[1] = new double[20];
		//for (int i = 0; i < 20; i++) {
		//	w[0][i] = rand()%10+20;
		//}
		//my_heap mh = my_heap(21, w);
		//for(int i=0;i<20;i++)
		//	mh.update(19-i,i*1.0);
		//heapnode hp;
		//for (int i = 0; i < 20; i++) {
		//	hp = mh.getTop();
		//	cout << hp.index << "  " << hp.data << endl;
		//	mh.rm(hp.index);
		//}



	return 0;

	//double a1 = 0.0067241780429055957;
	//double b1 = 0.0010883219570943996;
	//double b2 = 0.027208048927359983;
	//double a2 = 0.16810445107263988;
	//double s1 = capacity(a1, b1);
	//double s2 = capacity(a2, b2);
	//double s3 = capacity(a1 + a2, b1 + b2);
	//cout<< capacity(a1, b1)  <<endl;
	//cout << capacity(a2, b2) << endl;
	//cout << capacity(a1 + a2, b1 + b2) << endl;
	//cout << s1+s2-s3 << endl;
	//cout << deltaIf(a1,b1,a2,b2) << endl;
	//system("pause");
}