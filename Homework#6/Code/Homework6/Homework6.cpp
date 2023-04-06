#include <iostream>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <string>
#include <queue>
#include <stack>
#include <map>
#include <vector>
#include <set>
#include <sstream>
#include <ctime>
#include <climits>
#include <tuple>
#define N 
#define ll long long
#define endl "\n"
#define MID(a,b) (a+b)/2
#define f double
using namespace std;


f** InverseMatrix(f** m) {
	f d = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	f** im = (f**)malloc(sizeof(f*)*3);
	for (int i = 0; i < 3; i++){
		im[i] = (f*)malloc(sizeof(f) * 3);
	}
	if (d != 0.0)

	{
		/*im[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * id;
		im[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * id;
		im[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * id;
		im[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * id;
		im[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * id;
		im[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * id;
		im[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * id;
		im[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * id;
		im[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * id;*/
	
		im[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) /d;
		im[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) /d;
		im[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) /d;
		im[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) /d;
		im[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) /d;
		im[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) /d;
		im[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) /d;
		im[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) /d;
		im[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) /d;
	}
	return im;
}

int main()
{
	int T = 1;
	while (T <= 3) {
		char filename[100];
		sprintf(filename, "fitdata%d.dat", T);
		printf("%s\n", filename);
		FILE* F = fopen(filename, "r");
		int n = 0;

		f x[100];
		f y[100];
		f x_[100];
		f y_[100];

		f sumx2 = 0;
		f sumxy = 0;
		f sumx = 0;
		f sumy2 = 0;
		f sumy = 0;
		f sum1 = 0;
		f sumxx_ = 0;
		f sumyx_ = 0;
		f sumx_ = 0;

		f sumxy_ = 0;
		f sumyy_ = 0;
		f sumy_ = 0;
		
		while (1) {
			int eof = fscanf(F, "%lf %lf %lf %lf", &x[n], &y[n], &x_[n], &y_[n]);
			if (eof == EOF) break;
			n++;
		}


		for (int i = 0; i < n; i++) {
			sumx += x[i];
			sumxy += x[i] * y[i];
			sumx2 += x[i] * x[i];
			sumy2 += y[i] * y[i];
			sumy += y[i];
			sum1 += 1;
			sumxx_ += x[i] * x_[i];
			sumyx_ += y[i] * x_[i];
			sumx_ += x_[i];
			sumxy_ += x[i] * y_[i];
			sumyy_ += y[i] * y_[i];
			sumy_ += y_[i];
		}
		
		f ans[3] = { 0,0,0 };
		f left[3] = { sumxx_,sumyx_,sumx_ };
		f left2[3] = { sumxy_,sumyy_,sumy_ };
		f** right = (f**)malloc(sizeof(f*) * 3);
		for (int i = 0; i < 3; i++) {
			right[i] = (f*)malloc(sizeof(f) * 3);
		}
		f tmp[3][3] = { {sumx2,sumxy,sumx}, {sumxy,sumy2,sumy}, {sumx,sumy,sum1 } };
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				right[i][j] = tmp[i][j];
			}
		}
		right = InverseMatrix(right);
		
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				ans[i] += right[i][j] * left[j];
			}
		}
		printf("-------------For fitdatad%d.dat-------------\n", T);
		for (int i = 0; i < 3; i++) {
			printf("a%d : %.15lf\n", i + 1, ans[i]);
			ans[i] = 0;
		}

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				ans[i] += right[i][j] * left2[j];
			}
		}
		for (int i = 3; i < 6; i++) {
			printf("a%d : %.15lf\n", i+1, ans[i - 3]);
		}



		T++;
		printf("-------------------------------------------\n");
		printf("\n\n");
	}

	return 0;
}