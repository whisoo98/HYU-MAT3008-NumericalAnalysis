#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
//#include "nr.h"
//#include "nrutil.h"

#define NP 20
#define MP 20
#define MAXSTR 80
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define TINY 1.0e-20
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static float maxarg1, maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static int iminarg1, iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

void gaussj(float **a, int n, float **b, int m, int* gaussE);
void ludcmp(float **a, int n, int *indx, float *d, int* ludcmpE);
void svdcmp(float **a, int m, int n, float w[], float **v);
void mprove(float **a, float **alud, int n, int indx[], float b[], float x[]);
void lubksb(float **a, int n, int *indx, float b[]);
int *ivector(long nl, long nh);
float *vector(long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
float **matrix(long nrl, long nrh, long ncl, long nch);

void free_ivector(int *v, long nl, long nh);
void free_vector(float *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch);


float pythag(float a, float b); // pythagoras

int main(void) {
	/*
	Routine:
	1. Gauss-Jordan
	2. LU Decomposition
	=> 1 & 2 ways can't handle sigular matrix that does not have a inverse matrix.
	3. Iterative Improvement
	4. Singular Value Decomposition
	*/
	int i, j, k, l, m, n;
	float **a, **b, **lu, **svdu, **svdv;
	float **gaussai, **ludai, **svdai;
	float **gaussx, *ludx, *mpvx, **svdx;
	int *indx;
	float d, *col, *svdw, *tb;
	a = matrix(1, NP, 1, NP);
	b = matrix(1, NP, 1, NP);
	lu = matrix(1, NP, 1, NP);
	svdu = matrix(1, NP, 1, NP);
	svdv = matrix(1, NP, 1, NP);
	gaussai = matrix(1, NP, 1, NP);
	ludai = matrix(1, NP, 1, NP);
	svdai = matrix(1, NP, 1, NP);
	gaussx = matrix(1, NP, 1, NP);
	ludx = vector(1, NP);
	mpvx = vector(1, NP);
	svdx = matrix(1, NP, 1, NP);
	indx = ivector(1, NP);
	col = vector(1, NP);
	svdw = vector(1, NP);
	tb = vector(1, NP);
	FILE *fp;
	int T = 0;
	char seq[3][50] = { "1st","2nd","3rd" };
	while (T < 3) {
		int gaussE = 0, ludcmpE = 0;
		d = 0;
		T++;
		char filename[100];
		sprintf(filename, "lineq%d.dat", T);
		//printf("%s\n", filename);
		fp = fopen(filename, "r");
		for (i = 0; i < NP; i++) { // Matrix Initializing to zero
			for (j = 0; j < NP; j++) {
				a[i][j] = 0;
				b[i][j] = 0;
				lu[i][j] = 0;
				svdu[i][j] = 0;
				svdv[i][j] = 0;
				gaussai[i][j] = 0;
				ludai[i][j] = 0;
				svdai[i][j] = 0;
				gaussx[i][j] = 0;
				svdx[i][j] = 0;
			}
			indx[i] = 0;
			col[i] = 0;
			svdw[i] = 0;
			ludx[i] = 0;
			mpvx[i] = 0;
			tb[i] = 0;
		}
		char lineq[50];
		sprintf(lineq, "%s A Matrix of Linear Equation Ax=b", seq[T - 1]);
		printf("---------- %s ----------\n",lineq);
		printf("\n");

		fscanf(fp, "%d %d", &n, &m);
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= m; j++) {
				fscanf(fp, " %f", &a[i][j]); // A Matrix of Ax=b
				gaussai[i][j] = a[i][j]; // copy of a
				lu[i][j] = a[i][j]; // copy of a
				svdu[i][j] = a[i][j]; // copy of a
				printf("%7.3f ", a[i][j]);
			}
			printf("\n");
		}
		printf("\n");
		sprintf(lineq, "%s b Matrix of Linear Equation Ax=b", seq[T - 1]);
		printf("---------- %s ----------\n", lineq);
		printf("\n");

		for (i = 1; i <= m; i++) {
			fscanf(fp, " %f", &b[i][1]); // b Matrix of Ax=b
			gaussx[i][1] = b[i][1]; // gaussx is copy of b
			ludx[i] = b[i][1]; // ludx is copy of b
			tb[i] = b[i][1];
			printf("%7.3f ", b[i][1]);
		}
		printf("\n");
		printf("\n");


		gaussj(gaussai, n, gaussx, m, &gaussE); // Gauss-Jordan Elimination : gaussai will be an Inverse Matrix of a, x will be solution matrix of Ax=b

		ludcmp(lu, n, indx, &d, &ludcmpE); // LU Decompostion : lu will be a LU Matrix, indx will save biggest element idx
		if (ludcmpE == 0) {
			for (j = 1; j <= n; j++) d *= lu[j][j]; // Calculating the determinant of a matrix
			lubksb(lu, n, indx, ludx); // LU Backward substitution : ludx will be a solution matrix of Ax=b
			for (j = 1; j <= n; j++) { // Obatining inverse matrix using LU Decomposition
				for (i = 1; i <= n; i++) {
					col[i] = 0.0;
				}
				mpvx[j] = ludx[j];
				col[j] = 1.0;
				lubksb(lu, n, indx, col);

				for (i = 1; i < n; i++) ludai[i][j] = col[i];

			}
			mprove(a, lu, n, indx, tb, mpvx); // Iterative Improvement of a Solution
		}

		svdcmp(svdu, m, n, svdw, svdai); // svdu, svdw are matrix of SVD U & w, svdai is marix of SVD V(not a transpose)

		for (i = 1; i <= n; i++) {
			for (j = 1; j <= m; j++) {
				//if (fabs(svdw[j]) > 0.001) svdv[i][j] = svdai[i][j] / svdw[j]; // svdv is a matrix of V * w
				if ((svdw[j]) != 0.0) svdv[i][j] = svdai[i][j] / svdw[j]; // svdv is a matrix of V * w
				else svdv[i][j] = 0;
			}
		}

		for (i = 1; i <= n; i++) {
			for (j = 1; j <= m; j++) {
				svdai[i][j] = 0;
				for (k = 1; k <= m; k++) {
					svdai[i][j] += svdv[i][k] * svdu[j][k]; // svdai is inverse matrix of a using svd
				}
			}
		}

		for (i = 1; i <= n; i++) {
			for (k = 1; k <= m; k++) {
				svdx[i][1] += svdai[i][k] * b[k][1];
			}
		}


		printf("---------- About Matrix of Linear Equation %d ----------\n", T);
		printf("\n");
		

		if (gaussE == 0) {
			printf("---------- Determinant ----------\n");
			printf("\n");

			printf("%f \n", d);

			printf("\n");
			printf("---------- Solution of Gauss-Jordan Elimination ----------\n");
			printf("\n");

			for (i = 1; i <= m; i++) {
				printf("%10.6f ", gaussx[i][1]);
			}
			printf("\n");
			printf("\n");

			printf("---------- Inverse Matrix using Gauss-Jordan Elimination ----------\n");
			printf("\n");

			for (i = 1; i <= m; i++) {
				for (j = 1; j <= n; j++) printf("%10.6f ", gaussai[i][j]);
				printf("\n");
			}
		}
		else {
			printf("---------- Determinant ----------\n");
			printf("\n");
			printf("%f \n", 0);
			printf("\n");

			printf("!!!!! Gauss-Jordan Elimination CAN NOT HANDLE A Singular Matirx !!!!!\n");
			
		}
		printf("\n");


		if (ludcmpE == 0) {
			printf("---------- Solution of LU Decomposition ----------\n");
			printf("\n");

			for (i = 1; i <= m; i++) {
				printf("%10.6f ", ludx[i]);
			}
			printf("\n");
			printf("\n");

			printf("---------- Inverse Matrix using LU Decomposition ----------\n");
			printf("\n");

			for (i = 1; i <= m; i++) {
				for (j = 1; j <= n; j++) printf("%10.6f ", ludai[i][j]);
				printf("\n");
			}
			printf("\n");

			printf("---------- Solution of Iterative Improvement ----------\n");
			printf("\n");

			for (i = 1; i <= m; i++) {
				printf("%10.6f ", mpvx[i]);
			}
			printf("\n");
		}
		else {
			printf("!!!!! LU Decomposition CAN NOT HANDLE A Singular Matirx !!!!!\n");
		}
		printf("\n");
		printf("---------- Solution of Singular Value Decomposition ----------\n");
		printf("\n");

		for (i = 1; i <= m; i++) {
			printf("%10.6f ", svdx[i][1]);
		}
		printf("\n");
		printf("\n");

		printf("---------- Inverse Matrix using Singular Value Decomposition ----------\n");
		printf("\n");


		for (i = 1; i <= m; i++) {
			for (j = 1; j <= n; j++) printf("%10.6f ", svdai[i][j]);
			printf("\n");

		}
		printf("\n");
	}
	free_matrix(a, 1, NP, 1, NP);
	free_matrix(b, 1, NP, 1, NP);
	free_matrix(lu, 1, NP, 1, NP);
	free_matrix(svdu, 1, NP, 1, NP);
	free_matrix(svdv, 1, NP, 1, NP);
	free_matrix(gaussai, 1, NP, 1, NP);
	free_matrix(ludai, 1, NP, 1, NP);
	free_matrix(svdai, 1, NP, 1, NP);
	free_matrix(gaussx, 1, NP, 1, NP);
	free_vector(ludx, 1, NP);
	free_vector(mpvx, 1, NP);
	free_matrix(svdx, 1, NP, 1, NP);
	free_ivector(indx, 1, NP);
	free_vector(col, 1, NP);
	free_vector(svdw, 1, NP);
	free_vector(tb, 1, NP);
	fclose(fp);
	return 0;
}




void gaussj(float **a, int n, float **b, int m, int* gaussE)
{
	int *indxc, *indxr, *ipiv;
	int i, icol, irow, j, k, l, ll;
	float big, dum, pivinv, temp;

	indxc = ivector(1, n);
	indxr = ivector(1, n);
	ipiv = ivector(1, n);
	for (j = 1; j <= n; j++) ipiv[j] = 0;

	for (i = 1; i <= n; i++) {
		big = 0.0;
		for (j = 1; j <= n; j++) {
			if (ipiv[j] != 1) {
				for (k = 1; k <= n; k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) { // problem

							big = fabs(a[j][k]);
							irow = j;
							icol = k;

						}
					}
				}
			}
		}
		++(ipiv[icol]);

		if (irow != icol) {
			for (l = 1; l <= n; l++) SWAP(a[irow][l], a[icol][l])
				for (l = 1; l <= m; l++) SWAP(b[irow][l], b[icol][l])
		}

		indxr[i] = irow;
		indxc[i] = icol;
		if (a[icol][icol] == 0.0) {
			free_ivector(ipiv, 1, n);
			free_ivector(indxr, 1, n);
			free_ivector(indxc, 1, n);

			//printf("gaussj: CAN NOT HANDLE Singular Matrix\n");
			(*gaussE) = 1;
			return;
			//nrerror("gaussj: Singular Matrix");
		}
		pivinv = 1.0 / a[icol][icol];
		a[icol][icol] = 1.0;
		for (l = 1; l <= n; l++) a[icol][l] *= pivinv;
		for (l = 1; l <= m; l++) b[icol][l] *= pivinv;
		for (ll = 1; ll <= n; ll++)
			if (ll != icol) {
				dum = a[ll][icol];
				a[ll][icol] = 0.0;
				for (l = 1; l <= n; l++) a[ll][l] -= a[icol][l] * dum;
				for (l = 1; l <= m; l++) b[ll][l] -= b[icol][l] * dum;
			}
	}
	for (l = n; l >= 1; l--) {
		if (indxr[l] != indxc[l])
			for (k = 1; k <= n; k++)
				SWAP(a[k][indxr[l]], a[k][indxc[l]]);
	}
	free_ivector(ipiv, 1, n);
	free_ivector(indxr, 1, n);
	free_ivector(indxc, 1, n);
}
void ludcmp(float **a, int n, int *indx, float *d, int* ludcmpE)
{
	int i, imax, j, k;
	float big, dum, sum, temp;
	float *vv;
	vv = vector(1, n); 
	*d = 1.0; 
	for (i = 1; i <= n; i++) {
		big = 0.0;

		for (j = 1; j <= n; j++)
			if ((temp = fabs(a[i][j])) > big) big = temp;
		if (big == 0.0) {
			free_vector(vv, 1, n);
			(*ludcmpE) = 1;
			//printf("Singular matrix in routine ludcmp\n");
			return;
			//nrerror("Singular matrix in routine ludcmp"); // Exit
		}
		vv[i] = 1.0 / big; // divided by biggest number
	}

	for (j = 1; j <= n; j++) {
		for (i = 1; i < j; i++) {
			sum = a[i][j];
			for (k = 1; k < i; k++) sum -= a[i][k] * a[k][j];
			a[i][j] = sum; // calculating l[i][j]
		}
		big = 0.0;

		for (i = j; i <= n; i++) {
			sum = a[i][j];
			for (k = 1; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum; // calculating u[i][j]
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 1; k <= n; k++) {
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j] == 0.0) { // if a[j][j] is 0.0, it is sigular matrix.(even if it is caused by machine accuracy, it almost sigular one.) 
			//FOUNDIT
			(*ludcmpE) = 1;
			a[j][j] = TINY;
		}
		if (j != n) {
			dum = 1.0 / (a[j][j]);
			for (i = j + 1; i <= n; i++) a[i][j] *= dum;
		}
	}
	free_vector(vv, 1, n);

}
void svdcmp(float **a, int m, int n, float w[], float **v)
{
	float pythag(float a, float b);
	int flag, i, its, j, jj, k, l, nm;
	float anorm, c, f, g, h, s, scale, x, y, z, *rv1;

	rv1 = vector(1, n);
	g = scale = anorm = 0.0;
	for (i = 1; i <= n; i++) {
		l = i + 1;
		rv1[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m) {
			for (k = i; k <= m; k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k = i; k <= m; k++) {
					a[k][i] /= scale;
					s += a[k][i] * a[k][i];
				}
				f = a[i][i];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][i] = f - g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = i; k <= m; k++) s += a[k][i] * a[k][j];
					f = s / h;
					for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
				}
				for (k = i; k <= m; k++) a[k][i] *= scale;
			}
		}
		w[i] = scale * g;
		g = s = scale = 0.0;
		if (i <= m && i != n) {
			for (k = l; k <= n; k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k = l; k <= n; k++) {
					a[i][k] /= scale;
					s += a[i][k] * a[i][k];
				}
				f = a[i][l];
				g = -SIGN(sqrt(s), f);
				h = f * g - s;
				a[i][l] = f - g;
				for (k = l; k <= n; k++) rv1[k] = a[i][k] / h;
				for (j = l; j <= m; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[j][k] * a[i][k];
					for (k = l; k <= n; k++) a[j][k] += s * rv1[k];
				}
				for (k = l; k <= n; k++) a[i][k] *= scale;
			}
		}
		anorm = FMAX(anorm, (fabs(w[i]) + fabs(rv1[i])));
	}
	for (i = n; i >= 1; i--) {
		if (i < n) {
			if (g) {
				for (j = l; j <= n; j++)
					v[j][i] = (a[i][j] / a[i][l]) / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= n; k++) s += a[i][k] * v[k][j];
					for (k = l; k <= n; k++) v[k][j] += s * v[k][i];
				}
			}
			for (j = l; j <= n; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i] = 1.0;
		g = rv1[i];
		l = i;
	}
	for (i = IMIN(m, n); i >= 1; i--) {
		l = i + 1;
		g = w[i];
		for (j = l; j <= n; j++) a[i][j] = 0.0;
		if (g) {
			g = 1.0 / g;
			for (j = l; j <= n; j++) {
				for (s = 0.0, k = l; k <= m; k++) s += a[k][i] * a[k][j];
				f = (s / a[i][i])*g;
				for (k = i; k <= m; k++) a[k][j] += f * a[k][i];
			}
			for (j = i; j <= m; j++) a[j][i] *= g;
		}
		else for (j = i; j <= m; j++) a[j][i] = 0.0;
		++a[i][i];
	}
	for (k = n; k >= 1; k--) {
		for (its = 1; its <= 30; its++) {
			flag = 1;
			for (l = k; l >= 1; l--) {
				nm = l - 1;
				if ((float)(fabs(rv1[l]) + anorm) == anorm) {
					flag = 0;
					break;
				}
				if ((float)(fabs(w[nm]) + anorm) == anorm) break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i = l; i <= k; i++) {
					f = s * rv1[i];
					rv1[i] = c * rv1[i];
					if ((float)(fabs(f) + anorm) == anorm) break;
					g = w[i];
					h = pythag(f, g);
					w[i] = h;
					h = 1.0 / h;
					c = g * h;
					s = -f * h;
					for (j = 1; j <= m; j++) {
						y = a[j][nm];
						z = a[j][i];
						a[j][nm] = y * c + z * s;
						a[j][i] = z * c - y * s;
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j = 1; j <= n; j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) {
				free_vector(rv1, 1, n);

				printf("no convergence in 30 svdcmp iterations\n");
				return;
				//nrerror("no convergence in 30 svdcmp iterations");
			}
			x = w[l];
			nm = k - 1;
			y = w[nm];
			g = rv1[nm];
			h = rv1[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0*h*y);
			g = pythag(f, 1.0);
			f = ((x - z)*(x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;
			c = s = 1.0;
			for (j = l; j <= nm; j++) {
				i = j + 1;
				g = rv1[i];
				y = w[i];
				h = s * g;
				g = c * g;
				z = pythag(f, h);
				rv1[j] = z;
				c = f / z;
				s = h / z;
				f = x * c + g * s;
				g = g * c - x * s;
				h = y * s;
				y *= c;
				for (jj = 1; jj <= n; jj++) {
					x = v[jj][j];
					z = v[jj][i];
					v[jj][j] = x * c + z * s;
					v[jj][i] = z * c - x * s;
				}
				z = pythag(f, h);
				w[j] = z;
				if (z) {
					z = 1.0 / z;
					c = f * z;
					s = h * z;
				}
				f = c * g + s * y;
				x = c * y - s * g;
				for (jj = 1; jj <= m; jj++) {
					y = a[jj][j];
					z = a[jj][i];
					a[jj][j] = y * c + z * s;
					a[jj][i] = z * c - y * s;
				}
			}
			rv1[l] = 0.0;
			rv1[k] = f;
			w[k] = x;
		}
	}
	free_vector(rv1, 1, n);
}
void mprove(float **a, float **alud, int n, int indx[], float b[], float x[])
{
	void lubksb(float **a, int n, int *indx, float b[]);
	int j, i;
	double sdp;
	float *r;

	r = vector(1, n);
	for (i = 1; i <= n; i++) {
		sdp = -b[i];
		for (j = 1; j <= n; j++) sdp += a[i][j] * x[j];
		r[i] = sdp;
	}
	lubksb(alud, n, indx, r);
	for (i = 1; i <= n; i++) x[i] -= r[i];
	free_vector(r, 1, n);
}
void lubksb(float **a, int n, int *indx, float b[])
{
	int i, ii = 0, ip, j;
	float sum;

	for (i = 1; i <= n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = ii; j <= i - 1; j++) sum -= a[i][j] * b[j];
		else if (sum) ii = i;
		b[i] = sum;
	}
	for (i = n; i >= 1; i--) {
		sum = b[i];
		for (j = i + 1; j <= n; j++) sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
	}
}
float pythag(float a, float b) {
	return sqrtf(a*a + b*b);
}

int *ivector(long nl, long nh) {
	return (int*)calloc(sizeof(int), (nh - nl + 2));
}
float *vector(long nl, long nh) {
	return (float*)calloc(sizeof(float), (nh - nl + 2));
}
int **imatrix(long nrl, long nrh, long ncl, long nch) {
	long row = nrh - nrl + 2;
	long col = nch - ncl + 2;
	int** mat = (int**)malloc(sizeof(int*)*row);

	for (int i = 0; i < row; i++) {
		mat[i] = (int*)malloc(sizeof(int)*col);
	}
	return mat;
}
float **matrix(long nrl, long nrh, long ncl, long nch) {
	long row = nrh - nrl + 2;
	long col = nch - ncl + 2;
	float** mat = (float**)malloc(sizeof(float*)*row);

	for (int i = 0; i < row; i++) {
		mat[i] = (float*)malloc(sizeof(float)*col);
	}
	return mat;
}

void free_vector(float *v, long nl, long nh) {
	free(v);
}
void free_ivector(int *v, long nl, long nh) {
	free(v);
}
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch) {
	long row = nrh - nrl + 2;
	long col = nch - ncl + 2;
	for (int i = 0; i < row; i++) {
		free(m[i]);
	}
	free(m);
}
void free_matrix(float **m, long nrl, long nrh, long ncl, long nch) {
	long row = nrh - nrl + 2;
	long col = nch - ncl + 2;
	for (int i = 0; i < row; i++) {
		free(m[i]);
	}
	free(m);
}
