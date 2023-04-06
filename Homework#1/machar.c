#include <stdio.h>
#include <math.h>
//#define float double
#define CONVF(i) ((float)(i))
#define CONVD(i) ((double)(i))
/*
ibeta : the floating-point radix
it : the number of base-ibeta digis in the floating point mantissa
irnd : round style
ngrd : guard digits
machep : epsilon power
negep : neg epsilon power
iexp : exponent digits
minexp : min exponemt
maxexp : max exponent
eps : the smallest positive number that, added to 1.0, is not equal to 1.0
epsneg  the smallest positive number that, subtracted from 1.0, is not equal to 1.0
xmin : the smallest representable positive number
xman : the largest representable positive number
*/

void fmachar(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, float *eps, float *epsneg,
	float *xmin, float *xmax)
{
	int i, itemp, iz, j, k, mx, nxres;
	float a, b, beta, betah, betain, one, t, temp, temp1, tempa, two, y, z, zero;

	one = CONVF(1);
	two = one + one;
	zero = one - one;
	a = one;
	do {
		a += a;
		temp = a + one;
		temp1 = temp - a;
	} while (temp1 - one == zero);
	b = one;
	do {
		b += b;
		temp = a + b;
		itemp = (int)(temp - a);
	} while (itemp == 0);
	*ibeta = itemp;
	beta = CONVF(*ibeta);
	*it = 0;
	b = one;
	do {
		++(*it);
		b *= beta;
		temp = b + one;
		temp1 = temp - b;
	} while (temp1 - one == zero);
	*irnd = 0;
	betah = beta / two;
	temp = a + betah;
	if (temp - a != zero) *irnd = 1;
	tempa = a + beta;
	temp = tempa + betah;
	if (*irnd == 0 && temp - tempa != zero) *irnd = 2;
	*negep = (*it) + 3;
	betain = one / beta;
	a = one;
	for (i = 1; i <= (*negep); i++) a *= betain;
	b = a;
	for (;;) {
		temp = one - a;
		if (temp - one != zero) break;
		a *= beta;
		--(*negep);
	}
	*negep = -(*negep);
	*epsneg = a;
	*machep = -(*it) - 3;
	a = b;
	for (;;) {
		temp = one + a;
		if (temp - one != zero) break;
		a *= beta;
		++(*machep);
	}
	*eps = a;
	*ngrd = 0;
	temp = one + (*eps);
	if (*irnd == 0 && temp*one - one != zero) *ngrd = 1;
	i = 0;
	k = 1;
	z = betain;
	t = one + (*eps);
	nxres = 0;
	for (;;) {
		y = z;
		z = y * y;
		a = z * one;
		temp = z * t;
		if (a + a == zero || fabs(z) >= y) break;
		temp1 = temp * betain;
		if (temp1*beta == z) break;
		++i;
		k += k;
	}
	if (*ibeta != 10) {
		*iexp = i + 1;
		mx = k + k;
	}
	else {
		*iexp = 2;
		iz = (*ibeta);
		while (k >= iz) {
			iz *= *ibeta;
			++(*iexp);
		}
		mx = iz + iz - 1;
	}
	for (;;) {
		*xmin = y;
		y *= betain;
		a = y * one;
		temp = y * t;
		if (a + a != zero && fabs(y) < *xmin) {
			++k;
			temp1 = temp * betain;
			if (temp1*beta == y && temp != y) {
				nxres = 3;
				*xmin = y;
				break;
			}
		}
		else break;
	}
	*minexp = -k;
	if (mx <= k + k - 3 && *ibeta != 10) {
		mx += mx;
		++(*iexp);
	}
	*maxexp = mx + (*minexp);
	*irnd += nxres;
	if (*irnd >= 2) *maxexp -= 2;
	i = (*maxexp) + (*minexp);
	if (*ibeta == 2 && !i) --(*maxexp);
	if (i > 20) --(*maxexp);
	if (a != y) *maxexp -= 2;
	*xmax = one - (*epsneg);
	if ((*xmax)*one != *xmax) *xmax = one - beta * (*epsneg);
	*xmax /= (*xmin*beta*beta*beta);
	i = (*maxexp) + (*minexp) + 3;
	for (j = 1; j <= i; j++) {
		if (*ibeta == 2) *xmax += *xmax;
		else *xmax *= beta;
	}
}

void dmachar(int *ibeta, int *it, int *irnd, int *ngrd, int *machep, int *negep,
	int *iexp, int *minexp, int *maxexp, double *eps, double *epsneg,
	double *xmin, double *xmax)
{
	int i, itemp, iz, j, k, mx, nxres;
	double a, b, beta, betah, betain, one, t, temp, temp1, tempa, two, y, z, zero;

	one = CONVD(1);
	two = one + one;
	zero = one - one;
	a = one;
	do {
		a += a;
		temp = a + one;
		temp1 = temp - a;
	} while (temp1 - one == zero);
	b = one;
	do {
		b += b;
		temp = a + b;
		itemp = (int)(temp - a);
	} while (itemp == 0);
	*ibeta = itemp;
	beta = CONVD(*ibeta);
	*it = 0;
	b = one;
	do {
		++(*it);
		b *= beta;
		temp = b + one;
		temp1 = temp - b;
	} while (temp1 - one == zero);
	*irnd = 0;
	betah = beta / two;
	temp = a + betah;
	if (temp - a != zero) *irnd = 1;
	tempa = a + beta;
	temp = tempa + betah;
	if (*irnd == 0 && temp - tempa != zero) *irnd = 2;
	*negep = (*it) + 3;
	betain = one / beta;
	a = one;
	for (i = 1; i <= (*negep); i++) a *= betain;
	b = a;
	for (;;) {
		temp = one - a;
		if (temp - one != zero) break;
		a *= beta;
		--(*negep);
	}
	*negep = -(*negep);
	*epsneg = a;
	*machep = -(*it) - 3;
	a = b;
	for (;;) {
		temp = one + a;
		if (temp - one != zero) break;
		a *= beta;
		++(*machep);
	}
	*eps = a;
	*ngrd = 0;
	temp = one + (*eps);
	if (*irnd == 0 && temp*one - one != zero) *ngrd = 1;
	i = 0;
	k = 1;
	z = betain;
	t = one + (*eps);
	nxres = 0;
	for (;;) {
		y = z;
		z = y * y;
		a = z * one;
		temp = z * t;
		if (a + a == zero || fabs(z) >= y) break;
		temp1 = temp * betain;
		if (temp1*beta == z) break;
		++i;
		k += k;
	}
	if (*ibeta != 10) {
		*iexp = i + 1;
		mx = k + k;
	}
	else {
		*iexp = 2;
		iz = (*ibeta);
		while (k >= iz) {
			iz *= *ibeta;
			++(*iexp);
		}
		mx = iz + iz - 1;
	}
	for (;;) {
		*xmin = y;
		y *= betain;
		a = y * one;
		temp = y * t;
		if (a + a != zero && fabs(y) < *xmin) {
			++k;
			temp1 = temp * betain;
			if (temp1*beta == y && temp != y) {
				nxres = 3;
				*xmin = y;
				break;
			}
		}
		else break;
	}
	*minexp = -k;
	if (mx <= k + k - 3 && *ibeta != 10) {
		mx += mx;
		++(*iexp);
	}
	*maxexp = mx + (*minexp);
	*irnd += nxres;
	if (*irnd >= 2) *maxexp -= 2;
	i = (*maxexp) + (*minexp);
	if (*ibeta == 2 && !i) --(*maxexp);
	if (i > 20) --(*maxexp);
	if (a != y) *maxexp -= 2;
	*xmax = one - (*epsneg);
	if ((*xmax)*one != *xmax) *xmax = one - beta * (*epsneg);
	*xmax /= (*xmin*beta*beta*beta);
	i = (*maxexp) + (*minexp) + 3;
	for (j = 1; j <= i; j++) {
		if (*ibeta == 2) *xmax += *xmax;
		else *xmax *= beta;
	}
}


int main() {
	int ibeta, it, irnd, ngrd, machep, negep;
	int iexp, minexp, maxexp;
	float eps, epsneg, xmin, xmax;
	fmachar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &eps, &epsneg, &xmin, &xmax);
	printf("-----------------------------method1---------------------\n");
	printf("-----------------------------float---------------------\n");
	printf("ibeta : %d\n", ibeta);
	printf("it : %d\n", it);
	printf("irnd : %d\n", irnd);
	printf("ngrd : %d\n", ngrd);
	printf("machep : %d\n", machep);
	printf("negep : %d\n", negep);
	printf("iexp : %d\n", iexp);
	printf("minexp : %d\n", minexp);
	printf("maxexp : %d\n", maxexp);
	printf("eps : %.60f\n", eps); // This is the answer!
	printf("epsneg : %f\n", epsneg);
	printf("xmin : %f\n", xmin);
	printf("xmax : %f\n", xmax);
	printf("\n");
	printf("-----------------------------double---------------------\n");
	double deps, depsneg, dxmin, dxmax;
	dmachar(&ibeta, &it, &irnd, &ngrd, &machep, &negep, &iexp, &minexp, &maxexp, &deps, &depsneg, &dxmin, &dxmax);
	printf("ibeta : %d\n", ibeta);
	printf("it : %d\n", it);
	printf("irnd : %d\n", irnd);
	printf("ngrd : %d\n", ngrd);
	printf("machep : %d\n", machep);
	printf("negep : %d\n", negep);
	printf("iexp : %d\n", iexp);
	printf("minexp : %d\n", minexp);
	printf("maxexp : %d\n", maxexp);
	printf("eps : %.60f\n", deps); // This is the answer!
	printf("epsneg : %f\n", depsneg);
	printf("xmin : %f\n", dxmin);
	printf("xmax : %f\n", dxmax);
	printf("\n");
	printf("-----------------------------method2---------------------\n");
	int epsexp;
 	fget_eps(&epsexp);
	printf("-----------------------------float---------------------\n");

 	printf("epsexp : %d\n", epsexp);
 	float FMA = 2;
 	for (int i = 1; i <= epsexp; i++) {
 		FMA /= 2;
 	}
 	printf("eps : %.60f\n", FMA);
	dget_eps(&epsexp);
	double DMA = 2;
	for (int i = 1; i <= epsexp; i++) {
		DMA /= 2;
	}
	printf("-----------------------------double---------------------\n");
	printf("epsexp : %d\n", epsexp);
	printf("eps : %.60f\n", DMA);

	return 0;
}