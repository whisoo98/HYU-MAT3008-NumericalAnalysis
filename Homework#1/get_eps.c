#include <math.h>
#include <stdio.h>
//#define float double
void fget_eps(int* epsexp) {
	float one = ((float)(1));
	*epsexp = 0;
	float tmp = 1;
	float tmpplus = 1;
	while (one != tmp + tmpplus) {
		*epsexp = *epsexp + 1;
		tmpplus /= 2;
	}
}
void dget_eps(int* epsexp) {
	double one = ((double)(1));
	*epsexp = 0;
	double tmp = 1;
	double tmpplus = 1;
	while (one != tmp + tmpplus) {
		*epsexp = *epsexp + 1;
		tmpplus /= 2;
	}
}

 int main() {
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