
//#include <math.h>
#include <math.h>
#include <stdio.h>
#include <time.h>

float bessj1(float x) // Bessel1 Function
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = x * (72362614232.0 + y * (-7895059235.0 + y * (242396853.1
			+ y * (-2972611.439 + y * (15704.48260 + y * (-30.16036606))))));
		ans2 = 144725228442.0 + y * (2300535178.0 + y * (18583304.74
			+ y * (99447.43394 + y * (376.9991397 + y * 1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 2.356194491;
		ans1 = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4
			+ y * (0.2457520174e-5 + y * (-0.240337019e-6))));
		ans2 = 0.04687499995 + y * (-0.2002690873e-3
			+ y * (0.8449199096e-5 + y * (-0.88228987e-6
				+ y * 0.105787412e-6)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z * sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
	return ans;
}

float bessj0(float x) // Bessel0 Function
{
	float ax, z;
	double xx, y, ans, ans1, ans2;

	if ((ax = fabs(x)) < 8.0) {
		y = x * x;
		ans1 = 57568490574.0 + y * (-13362590354.0 + y * (651619640.7
			+ y * (-11214424.18 + y * (77392.33017 + y * (-184.9052456)))));
		ans2 = 57568490411.0 + y * (1029532985.0 + y * (9494680.718
			+ y * (59272.64853 + y * (267.8532712 + y * 1.0))));
		ans = ans1 / ans2;
	}
	else {
		z = 8.0 / ax;
		y = z * z;
		xx = ax - 0.785398164;
		ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
			+ y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
		ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
			+ y * (-0.6911147651e-5 + y * (0.7621095161e-6
				- y * 0.934945152e-7)));
		ans = sqrt(0.636619772 / ax)*(cos(xx)*ans1 - z * sin(xx)*ans2);
	}
	return ans;
}

float MyFunction(float x) // (x-2)(x-5)(x-8) Function
{
	float ans = (x - 2)*(x - 5)*(x - 8);
	return ans;
}

float DMyFunction(float x) // (x-2)(x-5)(x-8) Function
{
	float ans = (x - 5)*(x - 8) + (x-2)*(x-8) + (x-2)*(x-5);
	return ans;
}

#define JMAX 40
float rtbis(float(*func)(float), float x1, float x2, float xacc, int* cnt) // Bisection
{
	int j;
	float dx, f, fmid, xmid, rtb;

	f = (*func)(x1); // 초기 경계함수값
	fmid = (*func)(x2); // 초기 경계함수값
	if (f*fmid >= 0.0) { // 두 부호가 같으면 해가 존재 X
		printf("Root must be bracketed for bisection in rtbis\n");
	}
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		(*cnt)++; // 함수의 반복횟수를 확인하기 위한 변수
		fmid = (*func)(xmid = rtb + (dx *= 0.5)); // Bisection Part - 중간의 함수값을 확인
		if (fmid <= 0.0) rtb = xmid; // 함수값의 부호에 따라서 bracket을 줄여나감
		if (fabs(dx) < xacc || fmid == 0.0) return rtb; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
	}
	printf("Too many bisections in rtbis\n");
	return 0.0;
}
#undef JMAX

#define MAXIT 30
float rtflsp(float(*func)(float), float x1, float x2, float xacc, int* cnt) // Linear interpolation
{
	int j;
	float fl, fh, xl, xh, swap, dx, del, f, rtf;

	fl = (*func)(x1); // 왼쪽 함수값
	fh = (*func)(x2); // 오른쪽 함수값
	if (fl*fh > 0.0) printf("Root must be bracketed in rtflsp\n"); // 두 함수값의 부호가 같으면 해 존재하지 않음
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xl = x2;
		xh = x1;
		swap = fl;
		fl = fh;
		fh = swap;
	}
	dx = xh - xl; // 오차
	for (j = 1; j <= MAXIT; j++) {
		(*cnt)++; // 함수의 반복횟수를 확인하기 위한 변수
		rtf = xl + dx * fl / (fl - fh); // x_n+1
		f = (*func)(rtf); // f(x_n+1)
		if (f < 0.0) { // 새로운 x_n 설정
			del = xl - rtf; 
			xl = rtf;
			fl = f;
		}
		else {// 새로운 x_n 설정 -> 함수의 우측 수정
			del = xh - rtf;
			xh = rtf;
			fh = f;
		}
		dx = xh - xl;
		if (fabs(del) < xacc || f == 0.0) return rtf; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
	}
	printf("Maximum number of iterations exceeded in rtflsp\n");
	return 0.0;
}
#undef MAXIT

#define JMAX 20
float rtnewt(float(*func)(float), float(*funcd)(float), float x1, float x2,
	float xacc, int* cnt) // Newton-Raphson
{
	int j;
	float df, dx, f, rtn;
	int chk = 0;
	rtn = 0.5*(x1 + x2);
	for (j = 1; j <= JMAX; j++) {
		(*cnt)++; // 함수의 반복횟수를 확인하기 위한 변수
		f = (*func)(rtn); // f(x_n)
		df = -(*funcd)(rtn); // f'(x_n)

		dx = f / df;
		rtn -= dx; // x_n+1
		if ((x1 - rtn)*(rtn - x2) < 0.0) {// 정해진 구간 밖으로 벗어남
			printf("Jumped out of brackets in rtnewt\n");
		}
		if (fabs(dx) < xacc) {
			return rtn; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
		}
	}
	printf("Maximum number of iterations exceeded in rtnewt\n");
	return 0.0;
}
#undef JMAX

#define MAXIT 30
float rtsec(float(*func)(float), float x1, float x2, float xacc, int* cnt) // Secant method
{
	int j;
	float fl, f, dx, swap, xl, rts;

	fl = (*func)(x1); // f(p_n-1)의 초기값
	f = (*func)(x2); // f(p_n)의 초기값
	if (fabs(fl) < fabs(f)) {
		rts = x1;
		xl = x2;
		swap = fl;
		fl = f;
		f = swap;
	}
	else {
		xl = x1;
		rts = x2;
	}
	for (j = 1; j <= MAXIT; j++) {
		(*cnt)++; // 함수의 반복횟수를 확인하기 위한 변수
		dx = (xl - rts)*f / (f - fl);
		xl = rts;
		fl = f; // f(p_n-1)
		rts += dx; // p_n+1
		f = (*func)(rts); // f(p_n)
		if (fabs(dx) < xacc || f == 0.0) return rts; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
	}
	printf("Maximum number of iterations exceeded in rtsec\n");
	return 0.0;
}
#undef MAXIT

#define MAXIT 100
float rtsafe(float(*func)(float), float(*funcd)(float), float x1, float x2,
	float xacc,int* cnt) // Newton method with bracketing
{
	int j;
	float df, dx, dxold, f, fh, fl;
	float temp, xh, xl, rts;

	fl = (*func)(x1); // f(x1)
	df = -(*funcd)(x1); // f'(x1)
	fh = (*func)(x2); // f(x2)
	df = -(*funcd)(x2); // f'(x2)
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		printf("Root must be bracketed in rtsafe\n");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl = x1;
		xh = x2;
	}
	else {
		xh = x1;
		xl = x2;
	}
	rts = 0.5*(x1 + x2);
	dxold = fabs(x2 - x1);
	dx = dxold;
	f = (*func)(rts); // f(x_n)
	df = -(*funcd)(rts); // f'(x_n)
	for (j = 1; j <= MAXIT; j++) {
		(*cnt)++; // 함수의 반복횟수를 확인하기 위한 변수
		if ((((rts - xh)*df - f)*((rts - xl)*df - f) > 0.0) // Bisection Part - 중간의 함수값을 확인하여 bracket을 줄임
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold = dx;
			dx = 0.5*(xh - xl);
			rts = xl + dx;
			if (xl == rts) return rts; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
		}
		else {
			dxold = dx;
			dx = f / df;
			temp = rts;
			rts -= dx;
			if (temp == rts) return rts; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
		}
		if (fabs(dx) < xacc) return rts; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
		f = (*func)(rts);
		df = -(*funcd)(rts);
		if (f < 0.0) // Newton-Rappson Part - 함수값의 부호를 확인하여 x_n+1 갱신
			xl = rts;
		else
			xh = rts;
	}
	printf("Maximum number of iterations exceeded in rtsafe\n");
	return 0.0;
}
#undef MAXIT

#define MAXIT 30
float muller(float(*func)(float), float x1, float x2,
	float xacc, int* cnt) { // muller 

	float p0, p1, p2;
	float fpn2, fpn1, fpn0;
	float diff02, diff12, diff01;
	float a, b, c;
	float dx;


	p0 = x1;
	p1 = (x1 + x2) / 2;
	p2 = x2;

	fpn2 = (*func)(p2); // f(p2) -> f(p_n+1)
	fpn1 = (*func)(p1); // f(p1) -> f(p_n)
	fpn0 = (*func)(p0); // f(p0) -> f(p_n-1)

	diff02 = p0 - p2; // p0 - p2 -> p_n-1 - p_n+1
	diff12 = p1 - p2; // p1 - p2 -> p_n - p_n+1
	diff01 = p0 - p1; // p0 - p1 -> p_n-1 - p_n
	for (int j = 1; j <= MAXIT; j++) {
		float tmp;
		(*cnt)++; // 함수의 반복횟수를 확인하기 위한 변수
		c = fpn2; // c
		b = (diff02*diff02*(fpn1 - fpn2) - diff12 * diff12*(fpn0 - fpn2)) / (diff02*diff12*diff01); // b
		a = (diff12*(fpn0 - fpn2) - diff02 * (fpn1 - fpn2)) / (diff02*diff12*diff01); // a
		dx = 2 * c / (b + (b / fabs(b)*sqrtf(b*b - 4 * a*c))); // dx
		tmp = p2 - dx; // p_n+2

		p0 = p1;
		p1 = p2;
		p2 = tmp;

		fpn2 = (*func)(p2);
		fpn1 = (*func)(p1);
		fpn0 = (*func)(p0);

		diff02 = p0 - p2;
		diff12 = p1 - p2;
		diff01 = p0 - p1;
		// update

		if (fabs(dx) < xacc || fpn2 == 0.0) return p2; // 해를 찾았거나 오차가 xacc이하이면 함수 종료
	}
	printf("Maximum number of iterations exceeded in rtsec\n");

}
#undef MAXIT

void zbrak(float(*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb) // bracketing routine
{
	int nbb, i;
	float x, fp, fc, dx;	

	nbb = 0;
	dx = (x2 - x1) / n;
	fp = (*fx)(x = x1);
	for (i = 1; i <= n; i++) {
		fc = (*fx)(x += dx);
		if (fc*fp <= 0.0) {
			xb1[++nbb] = x - dx;
			xb2[nbb] = x;
			if (*nb == nbb) return;

		}
		fp = fc;
	}
	*nb = nbb;
}	




//2 3
//5 6
//8 9
int main() {
	float(*bessel0)(float);
	bessel0 = &bessj0;

	float(*bessel1)(float);
	bessel1 = &bessj1;

	float xacc = 0.000001;
	int nb=10;
	float arr[10] = { 0,0,0,0,0,0,0,0,0,0 };
	float brr[10] = { 0,0,0,0,0,0,0,0,0,0 };
	float crr[3];
	float drr[3];
	float ans;

	int cntbreak=9;
	zbrak(bessel0, 1.0, 10.0, cntbreak, arr, brr, &nb);
	
	
	for (int i = 0, j=0; i < 10; i++) {
		if (arr[i] * brr[i]) {
			crr[j] = arr[i];
			drr[j++] = brr[i];
			printf("Bracket of Root # %d : [%f , %f]\n", j, crr[j-1],drr[j-1]);
		}
	}
	printf("---------------------------------------------------------------\n");
	for (int i = 0; i < 3; i++) {
		int cnt=0;
		ans = rtbis(bessel0, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Bisection # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);
	}
	printf("---------------------------------------------------------------\n");
	
	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtflsp(bessel0, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Linear interpolation # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");
	
	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtsec(bessel0, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Secant # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");

	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtnewt(bessel0, bessel1, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Newton-Raphson # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");

	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtsafe(bessel0, bessel1, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Newton with bracketing # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");

	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = muller(bessel0, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Muller # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");
	printf("---------------------------------------------------------------\n");
	printf("-------------------For My New Function-------------------------\n");


	float(*my_func)(float);
	my_func = &MyFunction;
	
	float(*dmy_func)(float);
	dmy_func = &DMyFunction;

	
	cntbreak = 70;
	zbrak(my_func, 1.0, 10.0, cntbreak, arr, brr, &nb);


	for (int i = 0, j = 0; i < 10; i++) {
		if (arr[i] * brr[i]) {
			crr[j] = arr[i];
			drr[j++] = brr[i];
			printf("Bracket of Root # %d : [%f , %f]\n", j, crr[j - 1], drr[j - 1]);
		}
	}

	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtsafe(my_func, dmy_func, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Newton with bracketing # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");

	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtbis(my_func, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Bisection # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);
	}
	printf("---------------------------------------------------------------\n");

	for (int i = 0; i < 3; i++) {
		int cnt = 0;
		ans = rtnewt(my_func, dmy_func, crr[i], drr[i], xacc, &cnt);
		printf("Roots Using Newton-Raphson # %d : %f", i + 1, ans);
		printf("\tThe Number of Iteration %d\n", cnt);

	}
	printf("---------------------------------------------------------------\n");
	return 0;
}