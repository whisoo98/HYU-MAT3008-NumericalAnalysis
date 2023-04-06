
// Homework5View.cpp: CHomework5View 클래스의 구현
//

#include "pch.h"
#include "framework.h"
// SHARED_HANDLERS는 미리 보기, 축소판 그림 및 검색 필터 처리기를 구현하는 ATL 프로젝트에서 정의할 수 있으며
// 해당 프로젝트와 문서 코드를 공유하도록 해 줍니다.
#ifndef SHARED_HANDLERS
#include "Homework5.h"
#endif

#include "Homework5Doc.h"
#include "Homework5View.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#endif

//void polin2(int x1a[], int x2a[], long long **ya, int m, int n, float x1, float x2, long long *y, int *dy);
void polin2(int x1a[], int x2a[], COLORREF **ya, int m, int n, float x1, float x2, COLORREF *y, COLORREF *dy);
void Mypolin2(COLORREF **ya, int width, int height, float x, float y, COLORREF *ret);
// CHomework5View

IMPLEMENT_DYNCREATE(CHomework5View, CView)

BEGIN_MESSAGE_MAP(CHomework5View, CView)
	// 표준 인쇄 명령입니다.
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
END_MESSAGE_MAP()

// CHomework5View 생성/소멸

CHomework5View::CHomework5View() noexcept
{
	// TODO: 여기에 생성 코드를 추가합니다.

}

CHomework5View::~CHomework5View()
{
}

BOOL CHomework5View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: CREATESTRUCT cs를 수정하여 여기에서
	//  Window 클래스 또는 스타일을 수정합니다.

	return CView::PreCreateWindow(cs);
}

// CHomework5View 그리기

void CHomework5View::OnDraw(CDC* pDC)
{
	CHomework5Doc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	int NumArg;
	LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(),&NumArg);
	CString filename=argv[1];
	CImage image, image_;
	
	image.Load(filename); // just change extension to load jpg

	int n = image.GetWidth(); // Width of Pixel  //35
	int m = image.GetHeight(); // Height of Pixel  //39
	
	CString temp;

	int n_ = _wtoi(argv[2]); // Width
	int m_ = _wtoi(argv[3]); // Height
	
	float ratioX = (float)((float)n / (float)n_);
	float ratioY = (float)((float)m / (float)m_);
	
	int* x1a = (int*)malloc(sizeof(int)*(n + 1)); // what is in x1a? => list[0...n+1]
	int* x2a = (int*)malloc(sizeof(int)*(m + 1)); // what is in x2a? => list[0...m+1]
	
	COLORREF** ya = (COLORREF**)malloc((sizeof(COLORREF*)*(m + 1)));
	COLORREF** ty = (COLORREF**)malloc((sizeof(COLORREF*)*(m_ + 1)));
	for (int i = 0; i <= n; i++) {
		x1a[i] = i;
	}
	for (int i = 0; i <= m; i++) {
		x2a[i] = i;
		ya[i] = (COLORREF*)malloc((COLORREF)(sizeof(COLORREF)*(n + 1)));
	}

	for (int i = 0; i <= m_; i++) {
		ty[i] = (COLORREF*)malloc((COLORREF)(sizeof(COLORREF)*(n_ + 1)));
	}
	image_.Create(n_, m_, image.GetBPP());

	//for (int i = 1; i <= m; i++) for (int j = 1; j <= n; j++) ya[i][j] = image.GetPixel(x1a[j - 1], x2a[i - 1]);
	for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) {
		ya[i][j] = image.GetPixel(x1a[j], x2a[i]); 
		//rgb[i][j] = RGB(GetRValue(ya[i][j]), GetGValue(ya[i][j]), GetBValue(ya[i][j]));
		
	}
	
	//image_.BitBlt(dc.m_hDC, n + 100, 0); // express interpolating Image*/


	//for (int i = 0; i < m_; i++) for (int j = 0; j < n_; j++) ty[i][j] = (COLORREF)0;
	//pDC->TextOutW(300, 300, filename);
	image.Detach();

	// What I want to do ::
	// To get COLOR of pixel (i, j) on image_ that is interpolating image,
	// I make (x2, x1) using (i, j) that is to be interpolated point.
	// When I am doing interpolation using (x2, x1), then I can get interpolated COLOR that is on (i, j). 
	
	//for (int i = 1; i <= m_; i++) {
	for (int i = 0; i < m_; i++) {
		//for (int j = 1; j <= n_; j++) {
		for (int j = 0; j < n_; j++) {
			//float x1 = (float)(j+1) * ratioX;
			float x1 = (float)(j) * ratioX;
			//float x2 = (float)(i+1) * ratioY;
			float x2 = (float)(i) * ratioY;
			//COLORREF del;
			//polin2(x2a, x1a, ya, m, n, x2, x1, &ty[i][j], &del); 
			// This way is too slow because of it's time complexity
			// O(image_.width*image_.height) for Iterating Pixel of image
			// O(image_height) for polin2
			// O(image_width^2) for polint
			// Total Time Complexity is O(n^5)!!!
			// It means that if I put 50*50 image file, roughly said, then it'll be (50^5) operations.
			// It's 312'500'000 ...
			// So I decide to make own polin2
			// It is not an option. It is a must.
			Mypolin2(ya, n_, m_, x1, x2, &ty[i][j]);

			image_.SetPixel(j, i, ty[i][j]);
			//pDC->SetPixel(j, i, ty[i][j]);
			// I also want to override SetPixel function for more fast speed, but I'll put this aside.
		}
	}
	image_.BitBlt(pDC->m_hDC, 0, 0);
	
	
	int pos = filename.ReverseFind(_T('.'));
	CString extension;
	extension.Format(_T("Output%s"), filename.Right(filename.GetLength() - pos));

	image_.Save(extension);
	image_.Detach();
	// TODO: 여기에 원시 데이터에 대한 그리기 코드를 추가합니다.
}


// CHomework5View 인쇄

BOOL CHomework5View::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 기본적인 준비
	return DoPreparePrinting(pInfo);
}

void CHomework5View::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 인쇄하기 전에 추가 초기화 작업을 추가합니다.
}

void CHomework5View::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 인쇄 후 정리 작업을 추가합니다.
}


// CHomework5View 진단

#ifdef _DEBUG
void CHomework5View::AssertValid() const
{
	CView::AssertValid();
}

void CHomework5View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CHomework5Doc* CHomework5View::GetDocument() const // 디버그되지 않은 버전은 인라인으로 지정됩니다.
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CHomework5Doc)));
	return (CHomework5Doc*)m_pDocument;
}
#endif //_DEBUG


// CHomework5View 메시지 처리기
int *vector(int nl, int nh)
{
	int *v;

	v = (int *)malloc((unsigned int)((nh - nl + 1 + 1) * sizeof(int)));
	//	if (!v) nrerror("allocation failure in vector()");
	//return v - nl + 1;
	return v;
}

long long *Lvector(int nl, int nh)
{
	long long *v;

	v = (long long *)malloc((unsigned long long)((nh - nl + 1 + 1) * sizeof(long long)));
	//	if (!v) nrerror("allocation failure in vector()");
	//return v - nl + 1;
	return v;
}

COLORREF *Cvector(int nl, int nh)
{
	COLORREF *v;

	v = (COLORREF*)malloc(((nh - nl + 1 + 1) * sizeof(COLORREF)));
	//	if (!v) nrerror("allocation failure in vector()");
	//return v - nl + 1;
	return v;
}
void free_vector(int *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	//free((char*)(v + nl - 1));
	free(v);
}
void free_Lvector(long long *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	//free((char*)(v + nl - 1));
	free(v);
}

void free_Cvector(COLORREF *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	//free((char*)(v + nl - 1));
	free(v);
}

void polint(int xa[], COLORREF ya[], int n, float x, COLORREF *y, COLORREF *dy)
{
	int i, m, ns = 0;
	float dif, dift, ho, hp, den;
	COLORREF w; 
	COLORREF *c, *d;
	//dif = fabs(x - (float)xa[1]); //xa와 x의 차이 => 좌표 사이의 거리
	dif = fabs(x - (float)xa[0]); //xa와 x의 차이 => 좌표 사이의 거리
	c = Cvector(1, n);
	d = Cvector(1, n);
	//for (i = 1; i <= n; i++) { // base 1 !! Important !!
	for (i = 0; i < n; i++) { // base 1 !! Important !!
		if ((dift = fabs(x - (float)xa[i])) < dif) { // dif를 minmun diff로 만듦
			ns = i;
			dif = dift; 
		}
		c[i] = ya[i]; // ya저장 => pixel 값
		d[i] = ya[i]; // ya저장 => pixel 값
	}
	*y = ya[ns--]; //y의 값은 dif인 index의 값 & 해당 index를 1줄임
	//for (m = 1; m < n; m++) {
	for (m = 0; m < n-1; m++) {
		//for (i = 1; i <= n - m; i++) {
		for (i = 0; i < n - m; i++) {
			ho = (float)xa[i] - x;
			hp = (float)xa[i + m] - x;
			w = c[i + 1] - d[i];
			if ((den = ho - hp) == 0.0) {
				printf("Error in routine polint\n");
			}
			den = (float)w / den;
			d[i] = (COLORREF)(hp * den);
			c[i] = (COLORREF)(ho * den);
		}
		*y += (*dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
	}
	free_Cvector(d, 1, n);
	free_Cvector(c, 1, n);
}


void polin2(int x1a[], int x2a[], COLORREF **ya, int m, int n, float x1, float x2, COLORREF *y, COLORREF *dy)
{   //      array of indexes   |array of pixel |# rows & cols |interpolating points |return values
	void polint(int xa[], COLORREF ya[], int n, float x, COLORREF *y, COLORREF *dy);
	int j;
	COLORREF *ymtmp;

	//ymtmp = Lvector(1, m);
	ymtmp = Cvector(1, m);
	//for (j = 1; j <= m; j++) polint(x2a, ya[j], n, x2, &ymtmp[j], dy);
	for (j = 0; j < m; j++) polint(x2a, ya[j], n, x2, &ymtmp[j], dy);
	
	polint(x1a, ymtmp, m, x1, y, dy);
	free_Cvector(ymtmp, 1, m);
}

void Mypolin2(COLORREF **ya, int width, int height, float x, float y, COLORREF *ret) {
	int l, r, u, d;
	l = floor(x); // left of x
	u = floor(y); // up of y
	r = ceil(x);  // right of x
	d = ceil(y);  // down of y
	
	int m, n;
	m = abs(x - l); // difference between l and x cf) abs doesn't need
	n = abs(y - u); // difference between d and y cf) abs doesn't need
	
	int p1, p2, p3, p4;
	p1 = (1 - m)*(1 - n); // weight for ul
	p2 = (m)*(1 - n);     // weight for ur
	p3 = (1 - m)*(n);     // weight for dl
	p4 = (m)*(n);         // weight or dr
	
	*ret = 0;
	*ret += p1 * ya[u][l];
	*ret += p2 * ya[u][r];
	*ret += p3 * ya[d][l];
	*ret += p4 * ya[d][r];
}