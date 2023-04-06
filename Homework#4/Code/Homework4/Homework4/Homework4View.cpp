
// Homework4View.cpp: CHomework4View 클래스의 구현
//

#include "pch.h"
#include "framework.h"
// SHARED_HANDLERS는 미리 보기, 축소판 그림 및 검색 필터 처리기를 구현하는 ATL 프로젝트에서 정의할 수 있으며
// 해당 프로젝트와 문서 코드를 공유하도록 해 줍니다.
#ifndef SHARED_HANDLERS
#include "Homework4.h"
#endif

#include "Homework4Doc.h"
#include "Homework4View.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CHomework4View

IMPLEMENT_DYNCREATE(CHomework4View, CView)

BEGIN_MESSAGE_MAP(CHomework4View, CView)
	// 표준 인쇄 명령입니다.
	ON_COMMAND(ID_FILE_PRINT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_DIRECT, &CView::OnFilePrint)
	ON_COMMAND(ID_FILE_PRINT_PREVIEW, &CView::OnFilePrintPreview)
END_MESSAGE_MAP()

// CHomework4View 생성/소멸

CHomework4View::CHomework4View() noexcept
{
	// TODO: 여기에 생성 코드를 추가합니다.

}

CHomework4View::~CHomework4View()
{
}

BOOL CHomework4View::PreCreateWindow(CREATESTRUCT& cs)
{
	// TODO: CREATESTRUCT cs를 수정하여 여기에서
	//  Window 클래스 또는 스타일을 수정합니다.

	return CView::PreCreateWindow(cs);
}

// CHomework4View 그리기

void CHomework4View::OnDraw(CDC* pDC)
{
	CHomework4Doc* pDoc = GetDocument();
	ASSERT_VALID(pDoc);
	if (!pDoc)
		return;

	CPen penAxis(PS_DOT, 2, RGB(0, 0, 255));
	pDC->SelectObject(&penAxis);
	int x0 = 20;
	int y0_gauss =500;
	int y0_uniform = 850;
	pDC->MoveTo(x0, y0_uniform);
	pDC->LineTo(2000, y0_uniform);

	pDC->MoveTo(x0, y0_gauss);
	pDC->LineTo(2000, y0_gauss);

	pDC->MoveTo(x0, 20);
	pDC->LineTo(x0, y0_gauss);

	pDC->MoveTo(x0, y0_gauss+50);
	pDC->LineTo(x0, y0_uniform);


	//int Uniform;
	int gentimes[] = { 100,1000,10000,100000 }; // How many samples to generate?
	double resize[] = { 0.01,0.25,1,10 }; // To fit to the screen
	int k = 1; // To Select gentimes and resize

	std::random_device rd;
	std::mt19937 gen(rd()); // Gen a random number
	int st = -3, ed = 2; // range of Uniform Dist
	double mean = 0.5, stddev = 1.5; // Params of Gauss Dist
	double interval_Uniform = (ed - st) / 100.0; // Interval of Uniform Dist
	double interval_Gauss = ((mean + 4 * stddev) - (mean - 4 * stddev)) / 100.0; // Interval of Guass Dist, Range of Gauss Dist is from -4sigma to 4sigma, which covers almost 99.99%.
	std::uniform_real_distribution<> Uniform(st, ed); // Uniform Dist
	std::normal_distribution<> Gauss(mean, stddev); // Gauss Dist 
	std::vector<int> cnt_Uniform(105, 0); // Array of Uniform variables
	std::vector<int> cnt_Gauss(105, 0); // Array of Gauss variables


	for (int i = 1; i <= gentimes[k]; i++) { // Gen Gauss & Unifrom
		double uni = Uniform(gen);
		double gauss = Gauss(gen);
		int pos;
		for (int j = 1; j <= 100; j++) if (st + (j - 1)*interval_Uniform <= uni && uni < st + j * interval_Uniform) { pos = j; break; }
		cnt_Uniform[pos] += 1;

		for (int j = 1; j <= 100; j++) {
			if ((mean - 4 * stddev) + (j - 1)*interval_Gauss <= gauss && gauss < (mean - 4 * stddev) + j * interval_Gauss)
			{
				pos = j; break;
			}
		}
		cnt_Gauss[pos] += 1;
	}

	CClientDC dc(this); // To Draw on Screen
	dc.TextOutW(x0+10, 20 , _T("Histogram of Gauss Distribution"));
	dc.TextOutW(x0+10, y0_gauss + 50, _T("Histogram of Uniform Distribution"));
	CFont font;
	font.CreatePointFont(60,_T("Arial"));
	dc.SelectObject(font);
	CBrush brush;

	CString cstr;
	

	brush.CreateSolidBrush(RGB(255, 0, 0));
	dc.SelectObject(brush);

	int MAX = -1; // To MAXIMUM value
	int width = 18;
	//Gauss Dist
	for (int i = 1; i <= 100; i++) { // Print Histogram of Gauss
		//sine = (int)(50 * sin(	0.1*(double)k));
		int y = cnt_Gauss[i]; MAX = max(y, MAX);
		dc.Rectangle(x0 + (i * width), y0_gauss - 10- (cnt_Gauss[i]/resize[k]), x0 + (i + 1) * width - 2, y0_gauss - 10);
		if (i % 2) {
			cstr.Format(_T("%.2f"), (mean - 4 * stddev) + (i - 1)*interval_Gauss);
			dc.TextOutW(x0+(i*width), y0_gauss+5, cstr);

		}
		else {
			cstr.Format(_T("%.2f"), (mean - 4 * stddev) + (i - 1)*interval_Gauss);
			dc.TextOutW(x0 + (i * width), y0_gauss + 25, cstr);
		}

		//pDC->MoveTo(x0 + (i * 13), y0_gauss - 2 - cnt_Gauss[i]);
		//pDC->LineTo(x0 + (i * 13), y0_gauss- 2);
		//pDC->SetPixel(x0+(i+1), y + 100, RGB(255, 0, 0));
	}
	MAX = -1;
	brush.DeleteObject();

	brush.CreateSolidBrush(RGB(0, 255, 0));
	dc.SelectObject(brush);

	//Uniform Dist
	for (int i = 1; i <= 100; i++) { // Print Histogram of Uniform
		//sine = (int)(50 * sin(0.1*(double)k));
		int y = cnt_Uniform[i]; MAX = max(y, MAX);
		dc.Rectangle(x0 + (i * width), y0_uniform - 10 - (cnt_Uniform[i] / resize[k]), x0 + (i + 1) * width - 2, y0_uniform - 10);
		if (i % 2) {
			cstr.Format(_T("%.2f"), st + (i - 1)*interval_Uniform);
			dc.TextOutW(x0 + (i*width), y0_uniform + 5, cstr);

		}
		else {
			cstr.Format(_T("%.2f"), st + (i - 1)*interval_Uniform);
			dc.TextOutW(x0 + (i * width), y0_uniform + 25, cstr);
		}

	}
	brush.DeleteObject();

	// TODO: 여기에 원시 데이터에 대한 그리기 코드를 추가합니다.
}


// CHomework4View 인쇄

BOOL CHomework4View::OnPreparePrinting(CPrintInfo* pInfo)
{
	// 기본적인 준비
	return DoPreparePrinting(pInfo);
}

void CHomework4View::OnBeginPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 인쇄하기 전에 추가 초기화 작업을 추가합니다.
}

void CHomework4View::OnEndPrinting(CDC* /*pDC*/, CPrintInfo* /*pInfo*/)
{
	// TODO: 인쇄 후 정리 작업을 추가합니다.
}


// CHomework4View 진단

#ifdef _DEBUG
void CHomework4View::AssertValid() const
{
	CView::AssertValid();
}

void CHomework4View::Dump(CDumpContext& dc) const
{
	CView::Dump(dc);
}

CHomework4Doc* CHomework4View::GetDocument() const // 디버그되지 않은 버전은 인라인으로 지정됩니다.
{
	ASSERT(m_pDocument->IsKindOf(RUNTIME_CLASS(CHomework4Doc)));
	return (CHomework4Doc*)m_pDocument;
}
#endif //_DEBUG


// CHomework4View 메시지 처리기
