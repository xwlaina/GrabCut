// GrabCutDlg.cpp : 实现文件
//

#include "stdafx.h"
#include "GraphCut.h"
#include "GrabCutDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// CGrabCutDlg 对话框
IMPLEMENT_DYNAMIC(CGrabCutDlg, CDialog)

CGrabCutDlg::CGrabCutDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CGrabCutDlg::IDD, pParent)
{
	m_nNumBGMMs = 5;
	m_nNumFGMMs = 5;
	m_dLambda = 50.0;
	m_nNoiseRadius = 1;
	m_nTensorScale = 2;
    m_dTensorEsp = 2.0;
	m_dTensorTVPower = 0.6;
	m_nTensorNosteps = 2;
	m_dTensorSigma = 0.0;
	m_dTensorStepsize = 5000;
	m_nLowDim = 2;
}

CGrabCutDlg::~CGrabCutDlg()
{
}

void CGrabCutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Text(pDX, IDC_BGMMS_NUM, m_nNumBGMMs);
	DDX_Text(pDX, IDC_FGMMS_NUM, m_nNumFGMMs);
	DDX_Text(pDX, IDC_LAMBDA, m_dLambda);
	DDX_Text(pDX, IDC_NOISE_RADIUS, m_nNoiseRadius);
	DDX_Text(pDX, IDC_TENSOR_SCALE, m_nTensorScale);
	DDX_Text(pDX, IDC_TENSOR_ESP, m_dTensorEsp);
	DDX_Text(pDX, IDC_TENSOR_NOSTEPS, m_nTensorNosteps);
	DDX_Text(pDX, IDC_TENSOR_SIGMA, m_dTensorSigma);
	DDX_Text(pDX, IDC_TENSOR_STEPSIZE, m_dTensorStepsize);
	DDX_Text(pDX, IDC_TENSOR_TV_POWER, m_dTensorTVPower);
	DDX_Text(pDX, IDC_LOW_DIM, m_nLowDim);
}

BEGIN_MESSAGE_MAP(CGrabCutDlg, CDialog)
END_MESSAGE_MAP()

// CGrabCutDlg 消息处理程序
BOOL CGrabCutDlg::OnInitDialog() 
{
	CDialog::OnInitDialog();

	// TODO: Add extra initialization here
	CMainFrame* pFrame = (CMainFrame*) AfxGetApp()->GetMainWnd();
	Options *gOptions = pFrame->gOptions;
	m_dLambda = gOptions->m_grabcut.dLambda;
	m_nNumBGMMs = gOptions->m_grabcut.nNumBGMMs;
	m_nNumFGMMs = gOptions->m_grabcut.nNumFGMMs;
	m_nNoiseRadius = gOptions->m_grabcut.nNoiseRadius;
	m_nTensorScale=gOptions->m_texture.nTensorScale;
	m_nLowDim=gOptions->m_texture.nLowDim;
	m_dTensorEsp=gOptions->m_texture.dTensorEsp;
	m_nTensorNosteps=gOptions->m_texture.nTensorNosteps;
	m_dTensorSigma=gOptions->m_texture.dTensorSigma;
	m_dTensorStepsize=gOptions->m_texture.dTensorStepsize;
	m_dTensorTVPower=gOptions->m_texture.dTensorTVPower;
	UpdateData(FALSE);
	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}