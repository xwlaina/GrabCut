// OptionsDlg.cpp : ʵ���ļ�
//

#include "stdafx.h"
#include "GraphCut.h"
#include "OptionsDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

// COptionsDlg �Ի���

IMPLEMENT_DYNAMIC(COptionsDlg, CDialog)

COptionsDlg::COptionsDlg(CWnd* pParent /*=NULL*/)
	: CDialog(COptionsDlg::IDD, pParent)
{

}

COptionsDlg::~COptionsDlg()
{
}

void COptionsDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_OPTIONS_TAB, m_tab);
}


BEGIN_MESSAGE_MAP(COptionsDlg, CDialog)
END_MESSAGE_MAP()


// COptionsDlg ��Ϣ�������
BOOL COptionsDlg::OnInitDialog() 
{
	CDialog::OnInitDialog();
	// TODO: Add extra initialization here
	m_tab.InsertItem(0,_T("GrabCut"));
	m_grabcut.Create(IDD_GRABCUT,GetDlgItem(IDC_OPTIONS_TAB));
	//���IDC_TABTEST�ͻ�����С
	CRect rs;
	m_tab.GetClientRect(&rs);
	//�����ӶԻ����ڸ������е�λ��
	rs.top+=21; 
	rs.bottom-=4; 
	rs.left+=2; 
	rs.right-=4; 
	m_grabcut.MoveWindow(&rs);
	m_grabcut.ShowWindow(TRUE);
	//����Ĭ�ϵ�ѡ�
	m_tab.SetCurSel(0);
	return TRUE;  // return TRUE unless you set the focus to a control
	// EXCEPTION: OCX Property Pages should return FALSE
}

void COptionsDlg::OnOK() 
{
	// TODO: Add extra validation here
	CMainFrame* pFrame = (CMainFrame*) AfxGetApp()->GetMainWnd();
	Options *gOptions = pFrame->gOptions;
	m_grabcut.UpdateData(TRUE);
	gOptions->m_grabcut.dLambda = m_grabcut.m_dLambda;
	gOptions->m_grabcut.nNumBGMMs = m_grabcut.m_nNumBGMMs;
	gOptions->m_grabcut.nNumFGMMs = m_grabcut.m_nNumFGMMs;
	gOptions->m_grabcut.nNoiseRadius = m_grabcut.m_nNoiseRadius;
	gOptions->m_texture.nTensorScale = m_grabcut.m_nTensorScale;
	gOptions->m_texture.dTensorEsp = m_grabcut.m_dTensorEsp;
	gOptions->m_texture.dTensorTVPower = m_grabcut.m_dTensorTVPower;
	gOptions->m_texture.nTensorNosteps = m_grabcut.m_nTensorNosteps;
	gOptions->m_texture.dTensorSigma = m_grabcut.m_dTensorSigma;
	gOptions->m_texture.dTensorStepsize = m_grabcut.m_dTensorStepsize;
	gOptions->m_texture.nLowDim = m_grabcut.m_nLowDim;

	CDialog::OnOK();
}