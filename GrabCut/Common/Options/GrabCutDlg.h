#pragma once


// CGrabCutDlg 对话框

class CGrabCutDlg : public CDialog
{
	DECLARE_DYNAMIC(CGrabCutDlg)

public:
	CGrabCutDlg(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~CGrabCutDlg();

// 对话框数据
	enum { IDD = IDD_GRABCUT };
	UINT	m_nNumBGMMs;
	UINT	m_nNumFGMMs;
	double	m_dLambda;
	UINT	m_nNoiseRadius;
	UINT	m_nTensorScale;
	double	m_dTensorEsp;
	double	m_dTensorTVPower;
	UINT	m_nTensorNosteps;
	double	m_dTensorSigma;
	double	m_dTensorStepsize;
	UINT	m_nLowDim;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

protected:
	virtual BOOL OnInitDialog();

	DECLARE_MESSAGE_MAP()
};
