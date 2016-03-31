#pragma once

#include "GrabCutDlg.h"
// COptionsDlg 对话框

class COptionsDlg : public CDialog
{
	DECLARE_DYNAMIC(COptionsDlg)

public:
	COptionsDlg(CWnd* pParent = NULL);   // 标准构造函数
	virtual ~COptionsDlg();

// 对话框数据
	enum { IDD = IDD_OPTIONS };
	CTabCtrl	m_tab;
	CGrabCutDlg m_grabcut;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV 支持

protected:
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP()
};
