#pragma once

#include "GrabCutDlg.h"
// COptionsDlg �Ի���

class COptionsDlg : public CDialog
{
	DECLARE_DYNAMIC(COptionsDlg)

public:
	COptionsDlg(CWnd* pParent = NULL);   // ��׼���캯��
	virtual ~COptionsDlg();

// �Ի�������
	enum { IDD = IDD_OPTIONS };
	CTabCtrl	m_tab;
	CGrabCutDlg m_grabcut;

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV ֧��

protected:
	virtual BOOL OnInitDialog();
	virtual void OnOK();
	DECLARE_MESSAGE_MAP()
};
