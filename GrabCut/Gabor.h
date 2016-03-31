// Gabor.h: interface for the Gabor class.
/*
�޸���Ҧ��Gabor������룬
�ο����ף�
    Texture Features for Browsing and Retrieval of Image Data, by B.S.Manjunath and W.Y.Ma
�Լ���Ҧ������-��������������ң��ͼ���������
*/
//////////////////////////////////////////////////////////////////////



#if !defined(AFX_GABOR_H__896B99BC_E4E7_4068_8B90_D319214CA2C7__INCLUDED_)
#define AFX_GABOR_H__896B99BC_E4E7_4068_8B90_D319214CA2C7__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Helper.h"

class Gabor  
{
public:
	Image < CVector* >* GetFilters() const;
	double computeDistance2(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2) const;
	//void debug();
	void GetGaborFilter(CMatrix *Gr, CMatrix *Gi, int s, int n);
	void GaborFilteredImg();
	unsigned int Dim() const
	{
		return m_dim;
	}
	Gabor(const IplImage *cv_image);
	virtual ~Gabor();
private:
	unsigned int m_scale;//�߶ȸ�����S
	unsigned int m_orientation;//���������K
	unsigned int m_dim;
	double m_Uh;//highest spatial frequency
	double m_Ul;//lowest spatial frequency
	unsigned int m_flag;//if m_flag = 1, then remove the DC from the filter
	unsigned int m_side;//�˲���ϵ���뾶��filter dimension = (2*m_side+1)*(2*m_side+1)
	double m_alpha;
	CMatrix *m_image;//����ĻҶ�ͼ��
	unsigned int m_height, m_width;
	CHelper m_helper;
	Image < CVector* > *m_filters;
};

#endif // !defined(AFX_GABOR_H__896B99BC_E4E7_4068_8B90_D319214CA2C7__INCLUDED_)
