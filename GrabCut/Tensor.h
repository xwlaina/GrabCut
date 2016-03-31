// Tensor.h: interface for the Tensor class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_TENSOR_H)
#define _TENSOR_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Vector.h"
#include "Filter.h"
#include "Wavelet.h"
#include "Image.h"

class Tensor  
{
public:
	//////////////////////////////////////////////////////////////////////////
	//��߶ȷ����Խṹ�������캯������������
	Tensor(const IplImage *cv_image, BOOL isComputeGradient = FALSE);
	virtual ~Tensor();
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//���ܣ��õ����շ����Զ�߶Ƚṹ���������ս��
	//1������ͼ��Ķ�߶ȷ����Խṹ����
	Image < CVector* >* GetTensors() const;    
	//2�����ض�߶ȷ����Խṹ������ά��
	unsigned int Dim() const
	{
		return m_dim;
	}
	unsigned int GradDim() const
	{
		return m_grad_dim;
	}
	CMatrix ** GetGradient() const
	{
		return m_gradient;
	}
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//��ͼƬ����ʽ��ʾ�ṹ����
	void ShowTensorByColorImage();

	void ShowGradientImageOfAllScale();
	void test();
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//�������˲�
	//ֱ���ڹ��캯������ã�isComputeGradient = FALSE, isFiltering = TRUE 
	void nlST(BOOL isFiltering = TRUE);        	
	//////////////////////////////////////////////////////////////////////////

	double computeDistance2(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	double distance2_KL(CVector *tensors_1, CVector *tensors_2);
	///////////////////////////////////////////////////////////////////////////////////////	
	//�����ݶȣ�ֱ���ڹ��캯������ã�isComputeGradient = TRUE, isFiltering = FALSE or TRUE
	CMatrix ** CvtTensorToVectorBySVD(BOOL isFiltering = FALSE);
	///////////////////////////////////////////////////////////////////////////////////////
private:
	unsigned int m_w, m_h; //ͼ��Ŀ�Ⱥ͸߶�
	unsigned int m_levels; //�����߶�
	unsigned int m_dim;    //�����Զ�߶Ƚṹ����ά��
	CMatrix **m_tensor;             //�Ծ�����ʽ�洢�����Զ�߶Ƚṹ����,ÿһά��Ӧһ������
	Image < CVector* > *m_tensors;  //��������ʽ�洢�����Զ�߶Ƚṹ����

	CMatrix  **m_gradient;                   //m_gradient�洢ͼ���ݶȣ�����ÿһ��ͼ����ԣ��Ƕ�ÿһ���ص�Ľṹ��������SVD�ֽ�õ���
	unsigned int m_grad_dim;                 //�ݶ�ά��
	unsigned int m_axes_cnt;                 //m_axes_cnt=2,�������ά��
	IplImage * m_img;                        //����ԭͼ��
	Image < Color_RGB >  ** m_pImageTensorRGB;   //��ÿһ�߶ȵ�����ת��Ϊ��ɫͼ��洢����
};

#endif 
