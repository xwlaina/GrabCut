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
	//多尺度非线性结构张量构造函数与析构函数
	Tensor(const IplImage *cv_image, BOOL isComputeGradient = FALSE);
	virtual ~Tensor();
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//功能：得到最终非线性多尺度结构张量的最终结果
	//1：返回图像的多尺度非线性结构张量
	Image < CVector* >* GetTensors() const;    
	//2：返回多尺度非线性结构张量的维数
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
	//以图片的形式显示结构张量
	void ShowTensorByColorImage();

	void ShowGradientImageOfAllScale();
	void test();
	//////////////////////////////////////////////////////////////////////////

	//////////////////////////////////////////////////////////////////////////
	//非线性滤波
	//直接在构造函数后调用，isComputeGradient = FALSE, isFiltering = TRUE 
	void nlST(BOOL isFiltering = TRUE);        	
	//////////////////////////////////////////////////////////////////////////

	double computeDistance2(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	double distance2_KL(CVector *tensors_1, CVector *tensors_2);
	///////////////////////////////////////////////////////////////////////////////////////	
	//计算梯度，直接在构造函数后调用，isComputeGradient = TRUE, isFiltering = FALSE or TRUE
	CMatrix ** CvtTensorToVectorBySVD(BOOL isFiltering = FALSE);
	///////////////////////////////////////////////////////////////////////////////////////
private:
	unsigned int m_w, m_h; //图像的宽度和高度
	unsigned int m_levels; //张量尺度
	unsigned int m_dim;    //非线性多尺度结构张量维数
	CMatrix **m_tensor;             //以矩阵形式存储非线性多尺度结构张量,每一维对应一个矩阵
	Image < CVector* > *m_tensors;  //以向量形式存储非线性多尺度结构张量

	CMatrix  **m_gradient;                   //m_gradient存储图像梯度，对于每一幅图像而言，是对每一像素点的结构张量进行SVD分解得到的
	unsigned int m_grad_dim;                 //梯度维数
	unsigned int m_axes_cnt;                 //m_axes_cnt=2,坐标轴的维数
	IplImage * m_img;                        //保存原图像
	Image < Color_RGB >  ** m_pImageTensorRGB;   //将每一尺度的张量转化为彩色图像存储起来
};

#endif 
