// Tensor.cpp: implementation of the Tensor class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Tensor.h"
#include <opencv2\legacy\compat.hpp> 

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
Tensor::Tensor(const IplImage *cv_image, BOOL isComputeGradient)
{
	//保存原图像的副本
	m_img=cvCreateImage(cvSize(cv_image->width,cv_image->height),cv_image->depth,3);
	cvCopyImage(cv_image,m_img);

	//获取非线性多尺度结构张量的参数值
	m_levels = 2;

	ASSERT(m_levels > 0 );

	m_dim = m_levels * SiNGLE_TENSOR_DIM;    //SiNGLE_TENSOR_DIM单一张量

	//SiNGLE_TENSOR_DIM=n(n+1)/2;反解n=m_axes_cnt，m_axes_cnt为坐标抽的维数
	m_axes_cnt = (unsigned int)(sqrt(2 * SiNGLE_TENSOR_DIM + 0.25) - 0.5);   // 2

	m_grad_dim = m_levels * m_axes_cnt;      //m_grad_dim


	////////////////////////////////////////////////////////////////////////////
	//将多通道转化为单通道，默认为三个通道
	unsigned int x,y,i,n;
	m_w = cv_image->width;
	m_h = cv_image->height;
	IplImage *cv_channels[3];
	for (n = 0;n < 3;n++)
	{
		cv_channels[n] = cvCreateImage( cvGetSize(cv_image), cv_image->depth, 1 );
	}
	cvSplit(cv_image, cv_channels[0], cv_channels[1], cv_channels[2], NULL);


	////////////////////////////////////////////////////////////////////////////
	//初始化m_tensor，CMatrix(m_h,m_w)创建一个矩阵，其元素全为0
	m_tensor = new CMatrix *[m_dim];
	for (i=0;i<m_dim;i++)
	{   
		m_tensor[i] = new CMatrix(m_h,m_w);
	}

	////////////////////////////////////////////////////////////////////////////
	//将每一尺度的张量转化为彩色图像存储起来，申请空间	
	m_pImageTensorRGB=new Image<Color_RGB> *[m_levels];
	for (i=0;i<m_levels;i++)
	{   
		m_pImageTensorRGB[i] = new Image<Color_RGB> (m_w,m_h);
	}

	//初始化m_gradient
	if (isComputeGradient)
	{
		m_gradient = new CMatrix *[m_grad_dim];
		for (i=0;i<m_grad_dim;i++)
		{
			m_gradient[i] = new CMatrix(m_h,m_w);
		}
	}
	else
	{
		m_gradient = NULL;
	}


	//辅助矩阵
	CMatrix image(m_h, m_w);
	CMatrix dx(m_h,m_w);
	CMatrix dy(m_h,m_w);
	CMatrix dx2(m_h,m_w);
	CMatrix dy2(m_h,m_w);
	CMatrix dxdy(m_h,m_w);

	//利用固定数据创建一个矩阵
	CvMat cv_dx2 = cvMat(m_h, m_w, CV_64FC1, dx2.GetData());
	CvMat cv_dy2 = cvMat(m_h, m_w, CV_64FC1, dy2.GetData());
	CvMat cv_dxdy =cvMat(m_h, m_w, CV_64FC1, dxdy.GetData());


	//完成IplImage向CMatrix类型的转换，对每一个颜色通道分别进行处理
	for (n = 0;n <3;n++)	//n表示通道数，默认为3
	{  
		//将每一个通道的元素拷贝到image中
		for (y = 0; y < m_h; y++)
		{
			for (x = 0; x < m_w; x++)
			{
				uchar* dst = &CV_IMAGE_ELEM( cv_channels[n], uchar, y, x );
				image.SetElement(y, x, (double)(dst[0]));
			}
		}
		//计算每一个颜色通道的梯度(x方向,y方向)并分别赋给dx，dy
		image.centdiffX(dx);
		image.centdiffY(dy);

		//将dx，dy分别赋给cv_dx，cv_dy
		CvMat cv_dx = cvMat(m_h, m_w, CV_64FC1, dx.GetData());
		CvMat cv_dy = cvMat(m_h, m_w, CV_64FC1, dy.GetData());

		//初始化cv_tensor0，cv_tensor1，cv_tensor2，此时m_tensor[0],m_tensor[1],m_tensor[2]均初始化0
		CvMat cv_tensor0 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[0])->GetData());
		CvMat cv_tensor1 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[1])->GetData());
		CvMat cv_tensor2 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[2])->GetData());

		//计算图像的梯度，保存在cv_gradX，cv_gradY中，并赋值给m_gradient[0],m_gradient[1]
		if (isComputeGradient)
		{   
			//cv_gradX,cv_gradY初始化并计算
			CvMat cv_gradX = cvMat(m_h, m_w, CV_64FC1, (m_gradient[0])->GetData());
			CvMat cv_gradY = cvMat(m_h, m_w, CV_64FC1, (m_gradient[1])->GetData());
			cvAdd(&cv_gradX, &cv_dx, &cv_gradX);//对于三个通道进行累加
			cvAdd(&cv_gradY, &cv_dy, &cv_gradY);
		}

		//计算结构张量，cv_tensor0=dx*dx,cv_tensor1=dy*dy,cv_tensor2=dx*dy
		cvMul(&cv_dx, &cv_dx, &cv_dx2);
		cvAdd(&cv_tensor0, &cv_dx2, &cv_tensor0);
		cvMul(&cv_dy, &cv_dy, &cv_dy2);
		cvAdd(&cv_tensor1, &cv_dy2, &cv_tensor1);
		cvMul(&cv_dx, &cv_dy, &cv_dxdy);
		cvAdd(&cv_tensor2, &cv_dxdy, &cv_tensor2);


		//单尺度计算完毕，以下为多尺度非线性结构张量的计算方法
		if (m_levels > 1)
		{   
			unsigned int wavelet_levels = m_levels - 1;	//-1的原因是因为之前没有if (m_levels==1)的判断语句	
			double dMaxValue,dMinValue;
			cvMinMaxLoc(cv_channels[n], &dMinValue, &dMaxValue);//Finds global minimum, maximum 

			//将图像的像素值归一化到[0,1]
			Wavelet *wave = new Wavelet(&image, dMinValue, dMaxValue, wavelet_levels); //调用Wavelet的构造函数

			//新建WaveletDetailImages结构体的数组
			WaveletDetailImages *D_images = new WaveletDetailImages[wavelet_levels];

			for (i = 0; i < wavelet_levels; i++)
			{
				D_images[i].Detail_1 = new CMatrix(m_h, m_w);
				D_images[i].Detail_2 = new CMatrix(m_h, m_w);
			}

			wave->execute(D_images);//得到D(s,x),D(s,y)

			for (i = 0; i < wavelet_levels; i++)
			{   
				//默认多尺度结构张量的比例因子a=2
				double scale = pow((float)0.25, (int)(i + 1));              //见公式(2-15)
				CvMat cv_dx = cvMat(m_h, m_w, CV_64FC1, D_images[i].Detail_1->GetData());
				CvMat cv_dy = cvMat(m_h, m_w, CV_64FC1, D_images[i].Detail_2->GetData());
				CvMat cv_tensor0 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[(i+1) * SiNGLE_TENSOR_DIM])->GetData());
				CvMat cv_tensor1 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[(i+1) * SiNGLE_TENSOR_DIM + 1])->GetData());
				CvMat cv_tensor2 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[(i+1) * SiNGLE_TENSOR_DIM + 2])->GetData());
				//计算梯度
				if (isComputeGradient)
				{
					CvMat cv_gradX = cvMat(m_h, m_w, CV_64FC1, (m_gradient[(i+1) * m_axes_cnt])->GetData());
					CvMat cv_gradY = cvMat(m_h, m_w, CV_64FC1, (m_gradient[(i+1) * m_axes_cnt + 1])->GetData());
					cvAdd(&cv_gradX, &cv_dx, &cv_gradX);
					cvAdd(&cv_gradY, &cv_dy, &cv_gradY);
				}
				//计算张量
				cvMul(&cv_dx, &cv_dx, &cv_dx2, scale);
				cvAdd(&cv_tensor0, &cv_dx2, &cv_tensor0);
				cvMul(&cv_dy, &cv_dy, &cv_dy2, scale);
				cvAdd(&cv_tensor1, &cv_dy2, &cv_tensor1);
				cvMul(&cv_dx, &cv_dy, &cv_dxdy, scale);
				cvAdd(&cv_tensor2, &cv_dxdy, &cv_tensor2);
			}
			for (i = 0; i < wavelet_levels; i++)
			{
				delete D_images[i].Detail_1;
				delete D_images[i].Detail_2;
			}
			delete [] D_images;
			delete wave;
		}
		cvReleaseImage(&cv_channels[n]);
	}

	//将每一尺度的结构张量转换为彩色图像存储起来
	for (i=0;i<m_levels;i++)
	{
		for (y=0;y<m_h;y++)
		{
			for (x=0;x<m_w;x++)
			{
				(*m_pImageTensorRGB[i])(x,y).r=(m_tensor[i*SiNGLE_TENSOR_DIM])->GetElement(y,x);
				(*m_pImageTensorRGB[i])(x,y).g=(m_tensor[i*SiNGLE_TENSOR_DIM+1])->GetElement(y,x);
				(*m_pImageTensorRGB[i])(x,y).b=(m_tensor[i*SiNGLE_TENSOR_DIM+2])->GetElement(y,x);
			}
		}
	}
	m_tensors = NULL;	
}



Tensor::~Tensor()
{
	if (m_tensor != NULL)
	{
		for (int i=0;i<m_dim;i++)
		{
			if (m_tensor[i] != NULL)
			{
				delete m_tensor[i];
				m_tensor[i] = NULL;
			}
		}	
		delete [] m_tensor;
		m_tensor = NULL;
	}

	if (m_pImageTensorRGB != NULL)
	{
		for (int i=0;i<m_levels;i++)
		{
			if (m_pImageTensorRGB[i] != NULL)
			{
				delete m_pImageTensorRGB[i];
				m_pImageTensorRGB[i] = NULL;
			}
		}	
		delete [] m_pImageTensorRGB;
		m_pImageTensorRGB = NULL;
	}

	if (m_gradient != NULL)
	{
		for (int i=0;i<m_grad_dim;i++)
		{
			if (m_gradient[i] != NULL)
			{
				delete m_gradient[i];
				m_gradient[i] = NULL;
			}
		}	
		delete [] m_gradient;
		m_gradient = NULL;
	}


	if (m_tensors != NULL)
	{
		unsigned int x,y;
		for (y = 0; y < m_h; y++)
		{
			for (x = 0; x < m_w; x++)
			{
				if ((*m_tensors)(x,y) != NULL)
				{
					delete (*m_tensors)(x,y);
					(*m_tensors)(x,y) = NULL;
				}
			}
		}
		delete m_tensors;
		m_tensors = NULL;
	}
	if (m_img!=NULL)
	{
		cvReleaseImage(&m_img);
		m_img=NULL;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
//无论是否滤波，均将Tensor转化为向量m_tensors存储起来，并释放m_tensor所占内存
/////////////////////////////////////////////////////////////////////////////////////////
void Tensor::nlST(BOOL isFiltering)
{   
	//对于张量图像进行滤波
	
	if (isFiltering)
	{
		CFilter filter;
		filter.Diff(m_tensor, m_dim);
	}
	unsigned int x,y;
	int i;
	m_tensors = new Image < CVector* > (m_w, m_h);   //定义一个m_tensors的指向IplImage的指针

	for (y = 0; y < m_h; y++)
	{
		for (x = 0; x < m_w; x++)
		{
			(*m_tensors)(x,y) = new CVector(m_dim);
			for (i=0;i<m_dim;i++)
			{
				((*m_tensors)(x,y))->set(i, (m_tensor[i])->GetElement(y,x));
			}
		}
	}

	if (m_tensor != NULL)
	{
		for (int i=0;i<m_dim;i++)
		{
			if (m_tensor[i] != NULL)
			{
				delete m_tensor[i];
				m_tensor[i] = NULL;
			}
		}	
		delete [] m_tensor;
		m_tensor = NULL;
	}

}

////////////////////////////////////////////////////////////////////////////////////	
//显示每一尺度原图的梯度图像，梯度的计算公式为沿每一坐标抽的梯度的平均和的平方根
////////////////////////////////////////////////////////////////////////////////////
void Tensor::ShowGradientImageOfAllScale()
{
	unsigned int x,y;
	int i, k;
	double gradvalue;
	CString *pTitle=new CString [m_levels];


	IplImage ** grad= new IplImage *[m_levels];
	for (i=0;i<m_levels;i++)
	{
		grad[i] = cvCreateImage( cvSize( m_w, m_h ), 8, 1 );
		cvZero(grad[i]);
	}
	for (i=0;i<m_levels;i++)
	{
		for (y = 0;y < m_h;y++)
		{
			for (x = 0;x < m_w;x++)
			{
				uchar* dst = &CV_IMAGE_ELEM( grad[i], uchar, y, x);
				gradvalue = 0.0;
				for (k=0;k<m_axes_cnt;k++)
				{
					gradvalue += pow((m_gradient[i * m_axes_cnt + k])->GetElement(y, x), 2);
				}
				dst[0] = (uchar)(sqrt(gradvalue));
			}
		}

		pTitle[i].Format(_T("Image Gradient of Level %d"),i);
		cvNamedWindow( (char *)(LPCTSTR)pTitle[i], 1 );
		cvShowImage((char *)(LPCTSTR)pTitle[i], grad[i] );
	}
	for (i=0;i<m_levels;i++)
	{
		cvReleaseImage(&grad[i]);
		//cvDestroyWindow(pTitle[i]);
	}
}

////////////////////////////////////////////////////////////////////////////////////	
//以彩色图像显示每一尺度的张量信息
////////////////////////////////////////////////////////////////////////////////////
void Tensor::ShowTensorByColorImage()
{
	double ret_minr=0.0;
	double ret_maxr=0.0;
	double ret_ming=0.0;
	double ret_maxg=0.0;
	double ret_minb=0.0;
	double ret_maxb=0.0;
	int x,y,i;
	//纹理特征
	IplImage **pImg= new IplImage *[m_levels];
	for (i = 0;i < m_levels;i++)
	{
		pImg[i] = cvCreateImage( cvGetSize(m_img), m_img->depth, 3);
		cvZero(pImg[i]);
	}

	CString * ptitle=new CString [m_levels];

	for (i=0;i<m_levels;i++)
	{
		//找到每幅图像颜色通道的上限与下限值
		for (y=0; y<m_h;y++)
		{
			for (x=0;x<m_w;x++)
			{
				if((*m_pImageTensorRGB[i])(x,y).r>ret_maxr)
				{
					ret_maxr=(*m_pImageTensorRGB[i])(x,y).r;
				}
				if ((*m_pImageTensorRGB[i])(x,y).r<ret_minr)
				{
					ret_minr=(*m_pImageTensorRGB[i])(x,y).r;
				}

				if((*m_pImageTensorRGB[i])(x,y).g>ret_maxg)
				{
					ret_maxg=(*m_pImageTensorRGB[i])(x,y).g;
				}
				if ((*m_pImageTensorRGB[i])(x,y).g<ret_ming)
				{
					ret_ming=(*m_pImageTensorRGB[i])(x,y).g;
				}

				if((*m_pImageTensorRGB[i])(x,y).b>ret_maxb)
				{
					ret_maxb=(*m_pImageTensorRGB[i])(x,y).b;
				}
				if ((*m_pImageTensorRGB[i])(x,y).b<ret_minb)
				{
					ret_minb=(*m_pImageTensorRGB[i])(x,y).b;
				}

			}
		}
		uchar * dst=(uchar *)pImg[i]->imageData;
		for (y=0; y<m_h;y++)
		{
			for (x=0;x<m_w;x++)
			{
				int temp=y*(pImg[i]->widthStep)+3*x;
				dst[temp+2]=(uchar)(((*m_pImageTensorRGB[i])(x,y).r-ret_minr)/(ret_maxr-ret_minr)*256);
				dst[temp+1]=(uchar)(((*m_pImageTensorRGB[i])(x,y).g-ret_ming)/(ret_maxg-ret_ming)*256);
				dst[temp+0]=(uchar)(((*m_pImageTensorRGB[i])(x,y).b-ret_minb)/(ret_maxb-ret_minb)*256);
			}
		}
		ptitle[i].Format(_T("Image Texture of Level %d"),i);
		cvNamedWindow((char *)(LPCTSTR)ptitle[i],CV_WINDOW_AUTOSIZE);
		cvShowImage((char *)(LPCTSTR)ptitle[i],pImg[i]);
	}
	if (pImg != NULL)
	{
		for (i=0;i<m_levels;i++)
		{
			cvReleaseImage(&pImg[i]);
		}
		delete [] pImg;
	}
}

//获得图像的纹理特征
Image < CVector* >* Tensor::GetTensors() const
{
	return m_tensors;
}

////////////////////////////////////////////////////////////////////////////////////////
// 利用SVD对张量进行分解，提取最大特征值所对应的特征向量，存储在m_gradient中
////////////////////////////////////////////////////////////////////////////////////////
CMatrix ** Tensor::CvtTensorToVectorBySVD(BOOL isFiltering)
{
	//判断是否滤波
	if (isFiltering)
	{
		CFilter filter;
		filter.Diff(m_tensor, m_dim);
	}
	unsigned int x,y;
	int i,k;
	CMatrix *st = new CMatrix(m_axes_cnt, m_axes_cnt);
	CVector *eigenvalues = new CVector(m_axes_cnt);
	CMatrix *eigenvectors = new CMatrix(m_axes_cnt, m_axes_cnt);
	CvMat cv_st = cvMat(m_axes_cnt, m_axes_cnt, CV_TYPE , st->GetData());
	CvMat eval = cvMat(m_axes_cnt, 1, CV_TYPE, eigenvalues->addr());
	CvMat evec = cvMat(m_axes_cnt, m_axes_cnt, CV_TYPE, eigenvectors->GetData());

	double signFlag, sqrtLambda;
	for (y=0;y<m_h;y++)
	{
		for (x=0;x<m_w;x++)
		{
			for (i=0;i<m_levels;i++)
			{
				//将张量转换成矩阵存储起来
				st->SetElement(0, 0, (m_tensor[i * SiNGLE_TENSOR_DIM])->GetElement(y,x));
				st->SetElement(1, 1, (m_tensor[i * SiNGLE_TENSOR_DIM + 1])->GetElement(y,x));
				st->SetElement(0, 1, (m_tensor[i * SiNGLE_TENSOR_DIM + 2])->GetElement(y,x));
				st->SetElement(1, 0, (m_tensor[i * SiNGLE_TENSOR_DIM + 2])->GetElement(y,x));
				//SVD分解
				cvSVD(&cv_st, &eval, &evec);
				signFlag = 0.0;

				for (k=0;k<m_axes_cnt;k++)
				{
					signFlag += ((m_gradient[i * m_axes_cnt + k])->GetElement(y,x) * (eigenvectors->GetElement(k,0)));
				}
				//sqrtLambda 的初始化
				sqrtLambda = sqrt(eigenvalues->get(0));
				//判断向量的方向是否与之前的梯度方向一致
				if (signFlag >= 0)
				{
					for (k=0;k<m_axes_cnt;k++)
					{
						(m_gradient[i * m_axes_cnt + k])->SetElement(y, x, sqrtLambda * (eigenvectors->GetElement(k,0)));
					}
				} 
				else
				{
					for (k=0;k<m_axes_cnt;k++)
					{
						(m_gradient[i * m_axes_cnt + k])->SetElement(y, x, (- 1.0) * sqrtLambda * (eigenvectors->GetElement(k,0)));
					}
				}
			}
		}
	}

	delete st;
	delete eigenvalues;
	delete eigenvectors;

	if (m_tensor != NULL)
	{
		for (int i=0;i<m_dim;i++)
		{
			if (m_tensor[i] != NULL)
			{
				delete m_tensor[i];
				m_tensor[i] = NULL;
			}
		}	
		delete [] m_tensor;
		m_tensor = NULL;
	}
	return m_gradient;
}

////////////////////////////////////////////////////////////////////////////////////	
//显示梯度图像，梯度的计算公式为沿每一坐标抽的梯度的平均和的平方根
////////////////////////////////////////////////////////////////////////////////////
void Tensor::test()
{
	unsigned int x,y;
	CString title;
	int i, k;
	double gradvalue;
	uchar allgradvalue=0;

	IplImage* grad = cvCreateImage( cvSize( m_w, m_h ), 8, 1 );
	cvZero( grad );

	for (y = 0;y < m_h;y++)
	{
		for (x = 0;x < m_w;x++)
		{
			uchar* dst = &CV_IMAGE_ELEM( grad, uchar, y, x);
			gradvalue = 0.0;
			for (i=0;i<m_levels;i++)
			{
				for (k=0;k<m_axes_cnt;k++)
				{
					gradvalue += pow((m_gradient[i * m_axes_cnt + k])->GetElement(y, x), 2);
				}
				allgradvalue+=(uchar)(sqrt(gradvalue));
			}
			if (allgradvalue>255)
			{
				allgradvalue=255;
			}
			dst[0] =allgradvalue;
			allgradvalue=0;
		}
	}

	title.Format(_T("Image Gradient Sum of All Scales"));
	cvNamedWindow( (char *)(LPCTSTR)title, 1 );
	cvShowImage((char *)(LPCTSTR)title, grad);
	cvReleaseImage(&grad);
}


double Tensor::computeDistance2(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	CVector *tensors_1 = (*m_tensors)(x1,y1);
	CVector *tensors_2 = (*m_tensors)(x2,y2);
	return distance2_KL(tensors_1, tensors_2);
}


//计算两个张量之间的距离
double Tensor::distance2_KL(CVector *tensors_1,CVector *tensors_2)
{
	int levels = tensors_1->size() / SiNGLE_TENSOR_DIM;
	double w = 0.0;
	double value1, value2;
	double a_1, b_1, c_1, a_2, b_2, c_2;
	unsigned int pos, pos1, pos2;
	for (int k = 0; k < levels;k++)
	{
		pos = k * SiNGLE_TENSOR_DIM;
		pos1 = pos + 1;
		pos2 = pos + 2;
		a_1 = tensors_1->get(pos);
		b_1 = tensors_1->get(pos1);
		c_1 = tensors_1->get(pos2);
		a_2 = tensors_2->get(pos);
		b_2 = tensors_2->get(pos1);
		c_2 = tensors_2->get(pos2);
		ASSERT((a_1 * b_1 - c_1 * c_1) != 0 && (a_2 * b_2 - c_2 * c_2) != 0);
		value1 = (a_1 * b_2 + b_1 * a_2 - 2 * c_1 * c_2) / (a_1 * b_1 - c_1 * c_1);
		value2 = (a_2 * b_1 + b_2 * a_1 - 2 * c_2 * c_1) / (a_2 * b_2 - c_2 * c_2);
		w += ((value1 + value2 - 4) / 4);
	}
	return w;
}

