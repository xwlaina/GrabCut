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
	//����ԭͼ��ĸ���
	m_img=cvCreateImage(cvSize(cv_image->width,cv_image->height),cv_image->depth,3);
	cvCopyImage(cv_image,m_img);

	//��ȡ�����Զ�߶Ƚṹ�����Ĳ���ֵ
	m_levels = 2;

	ASSERT(m_levels > 0 );

	m_dim = m_levels * SiNGLE_TENSOR_DIM;    //SiNGLE_TENSOR_DIM��һ����

	//SiNGLE_TENSOR_DIM=n(n+1)/2;����n=m_axes_cnt��m_axes_cntΪ������ά��
	m_axes_cnt = (unsigned int)(sqrt(2 * SiNGLE_TENSOR_DIM + 0.25) - 0.5);   // 2

	m_grad_dim = m_levels * m_axes_cnt;      //m_grad_dim


	////////////////////////////////////////////////////////////////////////////
	//����ͨ��ת��Ϊ��ͨ����Ĭ��Ϊ����ͨ��
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
	//��ʼ��m_tensor��CMatrix(m_h,m_w)����һ��������Ԫ��ȫΪ0
	m_tensor = new CMatrix *[m_dim];
	for (i=0;i<m_dim;i++)
	{   
		m_tensor[i] = new CMatrix(m_h,m_w);
	}

	////////////////////////////////////////////////////////////////////////////
	//��ÿһ�߶ȵ�����ת��Ϊ��ɫͼ��洢����������ռ�	
	m_pImageTensorRGB=new Image<Color_RGB> *[m_levels];
	for (i=0;i<m_levels;i++)
	{   
		m_pImageTensorRGB[i] = new Image<Color_RGB> (m_w,m_h);
	}

	//��ʼ��m_gradient
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


	//��������
	CMatrix image(m_h, m_w);
	CMatrix dx(m_h,m_w);
	CMatrix dy(m_h,m_w);
	CMatrix dx2(m_h,m_w);
	CMatrix dy2(m_h,m_w);
	CMatrix dxdy(m_h,m_w);

	//���ù̶����ݴ���һ������
	CvMat cv_dx2 = cvMat(m_h, m_w, CV_64FC1, dx2.GetData());
	CvMat cv_dy2 = cvMat(m_h, m_w, CV_64FC1, dy2.GetData());
	CvMat cv_dxdy =cvMat(m_h, m_w, CV_64FC1, dxdy.GetData());


	//���IplImage��CMatrix���͵�ת������ÿһ����ɫͨ���ֱ���д���
	for (n = 0;n <3;n++)	//n��ʾͨ������Ĭ��Ϊ3
	{  
		//��ÿһ��ͨ����Ԫ�ؿ�����image��
		for (y = 0; y < m_h; y++)
		{
			for (x = 0; x < m_w; x++)
			{
				uchar* dst = &CV_IMAGE_ELEM( cv_channels[n], uchar, y, x );
				image.SetElement(y, x, (double)(dst[0]));
			}
		}
		//����ÿһ����ɫͨ�����ݶ�(x����,y����)���ֱ𸳸�dx��dy
		image.centdiffX(dx);
		image.centdiffY(dy);

		//��dx��dy�ֱ𸳸�cv_dx��cv_dy
		CvMat cv_dx = cvMat(m_h, m_w, CV_64FC1, dx.GetData());
		CvMat cv_dy = cvMat(m_h, m_w, CV_64FC1, dy.GetData());

		//��ʼ��cv_tensor0��cv_tensor1��cv_tensor2����ʱm_tensor[0],m_tensor[1],m_tensor[2]����ʼ��0
		CvMat cv_tensor0 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[0])->GetData());
		CvMat cv_tensor1 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[1])->GetData());
		CvMat cv_tensor2 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[2])->GetData());

		//����ͼ����ݶȣ�������cv_gradX��cv_gradY�У�����ֵ��m_gradient[0],m_gradient[1]
		if (isComputeGradient)
		{   
			//cv_gradX,cv_gradY��ʼ��������
			CvMat cv_gradX = cvMat(m_h, m_w, CV_64FC1, (m_gradient[0])->GetData());
			CvMat cv_gradY = cvMat(m_h, m_w, CV_64FC1, (m_gradient[1])->GetData());
			cvAdd(&cv_gradX, &cv_dx, &cv_gradX);//��������ͨ�������ۼ�
			cvAdd(&cv_gradY, &cv_dy, &cv_gradY);
		}

		//����ṹ������cv_tensor0=dx*dx,cv_tensor1=dy*dy,cv_tensor2=dx*dy
		cvMul(&cv_dx, &cv_dx, &cv_dx2);
		cvAdd(&cv_tensor0, &cv_dx2, &cv_tensor0);
		cvMul(&cv_dy, &cv_dy, &cv_dy2);
		cvAdd(&cv_tensor1, &cv_dy2, &cv_tensor1);
		cvMul(&cv_dx, &cv_dy, &cv_dxdy);
		cvAdd(&cv_tensor2, &cv_dxdy, &cv_tensor2);


		//���߶ȼ�����ϣ�����Ϊ��߶ȷ����Խṹ�����ļ��㷽��
		if (m_levels > 1)
		{   
			unsigned int wavelet_levels = m_levels - 1;	//-1��ԭ������Ϊ֮ǰû��if (m_levels==1)���ж����	
			double dMaxValue,dMinValue;
			cvMinMaxLoc(cv_channels[n], &dMinValue, &dMaxValue);//Finds global minimum, maximum 

			//��ͼ�������ֵ��һ����[0,1]
			Wavelet *wave = new Wavelet(&image, dMinValue, dMaxValue, wavelet_levels); //����Wavelet�Ĺ��캯��

			//�½�WaveletDetailImages�ṹ�������
			WaveletDetailImages *D_images = new WaveletDetailImages[wavelet_levels];

			for (i = 0; i < wavelet_levels; i++)
			{
				D_images[i].Detail_1 = new CMatrix(m_h, m_w);
				D_images[i].Detail_2 = new CMatrix(m_h, m_w);
			}

			wave->execute(D_images);//�õ�D(s,x),D(s,y)

			for (i = 0; i < wavelet_levels; i++)
			{   
				//Ĭ�϶�߶Ƚṹ�����ı�������a=2
				double scale = pow((float)0.25, (int)(i + 1));              //����ʽ(2-15)
				CvMat cv_dx = cvMat(m_h, m_w, CV_64FC1, D_images[i].Detail_1->GetData());
				CvMat cv_dy = cvMat(m_h, m_w, CV_64FC1, D_images[i].Detail_2->GetData());
				CvMat cv_tensor0 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[(i+1) * SiNGLE_TENSOR_DIM])->GetData());
				CvMat cv_tensor1 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[(i+1) * SiNGLE_TENSOR_DIM + 1])->GetData());
				CvMat cv_tensor2 = cvMat(m_h, m_w, CV_64FC1, (m_tensor[(i+1) * SiNGLE_TENSOR_DIM + 2])->GetData());
				//�����ݶ�
				if (isComputeGradient)
				{
					CvMat cv_gradX = cvMat(m_h, m_w, CV_64FC1, (m_gradient[(i+1) * m_axes_cnt])->GetData());
					CvMat cv_gradY = cvMat(m_h, m_w, CV_64FC1, (m_gradient[(i+1) * m_axes_cnt + 1])->GetData());
					cvAdd(&cv_gradX, &cv_dx, &cv_gradX);
					cvAdd(&cv_gradY, &cv_dy, &cv_gradY);
				}
				//��������
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

	//��ÿһ�߶ȵĽṹ����ת��Ϊ��ɫͼ��洢����
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
//�����Ƿ��˲�������Tensorת��Ϊ����m_tensors�洢���������ͷ�m_tensor��ռ�ڴ�
/////////////////////////////////////////////////////////////////////////////////////////
void Tensor::nlST(BOOL isFiltering)
{   
	//��������ͼ������˲�
	
	if (isFiltering)
	{
		CFilter filter;
		filter.Diff(m_tensor, m_dim);
	}
	unsigned int x,y;
	int i;
	m_tensors = new Image < CVector* > (m_w, m_h);   //����һ��m_tensors��ָ��IplImage��ָ��

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
//��ʾÿһ�߶�ԭͼ���ݶ�ͼ���ݶȵļ��㹫ʽΪ��ÿһ�������ݶȵ�ƽ���͵�ƽ����
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
//�Բ�ɫͼ����ʾÿһ�߶ȵ�������Ϣ
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
	//��������
	IplImage **pImg= new IplImage *[m_levels];
	for (i = 0;i < m_levels;i++)
	{
		pImg[i] = cvCreateImage( cvGetSize(m_img), m_img->depth, 3);
		cvZero(pImg[i]);
	}

	CString * ptitle=new CString [m_levels];

	for (i=0;i<m_levels;i++)
	{
		//�ҵ�ÿ��ͼ����ɫͨ��������������ֵ
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

//���ͼ�����������
Image < CVector* >* Tensor::GetTensors() const
{
	return m_tensors;
}

////////////////////////////////////////////////////////////////////////////////////////
// ����SVD���������зֽ⣬��ȡ�������ֵ����Ӧ�������������洢��m_gradient��
////////////////////////////////////////////////////////////////////////////////////////
CMatrix ** Tensor::CvtTensorToVectorBySVD(BOOL isFiltering)
{
	//�ж��Ƿ��˲�
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
				//������ת���ɾ���洢����
				st->SetElement(0, 0, (m_tensor[i * SiNGLE_TENSOR_DIM])->GetElement(y,x));
				st->SetElement(1, 1, (m_tensor[i * SiNGLE_TENSOR_DIM + 1])->GetElement(y,x));
				st->SetElement(0, 1, (m_tensor[i * SiNGLE_TENSOR_DIM + 2])->GetElement(y,x));
				st->SetElement(1, 0, (m_tensor[i * SiNGLE_TENSOR_DIM + 2])->GetElement(y,x));
				//SVD�ֽ�
				cvSVD(&cv_st, &eval, &evec);
				signFlag = 0.0;

				for (k=0;k<m_axes_cnt;k++)
				{
					signFlag += ((m_gradient[i * m_axes_cnt + k])->GetElement(y,x) * (eigenvectors->GetElement(k,0)));
				}
				//sqrtLambda �ĳ�ʼ��
				sqrtLambda = sqrt(eigenvalues->get(0));
				//�ж������ķ����Ƿ���֮ǰ���ݶȷ���һ��
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
//��ʾ�ݶ�ͼ���ݶȵļ��㹫ʽΪ��ÿһ�������ݶȵ�ƽ���͵�ƽ����
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


//������������֮��ľ���
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

