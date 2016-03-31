#include "stdafx.h"
#include "GraphCut.h"
#include "Helper.h"

//#include <cv.h>
//#include <cxcore.h>
//#include <highgui.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

CHelper::CHelper(void)
{
}


CHelper::~CHelper(void)
{
}

void CHelper::RGB2Lab(const Image < Color_RGB > * _rgb, Image < Color_Lab > * _lab)
{
	unsigned int x,y;
	unsigned int width, height;
	width = _rgb->width();
	height = _rgb->height();
	ASSERT (width == _lab->width() && height == _lab->height());
	IplImage *cv_image = cvCreateImage( cvSize(width,height), IPL_DEPTH_32F, 3 );
	for (y = 0;y < height;y++)
	{
		for (x = 0;x < width;x++)
		{
			float* dst = &CV_IMAGE_ELEM( cv_image, float, y, x*3 );
			dst[0] = (float)(((*_rgb)(x,y)).b);
			dst[1] = (float)(((*_rgb)(x,y)).g);
			dst[2] = (float)(((*_rgb)(x,y)).r);
		}
	}
	cvCvtColor( cv_image, cv_image, CV_BGR2Lab );   //Labģʽ�������ɫ����࣬������߼��豸�޹�
	for (y = 0;y < height;y++)
	{
		for (x = 0;x < width;x++)
		{
			float* dst = &CV_IMAGE_ELEM( cv_image, float, y, x*3 );
			((*_lab)(x,y)).L = (double)(dst[0]);
			((*_lab)(x,y)).a = (double)(dst[1]);
			((*_lab)(x,y)).b = (double)(dst[2]);
		}
	}
	cvReleaseImage(&cv_image);
}

/*************************************************************************
 *
 * �������ƣ�
 *   FFT()
 *
 * ����:
 *   std::complex<double> * TD	- ָ��ʱ�������ָ��
 *   std::complex<double> * FD	- ָ��Ƶ�������ָ��
 *   r						- 2��������������������2��r�η�=N
 *
 * ˵��:
 *   �ú�������ʵ�ֿ��ٸ���Ҷ�任��
 *
 ************************************************************************/
void CHelper::FFT(std::complex<double> * TD, std::complex<double> * FD, int r)
{
	// ����Ҷ�任����
	long	count;
	
	// ѭ������
	int		i,j,k;
	
	// �м����
	int		bfsize,p;
	
	// �Ƕ�
	double	angle;
	
	std::complex<double> *W,*X1,*X2,*X;
	
	// ���㸶��Ҷ�任����
	count = 1 << r;     //1����rλ����2��r�η�
	
	// ������������洢��
	W  = new std::complex<double>[count / 2];
	X1 = new std::complex<double>[count];
	X2 = new std::complex<double>[count];
	
	// �����Ȩϵ��
	for(i = 0; i < count / 2; i++)
	{
		angle = -i * PI * 2 / count;
		W[i] = std::complex<double> (cos(angle), sin(angle));
	}
	
	// ��ʱ���д��X1
	memcpy(X1, TD, sizeof(std::complex<double>) * count);
	
	// ���õ����㷨���п��ٸ���Ҷ�任
	for(k = 0; k < r; k++)  //FFT�ļ���
	{
		for(j = 0; j < 1 << k; j++)   //����
		{
			bfsize = 1 << (r-k);  //�����Ŀ��*2
			for(i = 0; i < bfsize / 2; i++)
			{
				p = j * bfsize;
				X2[i + p] = X1[i + p] + X1[i + p +bfsize/2];
				X2[i + p + bfsize / 2] = (X1[i + p] - X1[i + p + bfsize / 2]) * W[i * (1<<k)];
			}
		}
		X  = X1;
		X1 = X2;
		X2 = X;
	}
	
	// ��������(�������λ�����)
	for(j = 0; j < count; j++)
	{
		p = 0;
		for(i = 0; i < r; i++)
		{
			if (j&(1<<i))
			{
				p+=1<<(r-i-1);
			}
		}
		FD[j]=X1[p];
	}
	
	// �ͷ��ڴ�
	delete W;
	delete X1;
	delete X2;
}


/*************************************************************************
 *
 * �������ƣ�
 *   IFFT()
 *
 * ����:
 *   std::complex<double> * FD	- ָ��Ƶ��ֵ��ָ��
 *   std::complex<double> * TD	- ָ��ʱ��ֵ��ָ��
 *   r						- 2��������������������2��r�η�=N
 *
 * ˵��:
 *   �ú�������ʵ�ֿ��ٸ���Ҷ���任��
 *
 ************************************************************************/
void CHelper::IFFT(std::complex<double> * FD, std::complex<double> * TD, int r)
{
	// ����Ҷ�任����
	long	count;
	
	// ѭ������
	int		i;
	
	std::complex<double> *X;
	
	// ���㸶��Ҷ�任����
	count = 1 << r;   //1����rλ����2��r�η�
	
	// ������������洢��
	X = new std::complex<double>[count];
	
	// ��Ƶ���д��X
	memcpy(X, FD, sizeof(std::complex<double>) * count);

	// ����
	for(i = 0; i < count; i++)
	{
		X[i] = std::complex<double> (X[i].real(), -X[i].imag());
	}
	
	// ���ÿ��ٸ���Ҷ�任
	FFT(X, TD, r);
	
	// ��ʱ���Ĺ���
	for(i = 0; i < count; i++)
	{
		TD[i] = std::complex<double> (TD[i].real() / count, -TD[i].imag() / count);
	}
	
	// �ͷ��ڴ�
	delete X;
}

/*************************************************************************
 *
 * �������ƣ�
 *   CONV()
 *
 * ����:
 *   std::complex<double> * TD1	    - ָ��ʱ����������1��ָ��
 *   std::complex<double> * TD2	    - ָ��ʱ����������2��ָ��
 *   std::complex<double> * TDout	- ָ��ʱ���������ָ��
 *   M						    - ����1�ĳ���
 *   N                          - ����2�ĳ���
 *
 * ˵��:
 *   �ú�������FFTʵ�ֿ��پ����
 *
 ************************************************************************/
void CHelper::Convolution(std::complex<double> * TD1, std::complex<double> * TD2, std::complex<double> * TDout, int M, int N)
{
	// ����������
	int	count = M+N-1;

	// ����ʹ��FFT����count��չΪ2����
	int Lcount;
    int r=0;  // 2����������FFT����������2��r�η�=Lcount

	int temp;
	if (log((float)count)/log((float)2)-int(log((float)count)/log((float)2))==0)
      temp = log((float)count)/log((float)2);
	else
	  temp = log((float)count)/log((float)2)+1;
	r = temp;
	Lcount = 1<<r;	

	// ������������洢��
    std::complex<double> *X1, *X2, *FD1, *FD2, *FD12, *TD12;

	X1 = new std::complex<double>[Lcount];  //����������1
	X2 = new std::complex<double>[Lcount];  //����������2
	FD1 = new std::complex<double>[Lcount];   //����1�ĸ���Ҷ�任���
	FD2 = new std::complex<double>[Lcount];   //����2�ĸ���Ҷ�任���
	FD12 = new std::complex<double>[Lcount];   //����1,2��Ƶ����˽��
	TD12 = new std::complex<double>[Lcount];   //����1,2�ĸ���Ҷ���任���
	
    //�����в���ΪLcount����
	std::complex<double> *X, *Y;
	X = new std::complex<double>[M];  //��ʱ�洢��
	Y = new std::complex<double>[N];
	
	// ��ʱ���д��X,Y
	memcpy(X, TD1, sizeof(std::complex<double>) * M);
	memcpy(Y, TD2, sizeof(std::complex<double>) * N);

	// ѭ������
	int	i;

    for (i=0; i<M; i++)    //��������1����   
    {
        X1[i] = std::complex<double>(X[i].real(), X[i].imag());                                             
    }

    for (i=M; i<Lcount; i++)    //����1��0
    {
        X1[i] = std::complex<double>(0, 0);
    }

	for (i=0; i<N; i++)    //��������2����   
    {
        X2[i] = std::complex<double>(Y[i].real(), Y[i].imag());
    }

    for (i=N; i<Lcount; i++)    //����2��0
    {
        X2[i] = std::complex<double>(0, 0);
    }

    // �ͷ��ڴ�
	delete X;
	delete Y;

    //����1��FFT
	FFT(X1, FD1, r);

	//����2��FFT
    FFT(X2, FD2, r);

    //����1,2��Ƶ�����
    for (i=0; i<Lcount; i++)    //����1,2���
    {
        FD12[i] = std::complex<double>(FD1[i].real()*FD2[i].real()-FD1[i].imag()*FD2[i].imag(), FD1[i].real()*FD2[i].imag()+FD1[i].imag()*FD2[i].real());
    }

	//����1,2��Ƶ����˵�IFFT
    IFFT(FD12, TD12, r);

	//TD12�е�ǰM+N-1��Ϊ����������д��TDout
    memcpy(TDout, TD12, sizeof(std::complex<double>)*count);
	
	// �ͷ��ڴ�
	delete X1;
	delete X2;
	delete FD1;
	delete FD2;
	delete FD12;
    delete TD12;
}

//TDout��TD1�ĳ���һ����ȡM+N-1���м�M��
void CHelper::Convolution(double * TD1, double * TD2, double * TDout, int M, int N)
{
	std::complex<double>* c_TD1 = new std::complex<double>[M];
	std::complex<double>* c_TD2 = new std::complex<double>[N];
	std::complex<double>* c_TDout = new std::complex<double>[M+N-1];
	int i;
	for (i=0;i<M;i++)
	{
		c_TD1[i] = std::complex<double>(TD1[i], 0);
	}
	for (i=0;i<N;i++)
	{
		c_TD2[i] = std::complex<double>(TD2[i], 0);
	}
	Convolution(c_TD1, c_TD2, c_TDout, M, N);
	int pos = cvCeil( (N - 1) / 2.0 );
	for (i=0;i<M;i++)
	{
		TDout[i] = c_TDout[pos++].real();
	}
	delete c_TD1;
	delete c_TD2;
	delete c_TDout;
}
