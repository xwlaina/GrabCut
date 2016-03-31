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
	cvCvtColor( cv_image, cv_image, CV_BGR2Lab );   //Lab模式所定义的色彩最多，且与光线及设备无关
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
 * 函数名称：
 *   FFT()
 *
 * 参数:
 *   std::complex<double> * TD	- 指向时域数组的指针
 *   std::complex<double> * FD	- 指向频域数组的指针
 *   r						- 2的幂数，即迭代次数，2的r次方=N
 *
 * 说明:
 *   该函数用来实现快速付立叶变换。
 *
 ************************************************************************/
void CHelper::FFT(std::complex<double> * TD, std::complex<double> * FD, int r)
{
	// 付立叶变换点数
	long	count;
	
	// 循环变量
	int		i,j,k;
	
	// 中间变量
	int		bfsize,p;
	
	// 角度
	double	angle;
	
	std::complex<double> *W,*X1,*X2,*X;
	
	// 计算付立叶变换点数
	count = 1 << r;     //1左移r位，即2的r次方
	
	// 分配运算所需存储器
	W  = new std::complex<double>[count / 2];
	X1 = new std::complex<double>[count];
	X2 = new std::complex<double>[count];
	
	// 计算加权系数
	for(i = 0; i < count / 2; i++)
	{
		angle = -i * PI * 2 / count;
		W[i] = std::complex<double> (cos(angle), sin(angle));
	}
	
	// 将时域点写入X1
	memcpy(X1, TD, sizeof(std::complex<double>) * count);
	
	// 采用蝶形算法进行快速付立叶变换
	for(k = 0; k < r; k++)  //FFT的级数
	{
		for(j = 0; j < 1 << k; j++)   //组数
		{
			bfsize = 1 << (r-k);  //本级的跨度*2
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
	
	// 重新排序(输出是码位倒序的)
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
	
	// 释放内存
	delete W;
	delete X1;
	delete X2;
}


/*************************************************************************
 *
 * 函数名称：
 *   IFFT()
 *
 * 参数:
 *   std::complex<double> * FD	- 指向频域值的指针
 *   std::complex<double> * TD	- 指向时域值的指针
 *   r						- 2的幂数，即迭代次数，2的r次方=N
 *
 * 说明:
 *   该函数用来实现快速付立叶反变换。
 *
 ************************************************************************/
void CHelper::IFFT(std::complex<double> * FD, std::complex<double> * TD, int r)
{
	// 付立叶变换点数
	long	count;
	
	// 循环变量
	int		i;
	
	std::complex<double> *X;
	
	// 计算付立叶变换点数
	count = 1 << r;   //1左移r位，即2的r次方
	
	// 分配运算所需存储器
	X = new std::complex<double>[count];
	
	// 将频域点写入X
	memcpy(X, FD, sizeof(std::complex<double>) * count);

	// 求共轭
	for(i = 0; i < count; i++)
	{
		X[i] = std::complex<double> (X[i].real(), -X[i].imag());
	}
	
	// 调用快速付立叶变换
	FFT(X, TD, r);
	
	// 求时域点的共轭
	for(i = 0; i < count; i++)
	{
		TD[i] = std::complex<double> (TD[i].real() / count, -TD[i].imag() / count);
	}
	
	// 释放内存
	delete X;
}

/*************************************************************************
 *
 * 函数名称：
 *   CONV()
 *
 * 参数:
 *   std::complex<double> * TD1	    - 指向时域序列数组1的指针
 *   std::complex<double> * TD2	    - 指向时域序列数组2的指针
 *   std::complex<double> * TDout	- 指向时域结果数组的指针
 *   M						    - 序列1的长度
 *   N                          - 序列2的长度
 *
 * 说明:
 *   该函数利用FFT实现快速卷积。
 *
 ************************************************************************/
void CHelper::Convolution(std::complex<double> * TD1, std::complex<double> * TD2, std::complex<double> * TDout, int M, int N)
{
	// 卷积结果长度
	int	count = M+N-1;

	// 便于使用FFT，把count扩展为2的幂
	int Lcount;
    int r=0;  // 2的幂数，即FFT迭代次数，2的r次方=Lcount

	int temp;
	if (log((float)count)/log((float)2)-int(log((float)count)/log((float)2))==0)
      temp = log((float)count)/log((float)2);
	else
	  temp = log((float)count)/log((float)2)+1;
	r = temp;
	Lcount = 1<<r;	

	// 分配运算所需存储器
    std::complex<double> *X1, *X2, *FD1, *FD2, *FD12, *TD12;

	X1 = new std::complex<double>[Lcount];  //补齐后的序列1
	X2 = new std::complex<double>[Lcount];  //补齐后的序列2
	FD1 = new std::complex<double>[Lcount];   //序列1的傅立叶变换结果
	FD2 = new std::complex<double>[Lcount];   //序列2的傅立叶变换结果
	FD12 = new std::complex<double>[Lcount];   //序列1,2的频域相乘结果
	TD12 = new std::complex<double>[Lcount];   //序列1,2的傅立叶反变换结果
	
    //将序列补齐为Lcount长度
	std::complex<double> *X, *Y;
	X = new std::complex<double>[M];  //临时存储器
	Y = new std::complex<double>[N];
	
	// 将时域点写入X,Y
	memcpy(X, TD1, sizeof(std::complex<double>) * M);
	memcpy(Y, TD2, sizeof(std::complex<double>) * N);

	// 循环变量
	int	i;

    for (i=0; i<M; i++)    //拷贝序列1内容   
    {
        X1[i] = std::complex<double>(X[i].real(), X[i].imag());                                             
    }

    for (i=M; i<Lcount; i++)    //序列1补0
    {
        X1[i] = std::complex<double>(0, 0);
    }

	for (i=0; i<N; i++)    //拷贝序列2内容   
    {
        X2[i] = std::complex<double>(Y[i].real(), Y[i].imag());
    }

    for (i=N; i<Lcount; i++)    //序列2补0
    {
        X2[i] = std::complex<double>(0, 0);
    }

    // 释放内存
	delete X;
	delete Y;

    //序列1的FFT
	FFT(X1, FD1, r);

	//序列2的FFT
    FFT(X2, FD2, r);

    //序列1,2的频域相乘
    for (i=0; i<Lcount; i++)    //序列1,2相乘
    {
        FD12[i] = std::complex<double>(FD1[i].real()*FD2[i].real()-FD1[i].imag()*FD2[i].imag(), FD1[i].real()*FD2[i].imag()+FD1[i].imag()*FD2[i].real());
    }

	//序列1,2的频域相乘的IFFT
    IFFT(FD12, TD12, r);

	//TD12中的前M+N-1项为真正卷积结果写入TDout
    memcpy(TDout, TD12, sizeof(std::complex<double>)*count);
	
	// 释放内存
	delete X1;
	delete X2;
	delete FD1;
	delete FD2;
	delete FD12;
    delete TD12;
}

//TDout与TD1的长度一样，取M+N-1的中间M项
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
