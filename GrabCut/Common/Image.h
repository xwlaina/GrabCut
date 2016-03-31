/*
 * GrabCut implementation source code Copyright(c) 2005-2006 Justin Talbot
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 */

#ifndef IMAGE_H
#define IMAGE_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include "Global.h"
#include "Color_RGB.h"
// #include <vld.h>

// Images, really just a templatized 2D array. We use this for all the image variables.

template <class T> 
class Image
{

public:
        Image(unsigned int width, unsigned int height);
        ~Image();

        T *ptr()
        {
                return m_image;
        }

        T &operator()(int x, int y)
        {
                clampX(x);
                clampY(y);
                return m_image[y *m_width + x];   //返回图中某个位置...
        }
        const T &operator()(int x, int y)const
        {
                clampX(x);
                clampY(y);
                return m_image[y *m_width + x];
        }

        void fillRectangle(int x1, int y1, int x2, int y2, const T &t);
        void fill(const T &t);
		void fillCircle(int x, int y, int radius, const T &t);

        unsigned int width()const
        {
                return m_width;
        }
        unsigned int height()const
        {
                return m_height;
        }
private:

        void clampX(int &x)const
        {
                if (x < 0)
                {
                        x = 0;
                }
                if (x >= (int)m_width)
                {
                        x = m_width - 1;
                }
        }
        void clampY(int &y)const
        {
                if (y < 0)
                {
                        y = 0;
                }
                if (y >= (int)m_height)
                {
                        y = m_height - 1;
                }
        }
		
		friend void GetGradient(const Image < Color_RGB >  *image, FLOAT *deltar, FLOAT *deltasita);

        T *m_image;
        unsigned int m_width, m_height;
};


// Image member functions
template <class T> 
inline Image<T>::Image(unsigned int width, unsigned int height): m_width(width), m_height(height)
{
        m_image = new T[m_width *m_height];    //MonKey 有析构么？   
}

//-------------------------------------------------------------------------

template <class T> 
inline Image<T>::~Image()
{
        if (m_image!=NULL)
        {
                delete [] m_image;
				m_image = NULL;
        }
}

//-------------------------------------------------------------------------

template <class T> 
inline void Image<T>::fillRectangle(int x1, int y1, int x2, int y2, const T &t)
{
        clampX(x1);
        clampY(y1);
        clampX(x2);
        clampY(y2);

        if (y1> y2)
        {
                int t = y1;
                y1 = y2;
                y2 = t;
        }
        if (x1 > x2)
        {
                int t = x1;
                x1 = x2;
                x2 = t;
        }

        for (int i = y1; i <= y2; ++i)
        {
                for (int j = x1; j <= x2; ++j)
                {
                        m_image[i *m_width + j] = t;
                }
        }
}

//-------------------------------------------------------------------------

template <class T> 
inline void Image<T>::fill(const T &t)
{
        for (unsigned int i = 0; i<m_width *m_height; ++i)
        {
                m_image[i] = t;
        }
}

template <class T> 
inline void Image<T>::fillCircle(int x, int y, int radius, const T &t)
{
	//for (int i = x - radius; i <= x + radius; i++)
	//{
	//	for (int j = (int)(y - sqrt((pow((double)radius,2)-pow((double)(i-x),2))); j <= (int)(y + sqrt(pow((double)radius,2)-pow((double)(i-x),2))); j++)
	//	{
	//		if ((i >= 0) && (i < (int)m_width) && (j >= 0) && (j < (int)m_height))
	//		{
	//			m_image[j *m_width + i] = t;
	//		}
	//	}
	//}
	//MonKey 这里的逻辑判断顺序有问题
	//但是边界一般都会满足
	for (INT i=x-radius ; i<=x+radius ; i++)
	{
		for(INT j=y-radius ; j<=y+radius ; j++)
			if ( (i >= 0) && (i < (int)m_width) && 
				(j >= 0) && (j < (int)m_height)&&
				((i-x)*(i-x)+(j-y)*(j-y)<=radius*radius)
				)
			{
				m_image[j *m_width + i] = t;
			}
	}
}
//-------------------------------------------------------------------------

//得到输入图像的梯度
inline void GetGradient(const Image < Color_RGB >  *image, FLOAT *deltar, FLOAT *deltasita = NULL)
{
	//下面计算各像素在水平和垂直方向上的梯度,边缘点梯度计为0；
	INT* deltaxarr;
	INT* deltayarr;
	INT width = image->width();
	INT height = image->height();
	INT deltacount = width * height;
	deltaxarr = new INT[deltacount];
	deltayarr = new INT[deltacount];
	INT x,y;
	
    //暂不计算边缘点；
	for (y=1; y<height-1; y++)
	{
		for (x=1; x<width-1; x++)
		{
			INT deltaarrpos = y*width + x;//在梯度数组中的位置；
			//卷积计算；
			deltaxarr[deltaarrpos] = (INT) ( (
				((*image)(x+1,y-1)).g //右上
				+ ((*image)(x+1,y)).g //右
				+ ((*image)(x+1,y+1)).g //右下
				- ((*image)(x-1,y-1)).g //左上
				- ((*image)(x-1,y)).g //左
				- ((*image)(x-1,y+1)).g ) / 3 );//左下
			deltayarr[deltaarrpos] = (INT) ( ( 
				((*image)(x+1,y-1)).g //右上
				+ ((*image)(x,y-1)).g //上
				+ ((*image)(x-1,y-1)).g //左上
				- ((*image)(x-1,y+1)).g //左下
				- ((*image)(x,y+1)).g //下
				- ((*image)(x+1,y+1)).g) / 3 );//右下
		}
	}
	
	//边缘赋为其内侧点的值；
	for (y=0; y<height; y++)
	{
		INT x1 = 0;
		INT pos1 = y*width + x1;
		deltaxarr[pos1] = deltaxarr[pos1+1];
		deltayarr[pos1] = deltayarr[pos1+1];
		INT x2 = width-1;
		INT pos2 = y*width + x2;
		deltaxarr[pos2] = deltaxarr[pos2-1];
		deltayarr[pos2] = deltayarr[pos2-1];
	}
	for (x=0; x<width; x++)
	{
		INT y1 = 0;
		INT pos1 = x;
		INT inner = x + width;//下一行；
		deltaxarr[pos1] = deltaxarr[inner];
		deltayarr[pos1] = deltayarr[inner];
		INT y2 = height-1;
		INT pos2 = y2*width + x;
		inner = pos2 - width;//上一行；
		deltaxarr[pos2] = deltaxarr[inner];
		deltayarr[pos2] = deltayarr[inner];
	}
	
	
	for (y=0; y<height; y++)
	{
		for (x=0; x<width; x++)
		{
			INT temppos = y*width + x;
			if ( (deltaxarr[temppos])==0 )
			{
				if (deltayarr[temppos]!=0)
				{
					if (deltasita != NULL)
					{
						deltasita[temppos] = 0;//水平方向;
					}
					deltar[temppos] = (FLOAT) abs(deltayarr[temppos]);
				}else
				{
					if (deltasita != NULL)
					{
						deltasita[temppos] = -1;//无确定方向;
					}
					deltar[temppos] = (FLOAT) abs(deltayarr[temppos]);
				}
				continue;
			}
			if (deltasita != NULL)
			{
				deltasita[temppos] = (FLOAT) ( atan( 
					(FLOAT)deltayarr[temppos]
					/ (FLOAT)deltaxarr[temppos] ) + PI/2. );
			}
			deltar[temppos] = (FLOAT) sqrt((DOUBLE) 
				( deltayarr[temppos]*deltayarr[temppos]
				+ deltaxarr[temppos]*deltaxarr[temppos] ) );
		}
	}
	
	delete [] deltaxarr; deltaxarr = NULL; //删除水平和垂直梯度数组；
	delete [] deltayarr; deltayarr = NULL;
}

inline void GetGradient(const IplImage *cv_image)//, FLOAT *deltar)
{
	IplImage *cv_gray_image = cvCreateImage( cvGetSize(cv_image), cv_image->depth, 1 );
	cvCvtColor( cv_image, cv_gray_image, CV_BGR2GRAY );
	IplImage *cv_gradient_X = cvCreateImage(cvGetSize(cv_gray_image),IPL_DEPTH_32F, cv_gray_image->nChannels); 
	IplImage *cv_gradient_Y = cvCreateImage(cvGetSize(cv_gray_image),IPL_DEPTH_32F, cv_gray_image->nChannels); 
	IplImage *cv_result = cvCreateImage(cvGetSize(cv_gray_image),cv_gray_image->depth, cv_gray_image->nChannels); 
	cv_gradient_X->origin = cv_gray_image->origin;
	cv_gradient_Y->origin = cv_gray_image->origin;
	cvSobel(cv_gray_image,cv_gradient_X,1,0,3);
	cvSobel(cv_gray_image,cv_gradient_Y,0,1,3);
	cvMul(cv_gradient_X,cv_gradient_X,cv_gradient_X,1.0);
	cvMul(cv_gradient_Y,cv_gradient_Y,cv_gradient_Y,1.0);
	cvAdd(cv_gradient_X,cv_gradient_Y,cv_gradient_X,0); 
	cvPow(cv_gradient_X,cv_gradient_Y,0.5);
	cvConvert(cv_gradient_Y, cv_result);
	cvNamedWindow("Result", 1 );
	cvShowImage( "Result", cv_result);
	cvReleaseImage(&cv_gray_image);
	cvReleaseImage(&cv_gradient_X);
	cvReleaseImage(&cv_gradient_Y);
	cvReleaseImage(&cv_result);
}

#endif
