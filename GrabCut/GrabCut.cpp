// GrabCut.cpp : 定义控制台应用程序的入口点。
//

#include "stdafx.h"
#include <iostream>
#include <queue>
#include <opencv2\opencv.hpp> 
#include "CvvImage.h"
#include "Color_Lab.h"
#include "Color_RGB.h"
#include "Image.h"
#include "GrabCutTool.h"
#include "ImageProcessing.h"
using namespace std;

IplImage* saveimage(int m_w , int m_h ,const Image < SegmentationValue >  &hardSegmentation);
CvRect findrectpoint(const IplImage* im);
IplImage* graythresh(IplImage* image);

int _tmain(int argc, _TCHAR* argv[])
{ 
	queue<CvRect> rect;
	CFileFind finderrect;
	BOOL ans = finderrect.FindFile(_T("C:\\sources\\grabcut\\boundary_GT_rect\\*.*"));
	int lengthrect =0;
	while (ans)
	{
		ans = finderrect.FindNextFile();
		if (finderrect.IsDots())   continue;
		lengthrect++;
		CString filename = finderrect.GetFileName();
		CString LP = _T("C:\\sources\\grabcut\\boundary_GT_rect\\")+filename;
		CStringA strPathA(LP);
		const char *p = strPathA.GetString();
		//cvSaveImage("C:\\Users\\xu\\Desktop\\image\\2.jpg",cvLoadImage(p,-1));
		IplImage* image = graythresh(cvLoadImage(p,-1));
		//cvSaveImage("C:\\Users\\xu\\Desktop\\image\\myimg.jpg",image);
		CvRect rectpoint = findrectpoint(image);
		rect.push(rectpoint);
	}
	cout<<lengthrect<<endl;
	cout<<rect.size()<<endl;
	finderrect.Close();

	CFileFind finder;
	BOOL res = finder.FindFile(_T("C:\\sources\\grabcut\\data_GT\\*.*"));
	int length =0;
	while (res)
	{
		res = finder.FindNextFile();
		if (finder.IsDots())   continue;
		length++;
		CString filename = finder.GetFileName();
		CString LP = _T("C:\\sources\\grabcut\\data_GT\\")+filename;
		CImageProcessing m_ImageProcessing;
		m_ImageProcessing.LoadImage(LP);
		CvRect fl;
		if(!rect.empty())
		{
			fl = rect.front();
			rect.pop();
		}
		m_ImageProcessing.grabCut->ImageSegmentationByAGMMColor(fl.y,fl.x,(fl.y+fl.width-1),(fl.x+fl.height-1));
		//m_ImageProcessing.grabCut->ImageSegmentationByColor(fl.x,fl.y,(fl.x+fl.height-1),(fl.y+fl.width-1));
		IplImage* re = saveimage(m_ImageProcessing.GetWidth(),m_ImageProcessing.GetHeight(),*(m_ImageProcessing.grabCut->getSegmentationValue_AGMMcolor()));
		CString L1 = _T("C:\\sources\\grabcut\\data_segment(epsilon=0.9)\\")+filename;
		CStringA strPathA(L1);
		const char *p = strPathA.GetString();
		cvSaveImage(p, re);
	}
	finder.Close();
	cout<<length<<endl;
	system("pause");
	return 0;
}

IplImage* saveimage(int m_w , int m_h ,const Image < SegmentationValue >  &hardSegmentation)
{
	IplImage* src = cvCreateImage( cvSize(m_w,m_h), 8, 1 );
	cvZero(src);
	int i,j;
	LONG width,height;
	width = m_w;
	height = m_h;
	for ( i=0; i<height; i++)
	{
		for ( j=0; j<width; j++)
		{
			uchar* dst = &CV_IMAGE_ELEM( src, uchar, i, j );
			dst[0] = (uchar)((hardSegmentation(j,i)==SegmentationForeground)?255:0);
		}
	}
	return src;
}

IplImage* graythresh(IplImage* image)
{
	IplImage* src = cvCloneImage(image);
	int i,j;
	LONG width,height;
	width = image->width;
	height = image->height;
	for ( i=0; i<height; i++)
	{
		for ( j=0; j<width; j++)
		{
			uchar* dst = &CV_IMAGE_ELEM( src, uchar, i, j );
			if(dst[0] ==(uchar)128)
				dst[0] = (uchar)255;
			else
				dst[0] = (uchar)0;
		}
	}
	return src;
}

CvRect findrectpoint(const IplImage* im)
{
	CvRect res;
	CvPoint leftup,rightdown;
	int i,j;
	LONG width,height;
	width = im->width;
	height = im->height;
	int flag=0;
	for ( i=0; i<height; i++)
	{
		if(flag)
			break;
		for ( j=0; j<width; j++)
		{
			uchar* dst = &CV_IMAGE_ELEM( im, uchar, i, j );
			if(dst[0] != (uchar)0)
			{
				leftup.x = i;
				leftup.y = j;
				cout<<i<<j<<endl;
				flag =1;
				break;
			}
		}
	}
	flag =0;
	for ( i=height-1; i>=0; i--)
	{
		if(flag)
			break;
		for ( j=width-1; j>=0; j--)
		{
			uchar* dst = &CV_IMAGE_ELEM( im, uchar, i, j );
			if(dst[0] != (uchar)0)
			{
				rightdown.x = i;
				rightdown.y = j;
				cout<<i<<j<<endl;
				flag = 1;
				break;
			}
		}
	}
	res.x = leftup.x;
	res.y = leftup.y;
	res.height = rightdown.x - leftup.x +1 ;
	res.width = rightdown.y - leftup.y +1 ;
	return res;
}
