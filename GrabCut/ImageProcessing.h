// ImageProcessing.h: interface for the CImageProcessing class.
//
//////////////////////////////////////////////////////////////////////

#ifndef IMAGEPROCESSING_H
#define IMAGEPROCESSING_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#pragma warning(disable :4786)

#include "cv.h"
#include "highgui.h"
#include "GrabCutTool.h"
#include <vector>
#include <set>
#include "CvvImage.h"
#include "Color_Lab.h"
#include "Color_RGB.h"
#include "Image.h"

using namespace std;


class CImageProcessing  
{
public:
	CImageProcessing();
	virtual ~CImageProcessing();

public:
	BOOL isImageLoaded();
	LONG GetHeight();
	LONG GetWidth();
	Image < Color_RGB >  * GetImageData_RGB();
	Image < Color_Lab >  * GetImageData_Lab();
	IplImage * GetCVImage();
	void Destroy(void);
	BOOL SaveImage(LPCTSTR lpszPathName);
	BOOL LoadImage(LPCTSTR lpszPathName);
public:
	GrabCutTool *grabCut;   

		
private:
	BOOL  imageLoaded;
	LONG  imageWidth;
	LONG  imageHeight;
	Image < Color_RGB >  * imageData_RGB;
	Image < Color_Lab >  * imageData_Lab;
	CvvImage myImageObject;  

};


#endif 