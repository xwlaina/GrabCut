// ImageProcessing.cpp: implementation of the CImageProcessing class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "GraphCut.h"
#include "ImageProcessing.h"
#include "MainFrm.h"
#include "GraphCutView.h"
#include "Helper.h"
//#include "../Tensor.h"
// #include <vld.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
CImageProcessing::CImageProcessing()
{
	imageWidth = -1;
	imageHeight = -1;
	imageLoaded = FALSE;
	imageData_RGB = NULL;
	imageData_Lab = NULL;


	grabCut = NULL;

	imgMatting = NULL;

	m_mattingMask = NULL;
	
	//m_gabor = NULL;
}

CImageProcessing::~CImageProcessing()
{
	Destroy();
}

BOOL CImageProcessing::LoadImage(LPCTSTR lpszPathName)
{
	Destroy();
	m_szPathName = lpszPathName;

	CString szPathName=lpszPathName;	
	CStringA strPathA(szPathName);
	const char *p = strPathA.GetString();
	//MonKey VS2010 CSrting转char*

	if (!myImageObject.Load(p)) //cvvImage Load图像
	{                                             //MonKey  转型有问题
		//AfxMessageBox((LPCTSTR)cvErrorStr(cvGetErrStatus()));
		AfxMessageBox(_T("打开失败"));
		myImageObject.Destroy();
		return FALSE;
	}
	
	imageLoaded = TRUE;
	
	//以下将RGB数据存入数组以备处理；
	LONG width = myImageObject.Width();
	LONG height = myImageObject.Height();
	imageWidth = width;
	imageHeight = height; 
	if ( imageData_RGB != NULL )     //Image< Color_RGB > 与 Image< Color_Lab >  数据的清空
	{
		delete imageData_RGB;
		imageData_RGB = NULL;
	}
	if ( imageData_Lab != NULL )
	{
		delete imageData_Lab;
		imageData_Lab = NULL;
	}

	//以下将OpenCV中的imageData转化为width*height*3的矩阵，去掉因为widthStep而添加的0字节，存储顺序为BGR
	imageData_RGB = new Image < Color_RGB > (width, height);
	imageData_Lab = new Image < Color_Lab > (width, height);

	BYTE* cv_imageData = (BYTE*)(myImageObject.GetImage()->imageData);
	ASSERT(cv_imageData!=NULL);
    INT widthstep = myImageObject.GetImage()->widthStep;
	for (INT y=0; y<height; y++)
	{
		for (INT x=0; x<width; x++)
		{
			INT inarrpos = y*widthstep + x*3;
			(*imageData_RGB)(x, y) = Color_RGB(cv_imageData[inarrpos+2],cv_imageData[inarrpos+1],cv_imageData[inarrpos]);
		}
	}
	//CHelper helper;
	//helper.RGB2Lab(imageData_RGB, imageData_Lab);
	//MonKey 原RGB转Lab计算错误
	IplImage *img_BGR=cvLoadImage(p,1);
	IplImage *img_Lab=cvCreateImage(cvGetSize(img_BGR),IPL_DEPTH_8U,3);
	cvCvtColor( img_BGR, img_Lab, CV_BGR2Lab );
	//cvSaveImage("C:\\Users\\xu\\Desktop\\LAB\\1.bmp",img_Lab);
	BYTE* cv_imageLab = (BYTE*)(img_Lab->imageData);
	ASSERT(cv_imageData!=NULL);
	INT widthstepLab = img_Lab->widthStep;
	for (INT j=0; j<height; j++)
	{
		for (INT i=0; i<width; i++)
		{
			INT inarrpos = j*widthstep + i*3;
			(*imageData_Lab)(i, j) = Color_Lab(cv_imageData[inarrpos],cv_imageData[inarrpos+1],cv_imageData[inarrpos+2]);
		}
	}


//	m_gabor = new Gabor(GetCVImage());  //MonKey 纯颜色GC 没Gabor什么事
//	m_gabor->GaborFilteredImg();

	CMainFrame* pFrame = (CMainFrame*) AfxGetApp()->GetMainWnd();
	CGraphCutDoc* pDoc = ((CGraphCutView*)pFrame->GetPane(0,1))->GetDocument();
	ASSERT_VALID(pDoc);

	switch(pDoc->m_eCurCmd)
	{

	case NoneCmd:
		break;
	default:
		if (grabCut!=NULL)
		{
			delete grabCut;
			grabCut = NULL;
		}
		//Tensor *ts = new Tensor(GetCVImage(), TRUE);
		grabCut = new GrabCutTool(imageData_Lab,GetCVImage());  //GrabCut中有用到Gabor么？    //MonKey 有析构么？
		break;                                              //前面m_gabor 都注释掉了，这里传入构造函数的值NULL
                                                            //GrabCutTool有针对(XXX,null)的构造函数么？
	}

	m_mattingMask = new Image < unsigned int > (imageWidth, imageHeight);

	return TRUE;
}

BOOL CImageProcessing::SaveImage(LPCTSTR lpszPathName)
{
	CString szPathName=lpszPathName;  
	CStringA strPathA(szPathName);
	const char *p = strPathA.GetString();
	return myImageObject.Save(p);
}

void CImageProcessing::Destroy()
{
	if (imageData_RGB!=NULL)
	{
		delete imageData_RGB;
		imageData_RGB = NULL;
	}
	if (imageData_Lab!=NULL)
	{
		delete imageData_Lab;
		imageData_Lab = NULL;
	}

	myImageObject.Destroy();


	imageLoaded = FALSE;
	imageWidth = -1;
	imageHeight = -1;

	if (grabCut!=NULL)
	{
		delete grabCut;   //MonKey  
		grabCut = NULL;
	}

	if (m_mattingMask!=NULL)
	{
		delete m_mattingMask;
		m_mattingMask = NULL;
	}

	/*if (m_gabor != NULL)
	{
	delete m_gabor;
	m_gabor = NULL;
	}*/
}

void CImageProcessing::Show(HDC dc, int x, int y, int width, int height, int from_x, int from_y)
{
	myImageObject.Show(dc,x,y,width,height,from_x,from_y);
}

LONG CImageProcessing::GetWidth()
{
	return imageWidth;
}

Image < Color_RGB >  * CImageProcessing::GetImageData_RGB()
{
	return imageData_RGB;
}

Image < Color_Lab >  * CImageProcessing::GetImageData_Lab()
{
	return imageData_Lab;
}

LONG CImageProcessing::GetHeight()
{
	return imageHeight;
}

void CImageProcessing::DrawToHDC(HDC hDCDst, RECT *pDstRect)   //Show 与 DrawToHDC "重载关系"...
{
	myImageObject.DrawToHDC(hDCDst,pDstRect);
}

BOOL CImageProcessing::isImageLoaded()
{
	return imageLoaded;
}


extern int option;//定义在segmentMS头文件中；




void CImageProcessing::AddOutput(int flag, CString s)
{
	CMainFrame* pFrame = (CMainFrame*) AfxGetApp()->GetMainWnd();
	pFrame->myOutputBar.AddOutput(flag,s);
	pFrame->myOutputBar.Invalidate();
	CheckMessageQueue();
}

BOOL CImageProcessing::CheckMessageQueue()
{
	MSG msg;    
	while(PeekMessage(&msg,NULL,0,0,PM_REMOVE))
    { 
		if(msg.message==WM_QUIT) 
			return FALSE; 
		TranslateMessage(&msg); 
		DispatchMessage(&msg); 
    } 
	return TRUE;
}

void CImageProcessing::Display(const Image < SegmentationValue >  &hardSegmentation)
{
	IplImage* src = cvCreateImage( cvGetSize(myImageObject.GetImage()), 8, 1 );
	cvZero(src);
	int i,j;
	LONG width,height;
	width = myImageObject.Width();
	height = myImageObject.Height();
	for ( i=0; i<height; i++)
	{
		for ( j=0; j<width; j++)
		{
			uchar* dst = &CV_IMAGE_ELEM( src, uchar, i, j );
			dst[0] = (uchar)((hardSegmentation(j,i)==SegmentationForeground)?255:0);

			dst = &CV_IMAGE_ELEM( myImageObject.GetImage(), uchar, i, j*3 );
			Color_RGB mask;
			if (hardSegmentation(j,i)==SegmentationForeground)
			{
				mask = (*imageData_RGB)(j,i);
			} 
			else
			{
				mask = ((*imageData_RGB)(j,i))*0.5;
			}
			dst[0] = (uchar)(mask.b);
			dst[1] = (uchar)(mask.g);
			dst[2] = (uchar)(mask.r);
		}
	}

	CMainFrame* pFrame = (CMainFrame*) AfxGetApp()->GetMainWnd();
	CGraphCutView* pView = (CGraphCutView*)(pFrame->GetPane(0,1));
	if (pView->m_storage !=NULL)
	{
		cvReleaseMemStorage( &(pView->m_storage) );
		pView->m_storage = NULL;
	}
	pView->m_storage = cvCreateMemStorage(0);
	cvFindContours( src, pView->m_storage, &(pView->m_contour_d), sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_NONE );
	pView->Invalidate();
	pView->SetTimer(1,100,NULL);

	CvMemStorage* storage_tmp = cvCreateMemStorage(0);
	CvSeq* contour_tmp = NULL;
	cvFindContours( src, storage_tmp, &(contour_tmp), sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_NONE );
	IplImage* segment = cvCreateImage( cvGetSize(myImageObject.GetImage()), 8, 3 );
	segment =  cvCloneImage(myImageObject.GetImage());
	for( ; contour_tmp != 0; contour_tmp = contour_tmp->h_next )
	{
		//cvCvtColor(src,segment,CV_GRAY2BGR); 
		cvDrawContours(segment,contour_tmp,cvScalar(255,0,0),cvScalar(255,0,0),0,2,8); 
	}
	//cvSaveImage("C:\\Users\\xu\\Desktop\\bl result\\20.bmp",segment);
	cvReleaseImage(&src);
	cvReleaseImage(&segment);
}

//xw新建一个窗口显示纹理分割结果
void CImageProcessing::windowDisplay(const Image < SegmentationValue >  &hardSegmentation)
{
	IplImage* src = cvCreateImage( cvGetSize(myImageObject.GetImage()), 8, 1 );
	cvZero(src);
	int i,j;
	LONG width,height;
	width = myImageObject.Width();
	height = myImageObject.Height();
	for ( i=0; i<height; i++)
	{
		for ( j=0; j<width; j++)
		{
			uchar* dst = &CV_IMAGE_ELEM( src, uchar, i, j );
			dst[0] = (uchar)((hardSegmentation(j,i)==SegmentationForeground)?255:0);

			dst = &CV_IMAGE_ELEM( myImageObject.GetImage(), uchar, i, j*3 );
			Color_RGB mask;
			if (hardSegmentation(j,i)==SegmentationForeground)
			{
				mask = (*imageData_RGB)(j,i);
			} 
			else
			{
				mask = ((*imageData_RGB)(j,i))*0.5;
			}
			dst[0] = (uchar)(mask.b);
			dst[1] = (uchar)(mask.g);
			dst[2] = (uchar)(mask.r);
		}
	}
	CvMemStorage* tensor_storage = cvCreateMemStorage(0);
	CvSeq* tensor_contour = NULL;
	cvNamedWindow( "tensor segmentation", 1 );
	cvFindContours( src, tensor_storage, &(tensor_contour), sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_NONE );
	IplImage* segment = cvCreateImage( cvGetSize(myImageObject.GetImage()), 8, 3 );
	segment =  cvCloneImage(myImageObject.GetImage());
	for( ; tensor_contour != 0; tensor_contour = tensor_contour->h_next )
	{
		//cvCvtColor(src,segment,CV_GRAY2BGR); 
		cvDrawContours(segment,tensor_contour,cvScalar(255,0,0),cvScalar(0,0,255),0,2,8); 
	}
	cvShowImage("tensor segmentation", segment);
	/*
	int num_contour = tensor_contour->total; 
	for (int i=0; i<num_contour; i++)
	{
		CvPoint *p_point = (CvPoint*)cvGetSeqElem(tensor_contour, i);
		CvScalar s = cvGet2D(segment , p_point->x,p_point->y);
		s.val[0]=255;
		s.val[1]=0;
		s.val[2]=0;
		cvSet2D(segment,p_point->x,p_point->y,s);
	}
	cvNamedWindow( "tensor segmentation", 1 );
	cvShowImage("tensor segmentation", segment);
	*/
	tensor_contour = NULL;
	cvReleaseMemStorage( &tensor_storage ); 
	cvReleaseImage(&src);
	cvReleaseImage(&segment);
}

void CImageProcessing::updateMattingMask()
{
	CMainFrame* pFrame = (CMainFrame*) AfxGetApp()->GetMainWnd();
	CGraphCutDoc* pDoc = ((CGraphCutView*)pFrame->GetPane(0,1))->GetDocument();
	ASSERT_VALID(pDoc);
	for (unsigned int y = 0; y < imageHeight; ++y)
	{
		for (unsigned int x = 0; x < imageWidth; ++x)
		{
			switch(pDoc->m_eCurCmd)
			{
			case GrabCutCmd:
				(*m_mattingMask)(x, y) = ((*(grabCut->getSegmentationValue()))(x, y)==SegmentationForeground)?255:0;
				break;
			default:
				break;
			}	
		}
	}	
}

IplImage * CImageProcessing::GetCVImage()
{
	return myImageObject.GetImage();
}
