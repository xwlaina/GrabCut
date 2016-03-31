// Filter.h: interface for the CFilter class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_FILTER_H)
#define _FILTER_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "highgui.h"
#include "cv.h"
#include "Matrix.h"

//À©É¢½á¹¹
struct IntermediateData_Diffusivity
{
	CMatrix *diffusivity;
	CMatrix *dx;
	CMatrix *dy;
	CMatrix *smoothed_channel;
	IplImage *cv_channel;
	CvMat *cv_dx;
	CvMat *cv_dy;
	CvMat *cv_grad;
};

struct IntermediateData_AOS
{
	CMatrix *diffusivity;
	CMatrix *a_row;
	CMatrix *b_row;
	CMatrix *a_column;
	CMatrix *b_column;
	CMatrix *y_row;
	CMatrix *y_column;
};

class CFilter  
{
public:
	CFilter();
	virtual ~CFilter();
	void Diff(CMatrix **channels, unsigned int nochannels);
	void gaussDiff(CMatrix **channels, unsigned int nochannels); 
	void nlDiff(CMatrix **channels, unsigned int nochannels, int option = 0); //Nonlinear Diffusion 
private:
	void TVflow(CvMat *cv_grad);
	void AOS(CvMat *cv_grad);
	void computeDiffusivity(IntermediateData_Diffusivity &data, double sigma, int option);
	void AOS_scheme(IntermediateData_AOS &dada, double stepsize);

	CMatrix **m_channels;
	unsigned int m_nochannels;
	unsigned int m_w, m_h;
};

#endif 
