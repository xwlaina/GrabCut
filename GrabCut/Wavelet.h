// Wavelet.h: interface for the Wavelet class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(_WAVELET_H)
#define _WAVELET_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define FilterHLen 4
#define FilterGLen 2
#define FilterDLen 3

#include "Helper.h"
#include "Matrix.h"

struct WaveletDetailImages
{
	CMatrix *Detail_1;
	CMatrix *Detail_2;
};

class Wavelet  
{
public:
	void execute(WaveletDetailImages *results);
	Wavelet(const CMatrix *image, double dMinValue, double dMaxValue, unsigned int levels);
	virtual ~Wavelet();

private:
	void Convolution2(const CMatrix *image_in, CMatrix *image_out, double *columns_filter, int columns_filter_len, double *rows_filter, int rows_filter_len);

	int m_w, m_h;
	unsigned int m_levels;//变换层数
	double m_minValue, m_maxValue;
	CMatrix *m_image;

	CHelper m_helper;
	//常数
	double m_esp;
	unsigned int m_iterations;//迭代次数
	static double m_HFilter[FilterHLen];
	static double m_GFilter[FilterGLen];
	static double m_DFilter[FilterDLen];
};

#endif // !defined(AFX_WAVELET_H__31721814_7C7C_4671_A7AD_B0B6B9274710__INCLUDED_)
