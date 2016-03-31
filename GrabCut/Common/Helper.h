#pragma once

#include "GMM.h"
#include <assert.h>
#include "Matrix.h"
#include <float.h>




class CHelper
{
public:
	CHelper(void);
	~CHelper(void);
	void RGB2Lab(const Image < Color_RGB > * _rgb, Image < Color_Lab > * _lab);
	//利用FFT实现快速卷积
	void FFT(std::complex<double> * TD, std::complex<double> * FD, int r);
	void IFFT(std::complex<double> * FD, std::complex<double> * TD, int r);

	void Convolution(double * TD1, double * TD2, double * TDout, int M, int N);
	void Convolution(std::complex<double> * TD1, std::complex<double> * TD2, std::complex<double> * TDout, int M, int N);
};

