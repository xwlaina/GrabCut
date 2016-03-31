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
	double **dmatrix(int nrl, int nrh, int ncl, int nch);
	void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
	void four1(double *data, int nn, int isign);
	void four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign);
	double *dvector(int nl, int nh);
	void free_dvector(double *v, int nl, int nh);
	void RGB2Lab(const Image < Color_RGB > * _rgb, Image < Color_Lab > * _lab);
	//利用FFT实现快速卷积
	void FFT(std::complex<double> * TD, std::complex<double> * FD, int r);
	void IFFT(std::complex<double> * FD, std::complex<double> * TD, int r);


	void GetMagnitude(CMatrix *Output_mag, const CMatrix *Input_real, const CMatrix *Input_imag);
	void IFFT2(CMatrix *Output_real, CMatrix *Output_imag, const CMatrix *Input_real, const CMatrix *Input_imag);
	void FFT2(CMatrix *Output_real, CMatrix *Output_imag, const CMatrix *Input_real, const CMatrix *Input_imag);

	void Convolution(double * TD1, double * TD2, double * TDout, int M, int N);
	void Convolution(std::complex<double> * TD1, std::complex<double> * TD2, std::complex<double> * TDout, int M, int N);
};

