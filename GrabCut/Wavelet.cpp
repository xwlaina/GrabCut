// Wavelet.cpp: implementation of the Wavelet class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Wavelet.h"

//Use Matcom C++ library
//#include "matlib.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

double Wavelet::m_HFilter[FilterHLen] = {0.125,0.375,0.375,0.125};
double Wavelet::m_GFilter[FilterGLen] = {0.5,-0.5};
double Wavelet::m_DFilter[FilterDLen] = {1,0,0};


Wavelet::Wavelet(const CMatrix *image, double dMinValue, double dMaxValue, unsigned int levels)
{
	m_esp = 0.00001;
	m_iterations = 20;
	m_minValue = dMinValue;
	m_maxValue = dMaxValue;
	m_levels = levels;
	ASSERT(m_levels > 0);
	m_w = image->GetNumColumns();
	m_h = image->GetNumRows();
	m_image = new CMatrix(m_h, m_w);
	//归一化
	double dDiff = m_maxValue - m_minValue + m_esp;
	for (int y = 0 ; y < m_h; ++y)
	{
		for (int x = 0 ; x <  m_w; ++x)
		{
			m_image->SetElement(y, x, (image->GetElement(y, x) - m_minValue) / dDiff );
		}
	}
}

Wavelet::~Wavelet()
{
	if (m_image != NULL)
	{
		delete m_image;
		m_image = NULL;
	}
}

void Wavelet::Convolution2(const CMatrix *image_in, CMatrix *image_out, double *columns_filter, int columns_filter_len, double *rows_filter, int rows_filter_len)
{
	int width, height;
	width = image_in->GetNumColumns();
	height = image_in->GetNumRows();
	// 首先检查行列数是否相等
	ASSERT (width == image_out->GetNumColumns() && height == image_out->GetNumRows());
    
	
	CMatrix *image_tmp = new CMatrix(height,width);
	int i;
	double *DataRow = new double[width];
	CMatrix *Row_tmp = new CMatrix(1, width);
	double *DataCol = new double[height];
	CMatrix *Col_tmp = new CMatrix(height, 1);
	
	//按列卷积
	  for (i = 0; i < width; i++)
	  {
	  image_in->GetColVector(i, DataCol);
	  m_helper.Convolution(DataCol, columns_filter, Col_tmp->GetData(), height, columns_filter_len);
	  image_tmp->SetColVector(i, *Col_tmp);
	  }
	  //按行卷积
	  for (i = 0; i < height; i++)
	  {
	  image_tmp->GetRowVector(i, DataRow);
	  m_helper.Convolution(DataRow, rows_filter, Row_tmp->GetData(), width, rows_filter_len);
	  image_out->SetRowVector(i, *Row_tmp);
	  }
	  
		delete image_tmp;
		delete [] DataRow;
		delete Row_tmp;
		delete [] DataCol;
	delete Col_tmp;
	
	/*
	//matlab中矩阵按列存储，VC++中矩阵按行存储，需转换
	mwArray image_inM =  row2mat(height, width, image_in->GetData());

	mwArray colM(1, columns_filter_len, columns_filter);
	mwArray rowM(1, rows_filter_len, rows_filter);
	mwArray image_outM = transpose(conv2(colM, rowM, image_inM, "same"));
	memcpy(image_out->GetData(), mxGetPr(image_outM.GetData()), width * height * sizeof(double));
*/
}

void Wavelet::execute(WaveletDetailImages *results)
{
	//构建低分辨率的图像
	CMatrix **LowResolutionImages = new CMatrix *[m_levels];
	unsigned int n,j;
	for (j = 0; j < m_levels; j++)
	{
		LowResolutionImages[j] = new CMatrix(m_h, m_w);
	}
	Convolution2(m_image, LowResolutionImages[0], m_HFilter, FilterHLen, m_HFilter, FilterHLen);
	Convolution2(m_image, results[0].Detail_1, m_DFilter, FilterDLen, m_GFilter, FilterGLen);
	Convolution2(m_image, results[0].Detail_2, m_GFilter, FilterGLen, m_DFilter, FilterDLen);
	
	for (j = 1; j < m_levels; j++)
	{
		unsigned int scale =(unsigned int) pow((float)2,(int)j);
		unsigned int FilterHLen_j, FilterGLen_j;
		FilterHLen_j = scale * (FilterHLen - 1) + 1;
		FilterGLen_j = scale * (FilterGLen - 1) +1;
		double *HFilter_j = new double[FilterHLen_j];
		memset(HFilter_j, 0, sizeof(double) * FilterHLen_j);
		double *GFilter_j = new double[FilterGLen_j];
		memset(GFilter_j, 0, sizeof(double) * FilterGLen_j);
		for (n = 0; n < FilterHLen; n++)
		{
			HFilter_j[scale * n] = m_HFilter[n];
		}
		for (n = 0; n < FilterGLen; n++)
		{
			GFilter_j[scale * n] = m_GFilter[n];
		}

		Convolution2(LowResolutionImages[j - 1], LowResolutionImages[j], HFilter_j, FilterHLen_j, HFilter_j, FilterHLen_j);
		Convolution2(LowResolutionImages[j - 1], results[j].Detail_1, m_DFilter, FilterDLen, GFilter_j, FilterGLen_j);
		Convolution2(LowResolutionImages[j - 1], results[j].Detail_2, GFilter_j, FilterGLen_j, m_DFilter, FilterDLen);
	
		delete [] HFilter_j;
		delete [] GFilter_j;
	}

	//反归一化
	double dDiff = m_maxValue - m_minValue + m_esp;
	for (j = 0; j < m_levels; j++)
	{
		for (int y = 0 ; y < m_h; ++y)
		{
			for (int x = 0 ; x <  m_w; ++x)
			{
				results[j].Detail_1->SetElement(y, x, dDiff * results[j].Detail_1->GetElement(y, x) + m_minValue);
				results[j].Detail_2->SetElement(y, x, dDiff * results[j].Detail_2->GetElement(y, x) + m_minValue);
			}
		}
	}

	for (j = 0; j < m_levels; j++)
	{
		delete LowResolutionImages[j];
	}
	delete [] LowResolutionImages;
}