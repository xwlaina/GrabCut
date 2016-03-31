#include "stdafx.h"
#include "AGMM.h"


AGMM::AGMM(unsigned int K): m_K(K)
{
	m_Agaussians = new AGaussian[m_K];
}

//-------------------------------------------------------------------------
AGMM::~AGMM()
{
	if (NULL!=m_Agaussians)
	{
		delete [] m_Agaussians;
	}
}

//-------------------------------------------------------------------------
//计算非对称混合高斯模型的概率密度函数
double AGMM::p(const Color_Lab& c)
{
	double result = 0.0;

	if (NULL != m_Agaussians)
	{
		for (unsigned int i = 0; i < m_K; i++)
		{
			result += m_Agaussians[i].pi *p(i, c);
		}
	}

	return result;
}

//计算单个非对称高斯模型的概率密度函数
double AGMM::p(unsigned int i, const Color_Lab& c)
{
	double result = 0;

	if (m_Agaussians[i].pi > 0)
	{

		double c0 = c.L - m_Agaussians[i].mu.L;
		double c1 = c.a - m_Agaussians[i].mu.a;
		double c2 = c.b - m_Agaussians[i].mu.b;
		double d0,d1,d2;
		d0 = 0.0;
		d1 = 0.0;
		d2 = 0.0;
		if( c0>0 )
		{
			d0 = (double)((2.0 / (pow(2 * PI, 0.5) * (sqrt(m_Agaussians[i].leftcovariance[0]) + sqrt(m_Agaussians[i].rightcovariance[0])))) * exp( - 0.5 * pow(c0 , 2.0 ) / (m_Agaussians[i].leftcovariance[0])));	
		}else{
			d0 = (double)((2.0 / (pow(2 * PI, 0.5) * (sqrt(m_Agaussians[i].leftcovariance[0]) + sqrt(m_Agaussians[i].rightcovariance[0])))) * exp( - 0.5 * pow(c0 , 2.0 ) / (m_Agaussians[i].rightcovariance[0])));	
		}
		if( c1>0 )
		{
			d1 = (double)((2.0 / (pow(2 * PI, 0.5) * (sqrt(m_Agaussians[i].leftcovariance[1]) + sqrt(m_Agaussians[i].rightcovariance[1])))) * exp( - 0.5 * pow(c1 , 2.0 ) / (m_Agaussians[i].leftcovariance[1])));	
		}else{
			d1 = (double)((2.0 / (pow(2 * PI, 0.5) * (sqrt(m_Agaussians[i].leftcovariance[1]) + sqrt(m_Agaussians[i].rightcovariance[1])))) * exp( - 0.5 * pow(c1 , 2.0 ) / (m_Agaussians[i].rightcovariance[1])));	
		}
		if( c2>0 )
		{
			d2 = (double)((2.0 / (pow(2 * PI, 0.5) * (sqrt(m_Agaussians[i].leftcovariance[2]) + sqrt(m_Agaussians[i].rightcovariance[2])))) * exp( - 0.5 * pow(c2 , 2.0 ) / (m_Agaussians[i].leftcovariance[2])));	
		}else{
			d2 = (double)((2.0 / (pow(2 * PI, 0.5) * (sqrt(m_Agaussians[i].leftcovariance[2]) + sqrt(m_Agaussians[i].rightcovariance[2])))) * exp( - 0.5 * pow(c2 , 2.0 ) / (m_Agaussians[i].rightcovariance[2])));	
		}

		result = (double)(d0 * d1 * d2);

	}
	return result;
}

//-------------------------------------------------------------------------
//利用颜色特征图像和分割的中间结果图进行更新高斯混合的参数（利用OpenCV的kMeans函数）
void buildAGMMs_kmeans(AGMM &backgroundAGMM, AGMM &foregroundAGMM, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	unsigned int width, height;
	width  = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundAGMM.K();
	b_K = backgroundAGMM.K();

	AGaussianFitter *backFitters = new AGaussianFitter[b_K];
	AGaussianFitter *foreFitters = new AGaussianFitter[f_K];

	unsigned int foreCount = 0, backCount = 0;   //前景与背景的像素数
	unsigned int f_index = 0, b_index = 0;
	unsigned int x,y,i;

	//获得前景和背景的像素数
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				foreCount++;
			}
			else  if (hardSegmentation(x, y) == SegmentationBackground)
			{
				backCount++;
			}
		}
	}

	CvMat* f_Colors = cvCreateMat( foreCount, 1, CV_32FC3 );
	CvPoint3D32f* f_pt = (CvPoint3D32f*)f_Colors->data.fl;

	CvMat* b_Colors = cvCreateMat( backCount, 1, CV_32FC3 );
	CvPoint3D32f* b_pt = (CvPoint3D32f*)b_Colors->data.fl;

	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			const Color_Lab& c = image(x, y);
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				f_pt->x = (float)(c.L);
				f_pt->y = (float)(c.a);
				f_pt->z = (float)(c.b);
				f_pt++;
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				b_pt->x = (float)(c.L);
				b_pt->y = (float)(c.a);
				b_pt->z = (float)(c.b);
				b_pt++;
			}
		}
	}

	CvMat* f_clusters = cvCreateMat( foreCount, 1, CV_32SC1 );
	cvKMeans2( f_Colors, f_K, f_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, TERMI_CONST ));

	CvMat* b_clusters = cvCreateMat( backCount, 1, CV_32SC1 );
	cvKMeans2( b_Colors, b_K, b_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, TERMI_CONST ));

	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				foreFitters[f_clusters->data.i[f_index++]].add(image(x, y));
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				backFitters[b_clusters->data.i[b_index++]].add(image(x, y));
			}
		}
	}
	for (i = 0; i < b_K; i++)
	{
		backFitters[i].finalize((backgroundAGMM.getAGaussians())[i], backCount);
	}

	for (i = 0; i < f_K; i++)
	{
		foreFitters[i].finalize((foregroundAGMM.getAGaussians())[i], foreCount);
	}

	delete [] backFitters;
	delete [] foreFitters;
	cvReleaseMat( &f_clusters );
	cvReleaseMat( &f_Colors );
	cvReleaseMat( &b_clusters );
	cvReleaseMat( &b_Colors );
}

//-------------------------------------------------------------------------
//components返回值
void learnAGMMs(AGMM &backgroundAGMM, AGMM &foregroundAGMM, Image < unsigned int >  &components, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	unsigned int x, y;
	unsigned int i;

	unsigned int width, height;
	width  = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundAGMM.K();
	b_K = backgroundAGMM.K();

	// Step 4: Assign each pixel to the component which maximizes its probability，每一个像素属于前景与背景的那个高斯模型
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{
			const Color_Lab& c = image(x, y);

			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				int k = 0;
				double max = 0;

				for (i = 0; i < f_K; i++)
				{
					double p = foregroundAGMM.p(i, c);
					if (p > max)
					{
						k = i;
						max = p;
					}
				}

				components(x, y) = k;
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				int k = 0;
				double max = 0;

				for (i = 0; i < b_K; i++)
				{
					double p = backgroundAGMM.p(i, c);
					if (p > max)
					{
						k = i;
						max = p;
					}
				}

				components(x, y) = k;
			}
		}
	}
	// Step 5: Relearn GMMs from new component assignments
	// Set up Gaussian Fitters
	AGaussianFitter *backFitters = new AGaussianFitter[b_K];
	AGaussianFitter *foreFitters = new AGaussianFitter[f_K];
	unsigned int foreCount = 0, backCount = 0;
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{
			const Color_Lab& c = image(x, y);

			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				foreFitters[components(x, y)].add(c);
				foreCount++;
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				backFitters[components(x, y)].add(c);
				backCount++;
			}
		}
	}
	for (i = 0; i < b_K; i++)
	{
		backFitters[i].finalize((backgroundAGMM.getAGaussians())[i], backCount);
	}
	for (i = 0; i < f_K; i++)
	{
		foreFitters[i].finalize((foregroundAGMM.getAGaussians())[i], foreCount);
	}
	delete [] backFitters;
	delete [] foreFitters;
}

// GaussianFitter functions
AGaussianFitter::AGaussianFitter()
{
	clear();
	allcolor.reserve(100);
}

void AGaussianFitter::clear()
{
	s = Color_Lab();
	allcolor.clear();
	count = 0;
}

// Add a color sample
void AGaussianFitter::add(const Color_Lab& c)
{
	s.L += c.L;
	s.a += c.a;
	s.b += c.b;
	allcolor.push_back(c);	
	count++;
}

// Build the gaussian out of all the added colors
void AGaussianFitter::finalize(AGaussian &ag, unsigned int totalCount)const
{
	// Running into a singular covariance matrix is problematic. So we'll add a small epsilon
	// value to the diagonal elements to ensure a positive definite covariance matrix.
	if (count == 0)
	{
		ag.pi = 0;
	}
	else
	{
		// Compute mean of gaussian
		ag.mu.L = s.L / count;
		ag.mu.a = s.a / count;
		ag.mu.b = s.b / count;
		double leftsum0 = 0;
		double leftsum1 = 0;
		double leftsum2 = 0;
		double rightsum0 = 0;
		double rightsum1 = 0;
		double rightsum2 = 0;
		int leftcount0 = 0;
		int leftcount1 = 0;
		int leftcount2 = 0;
		int rightcount0 = 0;
		int rightcount1 = 0;
		int rightcount2 = 0;
		for(int i=0;i<allcolor.size();i++)
		{
			if(allcolor[i].L >= ag.mu.L)
			{
				rightsum0 += pow((allcolor[i].L - ag.mu.L) , 2.0 );
				rightcount0++;
			}else
			{
				leftsum0 += pow((allcolor[i].L - ag.mu.L) , 2.0 );
				leftcount0++;
			}
			if(allcolor[i].a >= ag.mu.a)
			{
				rightsum1 += pow((allcolor[i].a - ag.mu.a) , 2.0 );
				rightcount1++;
			}else
			{
				leftsum1 += pow((allcolor[i].a - ag.mu.a) , 2.0 );
				leftcount1++;
			}
			if(allcolor[i].b >= ag.mu.b)
			{
				rightsum2 += pow((allcolor[i].b - ag.mu.b) , 2.0 );
				rightcount2++;
			}else
			{
				leftsum2 += pow((allcolor[i].b - ag.mu.b) , 2.0 );
				leftcount2++;
			}
		}
		ag.leftcovariance[0] = (double) leftsum0 / (leftcount0 - 1) ;
		ag.leftcovariance[1] = (double) leftsum1 / (leftcount1 - 1) ;
		ag.leftcovariance[2] = (double) leftsum2 / (leftcount2 - 1) ;
		ag.rightcovariance[0] = (double) rightsum0 / (rightcount0 - 1) ;
		ag.rightcovariance[1] = (double) rightsum1 / (rightcount1 - 1) ;
		ag.rightcovariance[2] = (double) rightsum2 / (rightcount2 - 1) ;
		ag.pi = (double)count / totalCount;
	}
}

AGaussian * AGMM::getAGaussians()
{
	return m_Agaussians;
}