#include "stdafx.h"
#include "GMM_ex.h"
#include <vector>
#include <set>


Gaussian_ex::Gaussian_ex(unsigned int dim): m_dim(dim)
{
	unsigned int i;
	mu = new double[m_dim];
	memset(mu, 0, sizeof(double) * m_dim);
	covariance = new double *[m_dim];
	covariance[0] = new double[m_dim*m_dim];
	memset(covariance[0], 0, sizeof(double) * m_dim * m_dim);
	for (i = 1; i < m_dim; i++)
	{
		covariance[i] = covariance[i - 1] + m_dim;
	}
	determinant = 0.0;
	inverse = new double *[m_dim];
	inverse[0] = new double[m_dim*m_dim];
	memset(inverse[0], 0, sizeof(double) * m_dim * m_dim);
	for (i = 1; i < m_dim; i++)
	{
		inverse[i] = inverse[i - 1] + m_dim;
	}
	pi = 0.0;

	eigenvalues = new double[m_dim];
	memset(eigenvalues, 0, sizeof(double) * m_dim);
	eigenvectors = new double *[m_dim];
	eigenvectors[0] = new double[m_dim*m_dim];
	memset(eigenvectors[0], 0, sizeof(double) * m_dim * m_dim);
	for (i = 1; i < m_dim; i++)
	{
		eigenvectors[i] = eigenvectors[i - 1] + m_dim;
	}
}

Gaussian_ex::~Gaussian_ex()
{
	if (mu != NULL)
	{
		delete [] mu;
		mu = NULL;
	}
	if (covariance != NULL)
	{
		delete [] covariance[0];
		covariance[0] = NULL;
		delete [] covariance;
		covariance = NULL;
	}
	if (inverse != NULL)
	{
		delete [] inverse[0];
		inverse[0] = NULL;
		delete [] inverse;
		inverse = NULL;
	}
	if (eigenvalues != NULL)
	{
		delete [] eigenvalues;
		eigenvalues = NULL;
	}
	if (eigenvectors != NULL)
	{
		delete [] eigenvectors[0];
		eigenvectors[0] = NULL;
		delete [] eigenvectors;
		eigenvectors = NULL;
	}
}

GMM_ex::GMM_ex(unsigned int K, unsigned int dim): m_K(K), m_dim(dim)
{
	m_gaussians = (Gaussian_ex *)operator new(sizeof(Gaussian_ex) * m_K);

	for( unsigned int i=0; i<m_K; ++i )
	{
		new(&m_gaussians[i]) Gaussian_ex(m_dim);
	}

}

//-------------------------------------------------------------------------

GMM_ex::~GMM_ex()
{
	if (m_gaussians)
	{
		for( unsigned int i=0; i<m_K; ++i )
		{
			m_gaussians[m_K-1-i].~Gaussian_ex();
		}
		operator delete( m_gaussians );
	}
}

CString GMM_ex::getParameters(unsigned int component)
{
	CString strOut = _T(""), str, strWeight = _T(""), strMean = _T(""), strCovariance = _T("");
	if (m_gaussians)
	{
		unsigned int i,j,k;
		for (i = 0; i < m_K; i++)
		{
			str.Format(_T("W(%d) = %f;\r\n"), i+1, m_gaussians[i].pi);
			strWeight += str;

			str.Format(_T("Mu(:,%d,%d) = ["), i+1, component);
			strMean += str;
			for (j=0;j<m_dim;j++)
			{
				if (j == (m_dim - 1))
				{
					str.Format(_T("%f];\r\n"), m_gaussians[i].mu[j]);
				} 
				else
				{
					str.Format(_T("%f;"), m_gaussians[i].mu[j]);
				}
				strMean += str;
			}

			str.Format(_T("Cov(:,:,%d,%d) = ["), i+1, component);
			strCovariance += str;
			for (j=0;j<m_dim;j++)
			{
				for (k=0;k<m_dim;k++)
				{
					if (k == (m_dim - 1))
					{
						str.Format(_T("%f"), m_gaussians[i].covariance[j][k]);
					} 
					else
					{
						str.Format(_T("%f "), m_gaussians[i].covariance[j][k]);
					}
					strCovariance += str;
				}

				if (j == (m_dim - 1))
				{
					strCovariance += _T("];\r\n");
				} 
				else
				{
					strCovariance += _T(";");
				}
			}

		}	
		strOut = strWeight + strMean + strCovariance;
	}
	return strOut;
}

//-------------------------------------------------------------------------

double GMM_ex::p(const double *c) const
{
	double result = 0;

	if (m_gaussians)
	{
		for (unsigned int i = 0; i < m_K; i++)
		{
			result += m_gaussians[i].pi * p(i, c);
		}
	}

	return result;
}

//-------------------------------------------------------------------------

double GMM_ex::p(unsigned int i, const double *c) const
{
	int x,y;
	double *dif = new double[m_dim];
	double d, part, result = 0.0;

	if (m_gaussians[i].pi > 0)
	{
		if (m_gaussians[i].determinant > 0)
		{			
			for (x = 0;x<m_dim;x++)
			{
				dif[x] = c[x] - m_gaussians[i].mu[x];
			}

			d = 0.0;
			for (x = 0;x < m_dim;x++)
			{
				part = 0.0; 
				for (y = 0;y < m_dim;y++)
				{
					part += dif[y] * m_gaussians[i].inverse[y][x];
				}
				d += dif[x] * part;
			}

			result = (double)(1.0 / (pow(2 * PI, m_dim / 2.0) * sqrt(m_gaussians[i].determinant)) *exp( - 0.5 * d));
		}
	}

	delete [] dif;
	return result;
}

//-------------------------------------------------------------------------

void buildGMMs_kmeans(GMM_ex &backgroundGMM, GMM_ex &foregroundGMM, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	unsigned int f_dim, b_dim;
	f_dim = foregroundGMM.Dim();
	b_dim = backgroundGMM.Dim();

	unsigned int width, height;
	width = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundGMM.K();
	b_K = backgroundGMM.K();

	unsigned int x,y,i;
	GaussianFitter_ex **backFitters = new GaussianFitter_ex *[b_K];
	for (i=0;i<b_K;i++)
	{
		backFitters[i] = new GaussianFitter_ex(b_dim);
	}
	GaussianFitter_ex **foreFitters = new GaussianFitter_ex *[f_K];
	for (i=0;i<f_K;i++)
	{
		foreFitters[i] = new GaussianFitter_ex(f_dim);
	}

	unsigned int foreCount = 0, backCount = 0, f_index = 0, b_index = 0;


	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				foreCount++;
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				backCount++;
			}
		}
	}

	CvMat* f_features = cvCreateMat( foreCount, 1, CV_32FC(f_dim) );
	float* f_pt = (float*)f_features->data.fl;
	CvMat* b_features = cvCreateMat( backCount, 1, CV_32FC(b_dim) );
	float* b_pt = (float*)b_features->data.fl;

	CVector* pVec;
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			pVec = image(x,y);
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				for (i=0;i<f_dim;i++)
				{
					(*f_pt) = (float)(pVec->get(i));
					f_pt++;
				}
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				for (i=0;i<b_dim;i++)
				{
					(*b_pt) = (float)(pVec->get(i));
					b_pt++;
				}
			}
		}
	}
	CvMat* f_clusters = cvCreateMat( foreCount, 1, CV_32SC1 );
	cvKMeans2( f_features, f_K, f_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, TERMI_CONST ));
	CvMat* b_clusters = cvCreateMat( backCount, 1, CV_32SC1 );
	cvKMeans2( b_features, b_K, b_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, TERMI_CONST ));
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				(foreFitters[f_clusters->data.i[f_index++]])->add((image(x, y))->addr());
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				(backFitters[b_clusters->data.i[b_index++]])->add((image(x, y))->addr());
			}
		}
	}
	for (i = 0; i < b_K; i++)
	{
		(backFitters[i])->finalize((backgroundGMM.getGaussians())[i], backCount, false);
	}

	for (i = 0; i < f_K; i++)
	{
		(foreFitters[i])->finalize((foregroundGMM.getGaussians())[i], foreCount, false);
	}

	for (i=0;i<b_K;i++)
	{
		delete backFitters[i];
	}
	delete [] backFitters;
	for (i=0;i<f_K;i++)
	{
		delete foreFitters[i];
	}
	delete [] foreFitters;
	cvReleaseMat( &f_clusters );
	cvReleaseMat( &f_features );
	cvReleaseMat( &b_clusters );
	cvReleaseMat( &b_features );
}
//-------------------------------------------------------------------------

void buildGMMs(GMM_ex &backgroundGMM, GMM_ex &foregroundGMM, Image < unsigned int >  &components, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	// Step 3: Build GMMs using Orchard-Bouman clustering algorithm
	unsigned int f_dim, b_dim;
	f_dim = foregroundGMM.Dim();
	b_dim = backgroundGMM.Dim();

	unsigned int width, height;
	width = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundGMM.K();
	b_K = backgroundGMM.K();

	// Set up Gaussian Fitters
	unsigned int x,y,i,k;
	GaussianFitter_ex **backFitters = new GaussianFitter_ex *[b_K];
	for (i=0;i<b_K;i++)
	{
		backFitters[i] = new GaussianFitter_ex(b_dim);
	}
	GaussianFitter_ex **foreFitters = new GaussianFitter_ex *[f_K];
	for (i=0;i<f_K;i++)
	{
		foreFitters[i] = new GaussianFitter_ex(f_dim);
	}

	unsigned int foreCount = 0, backCount = 0;

	// Initialize the first foreground and background clusters
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{
			components(x, y) = 0;

			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				(foreFitters[0])->add((image(x, y))->addr());
				foreCount++;
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				(backFitters[0])->add((image(x, y))->addr());
				backCount++;
			}
		}
	}

	(backFitters[0])->finalize((backgroundGMM.getGaussians())[0], backCount, true);
	(foreFitters[0])->finalize((foregroundGMM.getGaussians())[0], foreCount, true);

	unsigned int nBack = 0, nFore = 0; // Which cluster will be split
	unsigned int maxK = b_K > f_K ? b_K: f_K;

	// Compute clusters
	for (i = 1; i < maxK; i++)
	{
		// Reset the fitters for the splitting clusters
		if (backFitters[nBack])
		{
			delete backFitters[nBack];
		}
		backFitters[nBack] = new GaussianFitter_ex(b_dim);
		if (foreFitters[nFore])
		{
			delete foreFitters[nFore];
		}
		foreFitters[nFore] = new GaussianFitter_ex(f_dim);

		// For brevity, get references to the splitting Gaussians
		Gaussian_ex &bg = (backgroundGMM.getGaussians())[nBack];
		Gaussian_ex &fg = (foregroundGMM.getGaussians())[nFore];

		// Compute splitting points
		double splitBack = 0.0, splitFore = 0.0;
		for (k = 0; k < f_dim; k++)
		{
			splitBack += bg.eigenvectors[k][0] * bg.mu[k];
			splitFore += fg.eigenvectors[k][0] * fg.mu[k];
		}

		// Split clusters nBack and nFore, place split portion into cluster i
		double* c;
		for (y = 0; y < height; ++y)
		{
			for (x = 0; x < width; ++x)
			{
				c = (image(x, y))->addr();

				// For each pixel
				if (i < f_K && hardSegmentation(x, y) == SegmentationForeground && components(x, y) == nFore)
				{
					double split = 0.0;
					for (k = 0; k < f_dim; k++)
					{
						split += fg.eigenvectors[k][0] * c[k];
					}			
					if (split > splitFore)
					{
						components(x, y) = i;
						(foreFitters[i])->add(c);
					}
					else
					{
						(foreFitters[nFore])->add(c);
					}
				}
				else if (i < b_K && hardSegmentation(x, y) == SegmentationBackground && components(x, y) == nBack)
				{
					double split = 0.0;
					for (k = 0; k < f_dim; k++)
					{
						split += bg.eigenvectors[k][0] * c[k];
					}
					if (split > splitBack)
					{
						components(x, y) = i;
						(backFitters[i])->add(c);
					}
					else
					{
						(backFitters[nBack])->add(c);
					}
				}
			}
		}


		// Compute new split Gaussians
		(backFitters[nBack])->finalize((backgroundGMM.getGaussians())[nBack], backCount, true);
		(foreFitters[nFore])->finalize((foregroundGMM.getGaussians())[nFore], foreCount, true);

		if (i < b_K)
		{
			(backFitters[i])->finalize((backgroundGMM.getGaussians())[i], backCount, true);
		}
		if (i < f_K)
		{
			(foreFitters[i])->finalize((foregroundGMM.getGaussians())[i], foreCount, true);
		}

		// Find clusters with highest eigenvalue
		nBack = 0;
		nFore = 0;

		for (unsigned int j = 0; j <= i; j++)
		{
			if (j < b_K && (backgroundGMM.getGaussians())[j].eigenvalues[0] > (backgroundGMM.getGaussians())[nBack].eigenvalues[0])
			{
				nBack = j;
			}

			if (j < f_K && (foregroundGMM.getGaussians())[j].eigenvalues[0] > (foregroundGMM.getGaussians())[nFore].eigenvalues[0])
			{
				nFore = j;
			}
		}
	}

	for (i=0;i<b_K;i++)
	{
		delete backFitters[i];
	}
	delete [] backFitters;
	for (i=0;i<f_K;i++)
	{
		delete foreFitters[i];
	}
	delete [] foreFitters;
}

//-------------------------------------------------------------------------

void learnGMMs(GMM_ex &backgroundGMM, GMM_ex &foregroundGMM, Image < unsigned int >  &components, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	unsigned int x, y;
	unsigned int i;

	unsigned int f_dim, b_dim;
	f_dim = foregroundGMM.Dim();
	b_dim = backgroundGMM.Dim();

	unsigned int width, height;
	width = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundGMM.K();
	b_K = backgroundGMM.K();

	// Step 4: Assign each pixel to the component which maximizes its probability
	double* c;
	double p, max;
	int k;
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{
			c = (image(x, y))->addr();

			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				k = 0;
				max = 0;

				for (i = 0; i < f_K; i++)
				{
					p = foregroundGMM.p(i, c);
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
				k = 0;
				max = 0;

				for (i = 0; i < b_K; i++)
				{
					p = backgroundGMM.p(i, c);
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
	GaussianFitter_ex **backFitters = new GaussianFitter_ex *[b_K];
	for (i=0;i<b_K;i++)
	{
		backFitters[i] = new GaussianFitter_ex(b_dim);
	}
	GaussianFitter_ex **foreFitters = new GaussianFitter_ex *[f_K];
	for (i=0;i<f_K;i++)
	{
		foreFitters[i] = new GaussianFitter_ex(f_dim);
	}

	unsigned int foreCount = 0, backCount = 0;

	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{
			c = (image(x, y))->addr();

			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				(foreFitters[components(x, y)])->add(c);
				foreCount++;
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				(backFitters[components(x, y)])->add(c);
				backCount++;
			}
		}
	}

	for (i = 0; i < b_K; i++)
	{
		(backFitters[i])->finalize((backgroundGMM.getGaussians())[i], backCount, false);
	}

	for (i = 0; i < f_K; i++)
	{
		(foreFitters[i])->finalize((foregroundGMM.getGaussians())[i], foreCount, false);
	}

	for (i=0;i<b_K;i++)
	{
		delete backFitters[i];
	}
	delete [] backFitters;
	for (i=0;i<f_K;i++)
	{
		delete foreFitters[i];
	}
	delete [] foreFitters;
}


// GaussianFitter functions
GaussianFitter_ex::GaussianFitter_ex(unsigned int dim): m_dim(dim)
{
	s = new double[m_dim];
	memset(s,0,sizeof(double)*m_dim);
	p = new double *[m_dim];
	p[0] = new double[m_dim*m_dim];
	memset(p[0], 0, sizeof(double) * m_dim * m_dim);
	for (unsigned int i = 1; i < m_dim; i++)
	{
		p[i] = p[i - 1] + m_dim;
	}
	count = 0;
}

GaussianFitter_ex::~GaussianFitter_ex()
{
	if (s != NULL)
	{
		delete [] s;
		s = NULL;
	}
	if (p != NULL)
	{
		delete [] p[0];
		p[0] = NULL;
		delete [] p;
		p = NULL;
	}
}

// Add a color sample, c µÄÎ¬ÊýÎªm_dim
void GaussianFitter_ex::add(double *c)
{
	int x,y;
	for (x = 0;x < m_dim;x++)
	{
		s[x] += c[x];
	}
	for (y = 0;y < m_dim;y++)
	{
		for (x = 0;x < m_dim;x++)
		{
			p[y][x] += c[y] * c[x]; 
		}
	}
	count++;
}

// Build the gaussian out of all the added colors
void GaussianFitter_ex::finalize(Gaussian_ex &g, unsigned int totalCount, bool computeEigens)const
{
	// Running into a singular covariance matrix is problematic. So we'll add a small epsilon
	// value to the diagonal elements to ensure a positive definite covariance matrix.
	const double Epsilon = (double)0.0001;

	if (count == 0)
	{
		g.pi = 0;
	}
	else
	{
		int x,y;
		// Compute mean of gaussian
		for (x = 0;x < m_dim;x++)
		{
			g.mu[x] = s[x] / count;
		}

		// Compute covariance matrix
		for (y = 0;y < m_dim;y++)
		{
			for (x = 0;x < m_dim;x++)
			{
				g.covariance[y][x] = (p[y][x] - count * g.mu[y] * g.mu[x]) / (count - 1);
				if (x == y)
				{
					g.covariance[y][x] += Epsilon;
				}
			}
		}

		// Compute determinant of covariance matrix
		// Compute inverse (cofactor matrix divided by determinant)
		CvMat covariance = cvMat(m_dim, m_dim, CV_TYPE, g.covariance[0]);
		CvMat inverse = cvMat(m_dim, m_dim, CV_TYPE, g.inverse[0]);
		g.determinant = cvInvert(&covariance, &inverse, CV_LU ); 

		// The weight of the gaussian is the fraction of the number of pixels in this Gaussian to the number of 
		// pixels in all the gaussians of this GMM.
		g.pi = (double)count / totalCount;

		if (computeEigens)
		{
			// Build OpenCV wrappers around our data.
			CvMat eval = cvMat(m_dim, 1, CV_TYPE, g.eigenvalues);
			CvMat evec = cvMat(m_dim, m_dim, CV_TYPE, g.eigenvectors[0]);

			// Compute eigenvalues and vectors using SVD
			cvSVD(&covariance, &eval, &evec);
		}
	}
}

Gaussian_ex * GMM_ex::getGaussians()
{
	return m_gaussians;
}
