#include "stdafx.h"
#include "GMM_Tensor.h"
#include <set>

Gaussian_Tensor::Gaussian_Tensor(unsigned int dim): m_dim(dim)
{
	unsigned int i;
	mu = new double[m_dim];
	memset(mu, 0, sizeof(double) * m_dim);
	invSqrtmMu = new double[m_dim];
	memset(invSqrtmMu, 0, sizeof(double) * m_dim);
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
	variance = 0.0;
	factor = 0.0;
}

Gaussian_Tensor::~Gaussian_Tensor()
{
	if (mu != NULL)
	{
		delete [] mu;
		mu = NULL;
	}
	if (invSqrtmMu != NULL)
	{
		delete [] invSqrtmMu;
		invSqrtmMu = NULL;
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
}

void Gaussian_Tensor::copy(const Gaussian_Tensor &gauss)
{
	ASSERT(m_dim == (gauss.Dim()));
	memcpy(mu, gauss.mu, sizeof(double) * m_dim);
	memcpy(invSqrtmMu, gauss.invSqrtmMu, sizeof(double) * m_dim);
	memcpy(covariance[0], gauss.covariance[0], sizeof(double) * m_dim * m_dim);
	determinant = gauss.determinant;
	memcpy(inverse[0], gauss.inverse[0], sizeof(double) * m_dim * m_dim);
	pi = gauss.pi;
	variance = gauss.variance;
	factor = gauss.factor;
}

//GMM_Tensor的构造函数
GMM_Tensor::GMM_Tensor( unsigned int K, unsigned int dim):  m_K(K), m_dim(dim)
{
	m_scaleCnt = dim / SiNGLE_TENSOR_DIM;
	unsigned int i;
	m_gaussians = new Gaussian_Tensor*[m_K];
	for (i = 0; i < m_K; ++i)
	{
		m_gaussians[i] = allocateGaussian();
	}
	m_unitedMean = new double[m_dim];
}

//-------------------------------------------------------------------------

GMM_Tensor::~GMM_Tensor()
{
	unsigned int i;
	if (m_gaussians != NULL)
	{
		for (i = 0; i < m_K; ++i)
		{
			deleteGaussian(m_gaussians[i]);
		}
		delete [] m_gaussians;
		m_gaussians = NULL;
	}

	if (m_unitedMean != NULL)
	{
		delete [] m_unitedMean;
		m_unitedMean = NULL;
	}
}

//-------------------------------------------------------------------------
Gaussian_Tensor* GMM_Tensor::allocateGaussian()
{
	Gaussian_Tensor* gauss = NULL;
	gauss = (Gaussian_Tensor*)operator new(sizeof(Gaussian_Tensor) * m_scaleCnt);
	for (unsigned int s = 0; s < m_scaleCnt; ++s)
	{
		new(&(gauss[s])) Gaussian_Tensor(SiNGLE_TENSOR_DIM);
	}
	return gauss;
}

//删除高斯张量
void GMM_Tensor::deleteGaussian(Gaussian_Tensor* gauss)
{
	if (gauss != NULL)
	{		
		for (unsigned int s = 0; s < m_scaleCnt; ++s)
		{
			gauss[m_scaleCnt-1-s].~Gaussian_Tensor();
		}
		operator delete( gauss );	
	}
}


//得到每一个高斯所对应的均值
const double *GMM_Tensor::getUnitedMean(const unsigned int& gaussIndex)
{
	if (m_gaussians != NULL)
	{
		Gaussian_Tensor* gauss = m_gaussians[gaussIndex];
		for (unsigned int k = 0; k < m_scaleCnt; ++k)
		{
			memcpy(m_unitedMean + (k * SiNGLE_TENSOR_DIM), gauss[k].mu, sizeof(double)*SiNGLE_TENSOR_DIM);
		}
		return m_unitedMean;
	}
	else
	{
		return NULL;
	}
}

//得到sample的概率密度函数
double GMM_Tensor::p(const double *sample)
{
	double result = 0.0;
	//m_gaussians：an array of K Gaussians 
	if (m_gaussians != NULL)
	{
		unsigned int i;
		for (i = 0; i < m_K; i++)
		{
			result += (m_gaussians[i][0].pi * p(i, sample));
		}
	}
	return result;
}

double GMM_Tensor::p(unsigned int i, const double *sample)
{//i表示第K个高斯分量
	double result = 1.0;
	unsigned int k;
	unsigned int x, y;
	double part, d;
	double *dif = new double[SiNGLE_TENSOR_DIM];

	for (k = 0; k < m_scaleCnt; k++)
	{
		if ( (m_gaussians[i][k].pi > 0) && (m_gaussians[i][k].determinant > 0) )
		{
			Vec(m_gaussians[i][k].mu, m_gaussians[i][k].invSqrtmMu, 1, sample + (k * SiNGLE_TENSOR_DIM), dif);
			//得到距离的均方导数的映射
			d = 0.0;
			for (x = 0;x < SiNGLE_TENSOR_DIM;x++)
			{
				part = 0.0; 
				for (y = 0;y < SiNGLE_TENSOR_DIM;y++)
				{
					part += dif[y] * m_gaussians[i][k].inverse[y][x];
				}
				d += dif[x] * part;
			}

			result *= (double)(1.0 / (pow(2 * PI, SiNGLE_TENSOR_DIM / 2.0) * sqrt(m_gaussians[i][k].determinant)) * exp( - 0.5 * d));			
		}
		else
		{
			delete [] dif;
			return 0.0;
		}
	}	
	delete [] dif;	
	return pow(result, 1.0 / m_scaleCnt);
}


//计算对角元素的最大值
double GMM_Tensor::maxDiag(const unsigned int dim, double** cov)
{
	double diag = cov[0][0];
	for (unsigned d = 1; d < dim; d++)
	{
		if (cov[d][d] > diag)
		{
			diag = cov[d][d]; 
		}
	}
	return diag;
}

//计算高斯结构张量的逆
double GMM_Tensor::computeInv(const unsigned int dim, double** cov, double** inv)
{
	memcpy(inv[0], cov[0], sizeof(double) * dim * dim);
	for (unsigned int d = 0; d < dim; d++)
	{
		inv[d][d] += REALMIN;
	}
	CvMat inverse = cvMat(dim, dim, CV_TYPE, inv[0]);
	return cvInvert(&inverse, &inverse, CV_LU );
}

//将一个高斯张量赋给另外一个高斯张量
void GMM_Tensor::copyGaussian(Gaussian_Tensor* dst, const Gaussian_Tensor* src)
{
	ASSERT((dst != NULL) && (src != NULL));
	for (unsigned int s = 0; s < m_scaleCnt; ++s)
	{
		(dst[s]).copy(src[s]);
	}
}

void GMM_Tensor::setGaussianMean(Gaussian_Tensor* gauss, const double *mu, BOOL isSetInvSqrtmMu)
{
	unsigned int s;
	double gPro[2][2];  //mat
	double gPro1[2][2]; //eval
	double gPro2[2][2]; //evec

	for (s = 0; s < m_scaleCnt; s++)
	{
		memcpy((gauss[s]).mu, mu + (s * SiNGLE_TENSOR_DIM), SiNGLE_TENSOR_DIM * sizeof(double));
	}

	if (isSetInvSqrtmMu)
	{
		unsigned int pos,pos1,pos2;
		for (s = 0; s < m_scaleCnt; ++s)
		{
			pos = s * SiNGLE_TENSOR_DIM;
			pos1 = pos + 1;
			pos2 = pos + 2;

			gPro[0][0] = mu[pos];
			gPro[1][1] = mu[pos1];
			gPro[0][1] = mu[pos2];
			gPro[1][0] = gPro[0][1];

			sqrtm2D(gPro, gPro1);
			invert2D(gPro1, gPro2);

			(gauss[s]).invSqrtmMu[0] = gPro2[0][0];
			(gauss[s]).invSqrtmMu[1] = gPro2[1][1];
			(gauss[s]).invSqrtmMu[2] = gPro2[0][1];
		}
	}
}


void KMeans_Tensor(const std::vector<CVector *>& samples, int cluster_count, unsigned int *labels, CvTermCriteria termcrit, CVector **initCenters)
{   //样本点数
	int pCnt = samples.size();
	if (pCnt < cluster_count)
	{
		return;
	}
	int k, index = 0;
	std::vector<CVector *>::const_iterator pIter = samples.begin();

	unsigned int dim = (*pIter)->size();
	unsigned int sCnt = dim / SiNGLE_TENSOR_DIM;

	CVector **formalClusters = new CVector *[cluster_count];//前次的聚类中心
	CVector **lastClusters = new CVector *[cluster_count];//本次的聚类中心

	for (k = 0; k < cluster_count; k++)
	{  
		formalClusters[k] = new CVector(dim);
		lastClusters[k] = new CVector(dim);
	}

	//初始的聚类中心(Arbitrarily assign a vector to each of the K clusters)
	if (initCenters != NULL)
	{
		for (k = 0; k < cluster_count; k++)
		{
			formalClusters[k]->set((initCenters[k])->addr());
		}
	} 
	else
	{
		srand(clock());
		int numLocalTries = 2 + (int)log((float)cluster_count);
		double currentPot = 0.0;
		double* closestDistSq = new double[pCnt];

		// Choose one random center and set the closestDistSq values
		index = GetRandom(pCnt);
		formalClusters[0]->set((samples[index])->addr());
		for (k = 0; k < pCnt; k++) 
		{
			closestDistSq[k] = distance2_KL(samples[k], samples[index]);
			currentPot += closestDistSq[k];
		}

		// Choose each center
		for (int centerCount = 1; centerCount < cluster_count; centerCount++) 
		{
			// Repeat several trials
			double bestNewPot = -1;
			int bestNewIndex;
			for (int localTrial = 0; localTrial < numLocalTries; localTrial++) 
			{
				// Choose our center - have to be slightly careful to return a valid answer even accounting
				// for possible rounding errors
				double randVal = (rand() / ((double)(RAND_MAX))) * currentPot;
				for (index = 0; index < (pCnt-1); index++) 
				{
					if (randVal <= closestDistSq[index])
						break;
					else
						randVal -= closestDistSq[index];
				}

				// Compute the new potential
				double newPot = 0.0;
				for (k = 0; k < pCnt; k++)
					newPot += __min( distance2_KL(samples[k], samples[index]), closestDistSq[k] );

				// Store the best result
				if ((bestNewPot < 0) || (newPot < bestNewPot)) 
				{
					bestNewPot = newPot;
					bestNewIndex = index;
				}
			}

			// Add the appropriate center
			formalClusters[centerCount]->set((samples[bestNewIndex])->addr());

			if (centerCount != (cluster_count - 1))
			{
				currentPot = bestNewPot;
				for (k = 0; k < pCnt; k++)
					closestDistSq[k] = __min( distance2_KL(samples[k], samples[bestNewIndex]), closestDistSq[k] );
			}
		}
		delete [] closestDistSq;
	}

	BOOL flag=TRUE;//迭代结束标志
	//迭代计算聚类中心
	while (flag) 
	{
		//根据前次的聚类中心，对各个采样点进行分类
		index = 0;
		for (pIter = samples.begin(); pIter != samples.end(); pIter++)
		{
			int mostSimIndex = 0;
			double mostSimValue = distance2_KL(*pIter, formalClusters[0]);
			for (k = 1; k < cluster_count; k++)
			{
				double curValue = distance2_KL(*pIter, formalClusters[k]);
				if (curValue < mostSimValue)
				{
					mostSimValue = curValue;
					mostSimIndex = k;
				}
			}
			labels[index++] = mostSimIndex;
		}

		//根据分类的象素点，重新计算聚类中心
		GaussianFitter_Tensor **Fitters = new GaussianFitter_Tensor *[cluster_count];	
		for (k = 0; k < cluster_count; k++)
		{
			Fitters[k] = new GaussianFitter_Tensor(sCnt);
		}
		GMM_Tensor *Gmms = new GMM_Tensor( cluster_count, dim);

		index = 0;
		for (pIter = samples.begin(); pIter != samples.end(); pIter++)
		{
			Fitters[labels[index++]]->add((*pIter)->addr());
		}
		for (k = 0; k < cluster_count; k++)
		{
			Fitters[k]->finalize( (Gmms->getGaussians())[k], pCnt, TRUE);
			lastClusters[k]->set(Gmms->getUnitedMean(k));
		}

		for (k = 0; k < cluster_count; k++)
		{
			delete Fitters[k];
		}
		delete [] Fitters;
		delete Gmms;

		//比较原聚类中心和新的聚类中心的误差和迭代次数
		if (termcrit.type & CV_TERMCRIT_ITER)
		{
			termcrit.max_iter--;
			if (termcrit.max_iter == 0)
			{
				break;
			}
		}
		BOOL unchanged=TRUE;
		if (termcrit.type & CV_TERMCRIT_EPS)
		{
			for (k = 0; k < cluster_count; k++)
			{
				if (distance2_KL(formalClusters[k], lastClusters[k]) > termcrit.epsilon)
				{
					unchanged = FALSE;
					break;
				}
			}
		}
		flag=!unchanged;

		//用新的聚类中心代替原聚类中心，再开始下次迭代
		if(flag)
		{
			for (k = 0; k < cluster_count; k++)
			{
				formalClusters[k]->set((lastClusters[k])->addr());
			}
		}
	}

	for (k = 0; k < cluster_count; k++)
	{
		delete formalClusters[k];
		delete lastClusters[k];
	}
	delete [] formalClusters;
	delete [] lastClusters;
}



void buildGMMs_kmeans(GMM_Tensor &backgroundGMM, GMM_Tensor &foregroundGMM, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation, BOOL isOnlyComputeMean)
{
	unsigned int f_sCnt, b_sCnt;
	f_sCnt = foregroundGMM.ScaleCnt();
	b_sCnt = backgroundGMM.ScaleCnt();

	unsigned int width, height;
	width = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundGMM.K();
	b_K = backgroundGMM.K();

	unsigned int x,y,i;
	GaussianFitter_Tensor **backFitters = new GaussianFitter_Tensor *[b_K];
	for (i=0;i<b_K;i++)
	{   //清除所有的样本点
		backFitters[i] = new GaussianFitter_Tensor(b_sCnt);
	}

	GaussianFitter_Tensor **foreFitters = new GaussianFitter_Tensor *[f_K];
	for (i=0;i<f_K;i++)
	{
		foreFitters[i] = new GaussianFitter_Tensor(f_sCnt);
	}

	std::vector<CVector *> f_features, b_features;
	f_features.clear();
	b_features.clear();
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				f_features.push_back(image(x,y));
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				b_features.push_back(image(x,y));
			}
		}
	}
	//labels数组
	unsigned int *f_clusters = new unsigned int[f_features.size()];
	unsigned int *b_clusters = new unsigned int[b_features.size()];

	//得到高斯模型的初始参数
	KMeans_Tensor(f_features, f_K, f_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, TERMI_CONST ), NULL);
	KMeans_Tensor(b_features, b_K, b_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, TERMI_CONST ), NULL);
	//重新扫描图像
	unsigned int f_index = 0, b_index = 0;
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{//将(x,y)划分成K类
				(foreFitters[f_clusters[f_index++]])->add((image(x, y))->addr());
			}
			else if (hardSegmentation(x, y) == SegmentationBackground)
			{
				(backFitters[b_clusters[b_index++]])->add((image(x, y))->addr());
			}
		}
	}
	//计算每一类高斯分布的参数
	for (i = 0; i < b_K; i++)
	{
		(backFitters[i])->finalize((backgroundGMM.getGaussians())[i], b_features.size(), isOnlyComputeMean);
	}

	for (i = 0; i < f_K; i++)
	{
		(foreFitters[i])->finalize((foregroundGMM.getGaussians())[i], f_features.size(), isOnlyComputeMean);
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

	delete [] f_clusters;
	delete [] b_clusters;
}





//-------------------------------------------------------------------------
//高斯参数学习
void learnGMMs(GMM_Tensor &backgroundGMM, GMM_Tensor &foregroundGMM, Image < unsigned int >  &components, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	unsigned int x, y;
	unsigned int i;

	unsigned int f_sCnt, b_sCnt;
	f_sCnt = foregroundGMM.ScaleCnt();
	b_sCnt = backgroundGMM.ScaleCnt();

	unsigned int width, height;
	width = image.width();
	height = image.height();

	unsigned int f_K, b_K;
	f_K = foregroundGMM.K();
	b_K = backgroundGMM.K();

	// Step 4: Assign each pixel to the component which maximizes its probability
	double* c;
	for (y = 0; y < height; ++y)
	{
		for (x = 0; x < width; ++x)
		{
			c = (image(x, y))->addr();

			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				int k = 0;
				double max = 0;

				for (i = 0; i < f_K; i++)
				{
					double p = foregroundGMM.p(i, c);
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
					double p = backgroundGMM.p(i, c);
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
	GaussianFitter_Tensor **backFitters = new GaussianFitter_Tensor *[b_K];
	for (i=0;i<b_K;i++)
	{
		backFitters[i] = new GaussianFitter_Tensor(b_sCnt);
	}
	GaussianFitter_Tensor **foreFitters = new GaussianFitter_Tensor *[f_K];
	for (i=0;i<f_K;i++)
	{
		foreFitters[i] = new GaussianFitter_Tensor(f_sCnt);
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
		(backFitters[i])->finalize((backgroundGMM.getGaussians())[i], backCount);
	}

	for (i = 0; i < f_K; i++)
	{
		(foreFitters[i])->finalize((foregroundGMM.getGaussians())[i], foreCount);
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
GaussianFitter_Tensor::GaussianFitter_Tensor(unsigned int scaleCnt): m_scaleCnt(scaleCnt)
{
	m_samples.clear();
}

GaussianFitter_Tensor::~GaussianFitter_Tensor()
{
	m_samples.swap(std::vector<double *> ());
	m_samples.clear();
}



void GaussianFitter_Tensor::add(double *sample)
{
	m_samples.push_back(sample);
}

void GaussianFitter_Tensor::clearSamples()
{
	m_samples.clear();
}

void GaussianFitter_Tensor::computeKLMean( Gaussian_Tensor *g, BOOL isOnlyComputeMean) const
{
	std::vector<double *>::const_iterator pIter;
	double* sam;
	unsigned int k,pos,pos1,pos2,count;
	count = m_samples.size();


	double gPro[2][2];  //mat
	double gPro1[2][2]; //evec
	double gPro2[2][2];

	double* A = new double[SiNGLE_TENSOR_DIM];
	double* B = new double[SiNGLE_TENSOR_DIM];
	double a, b, c, det;

	for (k = 0; k < m_scaleCnt; ++k)
	{
		pos = k * SiNGLE_TENSOR_DIM;
		pos1 = pos + 1;
		pos2 = pos + 2;
		memset(A,0,sizeof(double)*SiNGLE_TENSOR_DIM);
		memset(B,0,sizeof(double)*SiNGLE_TENSOR_DIM);
		for (pIter = m_samples.begin(); pIter != m_samples.end(); pIter++)
		{
			sam = *pIter;

			a = sam[pos];
			b = sam[pos1];
			c = sam[pos2];

			det = a * b - c * c;
			A[0] += a;
			B[0] += (b / det);
			A[1] += b;
			B[1] += (a / det);
			A[2] += c;
			B[2] -= (c / det);
		}

		//gPro = B
		gPro[0][0] = (B[0] / count);
		gPro[1][1] = (B[1] / count);
		gPro[0][1] = (B[2] / count);
		gPro[1][0] = gPro[0][1];

		sqrtm2D(gPro, gPro1); //gPro1 = sqrt(B)
		interMul2D(gPro1[0][0], gPro1[1][1], gPro1[0][1], A[0] / count, A[1] / count, A[2] / count, gPro2); //gPro2 = sqrt(B)*A*sqrt(B)
		sqrtm2D(gPro2, gPro1); //gPro1 =sqrt(sqrt(B)*A*sqrt(B))
		invert2D(gPro, gPro2); //gPro2 = inv(B) 
		sqrtm2D(gPro2, gPro); //gPro = sqrt(inv(B))
		interMul2D(gPro[0][0], gPro[1][1], gPro[0][1], gPro1[0][0], gPro1[1][1], gPro1[0][1], gPro2); //gPro2 = mu = sqrt(inv(B)) * sqrt(sqrt(B)*A*sqrt(B)) * sqrt(inv(B))
		if (isOnlyComputeMean == FALSE)
		{
			sqrtm2D(gPro2, gPro1); //gPro1 = sqrt(mu)
			invert2D(gPro1, gPro); //gPro = inv(sqrt(mu))
		}

		g[k].mu[0] = gPro2[0][0];
		g[k].mu[1] = gPro2[1][1];
		g[k].mu[2] = gPro2[0][1];
		if (isOnlyComputeMean == FALSE)
		{
			g[k].invSqrtmMu[0] = gPro[0][0];
			g[k].invSqrtmMu[1] = gPro[1][1];
			g[k].invSqrtmMu[2] = gPro[0][1];
		}

	}

	delete [] A;
	delete [] B;
}

//多尺度结构张量MSNST的梯度映射
void Vec(const double* mu, const double* invSqrtmMu, const unsigned int& scaleCnt, const double* sample, double* dif)
{
	unsigned int k, pos, pos1, pos2;
	double gPro[2][2];  //mat
	double gEigenvalues[2];
	double gPro2[2][2]; //evec
	for (k = 0; k < scaleCnt; k++)
	{
		pos = k * SiNGLE_TENSOR_DIM;
		pos1 = pos + 1;
		pos2 = pos + 2;

		interMul2D(invSqrtmMu[pos], invSqrtmMu[pos1], invSqrtmMu[pos2], sample[pos], sample[pos1], sample[pos2], gPro);
		JacobiEigenv2D(gPro, gEigenvalues, gPro2);
		diagMul2D(gPro2[0][0], gPro2[0][1], gPro2[1][0], gPro2[1][1], log(gEigenvalues[0]), log(gEigenvalues[1]), gPro);
		dif[pos] = gPro[0][0];
		dif[pos1] = SQRT_OF_2 * gPro[0][1];
		dif[pos2] = gPro[1][1];		
	}

}

// Build the gaussian out of all the added samples
void GaussianFitter_Tensor::finalize(Gaussian_Tensor *g, const unsigned int& totalCount, BOOL isOnlyComputeMean) const
{		
	const double Epsilon = (double)0.0001; // Running into a singular covariance matrix is problematic. So we'll add a small epsilon value to the diagonal elements to ensure a positive definite covariance matrix.

	unsigned int count = m_samples.size();

	unsigned int k, x, y;
	if (count == 0)
	{
		for (k = 0; k < m_scaleCnt; k++)
		{
			g[k].pi = 0;
		}
	}
	else
	{
		computeKLMean(g, isOnlyComputeMean);
		if (isOnlyComputeMean == FALSE)
		{
			//计算方差，见(3_14)
			std::vector<double *>::const_iterator pIter;

			unsigned int pos;
			double* dif = new double[SiNGLE_TENSOR_DIM];

			for (k = 0; k < m_scaleCnt; k++)
			{
				pos = k * SiNGLE_TENSOR_DIM;
				// Compute covariance matrix
				memset(g[k].covariance[0], 0, sizeof(double) * SiNGLE_TENSOR_DIM * SiNGLE_TENSOR_DIM);
				for (pIter = m_samples.begin(); pIter != m_samples.end(); pIter++)
				{   //梯度映射结果保存在dif中
					Vec(g[k].mu, g[k].invSqrtmMu, 1, (*pIter) + pos, dif);
					for (y = 0;y < SiNGLE_TENSOR_DIM;y++)
					{
						for (x = 0;x < SiNGLE_TENSOR_DIM;x++)
						{
							g[k].covariance[y][x] += (dif[y] * dif[x]);
						}
					}
				}

				for (y = 0;y < SiNGLE_TENSOR_DIM;y++)
				{
					for (x = 0;x < SiNGLE_TENSOR_DIM;x++)
					{
						g[k].covariance[y][x] *=  (1.0 / (count - 1));
						if (x == y)
						{
							g[k].covariance[y][x] += Epsilon;
						}
					}
				}

				CvMat covariance = cvMat(SiNGLE_TENSOR_DIM, SiNGLE_TENSOR_DIM, CV_TYPE, g[k].covariance[0]);
				CvMat inverse = cvMat(SiNGLE_TENSOR_DIM, SiNGLE_TENSOR_DIM, CV_TYPE, g[k].inverse[0]);
				g[k].determinant = cvInvert(&covariance, &inverse, CV_LU ); 
				// The weight of the gaussian is the fraction of the number of pixels in this Gaussian to the number of 
				// pixels in all the gaussians of this GMM.
				g[k].pi = (1.0 * count) / totalCount;
			}		

			delete [] dif;
		}
	} 		
}


//计算两个张量之间的距离
double distance2_KL(const CVector *tensors_1, const CVector *tensors_2)
{
	int levels = tensors_1->size() / SiNGLE_TENSOR_DIM;
	double w = 0.0;
	double value1, value2;
	double a_1, b_1, c_1, a_2, b_2, c_2;
	unsigned int pos, pos1, pos2;
	for (int k = 0; k < levels;k++)
	{
		pos = k * SiNGLE_TENSOR_DIM;
		pos1 = pos + 1;
		pos2 = pos + 2;
		a_1 = tensors_1->get(pos);
		b_1 = tensors_1->get(pos1);
		c_1 = tensors_1->get(pos2);
		a_2 = tensors_2->get(pos);
		b_2 = tensors_2->get(pos1);
		c_2 = tensors_2->get(pos2);
		ASSERT((a_1 * b_1 - c_1 * c_1) != 0 && (a_2 * b_2 - c_2 * c_2) != 0);
		value1 = (a_1 * b_2 + b_1 * a_2 - 2 * c_1 * c_2) / (a_1 * b_1 - c_1 * c_1);
		value2 = (a_2 * b_1 + b_2 * a_1 - 2 * c_2 * c_1) / (a_2 * b_2 - c_2 * c_2);
		w += ((value1 + value2 - 4) / 4);
	}
	return w;
}












