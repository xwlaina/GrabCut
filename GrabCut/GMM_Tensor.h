#ifndef GMM_TENSOR_H
#define GMM_TENSOR_H

#include "Image.h"
#include <vector>


class  Gaussian_Tensor
{
public:
	Gaussian_Tensor(unsigned int dim);
	~Gaussian_Tensor();

	unsigned int Dim() const
	{
		return m_dim;
	}
	void copy(const Gaussian_Tensor &gauss);

	double *mu;// mean of the gaussian: 0 Ix2, 1 Iy2, 2 IxIy; 3 Gx2, 4 Gy2, 5 GxGy;......
	double *invSqrtmMu;//均值平方根的逆
	double **covariance;   // covariance matrix of the gaussian
	double determinant;        // determinant of the covariance matrix
	double **inverse;      // inverse of the covariance matrix
	double pi;				 // weighting of this gaussian in the GMM.

	double variance; // variance like 1D situation, used for KL with                                  

	double factor; // factor = 1.0 / (pow(2 * PI, m_dim / 2.0) * sqrt(determinant)) or factor = 1.0 / (sqrt(2 * PI * variance))

private:
	unsigned int m_dim;
};


class GMM_Tensor
{
public:
	Gaussian_Tensor** getGaussians()
	{
		return m_gaussians;
	}
	GMM_Tensor( unsigned int K, unsigned int dim);  	// Initialize GMM with number of gaussians desired.
	~GMM_Tensor();


	unsigned int K()const
	{
		return m_K;
	}

	unsigned int Dim()const
	{
		return m_dim;
	}

	unsigned int ScaleCnt()
	{
		return m_scaleCnt;
	}

	// Returns the probability density of tensor vector c in this GMM,即得到某一个张量的概率密度函数
	double p(const double *sample);

	// Returns the probability density of tensor vector c in just Gaussian i
	double p(unsigned int i, const double *sample);

	//返回某一个高斯所对应的均值
	const double *getUnitedMean(const unsigned int& gaussIndex);

private:
	double *m_unitedMean; //Used for return temp united mean tensor of one of the K gaussians
	unsigned int m_K; // number of gaussians
	unsigned int m_scaleCnt;
	unsigned int m_dim; //m_scaleCnt*SiNGLE_TENSOR_DIM
	Gaussian_Tensor** m_gaussians; // an array of K CovGaussians or K*scaleCnt IndGaussians

	Gaussian_Tensor* allocateGaussian();
	void deleteGaussian(Gaussian_Tensor* gauss);
	//计算对角元素的最大值
	double maxDiag(const unsigned int dim, double** cov);

	// 计算协方差矩阵的逆
	double computeInv(const unsigned int dim, double** cov, double** inv);
	//拷贝高斯张量
	void copyGaussian(Gaussian_Tensor* dst, const Gaussian_Tensor* src);
	//设定高斯均值
	void setGaussianMean(Gaussian_Tensor* gauss, const double *mu, BOOL isSetInvSqrtmMu); //mu is mean of the gaussian: 0 Ix2, 1 Iy2, 2 IxIy; 3 Gx2, 4 Gy2, 5 GxGy;......

};

void KMeans_Tensor(const std::vector<CVector *>& samples, int cluster_count, unsigned int *labels, CvTermCriteria termcrit, CVector **initCenters = NULL);

void buildGMMs_kmeans(GMM_Tensor &backgroundGMM, GMM_Tensor &foregroundGMM, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation,  BOOL isOnlyComputeMean);

// Iteratively learn tensor GMMs using GrabCut updating algorithm
void learnGMMs(GMM_Tensor &backgroundGMM, GMM_Tensor &foregroundGMM, Image < unsigned int >  &components, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation);


void Vec(const double* mu, const double* invSqrtmMu, const unsigned int& scaleCnt, const double* sample, double* dif); // sample - mu

double distance2_KL(const CVector *tensors_1, const CVector *tensors_2);

// Helper class that fits a single Gaussian to color samples
class GaussianFitter_Tensor
{
public:
	GaussianFitter_Tensor(unsigned int scaleCnt);
	~GaussianFitter_Tensor();
	// Add a sample
	void add(double *sample);

	void clearSamples();

	// Build the gaussian out of all the added samples
	void finalize( Gaussian_Tensor *g, const unsigned int& totalCount, BOOL isOnlyComputeMean = FALSE) const;

private:
	unsigned int m_scaleCnt;
	std::vector<double *> m_samples; 
	void computeKLMean( Gaussian_Tensor *g, BOOL isOnlyComputeMean) const;
};



inline double tensorDistance2(const CVector *tensor_1, const CVector *tensor_2)
{
	return distance2_KL(tensor_1, tensor_2);
}

#endif //GMM_TENSOR_H
