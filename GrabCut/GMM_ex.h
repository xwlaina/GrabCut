#ifndef GMM_EX_H
#define GMM_EX_H

#include "Image.h"
#include "Matrix.h"
#include <vector>

class  Gaussian_ex
{
public:
	Gaussian_ex(unsigned int dim);
	~Gaussian_ex();

	double *mu;// mean of the gaussian
	double **covariance;   // covariance matrix of the gaussian
	double determinant;        // determinant of the covariance matrix
	double **inverse;      // inverse of the covariance matrix
	double pi;				 // weighting of this gaussian in the GMM.

	// These are only needed during Orchard and Bouman clustering.
	double *eigenvalues;     // eigenvalues of covariance matrix
	double **eigenvectors; // eigenvectors of   "          "

private:
	unsigned int m_dim;
};

class GMM_ex
{
public:
	Gaussian_ex * getGaussians();

	// Initialize GMM with number of gaussians desired.
	GMM_ex(unsigned int K, unsigned int dim);
	~GMM_ex();

	unsigned int K()const
	{
		return m_K;
	}

	unsigned int Dim()const
	{
		return m_dim;
	}

	// Returns the probability density of color c in this GMM
	double p(const double *c) const;

	// Returns the probability density of color c in just Gaussian k
	double p(unsigned int i, const double *c) const;

	CString getParameters(unsigned int component);

private:
	unsigned int m_K; // number of gaussians
	unsigned int m_dim;	
	Gaussian_ex *m_gaussians; // an array of K gaussians

	//	friend void buildGMMs_kmeans(GMM_gabor &backgroundGMM, GMM_gabor &foregroundGMM, const Image < CMatrix* >  &image, const Image < SegmentationValue >  &hardSegmentation);
	//	friend void buildGMMs(GMM_gabor &backgroundGMM, GMM_gabor &foregroundGMM, Image < unsigned int >  &components, const Image < CMatrix* >  &image, const Image < SegmentationValue >  &hardSegmentation);
	//	friend void learnGMMs(GMM_gabor &backgroundGMM, GMM_gabor &foregroundGMM, Image < unsigned int >  &components, const Image < CMatrix* >  &image, const Image < SegmentationValue >  &hardSegmentation);
};

void buildGMMs_kmeans(GMM_ex &backgroundGMM, GMM_ex &foregroundGMM, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation);

// Build the initial GMMs using the Orchard and Bouman color clustering algorithm
void buildGMMs(GMM_ex &backgroundGMM, GMM_ex &foregroundGMM, Image < unsigned int >  &components, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation);

// Iteratively learn GMMs using GrabCut updating algorithm
void learnGMMs(GMM_ex &backgroundGMM, GMM_ex &foregroundGMM, Image < unsigned int >  &components, const Image < CVector* >  &image, const Image < SegmentationValue >  &hardSegmentation);

// Helper class that fits a single Gaussian to color samples
class GaussianFitter_ex
{
public:
	GaussianFitter_ex(unsigned int dim);
	~GaussianFitter_ex();
	// Add a sample
	void add(double *c);

	// Build the gaussian out of all the added samples
	void finalize(Gaussian_ex &g, unsigned int totalCount, bool computeEigens = false)const;

private:
	unsigned int m_dim;
	double *s; // sum
	double **p; // matrix of products, some values are duplicated.

	unsigned int count; // count of samples added to the gaussian
};

#endif //GMM_EX_H
