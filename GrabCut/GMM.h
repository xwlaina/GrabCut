/*
 * GrabCut implementation source code Copyright(c) 2005-2006 Justin Talbot
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 */

#ifndef GMM_H
#define GMM_H

#include "Color_Lab.h"
#include "Image.h"
#include "Global.h"

struct Gaussian
{
        Color_Lab mu;                // mean of the gaussian
        double covariance[3][3];   // covariance matrix of the gaussian
        double determinant;        // determinant of the covariance matrix
        double inverse[3][3];      // inverse of the covariance matrix //逆矩阵
        double pi;				 // weighting of this gaussian in the GMM.

        // These are only needed during Orchard and Bouman clustering.
        double eigenvalues[3];     // eigenvalues of covariance matrix
        double eigenvectors[3][3]; // eigenvectors of   "          "
};

class GMM
{
public:
    	Gaussian * getGaussians();  //相当于一个数组，指针指向第一个高斯结构

        // Initialize GMM with number of gaussians desired.
        GMM(unsigned int K);
        ~GMM();

        unsigned int K()const   //这都什么鬼函数名...
        {
                return m_K;     //混合模型个数
        }

        // Returns the probability density of color c in this GMM
        double p(Color_Lab c);   //计算混合高斯模型的概率值

        // Returns the probability density of color c in just Gaussian k
        double p(unsigned int i, Color_Lab c);     //计算某一个高斯模型的值

private:

        unsigned int m_K; // number of gaussians
        Gaussian *m_gaussians; // an array of K gaussians

//		friend void buildGMMs_kmeans(GMM &backgroundGMM, GMM &foregroundGMM, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation);
        friend void buildGMMs(GMM &backgroundGMM, 
			                  GMM &foregroundGMM, 
							  Image < unsigned int >  &components, 
							  const Image < Color_Lab >  &image, 
							  const Image < SegmentationValue >  &hardSegmentation);
 
		friend void learnGMMs(GMM &backgroundGMM, 
			                  GMM &foregroundGMM, 
							  Image < unsigned int >  &components, 
							  const Image < Color_Lab >  &image, 
							  const Image < SegmentationValue >  &hardSegmentation);
};

//void buildGMMs_kmeans(GMM &backgroundGMM, GMM &foregroundGMM, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation);

// Build the initial GMMs using the Orchard and Bouman color clustering algorithm
void buildGMMs(GMM &backgroundGMM, 
			   GMM &foregroundGMM, 
			   Image < unsigned int >  &components, 
			   const Image < Color_Lab >  &image, 
			   const Image < SegmentationValue >  &hardSegmentation);

// Iteratively learn GMMs using GrabCut updating algorithm
void learnGMMs(GMM &backgroundGMM, 
			   GMM &foregroundGMM, 
			   Image < unsigned int >  &components, 
			   const Image < Color_Lab >  &image, 
			   const Image < SegmentationValue >  &hardSegmentation);


void buildGMMs_kmeans(GMM &backgroundGMM, GMM &foregroundGMM, 
	                  const Image < Color_Lab >  &image, 
					  const Image < SegmentationValue >  &hardSegmentation);



// Helper class that fits a single Gaussian to color samples
class GaussianFitter     //Fitter何用？
{
public:
        GaussianFitter();

        // Add a color sample
        void add(Color_Lab c);

        // Build the gaussian out of all the added color samples
        void finalize(Gaussian &g, unsigned int totalCount, bool computeEigens = false)const;

private:

        Color_Lab s; // sum of L, a, and b
        double p[3][3]; // matrix of products (i.e. L*L, L*a, L*b), some values are duplicated.  

        unsigned int count; // count of color samples added to the gaussian
};

#endif //GMM_H
