#ifndef AGMM_H
#define AGMM_H

#include "Color_Lab.h"
#include "Global.h"
#include "Image.h"
#include <vector>

struct AGaussian
{
	Color_Lab mu;                
	double leftcovariance[3];  
	double rightcovariance[3];
	double pi;	
};

class AGMM
{
public:
	AGaussian * getAGaussians();

	AGMM(unsigned int K);
	~AGMM();

	unsigned int K()const
	{
		return m_K;
	}

	double p(const Color_Lab& c);   
	double p(unsigned int i, const Color_Lab& c);   

private:

	unsigned int m_K;      
	AGaussian *m_Agaussians; // an array of m_K Agaussians	
};

//k-meansÀ„∑®
void buildAGMMs_kmeans(AGMM &backgroundAGMM, AGMM &foregroundAGMM, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation);

void learnAGMMs(AGMM &backgroundAGMM, AGMM &foregroundAGMM, Image < unsigned int >  &components, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation);

class AGaussianFitter
{
public:
	AGaussianFitter();
	void clear();
	void add(const Color_Lab & c);   // Add a color sample
	void finalize(AGaussian &ag, unsigned int totalCount)const;// Build the gaussian out of all the added color samples

private:
	Color_Lab s;            // sum of Lab
	vector <Color_Lab> allcolor;	
	unsigned int count; // count of color samples added to the gaussian
};

#endif