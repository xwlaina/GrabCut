/*
* GrabCutTool implementation source code Copyright(c) 2005-2006 Justin Talbot
*
* All Rights Reserved.
* For educational use only; commercial use expressly forbidden.
* NO WARRANTY, express or implied, for this software.
*/

#pragma once

#include "GMM.h"
#include "AGMM.h"
#include "Tensor.h"
#include "GMM_Tensor.h"
#include "Gabor.h"
#include "GMM_ex.h"
#include <vector>


class GrabCutTool
{
public:
	GrabCutTool(const Image < Color_Lab >  *image ,  const IplImage *cv_image);
	~GrabCutTool();

	void initialize_color(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void initialize_tensor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void initialize_twoGraph(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void initialize_balanceNode(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void initialize_AGMMcolor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void initialize_gabor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);

	void setTrimap(int x1, int y1, int x2, int y2, const TrimapValue &t);

	void fitGMMs_color();
	void fitGMMs_tensor();
	void fitGMMs_twoGraph();
	void fitGMMs_balanceNode();
	void fitGMMs_AGMMcolor();
	void fitGMMs_gabor();

	void refine_color();
	void refine_tensor();
	void refine_twoGraph();
	void refine_balanceNode();
	void refine_AGMMcolor();
	void refine_gabor();

	int refineOnce_color(BOOL isNoiseRemoval = FALSE);
	int refineOnce_tensor(BOOL isNoiseRemoval = FALSE);
	int refineOnce_twoGraph(BOOL isNoiseRemoval = FALSE);
	int refineOnce_balanceNode(BOOL isNoiseRemoval = FALSE);
	int refineOnce_AGMMcolor(BOOL isNoiseRemoval = FALSE);
	int refineOnce_gabor(BOOL isNoiseRemoval = FALSE);

	const Image < SegmentationValue > * getSegmentationValue();
	const Image < SegmentationValue > * getSegmentationValue_tensor();
	const Image < SegmentationValue > * getSegmentationValue_twoGraph();
	const Image < SegmentationValue > * getSegmentationValue_balanceNode();
	const Image < SegmentationValue > * getSegmentationValue_AGMMcolor();
	const Image < SegmentationValue > * getSegmentationValue_gabor();

public:
	void setTrimap(int x, int y, int radius, const TrimapValue &t);
	GMM *m_backgroundGMM_color,  *m_foregroundGMM_color;
	AGMM *m_backgroundAGMM_color, *m_foregroundAGMM_color;
	GMM_Tensor *m_backgroundGMM_tensor,  *m_foregroundGMM_tensor;         /*纹理特征的前背景的高斯混合模型*/
	GMM_ex *m_backgroundGMM_gabor,  *m_foregroundGMM_gabor;  //Used for Gabor


private:
	//BOOL m_bStop;
	const IplImage *cvimage;
	unsigned int m_w, m_h; // All the following Image<*> variables will be the same width and height.
	unsigned int m_left, m_top, m_right, m_bottom;
	unsigned int m_nNumFGMMs, m_nNumBGMMs;
	//long time;

	Image < FLOAT >  *m_gradient_mag;
	// Store them here so we don't have to keep asking for them.
	const Image < Color_Lab >  *m_Grabimg;    

	Tensor* m_tensor;               //Texture image
	Gabor* m_gabor;

	//const Gabor* m_gabor;
	Image < TrimapValue >  *m_trimap;
	Image < unsigned int >  *m_GMMcomponent;
	Image < SegmentationValue >  *m_hardSegmentation_color;
	Image < SegmentationValue >  *m_hardSegmentation_tensor;   //Used to store final segmentation result
	Image < SegmentationValue >  *m_hardSegmentation_twoGraph;   //Used to store final segmentation result
	Image < SegmentationValue >  *m_hardSegmentation_balanceNode;   //Used to store final segmentation result
	Image < SegmentationValue >  *m_hardSegmentation_AGMMcolor;   //Used to store final segmentation result
	Image < SegmentationValue >  *m_hardSegmentation_gabor;   //Used to store final segmentation result

	double CalDEdgeWeight(unsigned int x, unsigned int y);

	int updateHardSegmentation_color(BOOL isNoiseRemoval = FALSE); 
	int updateHardSegmentation_tensor(BOOL isNoiseRemoval = FALSE); 
	int updateHardSegmentation_twoGraph(BOOL isNoiseRemoval = FALSE); 
	int updateHardSegmentation_balanceNode(BOOL isNoiseRemoval = FALSE); 
	int updateHardSegmentation_AGMMcolor(BOOL isNoiseRemoval = FALSE); 
	int updateHardSegmentation_gabor(BOOL isNoiseRemoval = FALSE); 

	// Variables used in formulas from the paper.
	double m_lambda; // lambda = 50. This value was suggested the GrabCutTool paper.
	double m_beta_color; 
	double m_beta_tensor; 
	double m_beta_gabor; 
	double m_L; // L = a large value to force a pixel to be foreground or background
	double m_dNoiseConstant; 

	void compute_tensorBeta();
	void compute_colorBeta();
	void compute_gaborBeta();
	void computeL();

	// Precomputed N-link weights
	Image < NLinks >  *m_NLinks_color;
	Image < NLinks >  *m_NLinks_tensor;  
	Image < NLinks >  *m_NLinks_gabor;  

	void computeNLinks_color();
	double computeNLink_color(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	double computeNLink_tensor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);//Compute the N-Link between two point in color or texture image
	void computeNLinks_tensor();

	double computeNLink_gabor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);//Compute the N-Link between two point in color or texture image
	void computeNLinks_gabor();

	// Graph for Graphcut
	GraphType *m_graph;
	Image < GraphType::node_id >  *m_nodes;

	void initGraph_color(); 
	void initGraph_tensor(); 
	void initGraph_twoGraph(); 
	void initGraph_balanceNode(); 
	void initGraph_AGMMcolor(); 
	void initGraph_gabor(); 

	void noiseRemoval(Image < SegmentationValue >  &hardSegmentation);
	void computeGradient();

public:
	void ImageSegmentationByColor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void ImageSegmentationBytensor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void ImageSegmentationBytwoGraph(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void ImageSegmentationBybalanceNode(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void ImageSegmentationByAGMMColor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
	void ImageSegmentationByGabor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2);
};

