#include "stdafx.h"
#include "GrabCutTool.h"
#include "Vector.h"
#include "Helper.h"
#include "Gabor.h"
#include <fstream>
#include <iomanip>
#include <ctime>
// #include <vld.h>


enum eCommand { GrabCutCmd,GrabCutMSNSTCmd,GrabCutTwoGraphCmd,GrabCutBalanceNodeCmd,GrabCutGaborCmd,AGMMGrabCutCmd,GrabCutCosegCmd,NoneCmd};

GrabCutTool::GrabCutTool(const Image < Color_Lab>  *image,  const IplImage *cv_image)
{
	//m_bStop = FALSE;
	m_Grabimg = image;
	cvimage = cv_image;

	m_w = m_Grabimg->width();
	m_h = m_Grabimg->height();

	m_trimap = new Image < TrimapValue > (m_w, m_h);  
	m_trimap->fill(TrimapUnknown);

	m_GMMcomponent = new Image < unsigned int > (m_w, m_h); 

	m_hardSegmentation_color = new Image < SegmentationValue > (m_w, m_h);  
	m_hardSegmentation_color->fill(SegmentationBackground);

	m_hardSegmentation_tensor = new Image < SegmentationValue > (m_w, m_h);  
	m_hardSegmentation_tensor->fill(SegmentationBackground);

	m_hardSegmentation_twoGraph = new Image < SegmentationValue > (m_w, m_h);  
	m_hardSegmentation_twoGraph->fill(SegmentationBackground);

	m_hardSegmentation_balanceNode = new Image < SegmentationValue > (m_w, m_h);  
	m_hardSegmentation_balanceNode->fill(SegmentationBackground);

	m_hardSegmentation_AGMMcolor = new Image < SegmentationValue > (m_w, m_h);  
	m_hardSegmentation_AGMMcolor->fill(SegmentationBackground);

	m_hardSegmentation_gabor = new Image < SegmentationValue > (m_w, m_h);  
	m_hardSegmentation_gabor->fill(SegmentationBackground);
	//	m_softSegmentation = 0; // Not yet implemented


	m_dNoiseConstant =1; 	

	m_nNumFGMMs = 5;
	m_nNumBGMMs = 5;


	//set some constants
	m_lambda = 50;   //默认值 50

	m_tensor = NULL;
	m_gabor = NULL;

	m_foregroundGMM_color = NULL;
	m_backgroundGMM_color = NULL;

	m_foregroundAGMM_color = NULL;
	m_backgroundAGMM_color = NULL;

	m_foregroundGMM_tensor = NULL;
	m_backgroundGMM_tensor = NULL;

	m_backgroundGMM_gabor = NULL;
	m_foregroundGMM_gabor = NULL;

	m_NLinks_color = NULL;
	m_NLinks_tensor = NULL;
	m_NLinks_gabor = NULL;

	computeL();
	m_graph = 0;
	eCommand m_eCurCmd = AGMMGrabCutCmd;

	switch(m_eCurCmd)
	{
	case GrabCutCmd:
		compute_colorBeta();
		computeNLinks_color();	
		m_nodes = new Image < GraphType::node_id > (m_w, m_h);   
		break;
	case AGMMGrabCutCmd:
		m_gabor = new Gabor(cvimage); 
		m_gabor->GaborFilteredImg();
		compute_colorBeta();
		computeNLinks_color();	
		compute_gaborBeta();
		computeNLinks_gabor();	
		m_nodes = new Image < GraphType::node_id > (m_w, m_h);   
		break;
	case GrabCutMSNSTCmd:
		m_tensor = new Tensor(cvimage);
		m_tensor->nlST();	
		compute_tensorBeta();	
		computeNLinks_tensor();
		m_nodes = new Image < GraphType::node_id > (m_w, m_h);   
		break;
	case GrabCutTwoGraphCmd:
		//long start,end;
		//start = ::GetTickCount();
		m_gabor = new Gabor(cvimage); 
		m_gabor->GaborFilteredImg();
		//end = ::GetTickCount();
		//time = end - start;
		compute_gaborBeta();
		computeNLinks_gabor();	
		compute_colorBeta();
		computeNLinks_color();	
		m_nodes = new Image < GraphType::node_id > (2*m_w, m_h); 
		break;
	case GrabCutBalanceNodeCmd:
		m_gabor = new Gabor(cvimage); 
		m_gabor->GaborFilteredImg();
		compute_gaborBeta();
		computeNLinks_gabor();	
		compute_colorBeta();
		computeNLinks_color();	
		m_nodes = new Image < GraphType::node_id > (3*m_w, m_h); 
		break;
	case GrabCutGaborCmd:
		m_gabor = new Gabor(cvimage); 
		m_gabor->GaborFilteredImg();
		compute_gaborBeta();
		computeNLinks_gabor();	
		m_nodes = new Image < GraphType::node_id > (m_w, m_h);
        break;
	default:
		m_nodes = NULL;
		break;
	}
	//ofstream of;
	//of.open("C:\\Users\\xu\\Desktop\\time.txt",ios_base::app);
	//of<<setw(10)<<setiosflags(ios::left)<<"gabor建模:"<<time<<endl;
	//of.close();
}

GrabCutTool::~GrabCutTool(void)
{
	if (m_tensor != NULL)
	{
		delete m_tensor;
		m_tensor = NULL;
	}
	if (m_gabor != NULL)
	{
		delete m_gabor;
		m_gabor = NULL;
	}
	if (m_trimap != NULL)
	{
		delete m_trimap;
		m_trimap = NULL;
	}
	if (m_GMMcomponent != NULL)
	{
		delete m_GMMcomponent;
		m_GMMcomponent = NULL;
	}
	if (m_hardSegmentation_color != NULL)
	{
		delete m_hardSegmentation_color;
		m_hardSegmentation_color = NULL;
	}
	if (m_hardSegmentation_tensor != NULL)
	{
		delete m_hardSegmentation_tensor;
		m_hardSegmentation_tensor = NULL;
	}
	if (m_hardSegmentation_twoGraph != NULL)
	{
		delete m_hardSegmentation_twoGraph;
		m_hardSegmentation_twoGraph = NULL;
	}
	if (m_hardSegmentation_balanceNode != NULL)
	{
		delete m_hardSegmentation_balanceNode;
		m_hardSegmentation_balanceNode = NULL;
	}
	if (m_hardSegmentation_AGMMcolor != NULL)
	{
		delete m_hardSegmentation_AGMMcolor;
		m_hardSegmentation_AGMMcolor = NULL;
	}
	if (m_hardSegmentation_gabor != NULL)
	{
		delete m_hardSegmentation_gabor;
		m_hardSegmentation_gabor = NULL;
	}
	if (m_foregroundGMM_color != NULL)
	{
		delete m_foregroundGMM_color;
		m_foregroundGMM_color = NULL;
	}
	if (m_backgroundGMM_color != NULL)
	{
		delete m_backgroundGMM_color;
		m_backgroundGMM_color = NULL;
	}
	if (m_foregroundAGMM_color != NULL)
	{
		delete m_foregroundAGMM_color;
		m_foregroundAGMM_color = NULL;
	}
	if (m_backgroundAGMM_color != NULL)
	{
		delete m_backgroundAGMM_color;
		m_backgroundAGMM_color = NULL;
	}
	if (m_foregroundGMM_tensor != NULL)
	{
		delete m_foregroundGMM_tensor;
		m_foregroundGMM_tensor = NULL;
	}
	if (m_backgroundGMM_tensor != NULL)
	{
		delete m_backgroundGMM_tensor;
		m_backgroundGMM_tensor= NULL;
	}
	if (m_foregroundGMM_gabor != NULL)
	{
		delete m_foregroundGMM_gabor;
		m_foregroundGMM_gabor = NULL;
	}
	if (m_backgroundGMM_gabor != NULL)
	{
		delete m_backgroundGMM_gabor;
		m_backgroundGMM_gabor = NULL;
	}
	if (m_NLinks_color != NULL)
	{
		delete m_NLinks_color;
		m_NLinks_color = NULL;
	}
	if (m_NLinks_tensor != NULL)
	{
		delete m_NLinks_tensor;
		m_NLinks_tensor = NULL;
	}
	if (m_NLinks_gabor != NULL)
	{
		delete m_NLinks_gabor;
		m_NLinks_gabor = NULL;
	}
	if (m_nodes != NULL)
	{
		delete m_nodes;
		m_nodes = NULL;
	}
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}
}

void GrabCutTool::initialize_color(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	m_left = __min(x1,x2);
	m_top = __min(y1,y2);
	m_right = __max(x1,x2);
	m_bottom = __max(y1,y2);

	// Step 1: User creates inital Trimap with rectangle, Background outside, Unknown inside
	m_trimap->fill(TrimapBackground);
	m_trimap->fillRectangle(x1, y1, x2, y2, TrimapUnknown);

	// Step 2: Initial segmentation, Background where Trimap is Background, Foreground where Trimap is Unknown.
	m_hardSegmentation_color->fill(SegmentationBackground);
	m_hardSegmentation_color->fillRectangle(x1, y1, x2, y2, SegmentationForeground);
}

void GrabCutTool::initialize_tensor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	m_left = __min(x1,x2);
	m_top = __min(y1,y2);
	m_right = __max(x1,x2);
	m_bottom = __max(y1,y2);


	// Step 1: User creates inital Trimap with rectangle, Background outside, Unknown inside
	m_trimap->fill(TrimapBackground);
	m_trimap->fillRectangle(x1, y1, x2, y2, TrimapUnknown);

	// Step 2: Initial segmentation, Background where Trimap is Background, Foreground where Trimap is Unknown.
	m_hardSegmentation_tensor->fill(SegmentationBackground);
	m_hardSegmentation_tensor->fillRectangle(x1, y1, x2, y2, SegmentationForeground);
}

void GrabCutTool::initialize_twoGraph(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	m_left = __min(x1,x2);
	m_top = __min(y1,y2);
	m_right = __max(x1,x2);
	m_bottom = __max(y1,y2);

	m_trimap->fill(TrimapBackground);
	m_trimap->fillRectangle(x1, y1, x2, y2, TrimapUnknown);

	m_hardSegmentation_color->fill(SegmentationBackground);
	m_hardSegmentation_color->fillRectangle(x1, y1, x2, y2, SegmentationForeground);

	m_hardSegmentation_gabor->fill(SegmentationBackground);
	m_hardSegmentation_gabor->fillRectangle(x1, y1, x2, y2, SegmentationForeground);

	m_hardSegmentation_twoGraph->fill(SegmentationBackground);
	m_hardSegmentation_twoGraph->fillRectangle(x1, y1, x2, y2, SegmentationForeground);
}

void GrabCutTool::initialize_gabor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	m_left = __min(x1,x2);
	m_top = __min(y1,y2);
	m_right = __max(x1,x2);
	m_bottom = __max(y1,y2);

	m_trimap->fill(TrimapBackground);
	m_trimap->fillRectangle(x1, y1, x2, y2, TrimapUnknown);

	m_hardSegmentation_gabor->fill(SegmentationBackground);
	m_hardSegmentation_gabor->fillRectangle(x1, y1, x2, y2, SegmentationForeground);

}

void GrabCutTool::initialize_balanceNode(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
    m_left = __min(x1,x2);
	m_top = __min(y1,y2);
	m_right = __max(x1,x2);
	m_bottom = __max(y1,y2);

	m_trimap->fill(TrimapBackground);
	m_trimap->fillRectangle(x1, y1, x2, y2, TrimapUnknown);

	m_hardSegmentation_color->fill(SegmentationBackground);
	m_hardSegmentation_color->fillRectangle(x1, y1, x2, y2, SegmentationForeground);

	m_hardSegmentation_gabor->fill(SegmentationBackground);
	m_hardSegmentation_gabor->fillRectangle(x1, y1, x2, y2, SegmentationForeground);

	m_hardSegmentation_balanceNode->fill(SegmentationBackground);
	m_hardSegmentation_balanceNode->fillRectangle(x1, y1, x2, y2, SegmentationForeground);
}

void GrabCutTool::initialize_AGMMcolor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	m_left = __min(x1,x2);
	m_top = __min(y1,y2);
	m_right = __max(x1,x2);
	m_bottom = __max(y1,y2);

	// Step 1: User creates inital Trimap with rectangle, Background outside, Unknown inside
	m_trimap->fill(TrimapBackground);
	m_trimap->fillRectangle(x1, y1, x2, y2, TrimapUnknown);

	// Step 2: Initial segmentation, Background where Trimap is Background, Foreground where Trimap is Unknown.
	m_hardSegmentation_AGMMcolor->fill(SegmentationBackground);
	m_hardSegmentation_AGMMcolor->fillRectangle(x1, y1, x2, y2, SegmentationForeground);
}

void GrabCutTool::fitGMMs_color()
{
	// Step 3: Build GMMs using Orchard-Bouman clustering algorithm
	//初始GMM就算的不对！
	m_foregroundGMM_color = new GMM(m_nNumFGMMs);   //默认值 5    //MonKey 有析构么？
	m_backgroundGMM_color = new GMM(m_nNumBGMMs);   //       5        //MonKey 有析构么？
	buildGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_color);	
	//MonKey 试试buildGMMs_kmeans
	//buildGMMs_kmeans(*m_backgroundGMM, *m_foregroundGMM, *m_Grabimg, *m_hardSegmentation);
}

void GrabCutTool::fitGMMs_tensor()
{
	// Step 3: Build GMMs using Orchard-Bouman clustering algorithm
	m_foregroundGMM_tensor = new GMM_Tensor(m_nNumFGMMs, m_tensor->Dim());
	m_backgroundGMM_tensor = new GMM_Tensor(m_nNumBGMMs, m_tensor->Dim());
	buildGMMs_kmeans(*m_backgroundGMM_tensor, *m_foregroundGMM_tensor, *(m_tensor->GetTensors()), *m_hardSegmentation_tensor, FALSE);
}

void GrabCutTool::fitGMMs_twoGraph()
{
	// Step 3: Build GMMs using Orchard-Bouman clustering algorithm
	m_foregroundGMM_color = new GMM(m_nNumFGMMs);   //默认值 5    //MonKey 有析构么？
	m_backgroundGMM_color = new GMM(m_nNumBGMMs);   //       5        //MonKey 有析构么？
	buildGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_color);

	m_foregroundGMM_gabor = new GMM_ex(m_nNumFGMMs, m_gabor->Dim());
	m_backgroundGMM_gabor = new GMM_ex(m_nNumBGMMs, m_gabor->Dim());
	buildGMMs_kmeans(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *(m_gabor->GetFilters()),  *m_hardSegmentation_gabor);
}

void GrabCutTool::fitGMMs_gabor()
{
	// Step 3: Build GMMs using Orchard-Bouman clustering algorithm
	m_foregroundGMM_gabor = new GMM_ex(m_nNumFGMMs, m_gabor->Dim());
	m_backgroundGMM_gabor = new GMM_ex(m_nNumBGMMs, m_gabor->Dim());
	buildGMMs_kmeans(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *(m_gabor->GetFilters()),  *m_hardSegmentation_gabor);
}

void GrabCutTool::fitGMMs_balanceNode()
{
	// Step 3: Build GMMs using Orchard-Bouman clustering algorithm
	m_foregroundGMM_color = new GMM(m_nNumFGMMs);   //默认值 5    //MonKey 有析构么？
	m_backgroundGMM_color = new GMM(m_nNumBGMMs);   //       5        //MonKey 有析构么？
	buildGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_color);

	m_foregroundGMM_gabor = new GMM_ex(m_nNumFGMMs, m_gabor->Dim());
	m_backgroundGMM_gabor = new GMM_ex(m_nNumBGMMs, m_gabor->Dim());
	buildGMMs_kmeans(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *(m_gabor->GetFilters()),  *m_hardSegmentation_gabor);
}

void GrabCutTool::fitGMMs_AGMMcolor()
{
	//m_foregroundAGMM_color = new AGMM(m_nNumFGMMs); 
	//m_backgroundAGMM_color = new AGMM(m_nNumBGMMs);  
	//buildAGMMs_kmeans(*m_backgroundAGMM_color,  *m_foregroundAGMM_color,  *m_Grabimg,  *m_hardSegmentation_AGMMcolor);	
	m_foregroundGMM_color = new GMM(m_nNumFGMMs);   //默认值 5    //MonKey 有析构么？
	m_backgroundGMM_color = new GMM(m_nNumBGMMs);   //       5        //MonKey 有析构么？
	buildGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_AGMMcolor);

	m_foregroundGMM_gabor = new GMM_ex(m_nNumFGMMs, m_gabor->Dim());
	m_backgroundGMM_gabor = new GMM_ex(m_nNumBGMMs, m_gabor->Dim());
	buildGMMs_kmeans(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *(m_gabor->GetFilters()),  *m_hardSegmentation_AGMMcolor);
}

int GrabCutTool::refineOnce_color(BOOL isNoiseRemoval)
{
	double flow = 0;

	// Steps 4 and 5: Learn new GMMs from current segmentation
	learnGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_color);

	// Step 6: Run GraphCut and update segmentation
	initGraph_color();
	if (m_graph)
	{
		flow = m_graph->maxflow();
	}

	int changed = updateHardSegmentation_color(isNoiseRemoval);  //MonKey into this function leak problem !
	return changed;
}

int GrabCutTool::refineOnce_gabor(BOOL isNoiseRemoval)
{

	double flow = 0;

	learnGMMs(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *m_GMMcomponent,  *(m_gabor->GetFilters()),  *m_hardSegmentation_gabor);

	initGraph_gabor();
	if (m_graph)
	{
		flow = m_graph->maxflow();
	}

	int changed = updateHardSegmentation_gabor(isNoiseRemoval);  //MonKey into this function leak problem !

	return changed;
}

int GrabCutTool::refineOnce_tensor(BOOL isNoiseRemoval)
{

	double flow = 0;

	// Steps 4 and 5: Learn new GMMs from current segmentation
	learnGMMs(*m_backgroundGMM_tensor, *m_foregroundGMM_tensor, *m_GMMcomponent, *(m_tensor->GetTensors()), *m_hardSegmentation_tensor);

	// Step 6: Run GraphCut and update segmentation
	initGraph_tensor();
	if (m_graph)
	{
		flow = m_graph->maxflow();
	}

	int changed = updateHardSegmentation_tensor(isNoiseRemoval);  //MonKey into this function leak problem !


	return changed;
}


int GrabCutTool::refineOnce_twoGraph(BOOL isNoiseRemoval)
{

	double flow = 0;

	// Steps 4 and 5: Learn new GMMs from current segmentation
	learnGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_twoGraph);
	learnGMMs(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *m_GMMcomponent,  *(m_gabor->GetFilters()),  *m_hardSegmentation_twoGraph);

	// Step 6: Run GraphCut and update segmentation
	initGraph_twoGraph();
	if (m_graph)
	{
		//long start,end;
		//start = ::GetTickCount();
		flow = m_graph->maxflow();
		//end = ::GetTickCount();
		//ofstream of;
		//of.open("C:\\Users\\xu\\Desktop\\time.txt",ios_base::app);
		//of<<setw(10)<<setiosflags(ios::left)<<"maxflow:"<<end-start<<endl;
		//of.close();
	}

	int changed = updateHardSegmentation_twoGraph(isNoiseRemoval);  
	return changed;
}

int GrabCutTool::refineOnce_balanceNode(BOOL isNoiseRemoval)
{

	double flow = 0;

	// Steps 4 and 5: Learn new GMMs from current segmentation
	learnGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_balanceNode);
	learnGMMs(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *m_GMMcomponent,  *(m_gabor->GetFilters()),  *m_hardSegmentation_balanceNode);

	// Step 6: Run GraphCut and update segmentation
	initGraph_balanceNode();
	if (m_graph)
	{
		flow = m_graph->maxflow();
	}

	int changed = updateHardSegmentation_balanceNode(isNoiseRemoval);  //MonKey into this function leak problem !

	return changed;
}


int GrabCutTool::refineOnce_AGMMcolor(BOOL isNoiseRemoval)
{

	double flow = 0;

	// Steps 4 and 5: Learn new GMMs from current segmentation
	//learnAGMMs(*m_backgroundAGMM_color,  *m_foregroundAGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_AGMMcolor);
	learnGMMs(*m_backgroundGMM_color,  *m_foregroundGMM_color,  *m_GMMcomponent,  *m_Grabimg,  *m_hardSegmentation_AGMMcolor);
	learnGMMs(*m_backgroundGMM_gabor,  *m_foregroundGMM_gabor,  *m_GMMcomponent,  *(m_gabor->GetFilters()),  *m_hardSegmentation_AGMMcolor);

	// Step 6: Run GraphCut and update segmentation
	initGraph_AGMMcolor();
	if (m_graph)
	{
		flow = m_graph->maxflow();
	}

	int changed = updateHardSegmentation_AGMMcolor(isNoiseRemoval);  

	return changed;
}

void GrabCutTool::refine_color()
{
	int changed = m_w * m_h;
	//m_bStop = FALSE;
	while (changed > 10)
	{
		changed = refineOnce_color(TRUE);  //MonKey into this function leak problem !
	}
	//m_bStop = FALSE;
}

void GrabCutTool::refine_tensor()
{
	int changed = m_w * m_h;
	/*m_bStop = FALSE;*/
	while (changed > 10)
	{
		changed = refineOnce_tensor(TRUE);  //MonKey into this function leak problem !
	}
	//m_bStop = FALSE;
}

void GrabCutTool::refine_twoGraph()
{
	int changed = m_w * m_h;
	//m_bStop = FALSE;
	int count = 6;
	while (changed > 10 && count > 1)
	{
		changed = refineOnce_twoGraph(TRUE);  
		count--;
	}
	//m_bStop = FALSE;
}

void GrabCutTool::refine_balanceNode()
{
	int changed = m_w * m_h;
	/*m_bStop = FALSE;*/
	int count = 8 ;
	while (changed > 10 && count > 1)
	{
		changed = refineOnce_balanceNode(TRUE);  
		count--;
	}
	/*m_bStop = FALSE;*/
}

void GrabCutTool::refine_AGMMcolor()
 {
	int changed = m_w * m_h;
	//m_bStop = FALSE;
	while (changed > 10)
	{
		changed = refineOnce_AGMMcolor(TRUE);  //MonKey into this function leak problem !
	}
	//m_bStop = FALSE;
}

void GrabCutTool::refine_gabor()
{
	int changed = m_w * m_h;
	/*m_bStop = FALSE;*/
	int count = 5 ;
	while (changed > 10 && count > 1)
	{
		changed = refineOnce_gabor(TRUE);  //MonKey into this function leak problem !
		count--;
	}
	//m_bStop = FALSE;
}

int GrabCutTool::updateHardSegmentation_color(BOOL isNoiseRemoval)
{
	int changed = 0;
	unsigned int x,y;

	if (isNoiseRemoval == TRUE)
	{
		unsigned int rect_w,rect_h;
		rect_w = m_right - m_left + 1;
		rect_h = m_bottom - m_top + 1;
		SegmentationValue *old_seg = new SegmentationValue[rect_w * rect_h]; 
		unsigned int index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				old_seg[index++] = (*m_hardSegmentation_color)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationBackground;
					}
				}
			}
		}

		noiseRemoval(*m_hardSegmentation_color);

		index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				if (old_seg[index++] != (*m_hardSegmentation_color)(x, y))
				{
					changed++;
				}
			}
		}
		if (old_seg != NULL)
		{
			delete old_seg;
			old_seg = NULL;
		}
	} 
	else 
	{
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				SegmentationValue oldValue = (*m_hardSegmentation_color)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationBackground;
					}
				}

				if (oldValue != (*m_hardSegmentation_color)(x, y))
				{
					changed++;
				}
			}
		}
	}
	return changed;
}

int GrabCutTool::updateHardSegmentation_AGMMcolor(BOOL isNoiseRemoval)
{
	int changed = 0;
	unsigned int x,y;

	if (isNoiseRemoval == TRUE)
	{
		unsigned int rect_w,rect_h;
		rect_w = m_right - m_left + 1;
		rect_h = m_bottom - m_top + 1;
		SegmentationValue *old_seg = new SegmentationValue[rect_w * rect_h];  //MonKey new 没有对应的delete
		unsigned int index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				old_seg[index++] = (*m_hardSegmentation_AGMMcolor)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_AGMMcolor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_AGMMcolor)(x, y) = SegmentationBackground;
					}
				}
			}
		}

		noiseRemoval(*m_hardSegmentation_AGMMcolor);

		index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				if (old_seg[index++] != (*m_hardSegmentation_AGMMcolor)(x, y))
				{
					changed++;
				}
			}
		}
		if (old_seg != NULL)
		{
			delete old_seg;
			old_seg = NULL;
		}
	} 
	else 
	{
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				SegmentationValue oldValue = (*m_hardSegmentation_AGMMcolor)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_AGMMcolor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_AGMMcolor)(x, y) = SegmentationBackground;
					}
				}

				if (oldValue != (*m_hardSegmentation_AGMMcolor)(x, y))
				{
					changed++;
				}
			}
		}
	}
	return changed;
}

int GrabCutTool::updateHardSegmentation_tensor(BOOL isNoiseRemoval)
{
	int changed = 0;
	unsigned int x,y;

	if (isNoiseRemoval == TRUE)
	{
		unsigned int rect_w,rect_h;
		rect_w = m_right - m_left + 1;
		rect_h = m_bottom - m_top + 1;
		SegmentationValue *old_seg = new SegmentationValue[rect_w * rect_h]; 
		unsigned int index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				old_seg[index++] = (*m_hardSegmentation_tensor)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_tensor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_tensor)(x, y) = SegmentationBackground;
					}
				}
			}
		}

		noiseRemoval(*m_hardSegmentation_tensor);

		index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				if (old_seg[index++] != (*m_hardSegmentation_tensor)(x, y))
				{
					changed++;
				}
			}
		}
		if (old_seg != NULL)
		{
			delete old_seg;
			old_seg = NULL;
		}

	} 

	else 
	{
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				SegmentationValue oldValue = (*m_hardSegmentation_tensor)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_tensor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_tensor)(x, y) = SegmentationBackground;
					}
				}

				if (oldValue != (*m_hardSegmentation_tensor)(x, y))
				{
					changed++;
				}
			}
		}
	}

	return changed;
}

int GrabCutTool::updateHardSegmentation_twoGraph(BOOL isNoiseRemoval)
{
	int changed = 0;
	unsigned int x,y;

	if (isNoiseRemoval == TRUE)
	{
		unsigned int rect_w,rect_h;
		rect_w = m_right - m_left + 1;
		rect_h = m_bottom - m_top + 1;
		SegmentationValue *old_seg = new SegmentationValue[rect_w * rect_h]; 
		unsigned int index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				old_seg[index++] = (*m_hardSegmentation_twoGraph)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x/*+m_w*/, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_twoGraph)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_twoGraph)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x+m_w, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationBackground;
					}
				}
			}
		}

		noiseRemoval(*m_hardSegmentation_twoGraph);

		index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				if (old_seg[index++] != (*m_hardSegmentation_twoGraph)(x, y))
				{
					changed++;
				}
			}
		}
		if (old_seg != NULL)
		{
			delete old_seg;
			old_seg = NULL;
		}

	} 
	else 
	{
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				SegmentationValue oldValue = (*m_hardSegmentation_twoGraph)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x/*+m_w*/, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_twoGraph)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_twoGraph)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x+m_w, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationBackground;
					}

				}

				if (oldValue != (*m_hardSegmentation_tensor)(x, y))
				{
					changed++;
				}
			}
		}
	}
	return changed;
}

int GrabCutTool::updateHardSegmentation_balanceNode(BOOL isNoiseRemoval)
{
	int changed = 0;
	unsigned int x,y;

	if (isNoiseRemoval == TRUE)
	{
		unsigned int rect_w,rect_h;
		rect_w = m_right - m_left + 1;
		rect_h = m_bottom - m_top + 1;
		SegmentationValue *old_seg = new SegmentationValue[rect_w * rect_h]; 
		unsigned int index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				old_seg[index++] = (*m_hardSegmentation_balanceNode)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x+m_w, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_balanceNode)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_balanceNode)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x+2*m_w, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationBackground;
					}
				}
			}
		}

		noiseRemoval(*m_hardSegmentation_balanceNode);

		index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				if (old_seg[index++] != (*m_hardSegmentation_balanceNode)(x, y))
				{
					changed++;
				}
			}
		}
		if (old_seg != NULL)
		{
			delete old_seg;
			old_seg = NULL;
		}

	} 

	else 
	{
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				SegmentationValue oldValue = (*m_hardSegmentation_balanceNode)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x+m_w, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_balanceNode)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_balanceNode)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_color)(x, y) = SegmentationBackground;
					}

					if (m_graph->what_segment((*m_nodes)(x+2*m_w, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationBackground;
					}
				}

				if (oldValue != (*m_hardSegmentation_balanceNode)(x, y))
				{
					changed++;
				}
			}
		}
	}

	return changed;
}

int GrabCutTool::updateHardSegmentation_gabor(BOOL isNoiseRemoval)
{
	int changed = 0;
	unsigned int x,y;

	if (isNoiseRemoval == TRUE)
	{
		unsigned int rect_w,rect_h;
		rect_w = m_right - m_left + 1;
		rect_h = m_bottom - m_top + 1;
		SegmentationValue *old_seg = new SegmentationValue[rect_w * rect_h]; 
		unsigned int index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				old_seg[index++] = (*m_hardSegmentation_gabor)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationBackground;
					}
				}
			}
		}

		noiseRemoval(*m_hardSegmentation_gabor);

		index = 0;
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				if (old_seg[index++] != (*m_hardSegmentation_gabor)(x, y))
				{
					changed++;
				}
			}
		}
		if (old_seg != NULL)
		{
			delete old_seg;
			old_seg = NULL;
		}

	} 

	else 
	{
		for (y = m_top; y <= m_bottom; ++y)
		{
			for (x = m_left; x <= m_right; ++x)
			{
				SegmentationValue oldValue = (*m_hardSegmentation_gabor)(x, y);
				if ((*m_trimap)(x, y) == TrimapUnknown)
				{
					if (m_graph->what_segment((*m_nodes)(x, y)) == GraphType::SOURCE)
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationForeground;
					}
					else
					{
						(*m_hardSegmentation_gabor)(x, y) = SegmentationBackground;
					}
				}

				if (oldValue != (*m_hardSegmentation_gabor)(x, y))
				{
					changed++;
				}
			}
		}
	}

	return changed;
}

void GrabCutTool::setTrimap(int x1, int y1, int x2, int y2, const TrimapValue &t)
{
	(*m_trimap).fillRectangle(x1, y1, x2, y2, t);

	// Immediately set the segmentation as well so that the display will update.
	if (t == TrimapForeground)
	{
		(*m_hardSegmentation_twoGraph).fillRectangle(x1, y1, x2, y2, SegmentationForeground);
	}
	else if (t == TrimapBackground)
	{
		(*m_hardSegmentation_twoGraph).fillRectangle(x1, y1, x2, y2, SegmentationBackground);
	}

}

//private functions
void GrabCutTool::initGraph_color()
{
	unsigned int x, y, nx, ny;


	// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}

	// 只是在未知的区域构建图

	unsigned int nodesCnt = 0;
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				nodesCnt++;
		}
	}

	/*增加结点*/
	m_graph = new GraphType(nodesCnt,8*nodesCnt,error_function);
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				(*m_nodes)(x, y) = m_graph->add_node();
			}
		}
	}

	double backColor, foreColor,edgeweightColor;

	//构造传统的图结构
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				//注意计算方式
				foreColor =  - log(m_backgroundGMM_color->p((*m_Grabimg)(x, y)));
				backColor =  - log(m_foregroundGMM_color->p((*m_Grabimg)(x, y)));
				//FILE* pFile;
				//pFile=fopen("D:\\Program Files\\projects\\a.txt","a");   
				//fprintf(pFile,"fore: %f\n\n",foreColor);
				//fprintf(pFile,"back: %f\n\n",backColor);
				//fclose(pFile);
				if ((x > 0) && (y > 0))//upleft
				{
					nx = x - 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upleft;
					assert (edgeweightColor > 0 );
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						break;
					case TrimapForeground:
						foreColor   +=  edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}

				if (y > 0)//up
				{
					nx = x;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).up;
					assert (edgeweightColor > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}

				if ((x < (m_w - 1)) && (y > 0))//upright
				{
					nx = x + 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upright;   	
					assert (edgeweightColor > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}

				if (x < (m_w - 1))//right
				{
					nx = x + 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(x, y).right;

					assert (edgeweightColor > 0 );
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}
				if ((x < (m_w - 1)) && (y < (m_h - 1)))//downright
				{
					nx = x + 1;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).upleft;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}

				if (y < (m_h - 1))//down
				{
					nx = x;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).up;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}

				if ((x > 0) && (y < (m_h - 1)))//downleft
				{
					nx = x - 1;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).upright;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}

				if (x > 0)//left
				{
					nx = x - 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(nx, ny).right;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						break;
					}
				}
				m_graph->add_tweights((*m_nodes)(x, y), foreColor, backColor);
			}
		}
	}
}

void GrabCutTool::initGraph_tensor()
{
	unsigned int x, y, nx, ny;

	// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}

	// 只是在未知的区域构建图
	unsigned int nodesCnt = 0;
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				nodesCnt++;
		}
	}

	m_graph = new GraphType(nodesCnt,8*nodesCnt,error_function);
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				(*m_nodes)(x, y) = m_graph->add_node();
		}
	}

	//获得梯度图像
	Image < CVector* >* texImage = NULL;
	texImage = m_tensor->GetTensors();


	double back, fore, edgeweight;

	//构造图
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				//注意计算方式
				fore = - log(m_backgroundGMM_tensor->p(((*texImage)(x, y))->addr()));
				back = - log(m_foregroundGMM_tensor->p(((*texImage)(x, y))->addr()));

				if ((x > 0) && (y > 0))//upleft
				{
					nx = x - 1;
					ny = y - 1;
					edgeweight = (*m_NLinks_tensor)(x, y).upleft;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (y > 0)//up
				{
					nx = x;
					ny = y - 1;
					edgeweight = (*m_NLinks_tensor)(x, y).up;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}


				if ((x < (m_w - 1)) && (y > 0))//upright
				{
					nx = x + 1;
					ny = y - 1;
					edgeweight = (*m_NLinks_tensor)(x, y).upright;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (x < (m_w - 1))//right
				{
					nx = x + 1;
					ny = y;
					edgeweight = (*m_NLinks_tensor)(x, y).right;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if ((x < (m_w - 1)) && (y < (m_h - 1)))//downright
				{
					nx = x + 1;
					ny = y + 1;
					edgeweight = (*m_NLinks_tensor)(nx, ny).upleft;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (y < (m_h - 1))//down
				{
					nx = x;
					ny = y + 1;
					edgeweight = (*m_NLinks_tensor)(nx, ny).up;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if ((x > 0) && (y < (m_h - 1)))//downleft
				{
					nx = x - 1;
					ny = y + 1;
					edgeweight = (*m_NLinks_tensor)(nx, ny).upright;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (x > 0)//left
				{
					nx = x - 1;
					ny = y;
					edgeweight = (*m_NLinks_tensor)(nx, ny).right;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}
				m_graph->add_tweights((*m_nodes)(x, y), fore, back);
			}
		}
	}
}

void GrabCutTool::initGraph_twoGraph()
{
	unsigned int x, y, nx, ny;


	// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}

	// 只是在未知的区域构建图
	unsigned int nodesCnt = 0;
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				nodesCnt++;
		}
	}

	//新的颜色纹理合成方式,两层图
	//结点数：nodesNum=2*nodesCnt;
	//边数：  4*nodesCnt+4*nodesCnt+2*nodesCnt=10*nodesCnt;
	//增加结点
	m_graph = new GraphType(2*nodesCnt,10*nodesCnt,error_function);
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				(*m_nodes)(x, y)      = m_graph->add_node();
				(*m_nodes)(x+m_w,y)   = m_graph->add_node();
			}
		}
	}

	
	Image < CVector* >* texImage = NULL;
	texImage = m_gabor->GetFilters();

	double backColor, foreColor, edgeweightColor;
	double backTexture, foreTexture, edgeweightTexture;
	double dEdgeWeight;

	//构造图
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				//注意计算方式
				foreColor =  - log(m_backgroundGMM_color->p((*m_Grabimg)(x, y)));
				foreTexture = - log(m_backgroundGMM_gabor->p(((*texImage)(x, y))->addr()));

				backColor =  - log(m_foregroundGMM_color->p((*m_Grabimg)(x, y)));
				backTexture = - log(m_foregroundGMM_gabor->p(((*texImage)(x, y))->addr()));

				double lamda;
				if (abs(backTexture-foreTexture) > abs(backColor-foreColor))
				{
					lamda = 50;
				}
				else
				{
					lamda = 5;
				}
			/*	ofstream of;
				of.open("C:\\Users\\xu\\Desktop\\P.txt",ios_base::app);
				of<<setw(10)<<setiosflags(ios::left)<<"foreColor"<<foreColor<<endl;
				of<<setw(10)<<setiosflags(ios::left)<<"backColor"<<backColor<<endl;
				of<<setw(10)<<setiosflags(ios::left)<<"foreTexture"<<foreTexture<<endl;
				of<<setw(10)<<setiosflags(ios::left)<<"backTexture"<<backTexture<<endl;
				of.close();*/

				dEdgeWeight=CalDEdgeWeight(x,y);
				//增加惩罚边，双向的
				//if(dEdgeWeight <= 0)
				//{
				//	dEdgeWeight = dEdgeWeight+1;
				//}
				if(dEdgeWeight <= 0)
					dEdgeWeight = -dEdgeWeight;
				assert (dEdgeWeight > 0 );
				m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(x+m_w, y), lamda*dEdgeWeight, lamda*dEdgeWeight);

				if ((x > 0) && (y > 0))//upleft
				{
					nx = x - 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upleft;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).upleft;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+m_w, y),(*m_nodes)(nx+m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor   +=  edgeweightColor;
						foreTexture +=  edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (y > 0)//up
				{
					nx = x;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).up;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).up;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;


					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+m_w, y),(*m_nodes)(nx+m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}


				if ((x < (m_w - 1)) && (y > 0))//upright
				{
					nx = x + 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upright;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).upright;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+m_w, y),(*m_nodes)(nx+m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (x < (m_w - 1))//right
				{
					nx = x + 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(x, y).right;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).right;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+m_w, y),(*m_nodes)(nx+m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}

				}


				if ((x < (m_w - 1)) && (y < (m_h - 1)))//downright
				{
					nx = x + 1;
					ny = y + 1;

					edgeweightColor = (*m_NLinks_color)(nx, ny).upleft;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).upleft;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (y < (m_h - 1))//down
				{
					nx = x;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).up;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).up;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if ((x > 0) && (y < (m_h - 1)))//downleft
				{
					nx = x - 1;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).upright;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).upright;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (x > 0)//left
				{
					nx = x - 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(nx, ny).right;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).right;

					//edgeweightColor *= lamda1;
					//edgeweightTexture *= lamda2;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}
				m_graph->add_tweights((*m_nodes)(x, y), foreColor, backColor);
				m_graph->add_tweights((*m_nodes)(x+m_w, y), foreTexture, backTexture);
			}
		}
	}
}

void GrabCutTool::initGraph_balanceNode()
{
	unsigned int x, y, nx, ny;


	// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}

	// 只是在未知的区域构建图
	unsigned int nodesCnt = 0;
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				nodesCnt++;
		}
	}

	//新的颜色纹理合成方式,两层图
	//结点数：nodesNum=2*nodesCnt;
	//边数：  4*nodesCnt+4*nodesCnt+2*nodesCnt=10*nodesCnt;
	//增加结点
	m_graph = new GraphType(3*nodesCnt,10*nodesCnt,error_function);
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				(*m_nodes)(x, y)      = m_graph->add_node();
				(*m_nodes)(x+m_w,y)   = m_graph->add_node();
				(*m_nodes)(x+2*m_w,y)  = m_graph->add_node();
			}
		}
	}

	Image < CVector* >* texImage = NULL;
	texImage = m_gabor->GetFilters();


	double backColor, foreColor, edgeweightColor;
	double backTexture, foreTexture, edgeweightTexture;
	double dEdgeWeight;
	//构造图
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				//注意计算方式
				foreColor =  - log(m_backgroundGMM_color->p((*m_Grabimg)(x, y)));
				foreTexture = - log(m_backgroundGMM_gabor->p(((*texImage)(x, y))->addr()));

				backColor =  - log(m_foregroundGMM_color->p((*m_Grabimg)(x, y)));
				backTexture = - log(m_foregroundGMM_gabor->p(((*texImage)(x, y))->addr()));

				double lamda1,lamda2;
				if (abs(backTexture-foreTexture) > abs(backColor-foreColor))
				{
					lamda1 = 5;
					lamda2 = 50;
					//lamda1 = 50;
					//lamda2 = 5;
				}
				else
				{
					lamda1 = 50;
					lamda2 = 5;
				}

				dEdgeWeight=CalDEdgeWeight(x,y);
				//增加惩罚边，双向的
				/*cout<<dEdgeWeight<<endl;*/
				if(dEdgeWeight == 0 ) dEdgeWeight = dEdgeWeight + 0.001;
				assert(dEdgeWeight > 0);
				m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(x+m_w, y), lamda1*dEdgeWeight, lamda1*dEdgeWeight);
				m_graph->add_edge((*m_nodes)(x+m_w, y), (*m_nodes)(x+2*m_w, y), lamda2*dEdgeWeight, lamda2*dEdgeWeight);

				if ((x > 0) && (y > 0))//upleft
				{
					nx = x - 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upleft;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).upleft;
					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+2*m_w, y),(*m_nodes)(nx+2*m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor   +=  edgeweightColor;
						foreTexture +=  edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (y > 0)//up
				{
					nx = x;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).up;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).up;
					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+2*m_w, y),(*m_nodes)(nx+2*m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}


				if ((x < (m_w - 1)) && (y > 0))//upright
				{
					nx = x + 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upright;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).upright;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+2*m_w, y),(*m_nodes)(nx+2*m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (x < (m_w - 1))//right
				{
					nx = x + 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(x, y).right;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).right;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
						m_graph->add_edge((*m_nodes)(x+2*m_w, y),(*m_nodes)(nx+2*m_w, ny),edgeweightTexture,edgeweightTexture);
						break;
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}

				}


				if ((x < (m_w - 1)) && (y < (m_h - 1)))//downright
				{
					nx = x + 1;
					ny = y + 1;

					edgeweightColor = (*m_NLinks_color)(nx, ny).upleft;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).upleft;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (y < (m_h - 1))//down
				{
					nx = x;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).up;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).up;
					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if ((x > 0) && (y < (m_h - 1)))//downleft
				{
					nx = x - 1;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).upright;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).upright;

					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}

				if (x > 0)//left
				{
					nx = x - 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(nx, ny).right;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).right;
					assert (edgeweightColor > 0 );
					assert (edgeweightTexture > 0);
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						foreColor += edgeweightColor;
						foreTexture += edgeweightTexture;
						break;
					case TrimapBackground:
						backColor += edgeweightColor;
						backTexture += edgeweightTexture;
						break;
					}
				}
				m_graph->add_tweights((*m_nodes)(x, y), foreColor, backColor);
				m_graph->add_tweights((*m_nodes)(x+2*m_w, y), foreTexture, backTexture);
			}
		}
	}
}

void GrabCutTool::initGraph_AGMMcolor()
{
	unsigned int x, y, nx, ny;


	// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}

	// 只是在未知的区域构建图
	unsigned int nodesCnt = 0;
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				nodesCnt++;
		}
	}

	//新的颜色纹理合成方式,两层图
	//结点数：nodesNum=2*nodesCnt;
	//边数：  4*nodesCnt+4*nodesCnt+2*nodesCnt=10*nodesCnt;
	//增加结点
	m_graph = new GraphType(nodesCnt,8*nodesCnt,error_function);
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				(*m_nodes)(x, y)      = m_graph->add_node();
			}
		}
	}

	Image < CVector* >* texImage = NULL;
	texImage = m_gabor->GetFilters();


	double backColor, foreColor, edgeweightColor;
	double backTexture, foreTexture, edgeweightTexture;
	double fore, back, edgeweight;
	double epsilon = 0.9;
	//构造图
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				//注意计算方式
				foreColor =  - log(m_backgroundGMM_color->p((*m_Grabimg)(x, y)));
				foreTexture = - log(m_backgroundGMM_gabor->p(((*texImage)(x, y))->addr()));

				backColor =  - log(m_foregroundGMM_color->p((*m_Grabimg)(x, y)));
				backTexture = - log(m_foregroundGMM_gabor->p(((*texImage)(x, y))->addr()));
				
				fore = (1-epsilon)*foreColor + epsilon*foreTexture;
				back = (1-epsilon)*backColor + epsilon*backTexture;

				if ((x > 0) && (y > 0))//upleft
				{
					nx = x - 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upleft;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).upleft;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore   +=  edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (y > 0)//up
				{
					nx = x;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).up;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).up;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore   +=  edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}


				if ((x < (m_w - 1)) && (y > 0))//upright
				{
					nx = x + 1;
					ny = y - 1;
					edgeweightColor = (*m_NLinks_color)(x, y).upright;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).upright;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore   +=  edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (x < (m_w - 1))//right
				{
					nx = x + 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(x, y).right;
					edgeweightTexture =(*m_NLinks_gabor)(x,y).right;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore   +=  edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}

				}


				if ((x < (m_w - 1)) && (y < (m_h - 1)))//downright
				{
					nx = x + 1;
					ny = y + 1;

					edgeweightColor = (*m_NLinks_color)(nx, ny).upleft;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).upleft;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );
		
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (y < (m_h - 1))//down
				{
					nx = x;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).up;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).up;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}

				}

				if ((x > 0) && (y < (m_h - 1)))//downleft
				{
					nx = x - 1;
					ny = y + 1;
					edgeweightColor = (*m_NLinks_color)(nx, ny).upright;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).upright;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}

				}

				if (x > 0)//left
				{
					nx = x - 1;
					ny = y;
					edgeweightColor = (*m_NLinks_color)(nx, ny).right;
					edgeweightTexture =(*m_NLinks_gabor)(nx,ny).right;
					edgeweight = (1-epsilon)*edgeweightColor + epsilon*edgeweightTexture;
					assert (edgeweight > 0 );

					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}
				m_graph->add_tweights((*m_nodes)(x, y), fore, back);
			}
		}
	}
	//unsigned int x, y, nx, ny;


	//// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	//if (m_graph != NULL)
	//{
	//	delete m_graph;
	//	m_graph = NULL;
	//}

	//// 只是在未知的区域构建图

	//unsigned int nodesCnt = 0;
	//for (y = 0; y < m_h; ++y)
	//{
	//	for (x = 0; x < m_w; ++x)
	//	{
	//		if ((*m_trimap)(x, y) == TrimapUnknown)
	//			nodesCnt++;
	//	}
	//}

	///*增加结点*/
	//m_graph = new GraphType(nodesCnt,8*nodesCnt,error_function);
	//for (y = 0; y < m_h; ++y)
	//{
	//	for (x = 0; x < m_w; ++x)
	//	{
	//		if ((*m_trimap)(x, y) == TrimapUnknown)
	//		{
	//			(*m_nodes)(x, y) = m_graph->add_node();
	//		}
	//	}
	//}

	//double backColor, foreColor,edgeweightColor;

	////构造传统的图结构
	//for (y = 0; y < m_h; ++y)
	//{
	//	for (x = 0; x < m_w; ++x)
	//	{
	//		if ((*m_trimap)(x, y) == TrimapUnknown)
	//		{
	//			//注意计算方式
	//			foreColor =  - log(m_backgroundAGMM_color->p((*m_Grabimg)(x, y)));
	//			backColor =  - log(m_foregroundAGMM_color->p((*m_Grabimg)(x, y)));

	//			if ((x > 0) && (y > 0))//upleft
	//			{
	//				nx = x - 1;
	//				ny = y - 1;
	//				edgeweightColor = (*m_NLinks_color)(x, y).upleft;
	//				assert (edgeweightColor > 0 );

	//				switch((*m_trimap)(nx, ny))
	//				{
	//				case TrimapUnknown:
	//					m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
	//					break;
	//				case TrimapForeground:
	//					foreColor   +=  edgeweightColor;
	//					break;
	//				case TrimapBackground:
	//					backColor += edgeweightColor;
	//					break;
	//				}
	//			}

	//			if (y > 0)//up
	//			{
	//				nx = x;
	//				ny = y - 1;
	//				edgeweightColor = (*m_NLinks_color)(x, y).up;
	//				assert (edgeweightColor > 0 );

	//				switch((*m_trimap)(nx, ny))
	//				{
	//				case TrimapUnknown:
	//					m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
	//					break;
	//				case TrimapForeground:
	//					foreColor += edgeweightColor;
	//					break;
	//				case TrimapBackground:
	//					backColor += edgeweightColor;
	//					break;
	//				}
	//			}


	//			if ((x < (m_w - 1)) && (y > 0))//upright
	//			{
	//				nx = x + 1;
	//				ny = y - 1;
	//				edgeweightColor = (*m_NLinks_color)(x, y).upright;   	
	//				assert (edgeweightColor > 0 );

	//				switch((*m_trimap)(nx, ny))
	//				{
	//				case TrimapUnknown:
	//					m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
	//					break;
	//				case TrimapForeground:
	//					foreColor += edgeweightColor;
	//					break;
	//				case TrimapBackground:
	//					backColor += edgeweightColor;
	//					break;
	//				}
	//			}

	//			if (x < (m_w - 1))//right
	//			{
	//				nx = x + 1;
	//				ny = y;
	//				edgeweightColor = (*m_NLinks_color)(x, y).right;

	//				assert (edgeweightColor > 0 );
	//				switch((*m_trimap)(nx, ny))
	//				{
	//				case TrimapUnknown:
	//					m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweightColor, edgeweightColor);
	//					break;
	//				case TrimapForeground:
	//					foreColor += edgeweightColor;
	//					break;
	//				case TrimapBackground:
	//					backColor += edgeweightColor;
	//					break;
	//				}

	//			}
	//			m_graph->add_tweights((*m_nodes)(x, y), foreColor, backColor);
	//		}
	//	}
	//}
}

void GrabCutTool::initGraph_gabor()
{
	unsigned int x, y, nx, ny;

	// Set up the graph (it can only be used once, so we have to recreate it each time the graph is updated)
	if (m_graph != NULL)
	{
		delete m_graph;
		m_graph = NULL;
	}

	// 只是在未知的区域构建图
	unsigned int nodesCnt = 0;
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				nodesCnt++;
		}
	}

	m_graph = new GraphType(nodesCnt,8*nodesCnt,error_function);
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
				(*m_nodes)(x, y) = m_graph->add_node();
		}
	}

	//获得梯度图像
	Image < CVector* >* texImage = NULL;
	texImage = m_gabor->GetFilters();


	double back, fore, edgeweight;

	//构造图
	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if ((*m_trimap)(x, y) == TrimapUnknown)
			{
				//注意计算方式
				fore = - log(m_backgroundGMM_gabor->p(((*texImage)(x, y))->addr()));
				back = - log(m_foregroundGMM_gabor->p(((*texImage)(x, y))->addr()));

				if ((x > 0) && (y > 0))//upleft
				{
					nx = x - 1;
					ny = y - 1;
					edgeweight = (*m_NLinks_gabor)(x, y).upleft;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (y > 0)//up
				{
					nx = x;
					ny = y - 1;
					edgeweight = (*m_NLinks_gabor)(x, y).up;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}


				if ((x < (m_w - 1)) && (y > 0))//upright
				{
					nx = x + 1;
					ny = y - 1;
					edgeweight = (*m_NLinks_gabor)(x, y).upright;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (x < (m_w - 1))//right
				{
					nx = x + 1;
					ny = y;
					edgeweight = (*m_NLinks_gabor)(x, y).right;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapUnknown:
						m_graph->add_edge((*m_nodes)(x, y), (*m_nodes)(nx, ny), edgeweight, edgeweight);
						break;
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if ((x < (m_w - 1)) && (y < (m_h - 1)))//downright
				{
					nx = x + 1;
					ny = y + 1;
					edgeweight = (*m_NLinks_gabor)(nx, ny).upleft;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (y < (m_h - 1))//down
				{
					nx = x;
					ny = y + 1;
					edgeweight = (*m_NLinks_gabor)(nx, ny).up;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if ((x > 0) && (y < (m_h - 1)))//downleft
				{
					nx = x - 1;
					ny = y + 1;
					edgeweight = (*m_NLinks_gabor)(nx, ny).upright;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}

				if (x > 0)//left
				{
					nx = x - 1;
					ny = y;
					edgeweight = (*m_NLinks_gabor)(nx, ny).right;
					switch((*m_trimap)(nx, ny))
					{
					case TrimapForeground:
						fore += edgeweight;
						break;
					case TrimapBackground:
						back += edgeweight;
						break;
					}
				}
				m_graph->add_tweights((*m_nodes)(x, y), fore, back);
			}
		}
	}
}

void GrabCutTool::computeNLinks_color()
{
	m_NLinks_color = new Image < NLinks > (m_w, m_h);   
	for (unsigned int y = 0; y < m_h; ++y)
	{
		for (unsigned int x = 0; x < m_w; ++x)
		{
			if (x > 0 && y > 0)
			{
				(*m_NLinks_color)(x, y).upleft = computeNLink_color(x, y, x - 1, y - 1);

				//ofstream of;
				//of.open("C:\\Users\\xu\\Desktop\\NlinkColor.txt",ios_base::app);
				//of<<setw(10)<<setiosflags(ios::left)<<(*m_NLinks_color)(x, y).upleft<<endl;
				//of.close();
			}

			if (y > 0)
			{
				(*m_NLinks_color)(x, y).up = computeNLink_color(x, y, x, y - 1);
		/*		ofstream of;
				of.open("C:\\Users\\xu\\Desktop\\NlinkColor.txt",ios_base::app);
				of<<setw(10)<<setiosflags(ios::left)<<(*m_NLinks_color)(x, y).up<<endl;
				of.close();*/
			}

			if (x < m_w - 1 && y > 0)
			{
				(*m_NLinks_color)(x, y).upright = computeNLink_color(x, y, x + 1, y - 1);
			/*	ofstream of;
				of.open("C:\\Users\\xu\\Desktop\\NlinkColor.txt",ios_base::app);
				of<<setw(10)<<setiosflags(ios::left)<<(*m_NLinks_color)(x, y).upright<<endl;
				of.close();*/
			}

			if (x < m_w - 1)
			{
				(*m_NLinks_color)(x, y).right = computeNLink_color(x, y, x + 1, y);
			/*	ofstream of;
				of.open("C:\\Users\\xu\\Desktop\\NlinkColor.txt",ios_base::app);
				of<<setw(10)<<setiosflags(ios::left)<<(*m_NLinks_color)(x, y).right<<endl;
				of.close();*/
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void GrabCutTool::computeNLinks_tensor()
{
	m_NLinks_tensor = new Image < NLinks > (m_w, m_h);  

	unsigned int x,y;

	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if (x > 0 && y > 0)
			{
				(*m_NLinks_tensor)(x, y).upleft = computeNLink_tensor(x, y, x - 1, y - 1);
				//ofstream of;
				//of.open("C:\\Users\\xu\\Desktop\\NlinkTensor.txt",ios_base::app);
				//of<<setw(10)<<setiosflags(ios::left)<<(*m_NLinks_tensor)(x, y).upleft<<endl;
				//of.close();
			}
			if (y > 0)
			{
				(*m_NLinks_tensor)(x, y).up = computeNLink_tensor(x, y, x, y - 1);
			}

			if (x < m_w - 1 && y > 0)
			{
				(*m_NLinks_tensor)(x, y).upright = computeNLink_tensor(x, y, x + 1, y - 1);
			}

			if (x < m_w - 1)
			{
				(*m_NLinks_tensor)(x, y).right = computeNLink_tensor(x, y, x + 1, y);
			}
		}
	}
}

void GrabCutTool::computeNLinks_gabor()
{
	m_NLinks_gabor = new Image < NLinks > (m_w, m_h);  

	unsigned int x,y;

	for (y = 0; y < m_h; ++y)
	{
		for (x = 0; x < m_w; ++x)
		{
			if (x > 0 && y > 0)
			{
				(*m_NLinks_gabor)(x, y).upleft = computeNLink_gabor(x, y, x - 1, y - 1);
				//ofstream of;
				//of.open("C:\\Users\\xu\\Desktop\\NlinkGabor.txt",ios_base::app);
				//of<<setw(10)<<setiosflags(ios::left)<<(*m_NLinks_gabor)(x, y).upleft<<endl;
				//of.close();
			}
			if (y > 0)
			{
				(*m_NLinks_gabor)(x, y).up = computeNLink_gabor(x, y, x, y - 1);
			}

			if (x < m_w - 1 && y > 0)
			{
				(*m_NLinks_gabor)(x, y).upright = computeNLink_gabor(x, y, x + 1, y - 1);
			}

			if (x < m_w - 1)
			{
				(*m_NLinks_gabor)(x, y).right = computeNLink_gabor(x, y, x + 1, y);
			}
		}
	}
}

double GrabCutTool::computeNLink_color(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	return /*m_lambda*/ 10*exp( - m_beta_color * distance2((*m_Grabimg)(x1, y1), (*m_Grabimg)(x2, y2))) / distance(x1, y1, x2, y2);
}

double GrabCutTool::computeNLink_tensor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	return m_lambda * exp( - m_beta_tensor * m_tensor->computeDistance2( x1, y1, x2, y2)) / distance(x1, y1, x2, y2) + m_dNoiseConstant;
}

double GrabCutTool::computeNLink_gabor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	return /*m_lambda*/50 * exp( - m_beta_gabor * m_gabor->computeDistance2( x1, y1, x2, y2)) / distance(x1, y1, x2, y2) + m_dNoiseConstant;
}

void GrabCutTool::compute_tensorBeta()
{
	int edges = 0;

	double sum_texture = 0.0;

	for (unsigned int y = 0; y < m_h; ++y)
	{
		for (unsigned int x = 0; x < m_w; ++x)
		{
			if (x > 0 && y >0)
				// upleft
			{
				sum_texture += m_tensor->computeDistance2( x, y, x - 1, y - 1);
				edges++;
			}

			if (y > 0)
				// up
			{
				sum_texture += m_tensor->computeDistance2(x, y, x, y - 1);

				edges++;
			}

			if (x < m_w - 1 && y > 0)
				// upright
			{
				sum_texture += m_tensor->computeDistance2(x, y, x + 1, y - 1);

				edges++;
			}

			if (x < m_w - 1)
				// right
			{
				sum_texture += m_tensor->computeDistance2(x, y, x + 1, y);

				edges++;
			}
		}
	}

	m_beta_tensor = 0.0;
	m_beta_tensor = (double)(1.0 / (2.0 * sum_texture / edges));

}

void GrabCutTool::compute_gaborBeta()
{
	int edges = 0;

	double sum_texture = 0.0;

	for (unsigned int y = 0; y < m_h; ++y)
	{
		for (unsigned int x = 0; x < m_w; ++x)
		{
			if (x > 0 && y >0)
				// upleft
			{
				sum_texture += m_gabor->computeDistance2( x, y, x - 1, y - 1);
				edges++;
			}

			if (y > 0)
				// up
			{
				sum_texture += m_gabor->computeDistance2(x, y, x, y - 1);

				edges++;
			}

			if (x < m_w - 1 && y > 0)
				// upright
			{
				sum_texture += m_gabor->computeDistance2(x, y, x + 1, y - 1);

				edges++;
			}

			if (x < m_w - 1)
				// right
			{
				sum_texture += m_gabor->computeDistance2(x, y, x + 1, y);

				edges++;
			}
		}
	}

	m_beta_gabor = 0.0;
	m_beta_gabor = (double)(1.0 / (2.0 * sum_texture / edges));

}

void GrabCutTool::compute_colorBeta()
{

	double result = 0;
	int edges = 0;


	for (unsigned int y = 0; y < m_h; ++y)
	{
		for (unsigned int x = 0; x < m_w; ++x)
		{
			if (x > 0 && y >0)
				// upleft
			{
				result += distance2((*m_Grabimg)(x, y), (*m_Grabimg)(x - 1, y - 1));
				edges++;
			}

			if (y > 0)
				// up
			{
				result += distance2((*m_Grabimg)(x, y), (*m_Grabimg)(x, y - 1));

				edges++;
			}

			if (x < m_w - 1 && y > 0)
				// upright
			{
				result += distance2((*m_Grabimg)(x, y), (*m_Grabimg)(x + 1, y - 1));

				edges++;
			}

			if (x < m_w - 1)
				// right
			{
				result += distance2((*m_Grabimg)(x, y), (*m_Grabimg)(x + 1, y));
				edges++;
			}
		}
	}
	m_beta_color = 0.0;
	m_beta_color = (double)(1.0 / (2 *result / edges)); 
}

void GrabCutTool::computeL()
{
	m_L = 10000;//8 * m_lambda + 1;
}

const Image < SegmentationValue > * GrabCutTool::getSegmentationValue()
{
	return m_hardSegmentation_color;
}

const Image < SegmentationValue > * GrabCutTool::getSegmentationValue_tensor()
{
	return m_hardSegmentation_tensor;
}

const Image < SegmentationValue > * GrabCutTool::getSegmentationValue_gabor()
{
	return m_hardSegmentation_gabor;
}

const Image < SegmentationValue > * GrabCutTool::getSegmentationValue_twoGraph()
{
	return m_hardSegmentation_twoGraph;
}

const Image < SegmentationValue > * GrabCutTool::getSegmentationValue_balanceNode()
{
	return m_hardSegmentation_balanceNode;
}

const Image < SegmentationValue > * GrabCutTool::getSegmentationValue_AGMMcolor()
{
	return m_hardSegmentation_AGMMcolor;
}

//void GrabCutTool::stopRefine()
//{
//	m_bStop = TRUE;
//}

void GrabCutTool::setTrimap(int x, int y, int radius, const TrimapValue &t)
{
	(*m_trimap).fillCircle(x, y, radius, t);
	if (t == TrimapForeground)
	{
		(*m_hardSegmentation_twoGraph).fillCircle(x, y, radius, SegmentationForeground);
	}
	else if (t == TrimapBackground)
	{
		(*m_hardSegmentation_twoGraph).fillCircle(x, y, radius, SegmentationBackground);
	}
}

void GrabCutTool::noiseRemoval(Image < SegmentationValue >  &hardSegmentation)
{
	unsigned int x,y;
	unsigned int radius = 1;
	IplImage* segImg = cvCreateImage( cvSize( m_w, m_h ), 8, 1 );
	cvZero(segImg);

	for (y=0;y<m_h;y++)
	{
		for (x=0;x<m_w;x++)
		{
			uchar* dst = &CV_IMAGE_ELEM( segImg, uchar, y, x );
			dst[0] = (uchar)((((hardSegmentation)(x,y))==SegmentationForeground)?255:0);
		}
	}
	cvDilate(segImg, segImg, NULL, radius); 
	cvErode(segImg, segImg, NULL, radius);
	cvErode(segImg, segImg, NULL, radius);
	cvDilate(segImg, segImg, NULL, radius);

	for (y=0;y<m_h;y++)
	{
		for (x=0;x<m_w;x++)
		{
			uchar* dst = &CV_IMAGE_ELEM( segImg, uchar, y, x );
			if(dst[0] != 0)
				(hardSegmentation)(x,y) = SegmentationForeground;
			else
				(hardSegmentation)(x,y) = SegmentationBackground;
		}
	}

	cvReleaseImage(&segImg);
}

void GrabCutTool::computeGradient()
{
	m_gradient_mag = new Image < FLOAT > (m_w, m_h);
	//	GetGradient(m_image,m_gradient_mag->ptr());
}

void GrabCutTool::ImageSegmentationByColor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{

	initialize_color(x1,y1,x2,y2);
	fitGMMs_color();
	refine_color(); 
}

void GrabCutTool::ImageSegmentationBytensor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	initialize_tensor(x1,y1,x2,y2);
	fitGMMs_tensor();
	refine_tensor(); 
}

void GrabCutTool::ImageSegmentationBytwoGraph(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	initialize_twoGraph(x1,y1,x2,y2);
	fitGMMs_twoGraph();
	refine_twoGraph(); 
}

void GrabCutTool::ImageSegmentationBybalanceNode(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	initialize_balanceNode(x1,y1,x2,y2);
	fitGMMs_balanceNode();
	refine_balanceNode(); 
}

void GrabCutTool::ImageSegmentationByAGMMColor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{

	initialize_AGMMcolor(x1,y1,x2,y2);
	fitGMMs_AGMMcolor();
	refine_AGMMcolor(); 
}

void GrabCutTool::ImageSegmentationByGabor(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
	initialize_gabor(x1,y1,x2,y2);
	fitGMMs_gabor();
	refine_gabor();
}

double GrabCutTool::CalDEdgeWeight(unsigned int x, unsigned int y)
{
	double foreColor,backColor;
	double foreTexture,backTexture;
	double dEdgeWeight=0.0;

	Image < CVector* >* texImage = NULL;
	texImage = m_gabor->GetFilters();

	foreColor = m_foregroundGMM_color->p((*m_Grabimg)(x, y));
	backColor = m_backgroundGMM_color->p((*m_Grabimg)(x, y));
	foreColor = foreColor/(foreColor+backColor);

	foreTexture = m_foregroundGMM_gabor->p(((*texImage)(x, y))->addr());
	backTexture = m_backgroundGMM_gabor->p(((*texImage)(x, y))->addr());
	foreTexture = foreTexture/(foreTexture + backTexture);

	//	TRACE("(%f,%f)\n",foreColor,foreTexture);

	dEdgeWeight = 1-(fabs(foreColor-foreTexture));
	//	TRACE("dEdgeWeight=%f\n",dEdgeWeight);

	dEdgeWeight = 0.5*dEdgeWeight+0.1;
	
	//dEdgeWeight = 2*dEdgeWeight;
	//	TRACE("dEdgeWeight=%f\n",dEdgeWeight);

	//ofstream of;
	//of.open("C:\\Users\\xu\\Desktop\\dEdge.txt",ios_base::app);
	//of<<setw(10)<<setiosflags(ios::left)<<"dEdgeWeight:"<<foreColor<<endl;
	//of.close();
	//dEdgeWeight = 10*dEdgeWeight;
	//dEdgeWeight += 0.1;
	return dEdgeWeight;
}