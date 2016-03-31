/*
 * GrabCut implementation source code Copyright(c) 2005-2006 Justin Talbot
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 */
#include "stdafx.h"
#include "GMM.h"

#define CV_TYPE CV_64FC1

GMM::GMM(unsigned int K): m_K(K)
{
        m_gaussians = new Gaussian[m_K];  //Gaussian 数组     //MonKey 有析构么？
}

//-------------------------------------------------------------------------

GMM::~GMM()
{
        if (m_gaussians)
        {
                delete [] m_gaussians;
        }
}

//-------------------------------------------------------------------------

double GMM::p(Color_Lab c)
{
	double result = 0;
	
	if (m_gaussians)
	{
		for (unsigned int i = 0; i < m_K; i++)
		{
			result += m_gaussians[i].pi *p(i, c);
		}
	}
	
	return result;
}

//-------------------------------------------------------------------------

double GMM::p(unsigned int i, Color_Lab c)    //高斯概率密度函数在某点处的值
{
        double result = 0;

        if (m_gaussians[i].pi > 0)
        {
                if (m_gaussians[i].determinant > 0)
                {
                        double L = c.L - m_gaussians[i].mu.L;
                        double a = c.a - m_gaussians[i].mu.a;
                        double b = c.b - m_gaussians[i].mu.b;

                        double d = L *(L *m_gaussians[i].inverse[0][0] + 
							          a * m_gaussians[i].inverse[1][0] + 
									  b * m_gaussians[i].inverse[2][0]) + 

									  a *(L *m_gaussians[i].inverse[0][1] + 
									      a * m_gaussians[i].inverse[1][1] + 
									      b * m_gaussians[i].inverse[2][1]) + 

									  b *(L *m_gaussians[i].inverse[0][2] + 
									      a * m_gaussians[i].inverse[1][2] + 
									      b * m_gaussians[i].inverse[2][2]);    

                        result = (double)(1.0 / (sqrt(m_gaussians[i].determinant)) *exp( - 0.5 * d));
                }
        }
        return result;
}

//-------------------------------------------------------------------------

void buildGMMs_kmeans(GMM &backgroundGMM, GMM &foregroundGMM, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
	GaussianFitter *backFitters = new GaussianFitter[backgroundGMM.K()];
	GaussianFitter *foreFitters = new GaussianFitter[foregroundGMM.K()];

	unsigned int foreCount = 0, backCount = 0, f_index = 0, b_index = 0;
	unsigned int x,y,i;
	
	for (y = 0; y < image.height(); ++y)
	{
		for (x = 0; x < image.width(); ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				foreCount++;
			}
			else
			{
				backCount++;
			}
		}
	}

	CvMat* f_Colors = cvCreateMat( foreCount, 1, CV_32FC3 );
	CvPoint3D32f* f_pt = (CvPoint3D32f*)f_Colors->data.fl;

	CvMat* b_Colors = cvCreateMat( backCount, 1, CV_32FC3 );
	CvPoint3D32f* b_pt = (CvPoint3D32f*)b_Colors->data.fl;

	for (y = 0; y < image.height(); ++y)
	{
		for (x = 0; x < image.width(); ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				f_pt->x = (float)((image(x,y)).L);
				f_pt->y = (float)((image(x,y)).a);
				f_pt->z = (float)((image(x,y)).b);
				f_pt++;
			}
			else
			{
				b_pt->x = (float)((image(x,y)).L);
				b_pt->y = (float)((image(x,y)).a);
				b_pt->z = (float)((image(x,y)).b);
				b_pt++;
			}
		}
	}

	CvMat* f_clusters = cvCreateMat( foreCount, 1, CV_32SC1 );
	cvKMeans2( f_Colors, foregroundGMM.K(), f_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, 1.0 ));  //KMeans 聚类和下面的Orchard-Bouman 聚类有什么差别？
	
	CvMat* b_clusters = cvCreateMat( backCount, 1, CV_32SC1 );
	cvKMeans2( b_Colors, backgroundGMM.K(), b_clusters, cvTermCriteria( CV_TERMCRIT_EPS+CV_TERMCRIT_ITER, 10, 1.0 ));
	for (y = 0; y < image.height(); ++y)
	{
		for (x = 0; x < image.width(); ++x)
		{	
			if (hardSegmentation(x, y) == SegmentationForeground)
			{
				foreFitters[f_clusters->data.i[f_index++]].add(image(x, y));
			}
			else
			{
				backFitters[b_clusters->data.i[b_index++]].add(image(x, y));
			}
		}
	}


	/*
	for (i = 0; i < backgroundGMM.K(); i++)
	{
		backFitters[i].finalize(backgroundGMM.m_gaussians[i], backCount, false);
	}
	
	for (i = 0; i < foregroundGMM.K(); i++)
	{
		foreFitters[i].finalize(foregroundGMM.m_gaussians[i], foreCount, false);
	}
	*/

	for (i = 0; i <backgroundGMM.K(); i++)
	{
		backFitters[i].finalize((backgroundGMM.getGaussians())[i], backCount);
	}

	for (i = 0; i <foregroundGMM.K(); i++)
	{
		foreFitters[i].finalize((foregroundGMM.getGaussians())[i], foreCount);
	}
	
	delete [] backFitters;
	delete [] foreFitters;
	cvReleaseMat( &f_clusters );
	cvReleaseMat( &f_Colors );
	cvReleaseMat( &b_clusters );
	cvReleaseMat( &b_Colors );
}

//-------------------------------------------------------------------------

void buildGMMs(GMM &backgroundGMM, GMM &foregroundGMM, Image < unsigned int >  &components, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
        // Step 3: Build GMMs using Orchard-Bouman clustering algorithm

        // Set up Gaussian Fitters
        GaussianFitter *backFitters = new GaussianFitter[backgroundGMM.K()];   //5个fitter?   
        GaussianFitter *foreFitters = new GaussianFitter[foregroundGMM.K()];   // *XXXfitter 与  *m_gaussians 什么关系？   一一对应？
  
        unsigned int foreCount = 0, backCount = 0;

        // Initialize the first foreground and background clusters
		//扫描每个点 添加到 foreFitter OR backFitter
        for (unsigned int y = 0; y < image.height(); ++y)
        {
                for (unsigned int x = 0; x < image.width(); ++x)
                {
                        components(x, y) = 0;   //GMMcomponent 成员初始化

                        if (hardSegmentation(x, y) == SegmentationForeground)
                        {
                                foreFitters[0].add(image(x, y));
                                foreCount++;
                        }
                        else
                        {
                                backFitters[0].add(image(x, y));
                                backCount++;
                        }
                }
        }

        backFitters[0].finalize(backgroundGMM.m_gaussians[0], backCount, true);     //最初的finalize 要获取什么信息？
        foreFitters[0].finalize(foregroundGMM.m_gaussians[0], foreCount, true);     //获取最初的 分割点

        unsigned int nBack = 0, nFore = 0; // Which cluster will be split
        unsigned int maxK = backgroundGMM.K() > foregroundGMM.K() ? backgroundGMM.K(): foregroundGMM.K();
		
        // Compute clusters
        for (unsigned int i = 1; i < maxK; i++)
        {
                // Reset the fitters for the splitting clusters
                backFitters[nBack] = GaussianFitter();    
                foreFitters[nFore] = GaussianFitter();

                // For brevity, get references to the splitting Gaussians
                Gaussian &bg = backgroundGMM.m_gaussians[nBack];
                Gaussian &fg = foregroundGMM.m_gaussians[nFore];

                // Compute splitting points
                double splitBack = bg.eigenvectors[0][0] *bg.mu.L + bg.eigenvectors[1][0] *bg.mu.a + bg.eigenvectors[2][0] *bg.mu.b;
                double splitFore = fg.eigenvectors[0][0] *fg.mu.L + fg.eigenvectors[1][0] *fg.mu.a + fg.eigenvectors[2][0] *fg.mu.b;

                // Split clusters nBack and nFore, place split portion into cluster i
                for (unsigned int y = 0; y < image.height(); ++y)
                {
                        for (unsigned int x = 0; x < image.width(); ++x)
                        {
                                Color_Lab c = image(x, y);

                                // For each pixel
                                if (i < foregroundGMM.K() && hardSegmentation(x, y) == SegmentationForeground && components(x, y) == nFore)
                                {
                                        if (fg.eigenvectors[0][0] *c.L + fg.eigenvectors[1][0] *c.a + fg.eigenvectors[2][0] *c.b > splitFore)
                                        {
                                                components(x, y) = i;
                                                foreFitters[i].add(c);   //分到另一个 XXXFitter[i]中
                                        }
                                        else
                                        {
                                                foreFitters[nFore].add(c);
                                        }
                                }
                                else if (i < backgroundGMM.K() && hardSegmentation(x, y) == SegmentationBackground && components(x, y) == nBack)
                                {
                                        if (bg.eigenvectors[0][0] *c.L + bg.eigenvectors[1][0] *c.a + bg.eigenvectors[2][0] *c.b > splitBack)
                                        {
                                                components(x, y) = i;
                                                backFitters[i].add(c);
                                        }
                                        else
                                        {
                                                backFitters[nBack].add(c);
                                        }
                                }
                        }
                }


                // Compute new split Gaussians
                backFitters[nBack].finalize(backgroundGMM.m_gaussians[nBack], backCount, true);  //这里开始比例不是1了...
                foreFitters[nFore].finalize(foregroundGMM.m_gaussians[nFore], foreCount, true); 


                if (i < backgroundGMM.K())
                {
                        backFitters[i].finalize(backgroundGMM.m_gaussians[i], backCount, true);
                }
                if (i < foregroundGMM.K())
                {
                        foreFitters[i].finalize(foregroundGMM.m_gaussians[i], foreCount, true);
                }

                // Find clusters with highest eigenvalue
                nBack = 0;
                nFore = 0;

                for (unsigned int j = 0; j <= i; j++)
                {
                        if (j < backgroundGMM.K() && backgroundGMM.m_gaussians[j].eigenvalues[0] > backgroundGMM.m_gaussians[nBack].eigenvalues[0])
                        {
                                nBack = j;
                        }

                        if (j < foregroundGMM.K() && foregroundGMM.m_gaussians[j].eigenvalues[0] > foregroundGMM.m_gaussians[nFore].eigenvalues[0])
                        {
                                nFore = j;
                        }
                }
        }
        
        delete [] backFitters;
        delete [] foreFitters;
		//至此back and fore 各自5个高斯密度函数参数求得
}

//-------------------------------------------------------------------------

void learnGMMs(GMM &backgroundGMM, GMM &foregroundGMM, Image < unsigned int >  &components, const Image < Color_Lab >  &image, const Image < SegmentationValue >  &hardSegmentation)
{
        unsigned int x, y;
        unsigned int i;

        // Step 4: Assign each pixel to the component which maximizes its probability
        for (y = 0; y < image.height(); ++y)
        {
                for (x = 0; x < image.width(); ++x)
                {
                        Color_Lab c = image(x, y);

                        if (hardSegmentation(x, y) == SegmentationForeground)
                        {
                                int k = 0;
                                double max = 0;

                                for (i = 0; i < foregroundGMM.K(); i++)
                                {
                                        double p = foregroundGMM.p(i, c);
                                        if (p > max)
                                        {
                                                k = i;
                                                max = p;
                                        }
                                }

                                components(x, y) = k;   //把每个pixel分配到"概率最大"的m_gaussians[i] 的过程中 pixel(x,y) 对应的components(x,y) 标记值改变
                        }
                        else
                        {
                                int k = 0;
                                double max = 0;

                                for (i = 0; i < backgroundGMM.K(); i++)
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
        GaussianFitter *backFitters = new GaussianFitter[backgroundGMM.K()];     //重置GaussianFitter 数组
        GaussianFitter *foreFitters = new GaussianFitter[foregroundGMM.K()];

        unsigned int foreCount = 0, backCount = 0;

        for (y = 0; y < image.height(); ++y)
        {
                for (x = 0; x < image.width(); ++x)
                {
                        Color_Lab c = image(x, y);

                        if (hardSegmentation(x, y) == SegmentationForeground)
                        {
                                foreFitters[components(x, y)].add(c);
                                foreCount++;
                        }
                        else
                        {
                                backFitters[components(x, y)].add(c);
                                backCount++;
                        }
                }
        }

        for (i = 0; i < backgroundGMM.K(); i++)
        {
                backFitters[i].finalize(backgroundGMM.m_gaussians[i], backCount, false);
        }

        for (i = 0; i < foregroundGMM.K(); i++)
        {
                foreFitters[i].finalize(foregroundGMM.m_gaussians[i], foreCount, false);
        }

        delete [] backFitters;
        delete [] foreFitters;
}


// GaussianFitter functions
GaussianFitter::GaussianFitter()
{
        s = Color_Lab();

        p[0][0] = 0;   //这种初始化...醉了double p[3][3]={0};
        p[0][1] = 0;
        p[0][2] = 0;
        p[1][0] = 0;
        p[1][1] = 0;
        p[1][2] = 0;
        p[2][0] = 0;
        p[2][1] = 0;
        p[2][2] = 0;

        count = 0;
}

// Add a color sample
void GaussianFitter::add(Color_Lab c)
{
        s.L += c.L;
        s.a += c.a;
        s.b += c.b;

        p[0][0] += c.L * c.L;   //P矩阵何用？    // LL La Lb   对于每个Color_Lab c 累加3*3 个数据
        p[0][1] += c.L * c.a;                    // aL aa ab
        p[0][2] += c.L * c.b;   //    ↓         // bL ba bb
        p[1][0] += c.a * c.L;
        p[1][1] += c.a * c.a;   //用于下面计算协方差矩阵 (过渡数据)
        p[1][2] += c.a * c.b;
        p[2][0] += c.b * c.L;
        p[2][1] += c.b * c.a;
        p[2][2] += c.b * c.b;

        count++;
}

// Build the gaussian out of all the added colors
void GaussianFitter::finalize(Gaussian &g, unsigned int totalCount, bool computeEigens)const
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
                // Compute mean of gaussian
                g.mu.L = s.L / count;
                g.mu.a = s.a / count;
                g.mu.b = s.b / count;

                // Compute covariance matrix
				//
				g.covariance[0][0] = (p[0][0] - count*g.mu.L * g.mu.L) / (count - 1) + Epsilon;
                g.covariance[0][1] = (p[0][1] - count*g.mu.L * g.mu.a) / (count - 1);
                g.covariance[0][2] = (p[0][2] - count*g.mu.L * g.mu.b) / (count - 1);
                g.covariance[1][0] = (p[1][0] - count*g.mu.a * g.mu.L) / (count - 1);
                g.covariance[1][1] = (p[1][1] - count*g.mu.a * g.mu.a) / (count - 1) + Epsilon;
                g.covariance[1][2] = (p[1][2] - count*g.mu.a * g.mu.b) / (count - 1);
                g.covariance[2][0] = (p[2][0] - count*g.mu.b * g.mu.L) / (count - 1);
                g.covariance[2][1] = (p[2][1] - count*g.mu.b * g.mu.a) / (count - 1);
                g.covariance[2][2] = (p[2][2] - count*g.mu.b * g.mu.b) / (count - 1) + Epsilon;

				//协方差公式 系数为 1/(n-1)   ↑

                /*  
	            g.covariance[0][0] = p[0][0] / count - g.mu.L * g.mu.L + Epsilon;
                g.covariance[0][1] = p[0][1] / count - g.mu.L * g.mu.a;
                g.covariance[0][2] = p[0][2] / count - g.mu.L * g.mu.b;
                g.covariance[1][0] = p[1][0] / count - g.mu.a * g.mu.L;
                g.covariance[1][1] = p[1][1] / count - g.mu.a * g.mu.a + Epsilon;
                g.covariance[1][2] = p[1][2] / count - g.mu.a * g.mu.b;
                g.covariance[2][0] = p[2][0] / count - g.mu.b * g.mu.L;
                g.covariance[2][1] = p[2][1] / count - g.mu.b * g.mu.a;
                g.covariance[2][2] = p[2][2] / count - g.mu.b * g.mu.b + Epsilon;
	            */

                // Compute determinant of covariance matrix
                g.determinant = g.covariance[0][0]*(g.covariance[1][1] *g.covariance[2][2] - g.covariance[1][2] *g.covariance[2][1]) 
					          - g.covariance[0][1]*(g.covariance[1][0] *g.covariance[2][2] - g.covariance[1][2] *g.covariance[2][0]) 
							  + g.covariance[0][2]*(g.covariance[1][0] *g.covariance[2][1] - g.covariance[1][1] *g.covariance[2][0]);

                // Compute inverse (cofactor matrix divided by determinant)
                g.inverse[0][0] = (g.covariance[1][1] *g.covariance[2][2] - g.covariance[1][2] *g.covariance[2][1]) / g.determinant;
                g.inverse[1][0] =  - (g.covariance[1][0] *g.covariance[2][2] - g.covariance[1][2] *g.covariance[2][0]) / g.determinant;
                g.inverse[2][0] = (g.covariance[1][0] *g.covariance[2][1] - g.covariance[1][1] *g.covariance[2][0]) / g.determinant;
                g.inverse[0][1] =  - (g.covariance[0][1] *g.covariance[2][2] - g.covariance[0][2] *g.covariance[2][1]) / g.determinant;
                g.inverse[1][1] = (g.covariance[0][0] *g.covariance[2][2] - g.covariance[0][2] *g.covariance[2][0]) / g.determinant;
                g.inverse[2][1] =  - (g.covariance[0][0] *g.covariance[2][1] - g.covariance[0][1] *g.covariance[2][0]) / g.determinant;
                g.inverse[0][2] = (g.covariance[0][1] *g.covariance[1][2] - g.covariance[0][2] *g.covariance[1][1]) / g.determinant;
                g.inverse[1][2] =  - (g.covariance[0][0] *g.covariance[1][2] - g.covariance[0][2] *g.covariance[1][0]) / g.determinant;
                g.inverse[2][2] = (g.covariance[0][0] *g.covariance[1][1] - g.covariance[0][1] *g.covariance[1][0]) / g.determinant;

                // The weight of the gaussian is the fraction of the number of pixels in this Gaussian to the number of 
                // pixels in all the gaussians of this GMM.
                g.pi = (double)count / totalCount;			

                if (computeEigens)
                {
					// Build OpenCV wrappers around our data.
					CvMat mat = cvMat(3, 3, CV_TYPE, g.covariance);
					CvMat eval = cvMat(3, 1, CV_TYPE, g.eigenvalues);
					CvMat evec = cvMat(3, 3, CV_TYPE, g.eigenvectors);
					
					// Compute eigenvalues and vectors using SVD
					cvSVD(&mat, &eval, &evec);   //特征值 特征向量 计算出来是干嘛的...
                }
        }
}

Gaussian * GMM::getGaussians()
{
	return m_gaussians;
}
