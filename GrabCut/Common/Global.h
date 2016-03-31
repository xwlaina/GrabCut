#ifndef GLOBAL_H
#define GLOBAL_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#pragma warning(disable :4786)


#include "../stdafx.h"
#include <math.h>
#include "../Common\Maxflow\graph.h"
#include <opencv2\opencv.hpp> 
//#include "cv.h"
//#include "highgui.h"
#include <vector>
#include <set>
#include <map>
#include "../Vector.h"

using namespace std;

typedef Graph<double,double,double> GraphType;

typedef Vector<double> CVector; 

#define SiNGLE_TENSOR_DIM 3
#define TERMI_CONST 1.0

#define CV_TYPE CV_64FC1

#define REALMIN 2.2251e-308

#define SQRT_OF_2 1.4142135623730950488016887242097


struct	PointEx   
{   
    int x;
    int y; 
	PointEx(): x(0), y(0){}
	PointEx(int _x, int _y): x(_x), y(_y){}
	PointEx& operator = (const PointEx& RightSides)
	{
		x = RightSides.x;
		y = RightSides.y;
		return *this;
	}
	bool operator<(PointEx   const   &ref)const   
	{   
		return (x < ref.x || (x == ref.x && y < ref.y));   
	}   
	bool   operator==(PointEx   const   &ref)const   
	{   
		return   (x == ref.x && y == ref.y);   
	}   
}; 

#define loopi(X) for(int i = 0; i < X; i ++)
#define loopis(X,S) for(int i = 0; i < X; i += S)
#define loopj(X) for(int j = 0; j < X; j ++)
#define loopjs(X,S) for(int j = 0; j < X; j += S)
#define loopk(X) for(int k = 0; k < X; k ++)
#define loopks(X,S) for(int k = 0; k < X; k += S)

#define PI 3.14159265358979323846
#define E 2.71828183

#define NUM_ZERO 1.0e-10 

// User supplied Trimap values
enum TrimapValue
{
        TrimapUnknown, TrimapForeground, TrimapBackground
};

// hard segementation values
enum SegmentationValue
{
        SegmentationForeground, SegmentationBackground
};

enum FlagValue
{
	FlagForeground, FlagBackground, FlagMatting
};


// Storage for N-link weights, each pixel stores links to only four of its 8-neighborhood neighbors.
// This avoids duplication of links, while still allowing for relatively easy lookup.
struct NLinks
{
        double upleft;
        double up;
        double upright;
        double right;
		NLinks()
		{
			upleft = 0.0;
			up = 0.0;
			upright = 0.0;
			right = 0.0;
		}

};


////////////////////////////////////////////////////////////////////////
//将CVector转换成CvMat类型
inline CvMat * maketensortomat(const CVector *tensors)
{  
	int i=0;
	int j=0;
	//得到张量对应向量的维数
	int dim =tensors->size();
	//得到张量的维数
	int n= (unsigned int)(sqrt(2 *dim + 0.25) - 0.5);
	CvMat *mat=cvCreateMat(n,n,CV_32FC1);
	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			if (i>=j)
				cvmSet(mat,i,j,tensors->get((i*(i+1)/2+j)));
			else
				cvmSet(mat,i,j,tensors->get((j*(j+1)/2+i)));
	return  mat;
}

inline void sqrtm2D(double src[2][2], double dst[2][2])
{
	double delta = sqrt(src[0][0] * src[1][1] - src[0][1] * src[1][0]);
	double factor = 1.0 / sqrt(src[0][0] + src[1][1] + 2 * delta);
	dst[0][0] = factor * (delta + src[0][0]);
	dst[0][1] = factor * (src[0][1]);
	dst[1][0] = factor * (src[1][0]);
	dst[1][1] = factor * (delta + src[1][1]);
}

inline void invert2D(double src[2][2], double dst[2][2])
{
	double factor = 1.0 / (src[0][0] * src[1][1] - src[0][1] * src[1][0]);
	dst[0][0] = factor * (src[1][1]);
	dst[0][1] = factor * (- src[0][1]);
	dst[1][0] = factor * (- src[1][0]);
	dst[1][1] = factor * (src[0][0]);
}

inline void interMul2D(const double& x, const double& y, const double& z, const double& a, const double& b, const double& c, double pro[2][2])
{
	pro[0][0] = a*x*x + 2*c*x*z + b*z*z;
	pro[0][1] = a*x*z + c*z*z + c*x*y + b*y*z;
	pro[1][0] = pro[0][1];
	pro[1][1] = a*z*z + 2*c*y*z + b*y*y;
}

inline int GetRandom(int n) {
	int u = rand() * RAND_MAX + rand();
	return ((u % n) + n) % n;
}

inline void diagMul2D(const double& x, const double& y, const double& z, const double& w, const double& a, const double& b, double pro[2][2])
{
	pro[0][0] = a*x*x + b*y*y;
	pro[0][1] = a*x*z + b*y*w;
	pro[1][0] = pro[0][1];
	pro[1][1] = a*z*z + b*w*w;
}

inline void eye2D(double pro[2][2])
{
	pro[0][0] = 1.0;
	pro[0][1] = 0.0;
	pro[1][0] = 0.0;
	pro[1][1] = 1.0;
}

inline BOOL JacobiEigenv2D(double mat[2][2], double eval[2], double evec[2][2], unsigned int nMaxIt = 100, double eps = 0.000001)
{
	unsigned int i, cnt = 1;
	double fm, omega, x, y, cn, sn;
	eye2D(evec);
	while (TRUE)
	{
		if (fabs(mat[1][0]) < eps)
			break;

		if (cnt > nMaxIt)  
			return FALSE;

		cnt++;

		x = - mat[1][0]; 
		y = (mat[0][0] - mat[1][1]) / 2.0;
		omega = x / sqrt(x*x+y*y);
		if (y < 0.0) 
			omega = - omega;

		sn = omega / sqrt(2.0 * (1.0 + sqrt(1.0 - omega*omega)));
		cn = sqrt(1.0 - sn*sn);
		fm = mat[1][1];
		mat[1][1] = fm*cn*cn + (mat[0][0])*sn*sn + (mat[1][0])*omega;
		mat[0][0] = fm*sn*sn + (mat[0][0])*cn*cn - (mat[1][0])*omega;
		mat[1][0] = 0.0; 
		mat[0][1] = 0.0;

		for (i = 0; i < 2; i++)
		{ 
			fm = evec[i][1];
			evec[i][1] = fm*cn + (evec[i][0])*sn;
			evec[i][0] = - fm*sn + (evec[i][0])*cn;
		}
	}

	eval[0] = mat[0][0];
	eval[1] = mat[1][1];
	return TRUE;
}


// Helper function, finds distance between two pixels
inline double distance(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2)
{
        return sqrt((double)((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2)));
}

inline void error_function(char *error_string)
{
	AfxMessageBox((LPCTSTR)error_string);
}

//-------------------------------------------------------------------------

#endif
