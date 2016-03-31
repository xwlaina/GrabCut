// Gabor.cpp: implementation of the Gabor class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Gabor.h"


#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Gabor::Gabor(const IplImage *cv_image)
{
	m_scale = 2;
	m_orientation = 3;
	m_dim = m_scale * m_orientation;
	m_Ul = 0.1;
	m_Uh = 0.4;
	m_flag = 1;
	m_side = 60;
	m_alpha = pow((m_Uh/m_Ul), 1.0/(double)(m_scale-1));

	IplImage *cv_gray_image = cvCreateImage( cvGetSize(cv_image), cv_image->depth, 1 );
	cvCvtColor( cv_image, cv_gray_image, CV_BGR2GRAY );
	m_height = cv_gray_image->height;
	m_width = cv_gray_image->width;
	m_image = new CMatrix(m_height, m_width);
	int x,y;
	for (y = 0; y < m_height; y++)
	{
		for (x = 0; x < m_width; x++)
		{
			uchar* dst = &CV_IMAGE_ELEM( cv_gray_image, uchar, y, x );
			m_image->SetElement(y, x, (double)(dst[0]));
		}
	}
	cvReleaseImage(&cv_gray_image);
	m_filters = new Image < CVector* > (m_width, m_height);
	for (y = 0; y < m_height; y++)
	{
		for (x = 0; x < m_width; x++)
		{
			(*m_filters)(x,y) = new CVector(Dim());
		}
	}
	//	m_fPdf = new double[m_K];
	//	m_bPdf = new double[m_K];
}

Gabor::~Gabor()
{
	int x,y;
	if (m_image != NULL)
	{
		delete m_image;
		m_image = NULL;
	}
	if (m_filters != NULL)
	{
		for (y = 0; y < m_height; y++)
		{
			for (x = 0; x < m_width; x++)
			{
				if ((*m_filters)(x,y) != NULL)
				{
					delete (*m_filters)(x,y);
					(*m_filters)(x,y) = NULL;
				}
			}
		}
		delete m_filters;
		m_filters = NULL;
	}
}

//The GaborFilteredImg provides the outputs of the Gabor filter bank
void Gabor::GaborFilteredImg()
{
	int h, w, xs, ys, border, r1, r2, r3, r4, s, n;
	unsigned int i, j;
	CMatrix *IMG, *IMG_imag, *Gr, *Gi, *Tmp_1, *Tmp_2, *F_1, *F_2, *G_real, *G_imag, *F_real, *F_imag;
	CMatrix *m_output_R;//Gabor变换的实部
	CMatrix *m_output_I;//Gabor变换的虚部
	CMatrix *gabor_filters;//Gabor变换的幅值


	border = m_side;

	/* FFT2 */
	xs = (int) pow(2.0, ceil(log2((double)(m_height+2.0*border))));
	ys = (int) pow(2.0, ceil(log2((double)(m_width+2.0*border))));

	IMG = new CMatrix(xs, ys);

	r1 = m_width+border;
	r2 = m_width+2*border;
	for (h=0;h<border;h++) {
		for (w=0;w<border;w++)
			IMG->SetElement(h,w,m_image->GetElement(border-1-h,border-1-w));
		for (w=border;w<r1;w++)
			IMG->SetElement(h,w,m_image->GetElement(border-1-h,w-border));
		for (w=r1;w<r2;w++)
			IMG->SetElement(h,w,m_image->GetElement(border-1-h,2*m_width-w+border-1));
	}

	r1 = m_height+border;
	r2 = m_width+border;
	r3 = m_width+2*border;
	for (h=border;h<r1;h++) {
		for (w=0;w<border;w++)
			IMG->SetElement(h,w,m_image->GetElement(h-border,border-1-w));
		for (w=border;w<r2;w++)
			IMG->SetElement(h,w,m_image->GetElement(h-border,w-border));
		for (w=r2;w<r3;w++)
			IMG->SetElement(h,w,m_image->GetElement(h-border,2*m_width-w+border-1));
	}

	r1 = m_height+border;
	r2 = m_height+2*border;
	r3 = m_width+border;
	r4 = m_width+2*border;
	for (h=r1;h<r2;h++) {
		for (w=0;w<border;w++)
			IMG->SetElement(h,w,m_image->GetElement(2*m_height-h+border-1,border-1-w));
		for (w=border;w<r3;w++)
			IMG->SetElement(h,w,m_image->GetElement(2*m_height-h+border-1,w-border));
		for (w=r3;w<r4;w++)
			IMG->SetElement(h,w,m_image->GetElement(2*m_height-h+border-1,2*m_width-w+border-1));
	}

	F_real = new CMatrix(xs,ys);
	F_imag = new CMatrix(xs,ys);
	IMG_imag = new CMatrix(xs,ys);

	m_helper.FFT2(F_real, F_imag, IMG, IMG_imag);

	/* ----------- compute the Gabor filtered output ------------- */
	Gr = new CMatrix(2*m_side+1, 2*m_side+1);
	Gi = new CMatrix(2*m_side+1, 2*m_side+1);
	Tmp_1 = new CMatrix(xs,ys);
	Tmp_2 = new CMatrix(xs,ys);
	F_1 = new CMatrix(xs,ys);
	F_2 = new CMatrix(xs,ys);
	G_real = new CMatrix(xs,ys);
	G_imag = new CMatrix(xs,ys);
	m_output_R = new CMatrix((m_height)*m_scale, (m_width)*m_orientation);
	m_output_I = new CMatrix((m_height)*m_scale, (m_width)*m_orientation);
	gabor_filters = new CMatrix((m_height)*m_scale, (m_width)*m_orientation);

	for (s=0;s<m_scale;s++) 
	{
		for (n=0;n<m_orientation;n++) 
		{
			GetGaborFilter(Gr, Gi, s+1, n+1);
			F_1->Copy(0, 0, 0, 0, 2*m_side, 2*m_side, Gr);
			F_2->Copy(0, 0, 0, 0, 2*m_side, 2*m_side, Gi);
			m_helper.FFT2(G_real, G_imag, F_1, F_2);

			(*Tmp_1) = (*G_real) ^ (*F_real);
			(*Tmp_2) = (*G_imag) ^ (*F_imag);
			(*IMG) = (*Tmp_1) - (*Tmp_2);

			(*Tmp_1) = (*G_real) ^ (*F_imag);
			(*Tmp_2) = (*G_imag) ^ (*F_real);
			(*IMG_imag) = (*Tmp_1) + (*Tmp_2);

			m_helper.IFFT2(Tmp_1, Tmp_2, IMG, IMG_imag);

			m_output_R->Copy(s*m_height, n*m_width, 2*m_side, 2*m_side, m_height+2*m_side-1, m_width+2*m_side-1, Tmp_1);
			m_output_I->Copy(s*m_height, n*m_width, 2*m_side, 2*m_side, m_height+2*m_side-1, m_width+2*m_side-1, Tmp_2);
		}
	}
	m_helper.GetMagnitude(gabor_filters, m_output_R, m_output_I);
	for ( i=0; i<m_height; i++)
	{
		for ( j=0; j<m_width; j++)
		{
			for (s=0;s<m_scale;s++) 
			{
				for (n=0;n<m_orientation;n++)
				{
					((*m_filters)(j,i))->set(n+s*m_orientation, gabor_filters->GetElement(s*m_height + i,n*m_width + j));
				}
			}
		}
	}
	delete Gr;
	delete Gi;
	delete Tmp_1;
	delete Tmp_2;
	delete F_1;
	delete F_2;
	delete G_real;
	delete G_imag;
	delete F_real;
	delete F_imag;
	delete IMG;
	delete IMG_imag;
	delete m_output_R;
	delete m_output_I;
	delete gabor_filters;

}

/* ------------------------------------------------------------------------------------------------------
The function generates a Gabor filter with the selected index 's' and 'n' (scale and orientation, 
respectively) from a Gabor filter bank. The returned filter is stored in 'Gr' (real part) and 'Gi' (imaginary part).
--------------------------------------------------------------------------------------------------------*/
void Gabor::GetGaborFilter(CMatrix *Gr, CMatrix *Gi, int s, int n)
{
	double u0, z, Uvar, Vvar, Xvar, Yvar, X, Y, G, t1, t2, m;
	int x, y;

	u0 = m_Uh/pow(m_alpha, (double) m_scale - s);

	Uvar = (m_alpha-1.0)*u0/((m_alpha+1.0)*sqrt(2.0*log(2.0)));

	z = -2.0*log(2.0)*(Uvar*Uvar)/u0;
	Vvar = tan(PI/(2*m_orientation))*(u0+z)/sqrt(2.0*log(2.0)-z*z/(Uvar*Uvar));

	Xvar = 1.0/(2.0*PI*Uvar);
	Yvar = 1.0/(2.0*PI*Vvar);

	t1 = cos(PI/m_orientation*(n-1.0));
	t2 = sin(PI/m_orientation*(n-1.0));

	for (x=0;x<2*m_side+1;x++) {
		for (y=0;y<2*m_side+1;y++) {
			X = (double) (x-m_side)*t1+ (double) (y-m_side)*t2;
			Y = (double) -(x-m_side)*t2+ (double) (y-m_side)*t1;
			G = 1.0/(2.0*PI*Xvar*Yvar)*pow(m_alpha, (double) m_scale-s)*exp(-0.5*((X*X)/(Xvar*Xvar)+(Y*Y)/(Yvar*Yvar)));
			Gr->SetElement(x,y,(G*cos(2.0*PI*u0*X)));
			Gi->SetElement(x,y,(G*sin(2.0*PI*u0*X)));
		}
	}

	/* if flag = 1, then remove the DC from the filter */
	if (m_flag == 1) {
		m = 0;
		for (x=0;x<2*m_side+1;x++)
			for (y=0;y<2*m_side+1;y++)
				m += Gr->GetElement(x,y);

		m /= pow((double) 2.0*m_side+1, 2.0);

		for (x=0;x<2*m_side+1;x++)
			for (y=0;y<2*m_side+1;y++)
				Gr->SetElement(x,y,((Gr->GetElement(x,y)) - m));
	}	
}


//void Gabor::debug()
//{
//}


double Gabor::computeDistance2(unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2) const
{
	CVector *vector_1 = (*m_filters)(x1,y1);
	CVector *vector_2 = (*m_filters)(x2,y2);
	return distance2_Euclidean(vector_1, vector_2);
}


Image < CVector* >* Gabor::GetFilters() const
{
	return m_filters;
}