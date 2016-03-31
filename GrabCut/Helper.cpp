#include "stdafx.h"
#include "Helper.h"
#include <assert.h>

//#include <cv.h>
//#include <cxcore.h>
//#include <highgui.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

CHelper::CHelper(void)
{
}


CHelper::~CHelper(void)
{
}

double **CHelper::dmatrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;
	m = (double **) calloc((unsigned) (nrh-nrl+1), sizeof(double*));
	//if (!m) MessageBox(_T("allocation failure 1 in dmatrix()"));
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i] = (double *) calloc((unsigned) (nch-ncl+1), sizeof(double));
		//if (!m[i]) MessageBox(_T("allocation failure 2 in dmatrix()"));
		m[i] -= ncl;
	}
	return m;
}

void CHelper::free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void CHelper::four1(double *data, int nn, int isign)
{
	int n, mmax, m, j, istep, i;
	double wtemp, wr, wpr, wpi, wi, theta;
	double tempr, tempi;

	n = nn << 1;
	j = 1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP1(data[j],data[i]);
			SWAP1(data[j+1],data[i+1]);
		}
		m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax = 2;
	while (n > mmax) {
		istep = 2*mmax;
		theta = 6.28318530717959/(isign*mmax);
		wtemp = sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j = i+mmax;
				tempr = wr*data[j]-wi*data[j+1];
				tempi = wr*data[j+1]+wi*data[j];
				data[j] = data[i]-tempr;
				data[j+1] = data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp=wr)*wpr-wi*wpi+wr;
			wi = wi*wpr+wtemp*wpi+wi;
		}
		mmax = istep;
	}
}
void CHelper::four2(double **fftr, double **ffti, double **rdata, double **idata, int rs, int cs, int isign)
{
	double **T, *tmp1, *tmp2;
	int i, j;

	tmp1 = dvector(1,2*cs);
	tmp2 = dvector(1,2*rs);
	T = dmatrix(1,2*rs,1,cs);

	for (i=1;i<=rs;i++) {
		for (j=1;j<=cs;j++) {
			tmp1[j*2-1] = rdata[i][j];
			tmp1[j*2] = idata[i][j];
		}
		four1(tmp1, cs, isign);
		for (j=1;j<=cs;j++) {
			T[i*2-1][j] = tmp1[j*2-1];
			T[i*2][j] = tmp1[j*2];
		}
	}

	for (i=1;i<=cs;i++) {
		for (j=1;j<=rs;j++) {
			tmp2[j*2-1] = T[j*2-1][i];
			tmp2[j*2] = T[j*2][i];
		}
		four1(tmp2,rs,isign);
		for (j=1;j<=rs;j++) {
			fftr[j][i] = tmp2[j*2-1];
			ffti[j][i] = tmp2[j*2];
		}
	}
	free_dvector(tmp1, 1, 2*cs);
	free_dvector(tmp2, 1, 2*rs);
	free_dmatrix(T, 1, 2*rs, 1, cs); 
}
double *CHelper::dvector(int nl, int nh)
{
	double *v;

	v = (double *) calloc((unsigned) (nh-nl+1), sizeof(double));
	//if (!v) AfxMessageBox(_T("allocation failure in dvector()"));
	return v-nl;
}

void CHelper::free_dvector(double *v, int nl, int nh)
{
	free((char*) (v+nl));
}

void CHelper::RGB2Lab(const Image < Color_RGB > * _rgb, Image < Color_Lab > * _lab)
{
	unsigned int x,y;
	unsigned int width, height;
	width = _rgb->width();
	height = _rgb->height();
	assert (width == _lab->width() && height == _lab->height());
	IplImage *cv_image = cvCreateImage( cvSize(width,height), IPL_DEPTH_32F, 3 );
	for (y = 0;y < height;y++)
	{
		for (x = 0;x < width;x++)
		{
			float* dst = &CV_IMAGE_ELEM( cv_image, float, y, x*3 );
			dst[0] = (float)(((*_rgb)(x,y)).b);
			dst[1] = (float)(((*_rgb)(x,y)).g);
			dst[2] = (float)(((*_rgb)(x,y)).r);
		}
	}
	cvCvtColor( cv_image, cv_image, CV_BGR2Lab );   //Labģʽ�������ɫ����࣬������߼��豸�޹�
	for (y = 0;y < height;y++)
	{
		for (x = 0;x < width;x++)
		{
			float* dst = &CV_IMAGE_ELEM( cv_image, float, y, x*3 );
			((*_lab)(x,y)).L = (double)(dst[0]);
			((*_lab)(x,y)).a = (double)(dst[1]);
			((*_lab)(x,y)).b = (double)(dst[2]);
		}
	}
	cvReleaseImage(&cv_image);
}

/*************************************************************************
 *
 * �������ƣ�
 *   FFT()
 *
 * ����:
 *   std::complex<double> * TD	- ָ��ʱ�������ָ��
 *   std::complex<double> * FD	- ָ��Ƶ�������ָ��
 *   r						- 2��������������������2��r�η�=N
 *
 * ˵��:
 *   �ú�������ʵ�ֿ��ٸ���Ҷ�任��
 *
 ************************************************************************/
void CHelper::FFT(std::complex<double> * TD, std::complex<double> * FD, int r)
{
	// ����Ҷ�任����
	long	count;
	
	// ѭ������
	int		i,j,k;
	
	// �м����
	int		bfsize,p;
	
	// �Ƕ�
	double	angle;
	
	std::complex<double> *W,*X1,*X2,*X;
	
	// ���㸶��Ҷ�任����
	count = 1 << r;     //1����rλ����2��r�η�
	
	// ������������洢��
	W  = new std::complex<double>[count / 2];
	X1 = new std::complex<double>[count];
	X2 = new std::complex<double>[count];
	
	// �����Ȩϵ��
	for(i = 0; i < count / 2; i++)
	{
		angle = -i * PI * 2 / count;
		W[i] = std::complex<double> (cos(angle), sin(angle));
	}
	
	// ��ʱ���д��X1
	memcpy(X1, TD, sizeof(std::complex<double>) * count);
	
	// ���õ����㷨���п��ٸ���Ҷ�任
	for(k = 0; k < r; k++)  //FFT�ļ���
	{
		for(j = 0; j < 1 << k; j++)   //����
		{
			bfsize = 1 << (r-k);  //�����Ŀ��*2
			for(i = 0; i < bfsize / 2; i++)
			{
				p = j * bfsize;
				X2[i + p] = X1[i + p] + X1[i + p +bfsize/2];
				X2[i + p + bfsize / 2] = (X1[i + p] - X1[i + p + bfsize / 2]) * W[i * (1<<k)];
			}
		}
		X  = X1;
		X1 = X2;
		X2 = X;
	}
	
	// ��������(�������λ�����)
	for(j = 0; j < count; j++)
	{
		p = 0;
		for(i = 0; i < r; i++)
		{
			if (j&(1<<i))
			{
				p+=1<<(r-i-1);
			}
		}
		FD[j]=X1[p];
	}
	
	// �ͷ��ڴ�
	delete W;
	delete X1;
	delete X2;
}


/*************************************************************************
 *
 * �������ƣ�
 *   IFFT()
 *
 * ����:
 *   std::complex<double> * FD	- ָ��Ƶ��ֵ��ָ��
 *   std::complex<double> * TD	- ָ��ʱ��ֵ��ָ��
 *   r						- 2��������������������2��r�η�=N
 *
 * ˵��:
 *   �ú�������ʵ�ֿ��ٸ���Ҷ���任��
 *
 ************************************************************************/
void CHelper::IFFT(std::complex<double> * FD, std::complex<double> * TD, int r)
{
	// ����Ҷ�任����
	long	count;
	
	// ѭ������
	int		i;
	
	std::complex<double> *X;
	
	// ���㸶��Ҷ�任����
	count = 1 << r;   //1����rλ����2��r�η�
	
	// ������������洢��
	X = new std::complex<double>[count];
	
	// ��Ƶ���д��X
	memcpy(X, FD, sizeof(std::complex<double>) * count);

	// ����
	for(i = 0; i < count; i++)
	{
		X[i] = std::complex<double> (X[i].real(), -X[i].imag());
	}
	
	// ���ÿ��ٸ���Ҷ�任
	FFT(X, TD, r);
	
	// ��ʱ���Ĺ���
	for(i = 0; i < count; i++)
	{
		TD[i] = std::complex<double> (TD[i].real() / count, -TD[i].imag() / count);
	}
	
	// �ͷ��ڴ�
	delete X;
}

/*************************************************************************
 *
 * �������ƣ�
 *   CONV()
 *
 * ����:
 *   std::complex<double> * TD1	    - ָ��ʱ����������1��ָ��
 *   std::complex<double> * TD2	    - ָ��ʱ����������2��ָ��
 *   std::complex<double> * TDout	- ָ��ʱ���������ָ��
 *   M						    - ����1�ĳ���
 *   N                          - ����2�ĳ���
 *
 * ˵��:
 *   �ú�������FFTʵ�ֿ��پ����
 *
 ************************************************************************/
void CHelper::Convolution(std::complex<double> * TD1, std::complex<double> * TD2, std::complex<double> * TDout, int M, int N)
{
	// ����������
	int	count = M+N-1;

	// ����ʹ��FFT����count��չΪ2����
	int Lcount;
    int r=0;  // 2����������FFT����������2��r�η�=Lcount

	int temp;
	if (log((float)count)/log((float)2)-int(log((float)count)/log((float)2))==0)
      temp = log((float)count)/log((float)2);
	else
	  temp = log((float)count)/log((float)2)+1;
	r = temp;
	Lcount = 1<<r;	

	// ������������洢��
    std::complex<double> *X1, *X2, *FD1, *FD2, *FD12, *TD12;

	X1 = new std::complex<double>[Lcount];  //����������1
	X2 = new std::complex<double>[Lcount];  //����������2
	FD1 = new std::complex<double>[Lcount];   //����1�ĸ���Ҷ�任���
	FD2 = new std::complex<double>[Lcount];   //����2�ĸ���Ҷ�任���
	FD12 = new std::complex<double>[Lcount];   //����1,2��Ƶ����˽��
	TD12 = new std::complex<double>[Lcount];   //����1,2�ĸ���Ҷ���任���
	
    //�����в���ΪLcount����
	std::complex<double> *X, *Y;
	X = new std::complex<double>[M];  //��ʱ�洢��
	Y = new std::complex<double>[N];
	
	// ��ʱ���д��X,Y
	memcpy(X, TD1, sizeof(std::complex<double>) * M);
	memcpy(Y, TD2, sizeof(std::complex<double>) * N);

	// ѭ������
	int	i;

    for (i=0; i<M; i++)    //��������1����   
    {
        X1[i] = std::complex<double>(X[i].real(), X[i].imag());                                             
    }

    for (i=M; i<Lcount; i++)    //����1��0
    {
        X1[i] = std::complex<double>(0, 0);
    }

	for (i=0; i<N; i++)    //��������2����   
    {
        X2[i] = std::complex<double>(Y[i].real(), Y[i].imag());
    }

    for (i=N; i<Lcount; i++)    //����2��0
    {
        X2[i] = std::complex<double>(0, 0);
    }

    // �ͷ��ڴ�
	delete X;
	delete Y;

    //����1��FFT
	FFT(X1, FD1, r);

	//����2��FFT
    FFT(X2, FD2, r);

    //����1,2��Ƶ�����
    for (i=0; i<Lcount; i++)    //����1,2���
    {
        FD12[i] = std::complex<double>(FD1[i].real()*FD2[i].real()-FD1[i].imag()*FD2[i].imag(), FD1[i].real()*FD2[i].imag()+FD1[i].imag()*FD2[i].real());
    }

	//����1,2��Ƶ����˵�IFFT
    IFFT(FD12, TD12, r);

	//TD12�е�ǰM+N-1��Ϊ����������д��TDout
    memcpy(TDout, TD12, sizeof(std::complex<double>)*count);
	
	// �ͷ��ڴ�
	delete X1;
	delete X2;
	delete FD1;
	delete FD2;
	delete FD12;
    delete TD12;
}

//TDout��TD1�ĳ���һ����ȡM+N-1���м�M��
void CHelper::Convolution(double * TD1, double * TD2, double * TDout, int M, int N)
{
	std::complex<double>* c_TD1 = new std::complex<double>[M];
	std::complex<double>* c_TD2 = new std::complex<double>[N];
	std::complex<double>* c_TDout = new std::complex<double>[M+N-1];
	int i;
	for (i=0;i<M;i++)
	{
		c_TD1[i] = std::complex<double>(TD1[i], 0);
	}
	for (i=0;i<N;i++)
	{
		c_TD2[i] = std::complex<double>(TD2[i], 0);
	}
	Convolution(c_TD1, c_TD2, c_TDout, M, N);
	int pos = cvCeil( (N - 1) / 2.0 );
	for (i=0;i<M;i++)
	{
		TDout[i] = c_TDout[pos++].real();
	}
	delete c_TD1;
	delete c_TD2;
	delete c_TDout;
}

void CHelper::FFT2(CMatrix *Output_real, CMatrix *Output_imag, const CMatrix *Input_real, const CMatrix *Input_imag)
{
	int xs, ys, i, j;
	double **R, **I, **Fr, **Fi;

	xs = Input_real->GetNumRows();
	ys = Input_real->GetNumColumns();

	R  = dmatrix(1,xs,1,ys);
	I  = dmatrix(1,xs,1,ys);
	Fr = dmatrix(1,xs,1,ys);
	Fi = dmatrix(1,xs,1,ys);

	for (i=1;i<=Input_real->GetNumRows();i++) 
		for (j=1;j<=Input_real->GetNumColumns();j++) 
		{
			R[i][j] = Input_real->GetElement(i-1, j-1);
			I[i][j] = Input_imag->GetElement(i-1, j-1);
		}

		four2(Fr, Fi, R, I, xs, ys, 1);         /* 2-D FFT */

		for (i=1;i<=Input_real->GetNumRows();i++) 
			for (j=1;j<=Input_real->GetNumColumns();j++) 
			{
				Output_real->SetElement(i-1, j-1, Fr[i][j]);
				Output_imag->SetElement(i-1, j-1, Fi[i][j]);
			}

			free_dmatrix(R,1,xs,1,ys);
			free_dmatrix(I,1,xs,1,ys);   
			free_dmatrix(Fr,1,xs,1,ys);
			free_dmatrix(Fi,1,xs,1,ys);  
}


void CHelper::IFFT2(CMatrix *Output_real, CMatrix *Output_imag, const CMatrix *Input_real, const CMatrix *Input_imag)
{
	int xs, ys, i, j;
	double **R, **I, **Fr, **Fi, NN;

	xs = Input_real->GetNumRows();
	ys = Input_real->GetNumColumns();

	R  = dmatrix(1,xs,1,ys);
	I  = dmatrix(1,xs,1,ys);
	Fr = dmatrix(1,xs,1,ys);
	Fi = dmatrix(1,xs,1,ys);

	for (i=1;i<=Input_real->GetNumRows();i++) 
		for (j=1;j<=Input_real->GetNumColumns();j++) 
		{
			R[i][j] = Input_real->GetElement(i-1, j-1);
			I[i][j] = Input_imag->GetElement(i-1, j-1);
		}

		four2(Fr, Fi, R, I, xs, ys, -1);         /* 2-D IFFT */

		NN = (double) (xs*ys);

		for (i=1;i<=Input_real->GetNumRows();i++) 
			for (j=1;j<=Input_real->GetNumColumns();j++) 
			{
				Output_real->SetElement(i-1, j-1, Fr[i][j]/NN);
				Output_imag->SetElement(i-1, j-1, Fi[i][j]/NN);
			}

			free_dmatrix(R,1,xs,1,ys);
			free_dmatrix(I,1,xs,1,ys);   
			free_dmatrix(Fr,1,xs,1,ys);
			free_dmatrix(Fi,1,xs,1,ys); 
}

void CHelper::GetMagnitude(CMatrix *Output_mag, const CMatrix *Input_real, const CMatrix *Input_imag)
{
	// ���ȼ���������Ƿ����
	ASSERT (Input_real->GetNumColumns() == Input_imag->GetNumColumns() && Output_mag->GetNumColumns() == Input_imag->GetNumColumns() && Input_real->GetNumRows() == Input_imag->GetNumRows() && Output_mag->GetNumRows() == Input_imag->GetNumRows());

	for (int i = 0 ; i < Output_mag->GetNumRows() ; ++i)
	{
		for (int j = 0 ; j <  Output_mag->GetNumColumns(); ++j)
			Output_mag->SetElement(i, j, sqrt(((Input_real->GetElement(i,j))*(Input_real->GetElement(i,j)))+((Input_imag->GetElement(i,j))*(Input_imag->GetElement(i,j)))));
	}
}

