// Filter.cpp: implementation of the CFilter class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "Filter.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

CFilter::CFilter()
{   //m_channels�൱��m_tensor
	m_channels = NULL;
	// m_levels * SiNGLE_TENSOR_DIM;
	m_nochannels = 0;
	m_w = 0;
	m_h = 0;
}

CFilter::~CFilter()
{
	m_channels = NULL;
}

//����ʽ2-17��AOS��TVflow���ǹ����ݶȵĺ���
//cv_gradΪ�Ѿ��õ����ݶ�ֵͼ��
//�õ�g
void CFilter::TVflow(CvMat *cv_grad)
{
	double esp, TVPower;
	esp = 0.001;
	//dTensorTVPower�൱�ڹ�ʽ�е�P
	TVPower = 0.6;
	cvPow(cv_grad, cv_grad, TVPower);
	cvAddS(cv_grad, cvScalar(esp), cv_grad);
	cvDiv(NULL, cv_grad, cv_grad, 1.0);
}

//����һ�ֹ���g�ĺ���
void CFilter::AOS(CvMat *cv_grad)
{
	double esp;
	esp = 0.001;
	cvAddS(cv_grad, cvScalar(esp), cv_grad);
	cvConvertScale(cv_grad, cv_grad, 1.0 / 15.0);
	cvPow(cv_grad, cv_grad, 16);
	cvDiv(NULL, cv_grad, cv_grad, -4.2293);
	cvExp(cv_grad, cv_grad);
	cvSubRS(cv_grad, cvScalar(1), cv_grad);
}

//��������������������õ�g
void CFilter::computeDiffusivity(IntermediateData_Diffusivity &data, double sigma, int option)
{
	data.diffusivity->Zero();
	unsigned int x,y;
	//���ڲ�ͬ�ĳ߶Ⱥ������Ľṹ���д���
	for (unsigned int i=0;i<m_nochannels;i++)
	{
		if (sigma > 0)
		{
			for (y = 0;y < m_h;y++)
			{
				for (x = 0;x < m_w;x++)
				{
					float* dst = &CV_IMAGE_ELEM( data.cv_channel, float, y, x );
					dst[0] = (float)((m_channels[i])->GetElement(y,x));
				}
			}
			//GAUSSIAN�˲�
			cvSmooth(data.cv_channel, data.cv_channel, CV_GAUSSIAN, 0, 0, sigma);
			for (y = 0;y < m_h;y++)
			{
				for (x = 0;x < m_w;x++)
				{//��GAUSSIAN�˲�������ݱ�����smoothed_channel��ָ��ľ�����
					float* dst = &CV_IMAGE_ELEM( data.cv_channel, float, y, x );
					data.smoothed_channel->SetElement(y, x, (double)(dst[0]));
				}
			}
			//���˲�������ݽ���΢��
			//�ֱ���X,Y������ݶ�
			data.smoothed_channel->centdiffX(*(data.dx));
			data.smoothed_channel->centdiffY(*(data.dy));
		} 
		else
		{//���sigma<0,�Ͳ������˲���ֱ�Ӽ���X,Y������ݶ�
			(m_channels[i])->centdiffX(*(data.dx));
			(m_channels[i])->centdiffY(*(data.dy));
		}

		//��ʱֱ�ӵ��ýṹ���е�cv_dx��data.cv_dy��Ӧ�ô���һ����ֵ���ɣ�����
		//�����ݶ�
		//2*2�ṹ�����ĵ�һ��Ԫ��
		cvMul(data.cv_dx, data.cv_dx, data.cv_dx);
		//2*2�ṹ�����ĵ��ĸ�Ԫ��
		cvMul(data.cv_dy, data.cv_dy, data.cv_dy);
		//cv_grad��������ݶȵ�ƽ��
		cvAdd(data.cv_dx, data.cv_dy, data.cv_dx);
		cvAdd(data.cv_dx, data.cv_grad, data.cv_grad);
	}

	switch(option)
	{
	case 0:
		//�õ�g����ʽ(2_17)
		TVflow(data.cv_grad);
		break;
	case 1:
		AOS(data.cv_grad);
		break;
	case 2:
		break;
	default:
		break;
	}
}

//�Ѿ��õ�����ɢ�ṹg
void CFilter::AOS_scheme(IntermediateData_AOS &dada, double stepsize)
{
	unsigned int x,y,pos,pos_1;
	double *a, *b, *d;
	//Operating on Rows
	//dָ�����diffusivity�Ļ�����,diffusivityΪ��ɢ����
	//���õ�a=I-2*taoA
	//a,b,d��Ϊ���������е���ʽ��ϳ�һ��������
	d = dada.diffusivity->GetData();
	a = dada.a_row->GetData();
	b = dada.b_row->GetData();
	for (x = 0; x < m_w; x++)
	{
		for (y = 0; y < (m_h - 1); y++)
		{//�õ����ڴ洢�����λ��pos(x,y)
			pos = y * m_w + x;
			//pos_1(x,y+1)
			pos_1 = pos + m_w;
			//b
			b[pos] = d[pos] + d[pos_1];
		}
		dada.a_row->SetElement(0, x, dada.b_row->GetElement(0, x));
		dada.a_row->SetElement(m_h - 1, x, dada.b_row->GetElement(m_h - 2, x));
		//�õ�a���������ʽ
		for (y = 1; y < (m_h - 1); y++)
		{
			pos = y * m_w + x;
			pos_1 = pos - m_w;
			a[pos] = b[pos_1] + b[pos];
		}
		//b[pos]*(- stepsize)��Ϊ��ʲô
		for (y = 0; y < (m_h - 1); y++)
		{
			pos = y * m_w + x;
			b[pos] = b[pos] * (- stepsize);
		}
		//�õ�(2-32)�ұ߱��ʽ�ľ���
		for (y = 0; y < m_h; y++)
		{
			pos = y * m_w + x;
			a[pos] = 1 + stepsize * a[pos];
		}
	}

	//Operating on Columns
	a = dada.a_column->GetData();
	b = dada.b_column->GetData();
	for (y = 0; y < m_h; y++)
	{
		for (x = 0; x < (m_w - 1); x++)
		{
			pos = y * m_w + x;
			pos_1 = pos + 1;
			b[x * m_h + y] = d[pos] + d[pos_1];
		}
		dada.a_column->SetElement(0, y, dada.b_column->GetElement(0, y));
		dada.a_column->SetElement(m_w - 1, y, dada.b_column->GetElement(m_w - 2, y));
		for (x = 1; x < (m_w - 1); x++)
		{
			pos = x * m_h + y;
			pos_1 = pos - m_h;
			a[pos] = b[pos_1] + b[pos];
		}

		for (x = 0; x < (m_w - 1); x++)
		{
			pos = x * m_h + y;
			b[pos] = b[pos] * (- stepsize);
		}

		for (x = 0; x < m_w; x++)
		{
			pos = x * m_h + y;
			a[pos] = 1 + stepsize * a[pos];
		}
	}
	//�в�����Ϊ�����dada.a_row��dada.b_row
	CMatrix channel_trans(m_w, m_h);
	//thomas�㷨thomas(a,b,c,d)
	for (unsigned int i=0;i<m_nochannels;i++)
	{	 //�õ�(2-32)x����Ľ�
		dada.y_row->thomas(*(dada.a_row), *(dada.b_row), *(dada.b_row), *(m_channels[i]));	
		(m_channels[i])->Transpose(&channel_trans);
		//�õ�(2-32)y����Ľ�
		dada.y_column->thomas(*(dada.a_column), *(dada.b_column), *(dada.b_column), channel_trans);
		double *y_r, *y_c, *ch;
		y_r = dada.y_row->GetData();
		y_c = dada.y_column->GetData();
		ch = (m_channels[i])->GetData();
		for (y = 0;y < m_h; y++)
		{
			for (x = 0; x < m_w; x++)
			{
				pos = y * m_w + x;
				pos_1 = x * m_h + y;
				//y_r,��ʾʽ��2_32��ʽ�ұߵĵ�һ������Ӧ�Ľ�
				//y_c,��ʾʽ��2_32��ʽ�ұߵĵڶ�������Ӧ�Ľ�
				ch[pos] = 0.5 * (y_r[pos] + y_c[pos_1]);
			}
		}
	}
}


//
void CFilter::Diff(CMatrix **channels, unsigned int nochannels)
{   //�õ����������������
	m_channels = channels;
	m_nochannels = nochannels;
	m_w = (channels[0])->GetNumColumns();
	m_h = (channels[0])->GetNumRows();
	int nTensorNosteps = 2;
	int dTensorSigma = 0;
	if (nTensorNosteps > 0)
	{
		nlDiff(channels, nochannels);
	} 
	else if (dTensorSigma > 0)
	{
		gaussDiff(channels, nochannels);
	}

}


//��˹�˲�
void CFilter::gaussDiff(CMatrix **channels, unsigned int nochannels)
{
	double sigma = 0;
	if (sigma > 0)
	{
		IplImage *cv_channel = cvCreateImage( cvSize(m_w,m_h), IPL_DEPTH_32F, 1 );
		unsigned int x,y;
		for (unsigned int i=0;i<m_nochannels;i++)
		{
			for (y = 0;y < m_h;y++)
			{
				for (x = 0;x < m_w;x++)
				{
					float* dst = &CV_IMAGE_ELEM( cv_channel, float, y, x );
					dst[0] = (float)((m_channels[i])->GetElement(y,x));
				}
			}
			//��˹�˲�
			cvSmooth(cv_channel, cv_channel, CV_GAUSSIAN, 0, 0, sigma);
			for (y = 0;y < m_h;y++)
			{
				for (x = 0;x < m_w;x++)
				{
					float* dst = &CV_IMAGE_ELEM( cv_channel, float, y, x );
					(m_channels[i])->SetElement(y, x, (double)(dst[0]));
				}
			}
		}
		cvReleaseImage(&cv_channel);
	}
}

void CFilter::nlDiff(CMatrix **channels, unsigned int nochannels, int option)
{
	double sigma; //�ȶ�ͼ����и�˹ƽ������׼��
	double stepsize;
	unsigned int nosteps;
	sigma = 0;
	stepsize = 5000;
	nosteps = 2;
	//������ʼ��
	CMatrix *diffusivity = new CMatrix(m_h,m_w);

	IntermediateData_Diffusivity data_diff;
	data_diff.diffusivity = diffusivity;
	data_diff.dx = new CMatrix(m_h,m_w);
	data_diff.dy = new CMatrix(m_h,m_w);
	if (sigma > 0)
	{
		data_diff.smoothed_channel = new CMatrix(m_h,m_w);
		data_diff.cv_channel = cvCreateImage( cvSize(m_w,m_h), IPL_DEPTH_32F, 1 );
	} 
	else
	{
		data_diff.smoothed_channel = NULL;
		data_diff.cv_channel = NULL;
	}

	CvMat cv_dx = cvMat(m_h, m_w, CV_64FC1, data_diff.dx->GetData());
	CvMat cv_dy = cvMat(m_h, m_w, CV_64FC1, data_diff.dy->GetData());
	CvMat cv_grad = cvMat(m_h, m_w, CV_64FC1, data_diff.diffusivity->GetData());
	data_diff.cv_dx = &cv_dx;
	data_diff.cv_dy = &cv_dy;
	data_diff.cv_grad = &cv_grad;
	//һϵ�г�ʼ��
	IntermediateData_AOS data_aos;
	data_aos.diffusivity = diffusivity;
	data_aos.a_row = new CMatrix(m_h, m_w);
	data_aos.b_row = new CMatrix(m_h - 1, m_w);
	data_aos.a_column = new CMatrix(m_w, m_h);
	data_aos.b_column = new CMatrix(m_w - 1, m_h);
	data_aos.y_row = new CMatrix(m_h, m_w);
	data_aos.y_column = new CMatrix(m_w, m_h);


	for (unsigned int t = 0; t < nosteps; t++)
	{ //nosteps��������
		computeDiffusivity(data_diff, sigma, option);
		AOS_scheme(data_aos, stepsize);
	}


	delete diffusivity;

	delete data_diff.dx;
	delete data_diff.dy;

	if (data_diff.smoothed_channel != NULL)
	{
		delete data_diff.smoothed_channel;
		data_diff.smoothed_channel = NULL;
	}
	if (data_diff.cv_channel != NULL)
	{
		cvReleaseImage(&data_diff.cv_channel);
		data_diff.cv_channel = NULL;
	}

	delete data_aos.a_row;
	delete data_aos.b_row;
	delete data_aos.a_column;
	delete data_aos.b_column;
	delete data_aos.y_row;
	delete data_aos.y_column;
}
