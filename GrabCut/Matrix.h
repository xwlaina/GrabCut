// Matrix.h: interface for the CMatrix class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(MATRIX_H)
#define MATRIX_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <math.h>
#include <assert.h>

#define SWAP1(a,b) tempr=(a);(a)=(b);(b)=tempr

class CMatrix
{
public:
	void invert_sqrtm_2D(CMatrix* inv_s) const;
	CMatrix sqrtm2D() const;
	CMatrix invert2D() const;
	void Gram_Schmidt();
	double Delta_fabs(const CMatrix& other) const;
	double Max_fabs();
	void Fill(double value);
	BOOL Eig(double dblU[], double dblV[], CMatrix& mtxEigenVector, double eps = 0.000001);
	
	void thomas(const CMatrix& a, const CMatrix& b, const CMatrix& c, const CMatrix& d);
	void centdiffX(CMatrix& dx) const;
	void centdiffY(CMatrix& dy) const;
	BOOL Copy(int nRowTarget, int nColTarget, int nRowBegin, int nColBegin, int nRowEnd, int nColEnd, const CMatrix *value);
	double SquaredDistance();
	double Dot(CMatrix const &other);
	void Add(const CMatrix& other);
	void Sub(const CMatrix& other);
	void Mul(double value);
	void Div(double value);

	CMatrix();
	CMatrix(int nRows, int nCols);
    CMatrix(int nRows, int nCols, double value[]);
	CMatrix(int nSize);
	CMatrix(int nSize, double value[]);
	CMatrix(const CMatrix& other);
	~CMatrix();
	
	BOOL Zero();
	BOOL Init(int nRows, int nCols);
	BOOL MakeUnitMatrix(int nSize);
	BOOL FromString(CString s, const CString& sDelim /*= " "*/, BOOL bLineBreak /*= TRUE*/);
	CString ToString(const CString& sDelim = _T(" "), BOOL bLineBreak = TRUE, const CString& sFormat = _T("%f")) const;
	void ToString(CString& outS, const CString& sDelim = _T(" "), BOOL bLineBreak = TRUE, const CString& sFormat = _T("%f")) const;
	void SetData(double value[]);

	CString RowToString(int nRow, const CString& sDelim /*= " "*/) const;
	CString ColToString(int nCol, const CString& sDelim /*= " "*/) const;

	double GetElement(int nRow, int nCol) const;

	BOOL AddElement(int nRow, int nCol, double value);
	BOOL SubElement(int nRow, int nCol, double value);

	int	GetNumColumns() const;
	int	GetNumRows() const;

	//////////////////////////////////////////////////////////////////////
	// 获取矩阵的数据
	//
	// 参数：无
	//
	// 返回值：double型指针，指向矩阵各元素的数据缓冲区
	//////////////////////////////////////////////////////////////////////
	double* GetData() const
	{
		return m_pData;
	}

	int GetRowVector(int nRow, double* pVector) const;
	CMatrix GetRowVector(int nRow) const;
	void SetRowVector(int nRow, const CMatrix& other);
	int GetColVector(int nCol, double* pVector) const;
	void SetColVector(int nCol, const CMatrix& other);
	CMatrix& operator=(const CMatrix& other);
	BOOL operator==(const CMatrix& other) const;
	BOOL operator!=(const CMatrix& other) const;
	CMatrix	operator+(const CMatrix& other) const;
	CMatrix	operator-(const CMatrix& other) const;
	CMatrix	operator*(double value) const;
	CMatrix	operator/(double value) const;

	CMatrix	operator*(const CMatrix& other) const;
	CMatrix	operator^(const CMatrix& other) const;


	BOOL CMul(const CMatrix& AR, const CMatrix& AI, const CMatrix& BR, const CMatrix& BI, CMatrix& CR, CMatrix& CI) const;
	CMatrix Transpose() const;
	void Transpose(CMatrix* trans) const;
	BOOL InvertGaussJordan();
	BOOL InvertGaussJordan(CMatrix& mtxImag);
	BOOL InvertSsgj();
	BOOL InvertTrench();
	double DetGauss();
	int RankGauss();
	BOOL DetCholesky(double* dblDet);
	BOOL SplitLU(CMatrix& mtxL, CMatrix& mtxU);
	BOOL SplitQR(CMatrix& mtxQ);
	BOOL SplitUV(CMatrix& mtxU, CMatrix& mtxV, double eps = 0.000001);
	void ppp(double a[], double e[], double s[], double v[], int m, int n);
	void sss(double fg[2], double cs[2]);
	BOOL GInvertUV(CMatrix& mtxAP, CMatrix& mtxU, CMatrix& mtxV, double eps /*= 0.000001*/);
	BOOL MakeSymTri(CMatrix& mtxQ, CMatrix& mtxT, double dblB[], double dblC[]);
	BOOL SymTriEigenv(double dblB[], double dblC[], CMatrix& mtxQ, int nMaxIt /*= 60*/, double eps /*= 0.000001*/);
	void MakeHberg();
	BOOL HBergEigenv(double dblU[], double dblV[], int nMaxIt /*= 60*/, double eps /*= 0.000001*/);
	BOOL JacobiEigenv(double dblEigenValue[], CMatrix& mtxEigenVector, int nMaxIt = 60, double eps = 0.000001);
	BOOL JacobiEigenv2(double dblEigenValue[], CMatrix& mtxEigenVector, double eps = 0.000001);
	
	BOOL SetElement(int nRow, int nCol, double value);
	BOOL SetElements(int nStartRow, int nStartCol, int nRow, int nCol, CMatrix value);
	
private:
	int	m_nNumColumns;
	int	m_nNumRows;
	double*	m_pData;
};


//////////////////////////////////////////////////////////////////////
// 设置指定元素的值
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
// 3. double value - 指定元素的值
//
// 返回值：BOOL 型，说明设置是否成功
//////////////////////////////////////////////////////////////////////
inline BOOL CMatrix::SetElement(int nRow, int nCol, double value)
{
	if (nCol < 0 || nCol >= m_nNumColumns || nRow < 0 || nRow >= m_nNumRows)
		return FALSE;						// array bounds error
	if (m_pData == NULL)
		return FALSE;							// bad pointer error
	
	m_pData[nCol + nRow * m_nNumColumns] = value;

	return TRUE;
}

inline BOOL CMatrix::AddElement(int nRow, int nCol, double value)
{
	if (nCol < 0 || nCol >= m_nNumColumns || nRow < 0 || nRow >= m_nNumRows)
		return FALSE;						// array bounds error
	if (m_pData == NULL)
		return FALSE;							// bad pointer error
	
	m_pData[nCol + nRow * m_nNumColumns] += value;

	return TRUE;
}


inline BOOL CMatrix::SubElement(int nRow, int nCol, double value)
{
	if (nCol < 0 || nCol >= m_nNumColumns || nRow < 0 || nRow >= m_nNumRows)
		return FALSE;						// array bounds error
	if (m_pData == NULL)
		return FALSE;							// bad pointer error
	
	m_pData[nCol + nRow * m_nNumColumns] -= value;

	return TRUE;
}


//////////////////////////////////////////////////////////////////////
// 设置指定元素的值
//
// 参数：
// 1. int nRows - 指定的矩阵行数
// 2. int nCols - 指定的矩阵列数
//
// 返回值：double 型，指定元素的值
//////////////////////////////////////////////////////////////////////
inline double CMatrix::GetElement(int nRow, int nCol) const
{
	assert(nCol >= 0 && nCol < m_nNumColumns && nRow >= 0 && nRow < m_nNumRows); // array bounds error
	assert(m_pData);							// bad pointer error
	return m_pData[nCol + nRow * m_nNumColumns] ;
}

#endif // !defined(AFX_MATRIX_H__00A5FA3D_BEDB_4D66_91A5_0926A8AE045E__INCLUDED_)
