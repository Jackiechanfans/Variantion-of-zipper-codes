/**************************************************************************************

// FiniteField2.h: interface for the CFiniteField2 class.

**************************************************************************************/

#if !defined(AFX_FINITEFIELD2_H__CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_)
#define AFX_FINITEFIELD2_H__CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000




/*************************************************************************
Designed by Xiao Ma (maxiao@mail.sysu.edu.cn), Sun Yat-sen University.
This program can only be employed for academic research.
*************************************************************************/

class CFiniteField2  
{
public:
	CFiniteField2();
	virtual ~CFiniteField2();
	
	int m_len_element;
	int m_num_element;
	int *m_gen_poly;
	int *Exp;
	int *Log;
	int **Vec;

	int Add(int alpha, int beta);
	int Mult(int alpha, int beta);
	int Div(int alpha, int beta);
	int Pow(int alpha, int exp);

	void Malloc(int m);
	void Free();

	void VecAdd(int *alpha, int *beta, int *gamma, int len);
	void VecMult(int *alpha, int *beta, int *gamma, int len);
	void ScalarMult(int alpha, int *beta, int *gamma, int len);

	int PolyEvaluate(int x, int *poly, int deg);
	int PolyMult(int *A, int deg_A, int *B, int deg_B, int *C);
	int PolyAdd(int *A, int deg_A, int *B, int deg_B, int *C);
	int PolyDiv(int *A, int deg_A, int *B, int deg_B, int *C, int *D);
	void LagrangeInterpolation(int k, int *alpha, int *beta, int *poly);

private:
	static const int PrimitivePolynomial[17][6];
};

#endif // !defined(AFX_FINITEFIELD2_H__46CA0220_A3AA_432F_BBA0_5392F4A3E7D2__INCLUDED_)
