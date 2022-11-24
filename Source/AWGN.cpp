#include "Stdafx.h"
#include "AWGN.h"
#include "RandNum.h"

extern CWHRandNum rndGen1;


CAWGN::CAWGN(void)
{
}


CAWGN::~CAWGN(void)
{
}

void CAWGN::Malloc()
{
	m_dim_in_signal = 1;
	m_dim_out_signal = 1;

	return;
}

void CAWGN::Init(double var)
{
	m_var = var;
	m_sigma = sqrt(var);

	return;
}

void CAWGN::Free()
{
	return;
}

void CAWGN::InputOutput(double *xx, double *yy, int channeluse)
{
	int t;

	rndGen1.Normal(yy, channeluse);
	for (t = 0; t < channeluse; t++)
		yy[t] = xx[t] + m_sigma * yy[t];

	return;
}