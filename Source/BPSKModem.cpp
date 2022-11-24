#include "stdafx.h"
#include "BPSKModem.h"



CBPSKModem::CBPSKModem()
{
}

CBPSKModem::~CBPSKModem()
{
}

void CBPSKModem::Malloc() {

	m_num_signal = 2;
	m_dim_signal = 1;
	m_constellation = new double[m_num_signal];

	m_constellation[0] = +1;
	m_constellation[1] = -1;

	return;
}

void CBPSKModem::Free() {
	delete[]m_constellation;
	return;
}
void CBPSKModem::Modulation(int *uu, double *xx, int len_uu) {
	int i;
	for (i = 0; i < len_uu; ++i)
		xx[i] = m_constellation[uu[i]];
	return;
}
void CBPSKModem::Hard_DeModulation(double *yy, int *uu_hat, int len_yy) {
	int i;
	for (i = 0; i < len_yy; ++i) {
		if (yy[i] > 0)
			uu_hat[i] = 0;
		else
			uu_hat[i] = 1;
	}
	return;

}
void CBPSKModem::Soft_DeModulation(double *yy, double var, double *p0_cc, int len_yy) {
	int i;
	for (i = 0; i < len_yy; ++i) {//prob(+1|y)
		if (yy[i] > 0)
			p0_cc[i] = 1.0 / (1.0 + exp(-2.0 * yy[i] / var));
		else
			p0_cc[i] = 1.0 - 1.0 / (1.0 + exp(2.0 * yy[i] / var));
	}
}
