#pragma once
class CBPSKModem
{
public:
	CBPSKModem();
	~CBPSKModem();

	int m_num_signal;//星座点数
	int m_dim_signal;//星座维数

	void Malloc();//初始化参数，读取文件数据
	void Free();//释放空间
	void Modulation(int *uu, double *xx, int len_uu);//调制
	void Hard_DeModulation(double *yy, int *uu_hat, int len_yy);//解调，硬判决
	void Soft_DeModulation(double *yy, double var, double *p0_cc, int len_yy);//解调，软判决

private:
	double *m_constellation;//星座
};
