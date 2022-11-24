#pragma once
class CAWGN
{
public:
	CAWGN(void);
	~CAWGN(void);

	int m_dim_in_signal;//输入信道维数
	int m_dim_out_signal;//输出信道维数
	double m_var;//噪声方差
	double m_sigma;//噪声标准差


	void Malloc();
	void Init(double var);
	void Free();
	void InputOutput(double *xx, double *yy, int channeluse);

};

