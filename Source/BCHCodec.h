#pragma once

#include "FiniteField2.h"

/*BM: Berlekamp-Massey
GS: Guruswami-Sudan
TL: Tabel-looking*/
enum EBCHBasicDecoderType { BMBCH, GSBCH, TLBCH }; 

class CBCHCodec
{
public:
	CBCHCodec();
	~CBCHCodec();

	CFiniteField2 m_GF;

	int m_len_symbol;//The value of q for GF(2^q)
	int m_root0;//Primitive element
	int m_dist;//Designed distance
	int m_tmin;//Designed correction capacity

	int m_len_shorten;//The length to be shortened for BCH codes
	int m_codelen;//The code length for BCH codes
	int m_codedim;//The dimension for BCH codes
	int m_codechk;//The dimension for BCH codes
	double m_coderate;//The coderate for BCH codes
	
	int *m_syndrome;//The syndrome
	int *m_err;//The error pattern
	int *m_cc_hat;


	//Encoder
	int* m_rx;//用生成多项式编码时，用到的多项式余式
	int* m_generator;

	void Malloc(int code_no, char *file_name);
	void Free();

	//系统编码，信息位在次数高位
	void Encoder(int *uu, int *cc);
	bool Computing_Syndrome(int* rr, int* syndrome);//返回1表示detected_error
	bool Decoder_HIHO(int *zz, int *uu_hat, int *cc_hat, int *err);
	bool Decoder_HIHO(int *rr, int *uu_hat, int *cc_hat);
	bool Decoder_HIHO(int *rr, int *uu_hat);
	bool Decoder_HIHO_syndrome(int *syndrome, int *err);
	bool Decoder_HIHO_EraErr(int *rr, int *uu_hat, int *cc_hat);


	//系统编码，信息位在次数低位，缩短时，对应实际码字 u0c,缩短的是中间部分
	void Encoder_cs(int *uu, int *cc);
	bool Decoder_HIHO_cs(int *rr, int *uu_hat, int *cc_hat, int *err);


	//BMDecoder
	bool BMDecoder(int *rr, int *uu_hat, int *cc_hat, int *err);
	bool BMDecoder_syndrome(int *syndrome, int *err);

	//TLDecoder with t <= 6
	bool TLDecoder(int *rr, int *uu_hat, int *cc_hat, int *err);
	bool TLDecoder_syndrome(int *syndrome, int *err);

protected:
	EBCHBasicDecoderType m_BasicDecoderType;
	void calGenerator(int codelen, int root0, int dist, int*& generator, int& codechk, int eBCH_opt);

private:
	
	
	
	

	//Decoder
	
	char m_BasicDecoderStr[10];
	//for errors-and-erasures decoding
	int *m_rr_1;
	int *m_rr_2;
	int *m_uu_hat_1;
	int *m_uu_hat_2;
	int *m_cc_hat_1;
	int *m_cc_hat_2;
	
	//BMDecoder
	int *m_Lambda;
	int *m_B;
	int *m_T;

	//TLDecoder
	int m_Det[6];
	int m_sigma[6];
	int m_S[6];
};





