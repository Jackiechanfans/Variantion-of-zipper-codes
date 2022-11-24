#pragma once
#include "AWGN.h"
#include "BPSKModem.h"
#include "SourceSink.h"
#include "BCHCodec.h" 
#include "GIIBCHCodec.h"
#include "diag_zipper.h"
#include <vector>

class CBCH_BPSK
{
public:
	CBCH_BPSK();
	~CBCH_BPSK();


	CAWGN m_channel;
	CBPSKModem m_modem;
	CSourceSink m_sourcesink;
	CBCHCodec m_codec;
	CFiniteField2 m_GF;
	GIIBCHCodec m_giicodec;
	BSC_channel m_bsc_channel;
	Zipper_diag zipper;


	int** imf_bit;//imformation vector
	int* pos;

	int* r;
	int* L;

	int** m_B;
	int** m_Lambda;
	int*** syndrome;// syndrome
	int** err;//correcting error vector


	void StartSimulator_Initial_for_GII_and_BSC();
	void StartSimulator_Initial_for_zipper();
	void EndSimulator();
	void Simulator();
	void Encoding(int** code_bits);
	void Decoding_for_Zipper(long long int total_chunk);
	void Anchordecoding(int row_id, int i, int y, int last_row, int period_id, int* err_num, int score, vector<vector<vector<int>>>& AD_L, vector<vector<vector<int>>>& AD_epsilon);


};

