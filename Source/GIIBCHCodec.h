#include "BCHCodec.h"
#include <vector>

using namespace std;

class GIIBCHCodec: public CBCHCodec
{
public:
	GIIBCHCodec();
	~GIIBCHCodec();

	int num_of_dim; //The number of interleave for GII.
	int num_of_intl; //The number of nested interleave for GII.
	int* t; //The error correcting capability for different nested interleave.
	
	void Malloc();
	void Free_GII();
	void Encoder(int* uu, int* cc, int m_codechk, int m_codedim, int* m_generator);
	unsigned char Decoder(int* rr, int* err, int& r, int& L, int* m_B, int* m_Lambda, int& m_codelen, int* m_syndrome, int m_dist, int& L_1, int& deg_Lambda);
	void Computer_syd(int m_dist_0, int m_dist, int* syd0, int* syd1, int* I);
	int zipper_syd(int m_dist_0, int m_dist, int m_codelen, int* rr, int* m_syndrome);
	void LFSR_process(int m_dist_0, int m_dist, int* m_syndrome, int* m_Lambda, int deg_Lambda);
	void Updated_nestedsyn(int m_dist_0, int m_dist, int* m_syndrome, int* nested_syn, int l_val, int dim);
	void Computer_syn3(int m_dist_0, int m_dist, int* syd0, int* syd1, int* syd2, int** syn, vector<int> &I);
	void Computer_syn2(int m_dist_0, int m_dist, int* syd0, int* syd1, int** syn, vector<int> &I);
	int power(int i, int index);


	int* m_codedim;	//The length of different imformation vectors;
	int** gen; //The generators;
	int* dist;	//The distance for different generators;
	int* check_len;	//The length of different check vectors;
	int eBCH_option; //Extend BCH codes by one bit.	
	int** syndrome; //The syndromes.

};

class BSC_channel
{
public:
	BSC_channel();
	~BSC_channel();
	double cro_pro;//crossover_probability for BSC
	double max_err_num;//maximum_error_number

	void BSC_Channel(int* cc, int* rr, int start_len, int len, double pro);
};
