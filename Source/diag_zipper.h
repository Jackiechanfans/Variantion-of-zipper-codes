#include <vector>
using namespace std;

class Zipper_diag 
{

public:

	Zipper_diag();
	~Zipper_diag();

	int n; //Code length
	int m;	//The width of virtual buffer
	int w;//The size of tile
	int D;//The number of zipper codes
	int L;//L = m / w

	int multi_interleaver;// The flag of running MDC-TDZC ( multi-chain coupled TDZC )
	int option_to_running_anchor_decoding;//The flag of running anchor decoding (AD)

	int* AD_Threshold; // Threshold of AD
	int max_lookback; //The maimmum lookback for the encoding of zipper codes
	int buffer_len;//The decoding window size
	int extended_buffer_len;//The extended decoding window size
	int chunk_size;//Chunk size i.e, mu in the paper
	int max_iter;//The maximum iteration
	int*** buffer_sent;//The send buffer
	int*** buffer_recv;//The receive buffer

	//The interleaver map
	//Assume that there is a bit (i,j) located at the real buffer and the interleaver map is defined as phi, then
	//(map_row, map_col) = phi^(-1)(i, j).
	//If there is a bit (i,j) located at the virtual buffer, then
	//(i_map_row, i_map_col) = phi(i, j).
	int** map_row;
	int** map_col;
	int** i_map_row;
	int** i_map_col;
	int* t;//correcting capacity

	int ** status;//Status


	void start_zipper_simulator();
	void end_zipper_simulator();
};


