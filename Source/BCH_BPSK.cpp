///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Programmed by Xinwei Zhao (xinweizhao@stu2020.jnu.edu.cn), Jinan University. 
//This program can be employed for academic research only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "Stdafx.h"
#include "BCH_BPSK.h"
#include "RandNum.h"
#include "Utility.h"
#include <map>
#include <algorithm>

extern CWHRandNum rndGen1;



CBCH_BPSK::CBCH_BPSK()
{
}


CBCH_BPSK::~CBCH_BPSK()
{
}


void CBCH_BPSK::StartSimulator_Initial_for_GII_and_BSC()
{
	char file_name[80];
	char temp_str[80];
	FILE* fp;
	sprintf(file_name, "GIISetting.txt");
	printf("\nThe simulation is parameterized by  %s", file_name);
	if ((fp = fopen("GIISetting.txt", "r")) == NULL)
	{
		fprintf(stdout, "Cannot open GIISetting.txt file./n");
		exit(1);
	}
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_giicodec.eBCH_option); 
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &m_bsc_channel.cro_pro);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%lf", &m_bsc_channel.max_err_num);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_giicodec.num_of_dim);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_giicodec.num_of_intl);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_giicodec.m_len_symbol);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &m_giicodec.m_codelen);

	m_giicodec.m_len_shorten = (1 << m_giicodec.m_len_symbol)- 1 - m_giicodec.m_codelen;

	m_giicodec.t = new int[m_giicodec.num_of_intl+1];

	fscanf(fp, "%s", temp_str);
	for (int i = 0; i <= m_giicodec.num_of_intl; i++)
		fscanf(fp, "%d", &m_giicodec.t[i]); //The correcting capacity for different nested interleaves.

	fclose(fp);

	if ((fp = fopen("snrber.txt", "a+")) == NULL) {
		printf("\n Cannot open the snrber.txt file!!!\n");
		exit(1);
	}
	fprintf(fp, "\n\n%%%%The following results correspond to %s\n", file_name);
	fprintf(fp, "%%%%cro_pro =  %lf \n", m_bsc_channel.cro_pro);
	fprintf(fp, "%%%%The maxium number of error  = %d \n", m_bsc_channel.max_err_num);
	fclose(fp);

	m_giicodec.Malloc();		//Generate the generator of different BCH codes

	imf_bit = new int* [m_giicodec.num_of_dim];	
	for (int i = m_giicodec.num_of_intl; i < m_giicodec.num_of_dim; i++)
	{
		imf_bit[i] = new int[m_giicodec.m_codedim[0]];
	}
	for (int i = 0; i < m_giicodec.num_of_intl; i++)
	{
		imf_bit[i] = new int[m_giicodec.m_codedim[m_giicodec.num_of_intl - i]];
	}

	//Initialize some parameters of BM alrorithm
	m_B = new int* [m_giicodec.num_of_dim];		
	m_Lambda = new int* [m_giicodec.num_of_dim];	
	err = new int* [m_giicodec.num_of_dim];			//Error location vector
	r = new int[m_giicodec.num_of_dim];				
	L = new int[m_giicodec.num_of_dim];	

	for (int i = 0; i < m_giicodec.num_of_dim; i++)
	{
		m_B[i] = new int[m_giicodec.dist[m_giicodec.num_of_intl] + 1];
		m_Lambda[i] = new int[m_giicodec.dist[m_giicodec.num_of_intl] + 1];
		err[i] = new int[m_giicodec.m_codelen];
	}

}



void CBCH_BPSK::EndSimulator()
{
	m_giicodec.Free_GII();
	zipper.end_zipper_simulator();

	delete[] r;
	delete[] L;
	for (int i = 0; i < m_giicodec.num_of_dim; i++)
	{
		delete[] m_B[i];
		delete[] m_Lambda[i];
		delete[] err[i];
	}
	delete[] m_B;
	delete[] m_Lambda;
	delete[] err;

	for (int i = 0; i < m_giicodec.num_of_dim; i++)
	{
		delete[] imf_bit[i];

	}
	delete[] imf_bit;


	for (int i = 0; i < zipper.extended_buffer_len; i++)
	{
		for (int j = 0; j < m_giicodec.num_of_dim; j++)
		{
			delete[] syndrome[i][j];
		}
		delete[] syndrome[i];
	}
	delete[] syndrome;


}



void CBCH_BPSK::StartSimulator_Initial_for_zipper()
{
	zipper.n = m_giicodec.m_codelen;
	zipper.start_zipper_simulator(); //Set up zipper codes and anchor decoding

	syndrome = new int** [zipper.extended_buffer_len];	//syndrome
	for (int i = 0; i < zipper.extended_buffer_len; i++)
	{
		syndrome[i] = new int*[m_giicodec.num_of_dim];
		for (int j = 0; j < m_giicodec.num_of_dim; j++)
		{
			syndrome[i][j] = new int[m_giicodec.dist[m_giicodec.num_of_intl] + 1];
			for (int z = 0; z < m_giicodec.dist[m_giicodec.num_of_intl] + 1; z++)
				syndrome[i][j][z] = 0;
		}
	}	

}



void CBCH_BPSK::Simulator()
{
		
	StartSimulator_Initial_for_GII_and_BSC();	////Setup GII codes and BSC

	StartSimulator_Initial_for_zipper();	//Setup zipper codes and anchor decoder
	
	FILE* fp;
	static long long int total_chunk = 0; //Count total receive chunk
	double cro_pro = m_bsc_channel.cro_pro; //The crossover probability
	
		
		m_sourcesink.ClrCnt();
		while (m_sourcesink.m_num_err_bit < m_bsc_channel.max_err_num)
		{
			long long int chunk_row_start = total_chunk * zipper.chunk_size;
			int row_start = chunk_row_start % zipper.extended_buffer_len;
			int period_start = chunk_row_start % zipper.w;

			for (int i = 0; i < zipper.chunk_size; i++)		//Encode zipper codes
			{
				int row_id = (row_start + i) % zipper.extended_buffer_len;	
				int period_id = (period_start + i) % zipper.w;

				for (int y = m_giicodec.num_of_intl; y < m_giicodec.num_of_dim; y++)//Obtain random imformation bits
				{
					m_sourcesink.GetBitStr(imf_bit[y], m_giicodec.m_codedim[0] - zipper.m);	
				}
				for (int y = 0; y < m_giicodec.num_of_intl; y++)
				{
					m_sourcesink.GetBitStr(imf_bit[y], m_giicodec.m_codedim[m_giicodec.num_of_intl - y] - zipper.m);
				}

				for (int z = 0; z < zipper.D; z++)
				{
					int cnt = 0;
					int temp;
					if (zipper.multi_interleaver) // Couple in a cyclic manner
					{
						temp = z - 1;
						if (temp < 0)
							temp = zipper.D - 1;
					}
					else
					{
						temp = z;
					}

					for (int j = zipper.m; j < zipper.n; j++)			//Fill the entries of virtual buffers
					{
						int old_row = (row_id + zipper.map_row[period_id][j] + zipper.extended_buffer_len) % zipper.extended_buffer_len;
						if (z < m_giicodec.num_of_intl)
							imf_bit[z][zipper.m - m_giicodec.check_len[m_giicodec.num_of_intl - z] + cnt] = zipper.buffer_sent[old_row][temp][zipper.map_col[period_id][j]];
						else
							imf_bit[z][zipper.m - m_giicodec.check_len[0] + cnt] = zipper.buffer_sent[old_row][temp][zipper.map_col[period_id][j]];
						cnt++;
					}
				}

				Encoding(zipper.buffer_sent[row_id]); //Enchode for the current row



			}


			for (int i = 0; i < zipper.chunk_size; i++)	//Simulate BSC
			{
				int row_id = (row_start + i) % zipper.extended_buffer_len;
				int period_id = (period_start + i) % zipper.w;

				for (int z = 0; z < zipper.D; z++)
				{
					int temp;
					if (zipper.multi_interleaver)
					{
						temp = z - 1;
						if (temp < 0)
							temp = zipper.D - 1;
					}
					else
					{
						temp = z;
					}

					for (int j = zipper.m; j < zipper.n; j++)
					{
						int old_row = (row_id + zipper.map_row[period_id][j] + zipper.extended_buffer_len) % zipper.extended_buffer_len;
						zipper.buffer_recv[row_id][z][j] = zipper.buffer_recv[old_row][temp][zipper.map_col[period_id][j]];
					}

					m_bsc_channel.BSC_Channel(zipper.buffer_sent[row_id][z], zipper.buffer_recv[row_id][z], 0, zipper.m, cro_pro);

				}

				if (m_giicodec.num_of_dim == m_giicodec.num_of_intl)//Execute MDC-TDZC or zipper code
				{
					for (int z = 0; z < m_giicodec.num_of_dim; z++)//Compute syndrome for the receive chunk
					{
						
						zipper.status[row_id][z] = m_giicodec.zipper_syd(1, m_giicodec.dist[m_giicodec.num_of_dim - z], zipper.n, zipper.buffer_recv[row_id][z], syndrome[row_id][z]);
						if (zipper.status[row_id][z] == 3 && zipper.option_to_running_anchor_decoding == 1) // Set all status as 1 if running program of AD
							zipper.status[row_id][z] = 1;
					}
				}
				else //Execute GII-zipper code
				{
					if (zipper.option_to_running_anchor_decoding == 0) //Compute syndrome for the receive chunk
					{
						for (int z = 0; z < m_giicodec.num_of_dim; z++)
							zipper.status[row_id][z] = m_giicodec.zipper_syd(1, m_giicodec.dist[0], zipper.n, zipper.buffer_recv[row_id][z], syndrome[row_id][z]);
					}
					else
					{
						for (int z = 0; z < zipper.D; z++) // Set all status as 1 if running program of AD
							zipper.status[row_id][z] = 1;
					}

				}


			}

			
			Decoding_for_Zipper(total_chunk);//Decoding

			long long int last_row = total_chunk * zipper.chunk_size + zipper.chunk_size - 1;
			total_chunk++;	//The number of estimated chunk
			long long int oldest_row = last_row - zipper.buffer_len + 1;
			if (oldest_row < 0)
			{
				continue;
			}
			for (int i = 0; i < zipper.D; i++)//Count the number of errors
			{

				for (int j = 0; j < zipper.chunk_size; j++)
				{
					int row_id = (oldest_row + j) % zipper.extended_buffer_len;
					for (int z = 0; z < zipper.m; z++)
					{

						m_sourcesink.m_num_err_bit += zipper.buffer_sent[row_id][i][z] ^ zipper.buffer_recv[row_id][i][z];

					}
					m_sourcesink.m_num_tot_bit += zipper.m;
					m_sourcesink.m_num_tot_blk++;
				}

			}

			if (total_chunk % 50 == 0) {
				m_sourcesink.m_ber = m_sourcesink.m_num_err_bit / m_sourcesink.m_num_tot_bit;
				fprintf(stdout, "\ncro pro = %lf", cro_pro);
				m_sourcesink.PrintResult(stdout);


				if ((fp = fopen("temp_result.txt", "w+")) == NULL) {
					fprintf(stderr, "\n Cannot open the file!!!\n");
					exit(1);
				}
				fprintf(fp, "\ncro pro = %lf", cro_pro);
				m_sourcesink.PrintResult(fp);
				fclose(fp);
			}

		}

		m_sourcesink.m_ber = m_sourcesink.m_num_err_bit / m_sourcesink.m_num_tot_bit;
		printf("\n*** *** *** *** ***");
		m_sourcesink.PrintResult(stdout);
		printf("\n*** *** *** *** ***\n");


		if ((fp = fopen("snrber.txt", "a+")) == NULL) {
			printf("\n Cannot open the file!!!\n");
			exit(1);
		}

		fprintf(fp, "\n%lf %12.10lf", cro_pro, m_sourcesink.m_ber);
		fclose(fp);

	EndSimulator();	
	
}



void CBCH_BPSK::Decoding_for_Zipper(long long int total_chunk)
{
	int* err_num;//The number of errors
	int* deg_Lambda;// The degree of Lambda vector
	int count = 0;
	deg_Lambda = new int[zipper.D];
	err_num = new int[zipper.D];

	long long int last_row = total_chunk * zipper.chunk_size + zipper.chunk_size - 1;
	int row_end = last_row % zipper.extended_buffer_len;
	int period_end = last_row % zipper.w;


	 static vector<vector<vector<int>>> AD_L(zipper.extended_buffer_len, vector<vector<int>>(zipper.D, vector<int>())); // Conflict set for AD
	 static vector<vector<vector<int>>> AD_epsilon(zipper.extended_buffer_len, vector<vector<int>>(zipper.D, vector<int>())); // Correcting set for AD
	 static vector<int> AD_B;
	 for (long long int i = total_chunk * zipper.chunk_size; i <= last_row; i++)
	 {
		 for (int y = 0; y < zipper.D; y++)
		 {
			 int row_id = (i + zipper.extended_buffer_len) % zipper.extended_buffer_len;
			 if (!AD_L[row_id][y].empty())
				 AD_L[row_id][y].clear();
			 if (!AD_epsilon[row_id][y].empty())
				 AD_epsilon[row_id][y].clear();
		 }

	 }
	 
	int iter = 0;
	bool modified = true;
	while (modified && iter < zipper.max_iter)
	{

		modified = false;
		for (int i = 0; i < zipper.buffer_len; i++)//Iterative window decoding
		{
			int row_id = (row_end - i + zipper.extended_buffer_len) % zipper.extended_buffer_len;
			int period_id = (period_end - i % zipper.w + zipper.w) % zipper.w;

			vector<int> I; //The set of indices of sub-codewords declared failure
			map<int,int> I_t; //The set of indices of sub-codewords declared success
			vector<int> IIt;

			for (int y = 0; y < zipper.D; y++) {
				err_num[y] = 0;
			}

			bool luohy = false;		//If they are anchors, which means they don't need to be decoded again
			for (int y = 0; y < zipper.D; y++)
			{
				if (zipper.status[row_id][y] > 0)
				{
					luohy = true;
					break;
				}				
			}
			if (!luohy)
				continue;
		
			for (int y = 0; y < zipper.D; y++)
			{
				int temp, i_temp;

				if (zipper.multi_interleaver)
				{
					temp = y - 1;
					if (temp < 0)
						temp = zipper.D - 1;
					i_temp = (y + 1) % zipper.D;
				}
				else
				{
					temp = y;
					i_temp = y;
				}


				// If there exist sub-codewords with status 1
				if (zipper.status[row_id][y] == 1 || (zipper.status[row_id][y] == 2 && m_giicodec.num_of_dim != m_giicodec.num_of_intl))
				{

					if (m_giicodec.num_of_dim == m_giicodec.num_of_intl)//Execute MDC-TDZC or zipper code
					{

						for (int j = 0; j < m_giicodec.m_codelen; j++)
						{
							err[y][j] = 0;
						}

						zipper.status[row_id][y] = m_giicodec.zipper_syd(1, m_giicodec.dist[m_giicodec.num_of_dim - y], zipper.n, zipper.buffer_recv[row_id][y], syndrome[row_id][y]);

						if (zipper.status[row_id][y] == 1)
						{
							if (m_giicodec.t[0] == 0)
							{
								for (int j = 0; j < m_giicodec.dist[m_giicodec.num_of_intl]; j++) {
									m_Lambda[y][j] = 0;
									m_B[y][j] = 0;
								}
							}
							else
							{
								for (int j = 0; j < m_giicodec.dist[0]; j++) {
									m_Lambda[y][j] = 0;
									m_B[y][j] = 0;
								}
							}

							m_Lambda[y][0] = 1;
							m_B[y][0] = 1;
							r[y] = 0;
							L[y] = 0;

							//BM decoding
							zipper.status[row_id][y] = m_giicodec.Decoder(zipper.buffer_recv[row_id][y], err[y], r[y], L[y], m_B[y], m_Lambda[y], m_giicodec.m_codelen, syndrome[row_id][y], m_giicodec.dist[m_giicodec.num_of_dim - y], err_num[y], deg_Lambda[y]);

							if (m_giicodec.eBCH_option == 1) //If it is extended
							{
								if (zipper.status[row_id][y] == 1)
								{
									int XOR_temp = 0;
									for (int j = 0; j < m_giicodec.m_codelen; j++)
									{
										XOR_temp ^= zipper.buffer_recv[row_id][y][j];
									}
									if (XOR_temp != (deg_Lambda[y] % 2))
									{
										zipper.status[row_id][y] = 2;
									}
								}
							}
						}

						if (zipper.status[row_id][y] == 3)
						{
							zipper.status[row_id][y] = 1;
							err_num[y] = 0;
						}

					}
					else//Execute GII-zipper codes
					{
						for (int j = 0; j < m_giicodec.m_codelen; j++)
						{
							err[y][j] = 0;
						}

						zipper.status[row_id][y] = m_giicodec.zipper_syd(1, m_giicodec.dist[0], zipper.n, zipper.buffer_recv[row_id][y], syndrome[row_id][y]);

						if (zipper.status[row_id][y] == 1 || zipper.status[row_id][y] == 3)
						{
							for (int j = 0; j < m_giicodec.dist[m_giicodec.num_of_intl]; j++) {//BM相关参数初始化
								m_Lambda[y][j] = 0;
								m_B[y][j] = 0;
							}

							m_Lambda[y][0] = 1;
							m_B[y][0] = 1;
							r[y] = 0;
							L[y] = 0;
							//BM decoding
							zipper.status[row_id][y] = m_giicodec.Decoder(zipper.buffer_recv[row_id][y], err[y], r[y], L[y], m_B[y], m_Lambda[y], m_giicodec.m_codelen, syndrome[row_id][y], m_giicodec.dist[0], err_num[y], deg_Lambda[y]);
							if (m_giicodec.eBCH_option == 1)
							{
								if (zipper.status[row_id][y] == 1)
								{
									int XOR_temp = 0;
									for (int j = 0; j < m_giicodec.m_codelen; j++)
									{
										XOR_temp ^= zipper.buffer_recv[row_id][y][j];
									}
									if (XOR_temp != (deg_Lambda[y] % 2))
									{
										zipper.status[row_id][y] = 2;
										for (int j = 0; j < m_giicodec.m_codelen; j++)
										{
											err[y][j] = 0;
										}
									}
								}
							}

							if (zipper.status[row_id][y] == 1)// If decoding is successful, add its index into the set I_t
								I_t.insert(make_pair(y, deg_Lambda[y]));

						}

					}

				}

			}


			if (m_giicodec.num_of_dim == m_giicodec.num_of_intl)//Flipping errors if it is a MDC-TDZC or a zipper code
			{
				for (int y = 0; y < m_giicodec.num_of_dim; y++)
				{

					if (zipper.status[row_id][y] == 1)
					{

						int temp, i_temp;

						AD_B.clear(); //Clear backtracking set

						if (zipper.multi_interleaver)
						{
							temp = y - 1;
							if (temp < 0)
								temp = zipper.D - 1;
							i_temp = (y + 1) % zipper.D;
						}
						else
						{
							temp = y;
							i_temp = y;
						}

						if (err_num[y] != 0)
							pos = new int[err_num[y]];

						count = 0;
						for (int j = 0; j < err_num[y]; j++)
							pos[j] = 0;

						for (int j = 0; j < zipper.n; j++)
						{
							if (err[y][j])
							{
								pos[count] = j;
								count++;
							}
						}



						for (int j = 0; j < err_num[y]; j++)
						{
							if (pos[j] >= zipper.m)
							{
								long long int real_row = last_row - i + zipper.map_row[period_id][pos[j]];
								if (real_row < 0)
								{
									continue;
								}
								int real_row_id = real_row % zipper.extended_buffer_len;
								int real_pos = zipper.map_col[period_id][pos[j]];

								if (zipper.option_to_running_anchor_decoding == 1)//If it is anchor decoder
								{
									if (zipper.status[real_row_id][temp] == 0)		// Conflict with an anchor
									{
										if (AD_L[real_row_id][temp].size() >= zipper.AD_Threshold[temp])
										{
											AD_B.push_back(pos[j]);
										}
										else
										{
											zipper.status[row_id][y] = 3;
											AD_L[row_id][y].push_back(pos[j]);
											AD_L[real_row_id][temp].push_back(real_pos);

										}

									}
								}
								else
								{
									zipper.buffer_recv[real_row_id][temp][real_pos] ^= 0x01;
									zipper.status[real_row_id][temp] = 1;
								}
							}
							else
							{
								long long int virt_row = last_row - i + zipper.i_map_row[period_id][pos[j]];
								if (virt_row <= last_row)
								{
									int virt_row_id = virt_row % zipper.extended_buffer_len;
									int virt_pos = zipper.i_map_col[period_id][pos[j]];

									if (zipper.option_to_running_anchor_decoding == 1)
									{
										if (zipper.status[virt_row_id][i_temp] == 0) // Conflict with an anchor
										{
											if (AD_L[virt_row_id][i_temp].size() >= zipper.AD_Threshold[i_temp])
											{
												AD_B.push_back(pos[j]);
											}
											else
											{
												zipper.status[row_id][y] = 3;
												AD_L[row_id][y].push_back(pos[j]);
												AD_L[virt_row_id][i_temp].push_back(virt_pos);
											}
										}
									}
									else
									{
										zipper.buffer_recv[virt_row_id][i_temp][virt_pos] ^= 0x01;
										zipper.status[virt_row_id][i_temp] = 1;
									}

								}

							}
							if (zipper.option_to_running_anchor_decoding != 1)
								zipper.buffer_recv[row_id][y][pos[j]] ^= 0x01;

						}



						if (zipper.option_to_running_anchor_decoding == 1)
						{
							if (zipper.status[row_id][y] == 1)
							{

								for (int j = 0; j < err_num[y]; j++)    //If (row_id,y).status = 1 then do
								{
									long long int real_row;
									int real_row_id;
									int real_pos;
									int temp_a;
									if (pos[j] >= zipper.m)
									{
										real_row = last_row - i + zipper.map_row[period_id][pos[j]];
										real_row_id = real_row % zipper.extended_buffer_len;
										real_pos = zipper.map_col[period_id][pos[j]];
										if (real_row < 0)
										{
											continue;
										}
										temp_a = temp;
									}
									else
									{
										real_row = last_row - i + zipper.i_map_row[period_id][pos[j]];
										real_row_id = real_row % zipper.extended_buffer_len;
										real_pos = zipper.i_map_col[period_id][pos[j]];
										if (real_row > last_row)
										{
											zipper.buffer_recv[row_id][y][pos[j]] ^= 0x01; //Flip bit
											continue;
										}
										temp_a = i_temp;
									}
									if (!(zipper.status[row_id][y] == 0 && zipper.status[real_row_id][temp_a] == 0))
									{
										int temp_aa = temp_a;

										zipper.buffer_recv[real_row_id][temp_a][real_pos] ^= 0x01;
										zipper.buffer_recv[row_id][y][pos[j]] ^= 0x01;

										if (zipper.status[real_row_id][temp_a] == 2)
										{
											zipper.status[real_row_id][temp_a] = 1;
										}										
										else if (zipper.status[real_row_id][temp_a] == 3)
										{
											zipper.status[real_row_id][temp_a] = 1;

											for (vector<int>::iterator it = AD_L[real_row_id][temp_aa].begin(); it != AD_L[real_row_id][temp_aa].end(); it++)
											{
												int ad_period_id = (real_row % zipper.w + zipper.w) % zipper.w;
												if (*it >= zipper.m)
												{
													long long int ad_real_row = real_row + zipper.map_row[ad_period_id][*it];
													int ad_real_row_id = ad_real_row % zipper.extended_buffer_len;
													int ad_real_pos = zipper.map_col[ad_period_id][*it];
													if (ad_real_row < 0)
													{
														continue;
													}
													temp_a = temp_aa - 1;
													if (temp_a < 0)
														temp_a = zipper.D - 1;
													AD_L[ad_real_row_id][temp_a].erase(find(AD_L[ad_real_row_id][temp_a].begin(), AD_L[ad_real_row_id][temp_a].end(), ad_real_pos));
											
												}
												else
												{

													temp_a = (temp_aa + 1) % zipper.D;
													long long int ad_virt_row = real_row + zipper.i_map_row[ad_period_id][*it];
													int ad_virt_row_id = ad_virt_row % zipper.extended_buffer_len;
													int ad_virt_pos = zipper.i_map_col[ad_period_id][*it];
													if (ad_virt_row > last_row)
													{
														continue;
													}

													AD_L[ad_virt_row_id][temp_a].erase(find(AD_L[ad_virt_row_id][temp_a].begin(), AD_L[ad_virt_row_id][temp_a].end(), ad_virt_pos));

												}
									
											}
											AD_L[real_row_id][temp_aa].clear();
											
										}
									}

								}

								zipper.status[row_id][y] = 0; // Set anchor 
								
								if (!AD_epsilon[row_id][y].empty())
									AD_epsilon[row_id][y].clear();
								for (int jj = 0; jj < err_num[y]; jj++) // Save the decoding result of this anchor
								{
									AD_epsilon[row_id][y].push_back(pos[jj]);
								}

								if (!AD_B.empty())		
								{

									for (vector<int>::iterator iter1 = AD_B.begin(); iter1 != AD_B.end(); iter1++)
									{
										long long int ad_row;
										int ad_row_id, ad_temp, ad_period_id;
										if (*iter1 >= zipper.m)
										{
											ad_row = last_row - i + zipper.map_row[period_id][*iter1];
											ad_temp = temp;
										}
										else
										{
											ad_row = last_row - i + zipper.i_map_row[period_id][*iter1];
											ad_temp = i_temp;
										}
										ad_row_id = ad_row % zipper.extended_buffer_len;
										ad_period_id = (ad_row % zipper.w + zipper.w) % zipper.w;
										for (vector<int>::iterator iter2 = AD_L[ad_row_id][ad_temp].begin(); iter2 != AD_L[ad_row_id][ad_temp].end(); iter2++)
										{
											long long int aad_row;
											int aad_row_id, aad_pos, aad_temp;
											if (*iter2 >= zipper.m)
											{
												aad_row = ad_row + zipper.map_row[ad_period_id][*iter2];
												aad_temp = ad_temp - 1;
												if (aad_temp < 0)
													aad_temp = zipper.D - 1;
												aad_row_id = aad_row % zipper.extended_buffer_len;
												aad_pos = zipper.map_col[ad_period_id][*iter2];
											}
											else
											{
												aad_row = ad_row + zipper.i_map_row[ad_period_id][*iter2];
												aad_temp = (ad_temp + 1) % zipper.D;
												aad_row_id = aad_row % zipper.extended_buffer_len;
												aad_pos = zipper.i_map_col[ad_period_id][*iter2];
											}
											AD_L[aad_row_id][aad_temp].erase(find(AD_L[aad_row_id][aad_temp].begin(), AD_L[aad_row_id][aad_temp].end(), aad_pos));

											if (AD_L[aad_row_id][aad_temp].empty())
											{
												zipper.status[aad_row_id][aad_temp] = 1;
											}

										}
										AD_L[ad_row_id][ad_temp].clear();

										for (int j = 0; j < AD_epsilon[ad_row_id][ad_temp].size(); j++)//operating on errors location of an anchor
										{
											long long int real_row;
											int real_row_id;
											int real_pos;
											int aaad_temp;
											if (AD_epsilon[ad_row_id][ad_temp][j] >= zipper.m)
											{
												real_row = ad_row + zipper.map_row[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
												if (real_row < 0)
												{
													continue;
												}
												real_row_id = real_row % zipper.extended_buffer_len;
												real_pos = zipper.map_col[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
												aaad_temp = ad_temp - 1;
												if (aaad_temp < 0)
													aaad_temp = zipper.D - 1;
											}
											else
											{
												real_row = ad_row + zipper.i_map_row[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
												real_row_id = real_row % zipper.extended_buffer_len;
												real_pos = zipper.i_map_col[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
												if (real_row > last_row)
												{
													zipper.buffer_recv[ad_row_id][ad_temp][AD_epsilon[ad_row_id][ad_temp][j]] ^= 0x01;
													continue;
												}
												aaad_temp = (ad_temp + 1) % zipper.D;
											}
											if (!(zipper.status[ad_row_id][ad_temp] == 0 && zipper.status[real_row_id][aaad_temp] == 0))
											{

												zipper.buffer_recv[real_row_id][aaad_temp][real_pos] ^= 0x01;
												zipper.buffer_recv[ad_row_id][ad_temp][AD_epsilon[ad_row_id][ad_temp][j]] ^= 0x01;
												if (zipper.status[real_row_id][aaad_temp] == 2)
												{
													zipper.status[real_row_id][aaad_temp] = 1;
												}
												else if (zipper.status[real_row_id][aaad_temp] == 3)
												{
													zipper.status[real_row_id][aaad_temp] = 1;

													for (vector<int>::iterator it = AD_L[real_row_id][aaad_temp].begin(); it != AD_L[real_row_id][aaad_temp].end(); it++)
													{
														int aad_period_id = (real_row % zipper.w + zipper.w) % zipper.w;
														int Ad_temp;
														if (*it >= zipper.m)
														{
															long long int ad_real_row = real_row + zipper.map_row[aad_period_id][*it];
															int ad_real_row_id = ad_real_row % zipper.extended_buffer_len;
															int ad_real_pos = zipper.map_col[aad_period_id][*it];
															Ad_temp = aaad_temp - 1;
															if (Ad_temp < 0)
																Ad_temp = zipper.D - 1;
															AD_L[ad_real_row_id][Ad_temp].erase(find(AD_L[ad_real_row_id][Ad_temp].begin(), AD_L[ad_real_row_id][Ad_temp].end(), ad_real_pos));
													
														}
														else
														{
															Ad_temp = (aaad_temp + 1) % zipper.D;

															long long int ad_virt_row = real_row + zipper.i_map_row[aad_period_id][*it];
															int ad_virt_row_id = ad_virt_row % zipper.extended_buffer_len;
															int ad_virt_pos = zipper.i_map_col[aad_period_id][*it];
															AD_L[ad_virt_row_id][Ad_temp].erase(find(AD_L[ad_virt_row_id][Ad_temp].begin(), AD_L[ad_virt_row_id][Ad_temp].end(), ad_virt_pos));

														}
														
													}
													AD_L[real_row_id][aaad_temp].clear();

												}
											}

										}
										AD_epsilon[ad_row_id][ad_temp].clear();
										zipper.status[ad_row_id][ad_temp] = 3; //Frozen this anchor
									}


								}
								
								


							}
						}

						modified = true;
						if (err_num[y]!=0)
							delete[] pos;
						if (zipper.option_to_running_anchor_decoding != 1)
							zipper.status[row_id][y] = 3;
					}

				}
			}


			if (m_giicodec.num_of_dim != m_giicodec.num_of_intl)//Flipping errors if it is a GII-zipper codes
			{
				
				int* poly_nested0;
				int* poly_nested1;
				int* poly_nested2;
				int* nested_syn0;
				int* nested_syn1;
				int* nested_syn2;


				for (int y = 0; y < m_giicodec.num_of_dim; y++)
					if (zipper.status[row_id][y] == 2)
						I.push_back(y);


				if (I.size() > m_giicodec.num_of_intl)
				{
					for (int y = 0; y < m_giicodec.num_of_dim; y++)
					{
						if (zipper.status[row_id][y] == 2)
							zipper.status[row_id][y] = 4;
					}

				}

				bool luohy = false;
				for (int y = 0; y < m_giicodec.num_of_dim; y++) // If there exist sub-codewords with status 3 or 4, stop GII decoding
				{
					if (zipper.status[row_id][y] == 3 || zipper.status[row_id][y] == 4)
					{
						luohy = true;
						break;
					}
				}
				if (luohy)
				{
					for (int y = 0; y < m_giicodec.num_of_dim; y++)
					{
						if (zipper.status[row_id][y] == 2)
							zipper.status[row_id][y] = 4;
					}
					for (map<int, int>::iterator iter = I_t.begin(); iter != I_t.end(); iter++)
					{
						
						int y = iter->first;

						if (err_num[y] != 0)
							pos = new int[err_num[y]];

						count = 0;
						for (int j = 0; j < err_num[y]; j++)
							pos[j] = 0;

						for (int j = 0; j < zipper.n; j++)
						{
							if (err[y][j])
							{
								pos[count] = j;
								count++;
							}
						}

						Anchordecoding(row_id, i, y, last_row, period_id, err_num, -3, AD_L, AD_epsilon);

						modified = true;
						if (err_num[y] != 0)
							delete[] pos;
					}


					continue;
				}


				poly_nested0 = new int[m_giicodec.m_codelen];
				poly_nested1 = new int[m_giicodec.m_codelen + 3];
				poly_nested2 = new int[m_giicodec.m_codelen + 6];
				nested_syn0 = new int[m_giicodec.dist[3]];
				nested_syn1 = new int[m_giicodec.dist[3]];
				nested_syn2 = new int[m_giicodec.dist[3]];

				int round = 0;
				while (round < m_giicodec.num_of_intl) //Execute the second layer decoding
				{
					bool luoli = true;
					int ssstatus = -3 + round;

					if (I.size() == 0)
					{
						int* mixed_syd_0;
						mixed_syd_0 = new int[m_giicodec.dist[m_giicodec.num_of_intl]];

						for (int y = m_giicodec.dist[round]; y < m_giicodec.dist[3]; y++)
							mixed_syd_0[y] = nested_syn0[y];

						for (map<int, int>::iterator iter = I_t.begin(); iter != I_t.end(); iter++)
						{
							int y = iter->first;

							m_giicodec.LFSR_process(m_giicodec.dist[round], m_giicodec.dist[3], syndrome[row_id][y], m_Lambda[y], err_num[y]);
							m_giicodec.Updated_nestedsyn(m_giicodec.dist[round], m_giicodec.dist[3], syndrome[row_id][y], mixed_syd_0, 0, y);

						}

						for (int z = m_giicodec.dist[round]; z < m_giicodec.dist[3]; z++)
						{
							if (mixed_syd_0[z] != 0)
							{
								luoli = true;
								break;
							}
							else
							{
								luoli = false;
							}
						}
						delete[] mixed_syd_0;
					}

					if (!luoli) //If the higher syndromes are zeros, stop the GII decoding
					{
					
						for (map<int, int>::iterator iter = I_t.begin(); iter != I_t.end(); iter++)
						{

							int y = iter->first;

							if (err_num[y] != 0)
								pos = new int[err_num[y]];

							count = 0;
							for (int j = 0; j < err_num[y]; j++)
								pos[j] = 0;

							for (int j = 0; j < zipper.n; j++)
							{
								if (err[y][j])
								{
									pos[count] = j;
									count++;
								}
							}

							Anchordecoding(row_id, i, y, last_row, period_id, err_num, 0, AD_L, AD_epsilon);

							modified = true;
							if (err_num[y] != 0)
								delete[] pos;
						}
						break;
					}

					while (I.size() < m_giicodec.num_of_intl - round && I_t.size() != 0)
					{
						map<int, int>::iterator iter = I_t.begin();
						int max_value = iter->second;
						int dim_err = iter->first;
						for (iter = I_t.begin(); iter != I_t.end(); iter++)
						{
							if (iter->second > max_value)
							{
								max_value = iter->second;
								dim_err = iter->first;
							}
						}
						I.push_back(dim_err);
						I_t.erase(dim_err);
						if(round != 0)
							IIt.erase(find(IIt.begin(), IIt.end(), dim_err));
					}
					sort(I.begin(), I.end());

					if (I_t.size() != 0)
					{
						for (map<int, int>::iterator iter = I_t.begin(); iter != I_t.end(); iter++)
						{

							int y = iter->first;

							if (err_num[y] != 0)
								pos = new int[err_num[y]];

							count = 0;
							for (int j = 0; j < err_num[y]; j++)
								pos[j] = 0;

							for (int j = 0; j < zipper.n; j++)
							{
								if (err[y][j])
								{
									pos[count] = j;
									count++;
								}
							}

							Anchordecoding(row_id, i, y, last_row, period_id, err_num, ssstatus, AD_L, AD_epsilon);

							modified = true;
							if (err_num[y] != 0)
								delete[] pos;
						}
						I_t.clear();
					}

					for (int y = 0; y < m_giicodec.num_of_dim; y++) // If there exists sub-codewords with status 3, stop the GII decoding
					{
						if (zipper.status[row_id][y] == 3)
						{
							for (int y = 0; y < m_giicodec.num_of_dim; y++)
							{
								if (zipper.status[row_id][y] == 2)
									zipper.status[row_id][y] = 4;
							}
							luohy = true;
							break;
						}
					}
					if (luohy)
					{

						for (int y = 0; y < m_giicodec.num_of_dim; y++)
						{
							if (zipper.status[row_id][y] == 1)
							{

								if (err_num[y] != 0)
									pos = new int[err_num[y]];

								count = 0;
								for (int j = 0; j < err_num[y]; j++)
									pos[j] = 0;

								for (int j = 0; j < zipper.n; j++)
								{
									if (err[y][j])
									{
										pos[count] = j;
										count++;
									}
								}

								Anchordecoding(row_id, i, y, last_row, period_id, err_num, ssstatus, AD_L, AD_epsilon);

								modified = true;
								if (err_num[y] != 0)
									delete[] pos;
							}
						}
						break;
					}
					
					if (round == 0) //Compute higher nested syndromes of nested interleaves
					{
						if (I.size() == 3)
						{
							for (int j = 0; j < m_giicodec.m_codelen + 6; j++)
								poly_nested2[j] = 0;

							for (int j = 0; j < m_giicodec.m_codelen + 6; j++)
							{
								if (j < m_giicodec.m_codelen)
									poly_nested2[j] = zipper.buffer_recv[row_id][0][j];
								if (j >= 2 && j < m_giicodec.m_codelen + 2)
									poly_nested2[j] ^= zipper.buffer_recv[row_id][1][j - 2];
								if (j >= 4 && j < m_giicodec.m_codelen + 4)
									poly_nested2[j] ^= zipper.buffer_recv[row_id][2][j - 4];
								if (j >= 6)
									poly_nested2[j] ^= zipper.buffer_recv[row_id][3][j - 6];
							}

							m_giicodec.m_codelen += 6;
							m_giicodec.zipper_syd(m_giicodec.dist[0], m_giicodec.dist[1], m_giicodec.m_codelen, poly_nested2, nested_syn2);
							m_giicodec.m_codelen -= 6;
						}

						if (I.size() == 3 || I.size() == 2)
						{
							for (int j = 0; j < m_giicodec.m_codelen + 3; j++)
								poly_nested1[j] = 0;

							for (int j = 0; j < m_giicodec.m_codelen + 3; j++)
							{
								if (j < m_giicodec.m_codelen)
									poly_nested1[j] = zipper.buffer_recv[row_id][0][j];
								if (j >= 1 && j < m_giicodec.m_codelen + 1)
									poly_nested1[j] ^= zipper.buffer_recv[row_id][1][j - 1];
								if (j >= 2 && j < m_giicodec.m_codelen + 2)
									poly_nested1[j] ^= zipper.buffer_recv[row_id][2][j - 2];
								if (j >= 3)
									poly_nested1[j] ^= zipper.buffer_recv[row_id][3][j - 3];
							}

							m_giicodec.m_codelen += 3;
							m_giicodec.zipper_syd(m_giicodec.dist[0], m_giicodec.dist[2], m_giicodec.m_codelen, poly_nested1, nested_syn1);
							m_giicodec.m_codelen -= 3;
						}

						for (int j = 0; j < m_giicodec.m_codelen; j++)
							poly_nested0[j] = 0;

						for (int j = 0; j < m_giicodec.m_codelen; j++)
							for (int z = 0; z < m_giicodec.num_of_dim; z++)
								poly_nested0[j] ^= zipper.buffer_recv[row_id][z][j];

						m_giicodec.zipper_syd(m_giicodec.dist[0], m_giicodec.dist[3], zipper.n, poly_nested0, nested_syn0);

					}


					if (round != 0) //Update higher nested syndromes of nested interleaves
					{

						if (I.size() == 2)
						{
							for (int y = 0; y < IIt.size(); y++)
							{
								if (zipper.status[row_id][IIt[y]] <= 0)
								{
									m_giicodec.LFSR_process(m_giicodec.dist[round], m_giicodec.dist[2], syndrome[row_id][IIt[y]], m_Lambda[IIt[y]], err_num[IIt[y]]);
									m_giicodec.Updated_nestedsyn(m_giicodec.dist[round], m_giicodec.dist[2], syndrome[row_id][IIt[y]], nested_syn1, 1, IIt[y]);

									m_giicodec.LFSR_process(m_giicodec.dist[round], m_giicodec.dist[3], syndrome[row_id][IIt[y]], m_Lambda[IIt[y]], err_num[IIt[y]]);
									m_giicodec.Updated_nestedsyn(m_giicodec.dist[round], m_giicodec.dist[3], syndrome[row_id][IIt[y]], nested_syn0, 0, IIt[y]);

								}
									
							}
						}

						if (I.size() == 1)
						{
							for (int y = 0; y < IIt.size(); y++)
							{
								if (zipper.status[row_id][IIt[y]] <= 0)
								{

									m_giicodec.LFSR_process(m_giicodec.dist[round], m_giicodec.dist[3], syndrome[row_id][IIt[y]], m_Lambda[IIt[y]], err_num[IIt[y]]);
									m_giicodec.Updated_nestedsyn(m_giicodec.dist[round], m_giicodec.dist[3], syndrome[row_id][IIt[y]], nested_syn0, 0, IIt[y]);

								}

							}
						}

						IIt.clear();

					}

					//Derive higher-order syndromes for sub-codewords in I
					if (I.size() == 3)
					{	
						m_giicodec.Computer_syn3(m_giicodec.dist[0], m_giicodec.dist[1], nested_syn0, nested_syn1, nested_syn2, syndrome[row_id], I);
					}

					if (I.size() == 2)
					{
						m_giicodec.Computer_syn2(m_giicodec.dist[round], m_giicodec.dist[2], nested_syn0, nested_syn1, syndrome[row_id], I);
					}

					if (I.size() == 1)
					{

						for (int j = m_giicodec.dist[round]; j < m_giicodec.dist[3]; j++)
							syndrome[row_id][I[0]][j] = nested_syn0[j];
					}

					int rounder = 3 - I.size() + 1;//r'=v-|I|+1
					for (int z = 0; z < I.size(); z++)
					{
						int y = I[z];
						
						if (zipper.status[row_id][y] == 1)
						{
							for (int j = 0; j < zipper.n; j++)
								err[y][j] = 0;
						}
						//BM decoding
						zipper.status[row_id][y] = m_giicodec.Decoder(zipper.buffer_recv[row_id][y], err[y], r[y], L[y], m_B[y], m_Lambda[y], m_giicodec.m_codelen, syndrome[row_id][y], m_giicodec.dist[rounder], err_num[y], deg_Lambda[y]);
						if (m_giicodec.eBCH_option == 1)
						{
							if (zipper.status[row_id][y] == 1)
							{
								int XOR_temp = 0;
								for (int j = 0; j < m_giicodec.m_codelen; j++)
								{
									XOR_temp ^= zipper.buffer_recv[row_id][y][j];
								}
								if (XOR_temp != (deg_Lambda[y] % 2))
								{
									zipper.status[row_id][y] = 2;
									for (int j = 0; j < m_giicodec.m_codelen; j++)
									{
										err[y][j] = 0;
									}
								}
							}
						}

					}


					for (int y = 0; y < m_giicodec.num_of_dim; y++) //Update each set
					{
						if (zipper.status[row_id][y] == 1)
						{
							I.erase(find(I.begin(), I.end(), y));
							I_t.insert(make_pair(y, deg_Lambda[y]));	
							IIt.push_back(y);
						}
					}


					if (I_t.size() == 0 || (m_giicodec.t[rounder] == m_giicodec.t[3 - I.size() + 1] && I.size()!=0) || rounder == m_giicodec.num_of_intl )
					{
						int sttatus = -3 + rounder;

						for (int y = 0; y < m_giicodec.num_of_dim; y++)
						{
							if (zipper.status[row_id][y] == 2)
								zipper.status[row_id][y] = 4;
						}

						for (map<int, int>::iterator iter = I_t.begin(); iter != I_t.end(); iter++)
						{

							int y = iter->first;

							if (err_num[y] != 0)
								pos = new int[err_num[y]];

							count = 0;
							for (int j = 0; j < err_num[y]; j++)
								pos[j] = 0;

							for (int j = 0; j < zipper.n; j++)
							{
								if (err[y][j])
								{
									pos[count] = j;
									count++;
								}
							}

							Anchordecoding(row_id, i, y, last_row, period_id, err_num, sttatus, AD_L, AD_epsilon);

							modified = true;
							if (err_num[y] != 0)
								delete[] pos;
						}
						I_t.clear();


						break;
					}


					round = rounder; // r = r'


				}


				delete[] poly_nested0;
				delete[] poly_nested1;
				delete[] poly_nested2;
				delete[] nested_syn0;
				delete[] nested_syn1;
				delete[] nested_syn2;


				
			}


		
		}

		iter++;
	}
	delete[] err_num;
	delete[] deg_Lambda;


}


//Modified AD for GII-zipper codes
//score: The reliability score for the decoding result.
//AD_L: The conflict set
//AD_epsilon: The correcting set
void CBCH_BPSK::Anchordecoding(int row_id, int i, int y, int last_row, int period_id, int* err_num, int score, vector<vector<vector<int>>>& AD_L, vector<vector<vector<int>>>& AD_epsilon)
{
	int temp, i_temp;
	vector<int> AD_B; // Backtracking set for AD

	if (zipper.multi_interleaver)
	{
		temp = y - 1;
		if (temp < 0)
			temp = zipper.D - 1;
		i_temp = (y + 1) % zipper.D;
	}
	else
	{
		temp = y;
		i_temp = y;
	}


	for (int j = 0; j < err_num[y]; j++)
	{
		if (pos[j] >= zipper.m) //If it is located at virtual buffer
		{
			long long int real_row = last_row - i + zipper.map_row[period_id][pos[j]];
			if (real_row < 0)
			{
				continue;
			}
			int real_row_id = real_row % zipper.extended_buffer_len;
			int real_pos = zipper.map_col[period_id][pos[j]];

			if (zipper.option_to_running_anchor_decoding == 1)
			{
				if (zipper.status[real_row_id][temp] <= 0)
				{
					int AD_Thld;
					if (zipper.status[real_row_id][temp] == -3)
						AD_Thld = zipper.AD_Threshold[0];
					else if (zipper.status[real_row_id][temp] == -2)
						AD_Thld = zipper.AD_Threshold[1];
					else if (zipper.status[real_row_id][temp] == -1)
						AD_Thld = zipper.AD_Threshold[2];
					else
						AD_Thld = zipper.AD_Threshold[3];

					if (AD_L[real_row_id][temp].size() >= AD_Thld)
					{
						AD_B.push_back(pos[j]);
					}
					else
					{
						zipper.status[row_id][y] = 3;
						AD_L[row_id][y].push_back(pos[j]);
						AD_L[real_row_id][temp].push_back(real_pos);


					}

				}
			}
			else
			{
				zipper.buffer_recv[real_row_id][temp][real_pos] ^= 0x01;
				zipper.status[real_row_id][temp] = 1;
			}
		}
		else //If it is located at real buffer
		{
			long long int virt_row = last_row - i + zipper.i_map_row[period_id][pos[j]];
			if (virt_row <= last_row)
			{
				int virt_row_id = virt_row % zipper.extended_buffer_len;
				int virt_pos = zipper.i_map_col[period_id][pos[j]];

				if (zipper.option_to_running_anchor_decoding == 1)
				{
					if (zipper.status[virt_row_id][i_temp] <= 0)
					{
						int AD_Thld;
						if (zipper.status[virt_row_id][temp] == -3)
							AD_Thld = zipper.AD_Threshold[0];
						else if (zipper.status[virt_row_id][temp] == -2)
							AD_Thld = zipper.AD_Threshold[1];
						else if (zipper.status[virt_row_id][temp] == -1)
							AD_Thld = zipper.AD_Threshold[2];
						else
							AD_Thld = zipper.AD_Threshold[3];

						if (AD_L[virt_row_id][i_temp].size() >= AD_Thld)
						{
							AD_B.push_back(pos[j]);
						}
						else
						{
							zipper.status[row_id][y] = 3;
							AD_L[row_id][y].push_back(pos[j]);
							AD_L[virt_row_id][i_temp].push_back(virt_pos);


						}
					}
				}
				else
				{
					zipper.buffer_recv[virt_row_id][i_temp][virt_pos] ^= 0x01;
					zipper.status[virt_row_id][i_temp] = 1;
				}

			}

		}
		if (zipper.option_to_running_anchor_decoding != 1)
			zipper.buffer_recv[row_id][y][pos[j]] ^= 0x01;

	}



	if (zipper.option_to_running_anchor_decoding == 1)
	{
		if (zipper.status[row_id][y] == 1)
		{

			for (int j = 0; j < err_num[y]; j++)    //If (row_id,y).status = 1 then do:
			{
				long long int real_row;
				int real_row_id;
				int real_pos;
				int temp_a;
				if (pos[j] >= zipper.m)
				{
					real_row = last_row - i + zipper.map_row[period_id][pos[j]];
					real_row_id = real_row % zipper.extended_buffer_len;
					real_pos = zipper.map_col[period_id][pos[j]];
					if (real_row < 0)
					{
						continue;
					}
					temp_a = temp;
				}
				else
				{
					real_row = last_row - i + zipper.i_map_row[period_id][pos[j]];
					real_row_id = real_row % zipper.extended_buffer_len;
					real_pos = zipper.i_map_col[period_id][pos[j]];
					if (real_row > last_row)
					{
						zipper.buffer_recv[row_id][y][pos[j]] ^= 0x01;
						continue;
					}
					temp_a = i_temp;
				}
				if (!(zipper.status[row_id][y] <= 0 && zipper.status[real_row_id][temp_a] <= 0))
				{
					int temp_aa = temp_a;

					zipper.buffer_recv[real_row_id][temp_a][real_pos] ^= 0x01;
					zipper.buffer_recv[row_id][y][pos[j]] ^= 0x01;


					if (zipper.status[real_row_id][temp_a] == 2 || zipper.status[real_row_id][temp_a] == 4)
					{
						zipper.status[real_row_id][temp_a] = 1;
					}
					else if (zipper.status[real_row_id][temp_a] == 3)
					{
						zipper.status[real_row_id][temp_a] = 1;

						for (vector<int>::iterator it = AD_L[real_row_id][temp_aa].begin(); it != AD_L[real_row_id][temp_aa].end(); it++)
						{
							int ad_period_id = (real_row % zipper.w + zipper.w) % zipper.w;
							if (*it >= zipper.m)
							{
								long long int ad_real_row = real_row + zipper.map_row[ad_period_id][*it];
								int ad_real_row_id = ad_real_row % zipper.extended_buffer_len;
								int ad_real_pos = zipper.map_col[ad_period_id][*it];
								if (ad_real_row < 0)
								{
									continue;
								}
								if (zipper.multi_interleaver)
								{
									temp_a = temp_aa - 1;
									if (temp_a < 0)
										temp_a = zipper.D - 1;
								}
								else
								{
									temp_a = temp_aa;
								}

								AD_L[ad_real_row_id][temp_a].erase(find(AD_L[ad_real_row_id][temp_a].begin(), AD_L[ad_real_row_id][temp_a].end(), ad_real_pos));


							}
							else
							{
								if (zipper.multi_interleaver)
									temp_a = (temp_aa + 1) % zipper.D;
								else
									temp_a = temp_aa;
								long long int ad_virt_row = real_row + zipper.i_map_row[ad_period_id][*it];
								int ad_virt_row_id = ad_virt_row % zipper.extended_buffer_len;
								int ad_virt_pos = zipper.i_map_col[ad_period_id][*it];
								if (ad_virt_row > last_row)
								{
									continue;
								}

								AD_L[ad_virt_row_id][temp_a].erase(find(AD_L[ad_virt_row_id][temp_a].begin(), AD_L[ad_virt_row_id][temp_a].end(), ad_virt_pos));



							}

						}
						AD_L[real_row_id][temp_aa].clear();

					}
					for (int j = 0; j < zipper.D; j++)
					{
						if (zipper.status[real_row_id][j] == 4)
							zipper.status[real_row_id][j] = 1;
					}
				}

			}

			zipper.status[row_id][y] = score; //Set it as anchor with status "score"

			if (!AD_epsilon[row_id][y].empty())
				AD_epsilon[row_id][y].clear();
			for (int jj = 0; jj < err_num[y]; jj++)
			{
				AD_epsilon[row_id][y].push_back(pos[jj]);
			}


			if (!AD_B.empty()) // Backtracking anchors
			{

				for (vector<int>::iterator iter1 = AD_B.begin(); iter1 != AD_B.end(); iter1++)
				{
					long long int ad_row;
					int ad_row_id, ad_temp, ad_period_id;
					if (*iter1 >= zipper.m)
					{
						ad_row = last_row - i + zipper.map_row[period_id][*iter1];
						ad_temp = temp;
					}
					else
					{
						ad_row = last_row - i + zipper.i_map_row[period_id][*iter1];
						ad_temp = i_temp;
					}
					ad_row_id = ad_row % zipper.extended_buffer_len;
					ad_period_id = (ad_row % zipper.w + zipper.w) % zipper.w;
					for (vector<int>::iterator iter2 = AD_L[ad_row_id][ad_temp].begin(); iter2 != AD_L[ad_row_id][ad_temp].end(); iter2++)
					{
						long long int aad_row;
						int aad_row_id, aad_pos, aad_temp;
						if (*iter2 >= zipper.m)
						{
							aad_row = ad_row + zipper.map_row[ad_period_id][*iter2];
							if (zipper.multi_interleaver)
							{
								aad_temp = ad_temp - 1;
								if (aad_temp < 0)
									aad_temp = zipper.D - 1;
							}
							else
							{
								aad_temp = ad_temp;
							}

							aad_row_id = aad_row % zipper.extended_buffer_len;
							aad_pos = zipper.map_col[ad_period_id][*iter2];
						}
						else
						{
							aad_row = ad_row + zipper.i_map_row[ad_period_id][*iter2];
							if (zipper.multi_interleaver)
							{
								aad_temp = (ad_temp + 1) % zipper.D;
							}
							else
							{
								aad_temp = ad_temp;
							}

							aad_row_id = aad_row % zipper.extended_buffer_len;
							aad_pos = zipper.i_map_col[ad_period_id][*iter2];
						}



						AD_L[aad_row_id][aad_temp].erase(find(AD_L[aad_row_id][aad_temp].begin(), AD_L[aad_row_id][aad_temp].end(), aad_pos));


						if (AD_L[aad_row_id][aad_temp].empty())
						{
							zipper.status[aad_row_id][aad_temp] = 1;
							for (int j = 0; j < zipper.D; j++)  //The differences between the modified AD and the original AD
							{
								if (zipper.status[aad_row_id][j] == 4)
									zipper.status[aad_row_id][j] = 1;
							}
						}

					}
					AD_L[ad_row_id][ad_temp].clear();

					for (int j = 0; j < AD_epsilon[ad_row_id][ad_temp].size(); j++) 
					{
						long long int real_row;
						int real_row_id;
						int real_pos;
						int aaad_temp;
						if (AD_epsilon[ad_row_id][ad_temp][j] >= zipper.m)
						{
							real_row = ad_row + zipper.map_row[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
							if (real_row < 0)
							{
								continue;
							}
							real_row_id = real_row % zipper.extended_buffer_len;
							real_pos = zipper.map_col[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
							if (zipper.multi_interleaver)
							{
								aaad_temp = ad_temp - 1;
								if (aaad_temp < 0)
									aaad_temp = zipper.D - 1;
							}
							else
							{
								aaad_temp = ad_temp;
							}

						}
						else
						{
							real_row = ad_row + zipper.i_map_row[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
							real_row_id = real_row % zipper.extended_buffer_len;
							real_pos = zipper.i_map_col[ad_period_id][AD_epsilon[ad_row_id][ad_temp][j]];
							if (real_row > last_row)
							{
								zipper.buffer_recv[ad_row_id][ad_temp][AD_epsilon[ad_row_id][ad_temp][j]] ^= 0x01;
								continue;
							}
							if (zipper.multi_interleaver)
								aaad_temp = (ad_temp + 1) % zipper.D;
							else
								aaad_temp = ad_temp;
						}
						if (!(zipper.status[ad_row_id][ad_temp] <= 0 && zipper.status[real_row_id][aaad_temp] <= 0))
						{

							zipper.buffer_recv[real_row_id][aaad_temp][real_pos] ^= 0x01;
							zipper.buffer_recv[ad_row_id][ad_temp][AD_epsilon[ad_row_id][ad_temp][j]] ^= 0x01;


							if (zipper.status[real_row_id][aaad_temp] == 2 || zipper.status[real_row_id][aaad_temp] == 4)
							{
								zipper.status[real_row_id][aaad_temp] = 1;
							}
							else if (zipper.status[real_row_id][aaad_temp] == 3)
							{
								zipper.status[real_row_id][aaad_temp] = 1;

								for (vector<int>::iterator it = AD_L[real_row_id][aaad_temp].begin(); it != AD_L[real_row_id][aaad_temp].end(); it++)
								{
									int aad_period_id = (real_row % zipper.w + zipper.w) % zipper.w;
									int Ad_temp;
									if (*it >= zipper.m)
									{
										long long int ad_real_row = real_row + zipper.map_row[aad_period_id][*it];
										int ad_real_row_id = ad_real_row % zipper.extended_buffer_len;
										int ad_real_pos = zipper.map_col[aad_period_id][*it];
										if (zipper.multi_interleaver)
										{
											Ad_temp = aaad_temp - 1;
											if (Ad_temp < 0)
												Ad_temp = zipper.D - 1;
										}
										else
										{
											Ad_temp = aaad_temp;
										}

										AD_L[ad_real_row_id][Ad_temp].erase(find(AD_L[ad_real_row_id][Ad_temp].begin(), AD_L[ad_real_row_id][Ad_temp].end(), ad_real_pos));

									}
									else
									{
										if (zipper.multi_interleaver)
											Ad_temp = (aaad_temp + 1) % zipper.D;
										else
											Ad_temp = aaad_temp;

										long long int ad_virt_row = real_row + zipper.i_map_row[aad_period_id][*it];
										int ad_virt_row_id = ad_virt_row % zipper.extended_buffer_len;
										int ad_virt_pos = zipper.i_map_col[aad_period_id][*it];
										AD_L[ad_virt_row_id][Ad_temp].erase(find(AD_L[ad_virt_row_id][Ad_temp].begin(), AD_L[ad_virt_row_id][Ad_temp].end(), ad_virt_pos));

									}

								}
								AD_L[real_row_id][aaad_temp].clear();

							}
							for (int j = 0; j < zipper.D; j++)
							{
								if (zipper.status[real_row_id][j] == 4)//The differences between the modified AD and the original AD
									zipper.status[real_row_id][j] = 1;
							}

						}

					}
					AD_epsilon[ad_row_id][ad_temp].clear();

					for (int j = 0; j < zipper.D; j++)
					{
						if (zipper.status[ad_row_id][j] == 4)//The differences between the modified AD and the original AD
							zipper.status[ad_row_id][j] = 1;
					}
					zipper.status[ad_row_id][ad_temp] = 3; // Frozen the anchor


				}


			}


		}

	}
}


void CBCH_BPSK::Encoding(int **code_bits)
{


	for (int i = m_giicodec.num_of_intl; i < m_giicodec.num_of_dim; i++)//Encode for a_v(x),a_(v+1)(x),...,a_(m-1)(x) imformation polynomials using generator of C_0
		m_giicodec.Encoder(imf_bit[i], code_bits[i], (m_giicodec.m_codelen - m_giicodec.m_codedim[0]), m_giicodec.m_codedim[0], m_giicodec.gen[0]);

	if (m_giicodec.num_of_intl == m_giicodec.num_of_dim)//This case is served as MDC-TDZC or zipper code.
	{
		for (int i = 0; i < m_giicodec.num_of_intl; i++)
		{
				m_giicodec.Encoder(imf_bit[i], code_bits[i], m_giicodec.check_len[m_giicodec.num_of_intl - i], m_giicodec.m_codedim[m_giicodec.num_of_intl - i], m_giicodec.gen[m_giicodec.num_of_intl - i]);
		}
			
	}

	if (m_giicodec.num_of_intl == 3 && m_giicodec.num_of_dim == 4) //This case is served as GII(4,3)-TDZC code
	{
		int* mixed_imfbits;
		int* mixed_chkbits;
		int code_chk;
		int* temp_codebits;
		int* temp_codebits2;


		code_chk = m_giicodec.m_codelen - m_giicodec.m_codedim[1];	//The length of check bits for g_1(x);

		mixed_imfbits = new int[m_giicodec.m_codedim[1] + 2];			//Temp for nested encoding
		mixed_chkbits = new int[code_chk];
		temp_codebits = new int[m_giicodec.m_codelen + 2];


		for (int i = 0; i < m_giicodec.m_codedim[1] + 2; i++)	//Initial
			mixed_imfbits[i] = 0;

		temp_codebits[0] = 0;								
		temp_codebits[1] = 0;
		for (int i = 0; i < m_giicodec.m_codelen; i++) // (x^2 + x + 1)*c_3(x)
		{
			temp_codebits[i + 2] = code_bits[3][i];
			temp_codebits[i + 1] ^= code_bits[3][i];
			temp_codebits[i] = code_bits[3][i] ^ temp_codebits[i];
		}

		int* temp_code_bits;
		temp_code_bits = new int[m_giicodec.m_codelen + 2];

		for (int i = code_chk; i < m_giicodec.m_codelen + 2; i++)
			mixed_imfbits[i - code_chk] = temp_codebits[i];
		for (int i = 0; i < m_giicodec.m_codedim[1]; i++)
			mixed_imfbits[i] = mixed_imfbits[i] ^ imf_bit[2][i];

		m_giicodec.m_codelen += 2;

		m_giicodec.Encoder(mixed_imfbits, temp_code_bits, code_chk, m_giicodec.m_codedim[1]+2, m_giicodec.gen[1]);	// for p^(*)(x)

		m_giicodec.m_codelen -= 2;

		for (int i = 0; i < code_chk; i++)		//p(x) = p^(*)(x)+L(f_1(x))
		{
			mixed_chkbits[i] = temp_code_bits[i];
			mixed_chkbits[i] = mixed_chkbits[i] ^ temp_codebits[i];

		}
		for (int i = 0; i < code_chk; i++)		//c_2(x) = [p(x), a_1(x)]
			code_bits[2][i] = mixed_chkbits[i];
		for (int i = code_chk; i < m_giicodec.m_codelen; i++)
			code_bits[2][i] = imf_bit[2][i - code_chk];

		delete[] mixed_imfbits;
		delete[] mixed_chkbits;
		delete[] temp_code_bits;


		code_chk = m_giicodec.m_codelen - m_giicodec.m_codedim[2];	//The length of check bits for g_2(x);


		mixed_imfbits = new int[m_giicodec.m_codedim[2] + 2];			//Temp for nested encoding
		mixed_chkbits = new int[code_chk];
		temp_codebits2 = new int[m_giicodec.m_codelen + 2];

		for (int i = 0; i < m_giicodec.m_codedim[2] + 2; i++)			//Initial
			mixed_imfbits[i] = 0;

		temp_codebits2[0] = 0;								// (x+1)*c_2(x)
		for (int i = 0; i < m_giicodec.m_codelen; i++)
		{
			temp_codebits2[i + 1] = code_bits[2][i];
			temp_codebits2[i] = code_bits[2][i] ^ temp_codebits2[i];
		}

		temp_code_bits = new int[m_giicodec.m_codelen + 2];

		for (int i = code_chk; i < m_giicodec.m_codelen + 2; i++)
			mixed_imfbits[i - code_chk] = temp_codebits[i];

		for (int i = code_chk; i < m_giicodec.m_codelen + 1; i++)
			mixed_imfbits[i - code_chk] ^= temp_codebits2[i];

		for (int i = 0; i < m_giicodec.m_codedim[2]; i++)
			mixed_imfbits[i] = mixed_imfbits[i] ^ imf_bit[1][i];

		m_giicodec.m_codelen += 2;

		m_giicodec.Encoder(mixed_imfbits, temp_code_bits, code_chk, m_giicodec.m_codedim[2]+2, m_giicodec.gen[2]);	// for p^(*)(x)

		m_giicodec.m_codelen -= 2;

		for (int i = 0; i < code_chk; i++)		//p(x) = p^(*)(x)+L(f_2(x))
		{
			mixed_chkbits[i] = temp_code_bits[i];
			mixed_chkbits[i] = mixed_chkbits[i] ^ temp_codebits[i];
			mixed_chkbits[i] ^= temp_codebits2[i];

		}

		for (int i = 0; i < code_chk; i++)		//c_1(x) = [p(x), a_1(x)]
			code_bits[1][i] = mixed_chkbits[i];
		for (int i = code_chk; i < m_giicodec.m_codelen; i++)
			code_bits[1][i] = imf_bit[1][i - code_chk];


		delete[] mixed_imfbits;
		delete[] mixed_chkbits;
		delete[] temp_code_bits;
		delete[] temp_codebits;
		delete[] temp_codebits2;


		code_chk = m_giicodec.m_codelen - m_giicodec.m_codedim[3];	//The length of check bits for g_3(x);

		mixed_imfbits = new int[m_giicodec.m_codedim[3]];			//Temp for nested encoding
		mixed_chkbits = new int[code_chk];

		for (int i = 0; i < m_giicodec.m_codedim[3]; i++)			//Initial
			mixed_imfbits[i] = 0;


		for (int i = code_chk; i < m_giicodec.m_codelen; i++)		//Truncate f_0(x) = c_1(x) + c_2(x) + ...
			for (int j = 1; j < m_giicodec.num_of_dim; j++)
				mixed_imfbits[i - code_chk] = code_bits[j][i] ^ mixed_imfbits[i - code_chk];

		for (int i = 0; i < m_giicodec.m_codedim[3]; i++)			
			mixed_imfbits[i] = mixed_imfbits[i] ^ imf_bit[0][i];

		m_giicodec.Encoder(mixed_imfbits, code_bits[0], code_chk, m_giicodec.m_codedim[3], m_giicodec.gen[3]);	// for p^(*)(x)

		for (int i = 0; i < code_chk; i++)		//p(x) = p^(*)(x)+L(f_0(x))
		{
			mixed_chkbits[i] = code_bits[0][i];
			for (int j = 1; j < m_giicodec.num_of_dim; j++)
			{
				mixed_chkbits[i] = mixed_chkbits[i] ^ code_bits[j][i];
			}
		}
		for (int i = 0; i < code_chk; i++)		//c_0(x) = [p(x), a_0(x)]
			code_bits[0][i] = mixed_chkbits[i];
		for (int i = code_chk; i < m_giicodec.m_codelen; i++)
			code_bits[0][i] = imf_bit[0][i - code_chk];

		delete[] mixed_imfbits;
		delete[] mixed_chkbits;

	}

}




