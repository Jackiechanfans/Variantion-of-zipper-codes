#include "diag_zipper.h"
#include "Stdafx.h"


Zipper_diag::Zipper_diag()
{
}

Zipper_diag::~Zipper_diag()
{
}

void Zipper_diag::end_zipper_simulator()
{
	for (int i = 0; i < extended_buffer_len; i++)
		delete[] status[i];
	delete[] status;
	for (int i = 0; i < extended_buffer_len; i++)
	{
		for (int j = 0; j < D; j++)
		{
			delete[] buffer_sent[i][j];
			delete[] buffer_recv[i][j];
		}
		delete[] buffer_sent[i];
		delete[] buffer_recv[i];
	}
	delete[] buffer_sent;
	delete[] buffer_recv;
	delete[] t;
	delete[] AD_Threshold;
}

void Zipper_diag::start_zipper_simulator()
{
	char file_name[80];
	char temp_str[80];

	int delayed = 0; // Delay  zipper codes



	FILE* fp;
	sprintf(file_name, "Zippersetting.txt");
	printf("\nThe simulation is parameterized by  %s", file_name);
	if ((fp = fopen("Zippersetting.txt", "r")) == NULL)
	{
		fprintf(stdout, "Cannot open Zippersetting.txt file./n");
		exit(1);
	}
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &option_to_running_anchor_decoding);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &multi_interleaver);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &delayed);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &D);
	fscanf(fp, "%s", temp_str);
	AD_Threshold = new int[D];
	for (int i = 0; i < D; i++)
		fscanf(fp, "%d", &AD_Threshold[i]);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &w);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &chunk_size);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &buffer_len);
	fscanf(fp, "%s", temp_str);
	fscanf(fp, "%d", &max_iter);		

	fclose(fp);

	m = n / 2;
	L = m / w;
	max_lookback = m + delayed * w;
	extended_buffer_len = buffer_len + 2 * (chunk_size + max_lookback);
	
	status = new int* [extended_buffer_len];
	buffer_sent = new int** [extended_buffer_len];
	buffer_recv = new int** [extended_buffer_len];

	for (int i = 0; i < extended_buffer_len; i++)
	{
		buffer_sent[i] = new int*[D];
		buffer_recv[i] = new int*[D];

		for(int j = 0; j< D; j++)
		{ 
			buffer_sent[i][j] = new int[n];
			buffer_recv[i][j] = new int[n];
			for (int z = 0; z < n; z++)
			{
				buffer_sent[i][j][z] = 0;
				buffer_recv[i][j][z] = 0;
			}
		}
	}
	for (int j = 0; j < extended_buffer_len; j++)
	{ 
		status[j] = new int[D];
		for (int i = 0; i < D; i++)
			status[j][i] = 1;
	}
			
	map_row = new int*[w];
	map_col = new int* [w];
	i_map_row = new int* [w];
	i_map_col = new int* [w];

	for (int i = 0; i < w; i++)
	{
		map_row[i] = new int[n];
		map_col[i] = new int[n];
		i_map_row[i] = new int[m];
		i_map_col[i] = new int[m];
	}

	//The mapping rule for a TDZC
	for (int p = 0; p < w; p++)
	{
		for (int s = L; s < 2 * L; s++)
		{
			for (int c = 0; c < w; c++)
			{
				map_row[p][w * s + c] = w * (0 - s - delayed + L - 1) + c - p;
				map_col[p][w * s + c] = w * (s - L) + p;


			}
		}
	}
	for (int p = 0; p < w; p++)
	{
		for (int s = 0; s < L; s++)
		{
			for (int c = 0; c < w; c++)
			{
				i_map_row[p][w * s + c] = w * (s + 1 + delayed) + c - p;
				i_map_col[p][w * s + c] = w * (s + L) + p;


			}
		}
	}
	
}



