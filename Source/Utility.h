///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Designed by Xiao Ma (maxiao@mail.sysu.edu.cn), Sun Yat-sen University.
//This program can be employed for academic research only.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#pragma once

#ifndef _UTILITY_H
#define _UTILITY_H

#define PI 3.14159265358979
#define SMALLPROB 1.0e-12
#define ERASE -1

void ProbClip(double *xx, int len_xx);

void Dec2Bin(int d, int *b, int len_b);
void SeqDec2Bin(int *bin, int *dec, int len_dec, int len_symbol);

void SeqBin2Dec(int *bin, int *dec, int len_dec, int len_symbol);

int BitDotProd(int a, int b, int len);

int Systemizer(int num_row, int num_col, int **H, int **sysH, int *pai);

int min(int x, int y);
int max(int x, int y);

void BubbleSort(double *value, int *index, int len);

//
#endif