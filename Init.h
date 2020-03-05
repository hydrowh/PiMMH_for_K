#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iterator>     // std::back_inserter

using namespace std;

int I;
int J;
int K;
int T;
int Node_I; 
int Node_J;
int dim_theta;
int dim_x;
int dim_y;

int NI,NJ;
int ON;
int i, j, k, t;
int Burn_in, Out_i;
double M;
double *X;//[I / Node_I][J / Node_J][T][dim_x];
double *Wt_ksum;//[I / Node_I][J / Node_J][T];
double *theta;//[I / Node_I][J / Node_J][dim_theta];
vector<vector<double> > Y;

FILE* input1, input2, Output;

int Init(int rank, int numprocs) {
	NI = 0;
	NJ = 0;
	int rand_seed = 0;
	M = 1.0;
	fstream file;
	ifstream ifs("setting.csv");
	if (!ifs) {
		cout << "設定ファイル入力エラー";
		return 1;
	}
	else {
		cout << "設定ファイル入力OK" << "\n";
	}

	int num, count;
	string str;
	num = 0;
	while (getline(ifs, str)) {
		if (num == 1) {
			string token;
			istringstream stream(str);
			count = 0;
			while (getline(stream, token, ',')) {
				if (count == 0) {
					istringstream iss(token);
					iss >> I;
				}
				else if (count == 1) {
					istringstream iss(token);
					iss >> J;
				}
				else if (count == 2) {
					istringstream iss(token);
					iss >> K;
				}
				else if (count == 3) {
					istringstream iss(token);
					iss >> T;
				}
				else if (count == 4) {
					istringstream iss(token);
					iss >> dim_x;
				}
				else if (count == 5) {
					istringstream iss(token);
					iss >> dim_y;
				}
				else if (count == 6) {
					istringstream iss(token);
					iss >> dim_theta;
				}
				else if (count == 7) {
					istringstream iss(token);
					iss >> Node_I;
				}
				else if (count == 8) {
					istringstream iss(token);
					iss >> Node_J;
				}
				count++;
			}
		}
		if (num == 3) {
			string token;
			istringstream stream(str);
			count = 0;
			while (getline(stream, token, ',')) {
				if (count == 0) {
					istringstream iss(token);
					iss >> rand_seed;
				}
				else if (count == 1) {
					istringstream iss(token);
					iss >> M;
				}
				else if (count == 2) {
					istringstream iss(token);
					iss >> ON;
				}
				else if (count == 3) {
					istringstream iss(token);
					iss >> Burn_in;
					if (Burn_in < 1)Burn_in = 1;
				}
				else if (count == 4) {
					istringstream iss(token);
					iss >> Out_i;
					if (Out_i < 1)Out_i = 1;
				}
				count++;
			}
		}
		num += 1;
	}
	ifs.close();

	NI = I / Node_I;
	NJ = J / Node_J;
	if (rank == 0) {
		if (I == 0 || J == 0 || K == 0 || T == 0 ||
			dim_theta == 0 || dim_x == 0 || dim_y == 0 ||
			I % Node_I != 0 || J % Node_J != 0 || numprocs != Node_I*Node_J) {
			cout << "入力ファイルエラー";
			return 1;
		}
		else {
			cout << "入力OK" << "\n";
		}
	}
	init_genrand(rand_seed + rank);//rand seed

	X = (double*)malloc(NI * NJ * T * dim_x * sizeof(double));
	Wt_ksum = (double*)malloc(NI * NJ * T * sizeof(double));
	theta = (double*)malloc(NI * NJ * dim_theta * sizeof(double));

	Y.resize(T);
	for (t = 0; t < T; t++) {
		Y[t].resize(dim_y);
	}

	ifstream ifsy("y_data.csv");
	if (!ifsy) {
		cout << "Y入力エラー";
		return 1;
	}
	else {
		cout << "Y入力OK" << "\n";
	}
	int dim, time;
	dim = 0;
	while (getline(ifsy, str)) {
		string token;
		istringstream stream(str);
		time = 0;
		while (getline(stream, token, ',') && time < T) {
			istringstream iss(token);
			double tmp = 0.0;
			iss >> tmp;
			Y[time][dim] = tmp;
			time++;
		}
		dim++;
	}
	ifsy.close();
	cout << "Init OK" << "\n";
	return 0;
}
