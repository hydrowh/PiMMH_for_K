#include "mpi.h"
#include "Model.h"

int main(int argc, char *argv[]) {

	int rank, numprocs;
	int root;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	int akt;
	int i, j, k, l, t;
	double tmp;
	double p1, p2, p3;
	vector<vector<vector<double> > > Xtk;
	vector<vector<double> > X_akt;
	vector<vector<double> > X_t1;
	vector<vector<double> > Wtk;
	vector<double> theta_t;
	vector<vector<int> > AKT;
	root = 0;
	Init(rank, numprocs);
	
	AKT.resize(T - 1);
	for (t = 0; t < T; t++) {
		AKT[t].resize(K);
	}

	Xtk.resize(T);
	for (t = 0; t < T; t++) {
		Xtk[t].resize(K);
		for (k = 0; k < K; k++) {
			Xtk[t][k].resize(dim_x);
		}
	}
	
	X_akt.resize(K);
	for (k = 0; k < K; k++) {
		X_akt[k].resize(dim_x);
	}
	X_t1.resize(K);
	for (k = 0; k < K; k++) {
		X_t1[k].resize(dim_x);
	}
	
	Wtk.resize(T);
	for (t = 0; t < T; t++) {
		Wtk[t].resize(K);
	}

	theta_t.resize(dim_theta);

	for (i = 0; i < NI; i++) {
		for (j = 0; j < NJ; j++) {
			/////////////////PMMH method/////////////////////
			Q_theta_d(theta_t);//Algorithm8 line 3
			for (l = 0; l < dim_theta; l++) {
				theta[i * NJ * dim_theta + j * dim_theta + l] = theta_t[l];
			}
			for (k = 0; k < K; k++) {
				Xtk_init(Xtk[0][k], Y[0], theta_t);//Algorithm1 line 4
				Wtk[0][k] = Wtk_init(Xtk[0][k], Y[0], theta_t);//Algorithm1 line 5
			}

			for (t = 0; t < T - 1; t++) {
				for (k = 0; k < K; k++) {
					AKT[t][k] = Draw_akt(Wtk[t]);//Algorithm1 line 8
					for (l = 0; l < dim_x; l++) {
						X_akt[k][l] = Xtk[t][AKT[t][k]][l];
					}
					rand_Xt1__yt1_xakt_theta(X_t1[k], Y[t + 1], X_akt[k], theta_t);// Algorithm1 line 9

					for (l = 0; l < dim_x; l++) {
						Xtk[t + 1][k][l] = X_t1[k][l];
					}
					// Algorithm1 line 10
					p1 = P_yt1__xt1_theta(Y[t + 1], X_t1[k], theta_t);
					p2 = P_xt1__xakt_theta(X_t1[k], X_akt[k], theta_t);
					p3 = Q_xt1__yt1_xakt_theta(X_t1[k], Y[t + 1], X_akt[k], theta_t);
					Wtk[t + 1][k] = (p1 * p2 / p3);
				}
			}

			akt = Draw_akt(Wtk[T - 1]);// Algorithm1 line 11 first half
			for (l = 0; l < dim_x; l++) {// Algorithm1 line 11 latter half
				X[i * NJ * T * dim_x + j * T * dim_x + (T - 1) * dim_x + l] = Xtk[T - 1][akt][l];
			}
			for (t = T - 2; t >= 0; t--) {//Algorithm1 line 12~13
				akt = AKT[t][akt];
				for (l = 0; l < dim_x; l++) {
					X[i * NJ * T * dim_x + j * T * dim_x + t * dim_x + l] = Xtk[t][akt][l];
				}
			}

			for (t = 0; t < T; t++) {
				tmp = 0.0;
				for (k = 0; k < K; k++) {
					tmp += Wtk[t][k];
				}
				Wt_ksum[i * NJ * T + j * T + t] = tmp;
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (rank != root) {
		MPI_Send(&theta[0],   NI*NJ*dim_theta, MPI_DOUBLE, root, 0, MPI_COMM_WORLD);
		MPI_Send(&Wt_ksum[0], NI*NJ*T,         MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
		MPI_Send(&X[0],       NI*NJ*T*dim_x,   MPI_DOUBLE, root, 2, MPI_COMM_WORLD);
	}
	else if (rank == root) {
		////////// Algorithm 7////////////////
		double p_t, q_t, RJ_O;
		int ni, nj;
		int i_p, it, jt;
		char filename[100];

		vector<vector<double> > theta_I;// [I][dim_theta]
		vector<vector<double> > theta_m;// [I][dim_theta]
		vector<vector<vector<double> > > X_I;// [I][T][dim_x]
		vector<vector<vector<double> > >  X_m;// [I][T][dim_x]
		vector<vector<double> > Wt_ksum_I;// [I][T]
		vector<double> rj_ave;// [I]
		vector<double> rj_real;// [I]

		theta_I.resize(I);
		theta_m.resize(I);
		X_I.resize(I);
		X_m.resize(I);
		Wt_ksum_I.resize(I);
		rj_ave.resize(I);
		rj_real.resize(I);
		for (i = 0; i<I; i++) {
			theta_I[i].resize(dim_theta);
			theta_m[i].resize(dim_theta);
			X_I[i].resize(T);
			X_m[i].resize(T);
			Wt_ksum_I[i].resize(T);
			for (t = 0; t < T; t++) {
				X_I[i][t].resize(dim_x);
				X_m[i][t].resize(dim_x);
			}
			rj_ave[i] = 0.0;
			rj_real[i] = 0.0;
		}
		
		fstream file_r, file_th[100], file_xi[100];
		sprintf(filename, "./rj_result.csv");
		file_r.open(filename, ios::out);
		for (k = 0; k < dim_theta; k++) {
			sprintf(filename, "./theta_result_dim%d.csv", k + 1);
			file_th[k].open(filename, ios::out);
		}
		for (k = 0; k < dim_x; k++) {
			sprintf(filename, "./xi_result_dim%d.csv", k + 1);
			file_xi[k].open(filename, ios::out);
		}

		for (i = 0; i<I; i++) {
			if (Out_i == 1 || (Out_i > 1 && (i + 1) % Out_i == 0) || i == 0 || i == I - 1) {
				file_r << ",i=" << i + 1;
				for (k = 0; k < dim_theta; k++) {
					file_th[k] << ",i=" << i + 1;
				}
			}
		}
		file_r << "\n";
		for (k = 0; k < dim_theta; k++) {
			file_th[k] << "\n";
		}
		
		gmp_randstate_t state;
		mpf_t W_ij, Hp, Hp_t, WP, WP_t, DBLM;
		mpf_t RJ, RP, TMPW1, TMPW2, MK_inv_gmp, I_inv_gmp, P_gmp, Q_gmp;
		mpf_t RJ_ave;
		mpf_t *WIJ = (mpf_t*)malloc(sizeof(mpf_t)*I);
		mpf_t *WJ_1 = (mpf_t*)malloc(sizeof(mpf_t)*I);

		mpf_init(W_ij);
		mpf_init(Hp);
		mpf_init(Hp_t);
		mpf_init(WP);
		mpf_init(WP_t);
		mpf_init(DBLM);
		mpf_init(RJ);
		mpf_init(RP);
		mpf_init(TMPW1);
		mpf_init(TMPW2);
		mpf_init(MK_inv_gmp);
		mpf_init(I_inv_gmp);
		mpf_init(P_gmp);
		mpf_init(Q_gmp);
		mpf_init(RJ_ave);
		for (i = 0; i<I; i++) {
			mpf_init(WIJ[i]);
			mpf_init(WJ_1[i]);
		}
		
		printf("Alg7 start\n");
		mpf_set_d(TMPW1, M);
		mpf_set_d(TMPW2, (double)K);
		mpf_set_d(DBLM, DBL_MIN);
		mpf_div(MK_inv_gmp, TMPW1, TMPW2);

		mpf_set_d(RJ_ave, 0.0);
		gmp_randinit_mt(state);
		for (nj = 0; nj < Node_J; nj++) {
			for (j = 0; j < NJ; j++) {
				jt = j + nj * NJ;
				mpf_set_d(WP_t, 0.0);
				for (ni = 0; ni < Node_I; ni++) {
					printf("nj=%d ni=%d completed !\n", nj, ni);
					if (j == 0) {
						if ((Node_J > 1 && nj > 0) || (Node_J == 1)) {
							if ((Node_I > 1 && ni > 0) || (Node_I == 1)) {
								MPI_Recv(&theta[0],   NI*NJ*dim_theta, MPI_DOUBLE, ni * Node_J + nj, 0, MPI_COMM_WORLD, &status);
								MPI_Recv(&Wt_ksum[0], NI*NJ*T,         MPI_DOUBLE, ni * Node_J + nj, 1, MPI_COMM_WORLD, &status);
								MPI_Recv(&X[0],       NI*NJ*T*dim_x,   MPI_DOUBLE, ni * Node_J + nj, 2, MPI_COMM_WORLD, &status);
							}
						}
					}

					for (i = 0; i < NI; i++) {
						it = i + ni * NI;

						for (l = 0; l < dim_theta; l++) {
							theta_I[it][l] = theta[i * NJ * dim_theta + j * dim_theta + l];
						}
						for (t = 0; t < T; t++) {
							Wt_ksum_I[it][t] = Wt_ksum[i * NJ * T + j * T + t];
						}
						for (t = 0; t < T; t++) {
						    for (l = 0; l < dim_x; l++) {
								X_I[it][t][l] = X[i * NJ * T * dim_x + j * T * dim_x + t * dim_x + l];
							}
						}
						//h* 
						mpf_set_d(Hp, 1.0);
						for (t = 0; t < T; t++) {
							mpf_set_d(TMPW1, Wt_ksum_I[it][t]);
							mpf_mul(TMPW2, TMPW1, MK_inv_gmp);
							mpf_set(Hp_t, Hp);
							mpf_mul(Hp, TMPW2, Hp_t);
						}

						p_t = P_theta_p(theta_I[it]);
						if (isnan(p_t)) {
							p_t = DBL_MIN;
						}
						mpf_set_d(P_gmp, p_t);
						mpf_set(Hp_t, Hp);
						mpf_mul(Hp, Hp_t, P_gmp);
						//h*

						//Wij
						q_t = Q_theta_p(theta_I[it]);
						if (isnan(q_t)) {
							q_t = DBL_MIN;
						}
						mpf_set_d(Q_gmp, q_t);
						mpf_div(W_ij, Hp, Q_gmp);
						//Wij

						mpf_set(WIJ[it], W_ij);
						mpf_set(TMPW1, WP_t);
						mpf_add(WP_t, TMPW1, W_ij);

						mpf_urandomb(RP, state, ON);
						i_p = D_gmp(WIJ, RP, (it + 1));
						mpf_set_d(TMPW2, (double)(it + 1));
						mpf_div(WP, WP_t, TMPW2);

						if (jt == 0) {
							mpf_set_d(RJ, 1.0);
							mpf_set(WJ_1[it], WP);
							for (k = 0; k < dim_theta; k++) {
								theta_m[it][k] = theta_I[i_p][k];
							}

							for (k = 0; k < dim_x; k++) {
								for (t = 0; t < T; t++) {
									X_m[it][t][k] = X_I[i_p][t][k];
								}
							}
						}
						else {
							mpf_set_d(RJ, 1.0);
							if (0 < mpf_cmp_d(WJ_1[it], 0.0)) {
								mpf_div(RJ, WP, WJ_1[it]);
							}
							if (0 < mpf_cmp_d(RJ, 1.0)) {
								mpf_set_d(RJ, 1.0);
							}
							mpf_urandomb(RP, state, ON);
							if (mpf_cmp(RP, RJ) <= 0) {
								mpf_set(WJ_1[it], WP);
								if (jt >= Burn_in) {
									rj_real[it] += 1.0;
								}
								for (k = 0; k < dim_theta; k++) {
									theta_m[it][k] = theta_I[i_p][k];
								}
								for (k = 0; k < dim_x; k++) {
									for (t = 0; t < T; t++) {
										X_m[it][t][k] = X_I[i_p][t][k];
									}
								}
							}
						}

						if (mpf_cmp(DBLM, RJ) <= 0) {
							RJ_O = mpf_get_d(RJ);
						}
						else {
							RJ_O = DBL_MIN;
						}
						if (jt >= Burn_in) {
							rj_ave[it] += RJ_O;
						}

						if (it == 0) {
							file_r << "j=" << (jt + 1);
							for (k = 0; k < dim_theta; k++) {
								file_th[k] << "j=" << (jt + 1);
							}
						}
						if (it == 0 || (it + 1) % Out_i == 0 || it == I - 1) {
							file_r << "," << RJ_O;
							for (k = 0; k < dim_theta; k++) {
								file_th[k] << "," << theta_m[it][k];
							}
						}
						if (it == I - 1) {
							file_r << "\n";
							for (k = 0; k < dim_theta; k++) {
								file_th[k] << "\n";
							}
						}

					}//i
				}//ni
				for (k = 0; k < dim_x; k++) {
					file_xi[k] << "j=" << (jt + 1);
					for (t = 0; t < T; t++) {
						file_xi[k] << ", " << X_m[I - 1][t][k];
					}
					file_xi[k] << "\n";
				}
			}//j
		}//nj

		file_r << "Average acceptance probability";
		for (i = 0; i<I; i++) {
			if (Out_i == 1 || (Out_i > 1 && (i + 1) % Out_i == 0) || i == 0 || i == I - 1) {
				file_r << "," << (rj_ave[i] / (double)(J - Burn_in));
			}
		}
		file_r << "\n";
		file_r << "Realized overall acceptance rate";
		for (i = 0; i<I; i++) {
			if (Out_i == 1 || (Out_i > 1 && (i + 1) % Out_i == 0) || i == 0 || i == I - 1) {
				file_r << "," << (rj_real[i] / (double)(J - Burn_in));
			}
		}
		file_r << "\n";
		file_r.close();
		for (k = 0; k < dim_theta; k++) {
			file_th[k].close();
		}
		for (k = 0; k < dim_x; k++) {
			file_xi[k].close();
		}

		mpf_clear(W_ij);
		mpf_clear(Hp);
		mpf_clear(Hp_t);
		mpf_clear(WP);
		mpf_clear(WP_t);
		mpf_clear(DBLM);
		mpf_clear(RJ);
		mpf_clear(RP);
		mpf_clear(TMPW1);
		mpf_clear(TMPW2);
		mpf_clear(MK_inv_gmp);
		mpf_clear(I_inv_gmp);
		mpf_clear(P_gmp);
		mpf_clear(Q_gmp);
		for (i = 0; i<I; i++) {
			mpf_clear(WIJ[i]);
			mpf_clear(WJ_1[i]);
		}
		mpf_clear(RJ_ave);
		printf("All completed\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
