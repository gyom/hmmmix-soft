#ifndef HMMMIX_COMMON
#define HMMMIX_COMMON

#include <math.h>
/*#include "mex.h" probably not necessary, and messes up my Python */
#include <string.h>

void fwd_back(double * OUT_logSmoothedMarginals, double * OUT_logTwoSliceMarginals, double * logPrior, double * logTransmat, int transmatStride, double * logObslik, int K, int T, int twoSliceOption);

double hmmmix_p_Ypt_given_Mpt_Gp(double * Ypt, int Mgt, double * alpha_K, double * mu_Kp, double * nu_K, double * lambda_Kp, int K);

double hmmmix_Etgk(double * C, double * Y_Pt, int Mgt, double * alpha_K, double * mu_KP, double * nu_K, double * lambda_KP, int P, int G, int g, int K);

void transition_matrix_MLE_from_twoslice_marginals(double * A_KKg, double * xi_KKTg, int K, int Tminus1);
void transition_matrix_MLE_from_twoslice_marginals_pseudocounts_normalized(double * A_KKg, double * xi_KKTg, int K, int Tminus1, double * xi_KK_pseudocounts);


double chain_negKL_from_twoslice_marginals(double * Ag, double * eta_KKTg, double * hM_K1g, double * chainsPrior, int K, int Tminus1);

void normalizeRows(double * A, int M, int N, double padding);
void normalizeColumns(double * A, int M, int N, double padding);


void sampleFromSingleHiddenChain(int * samples_TN, double * smoothedFirstStates_K, double * nxi_KKTm1, int K, int T, int N, double * uniformRandomValues_TN);
int sampleVectorForIndex(double * V, int K, double u, int stride);

void normalize_xi_KKBG(double * xi_KKBG, int K, int B, int G);
void findFirstPositiveIndexInEveryRow(int * returnedIndices, double * E, int N, int D);

double vectorMinimum(double *A, int L);
double vectorMaximum(double *A, int L);

#endif

