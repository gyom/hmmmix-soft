#include "hmmmix_common.h"

/* using the full C instead of the partial C like in Matlab */
double hmmmix_Etgk(double * C, double * Y_Pt, int Mgt, double * alpha_K, double * mu_KP,\
		   double * nu_K, double * lambda_KP, int P, int G, int g, int K) {

    double Etgk = 0;
    int p;
    double safety;

    for (p=0;p<P;++p) {
	safety = hmmmix_p_Ypt_given_Mpt_Gp(Y_Pt + p, Mgt, alpha_K, mu_KP + p*K, nu_K, lambda_KP + p*K, K);
	if (safety > 0)
	    Etgk = Etgk + C[g+p*G] * log(safety );
	else {
	    Etgk = -1000;
	    /*mexPrintf("This should never happen. You have a problem in the design. Etgk = -1000.\n");*/
	}
    }
    return Etgk;
}

double hmmmix_p_Ypt_given_Mpt_Gp(double * Ypt, int Mgt, double * alpha_K, double * mu_Kp, double * nu_K, double * lambda_Kp, int K) {

    double p_Ypt_given_Mgt_Gp = 0;
    double contribution;
    int z;

    for (z=0;z<K;++z) {
	if (Mgt == z)
	    contribution = alpha_K[Mgt] / (alpha_K[Mgt] + 2);
	else
	    contribution = 1 / (alpha_K[Mgt] + 2);

	/* Omitting a multiplicative constant from the tpdf that depends on nu.
     * This was fine when all the patients had the same value of nu, but in the most
     * recent version we have that nu(1:K, 1:P) can vary. In that case, it would be wise
     * to put compute the values of gamma((nu+1)/2) and gamma(nu/2) somewhere in there.
     */
	p_Ypt_given_Mgt_Gp += \
	    contribution * sqrt(lambda_Kp[z])*pow(1 + pow((Ypt[0] - mu_Kp[z]),2)*lambda_Kp[z]/nu_K[z], -0.5*(nu_K[z] + 1));
    }
    
    return p_Ypt_given_Mgt_Gp;
}

/* For bad design reasons, I wind up usually passing a K-K-T array when I just
 * want to use the first T-1 matrices. I'll stop this madness here and the following
 * function will use all the matrices up to the last index. It's left to the user to
 * supply the -1 offset.
 */
void transition_matrix_MLE_from_twoslice_marginals(double * A_KKg, double * xi_KKTg, int K, int Tminus1) {

    int t,i,j;

    for(i=0;i<K;++i) {
    	for(j=0;j<K;++j) {
            A_KKg[i+j*K] = 0;
                for(t=0;t<Tminus1;++t) {
            		A_KKg[i+j*K] += xi_KKTg[i+j*K+t*K*K];
           	    }
        	    A_KKg[i+j*K] /= Tminus1;
    	}
    }
    return;
}

/* Sohrab's implementation has pseudocounts of 1 for transitions to different
 * states and 100 for same states. I'm having some troubles beating his
 * implementation so I'll try to add this.
 */
void transition_matrix_MLE_from_twoslice_marginals_pseudocounts_normalized(double * A_KKg, double * xi_KKTg, int K, int Tminus1, double * xi_KK_pseudocounts) {

    int t,i,j;
    
    for(i=0;i<K;++i) {
    	for(j=0;j<K;++j) {
            A_KKg[i+j*K] = 0;
            for(t=0;t<Tminus1;++t) {
         		A_KKg[i+j*K] += xi_KKTg[i+j*K+t*K*K];
      	    }
      	    A_KKg[i+j*K] += xi_KK_pseudocounts[i+K*j];
    	}
    }
    
    normalizeRows( A_KKg, K, K, 0);
    
    return;
}



/* We're assuming that the rows of eta and A are normalized. If they're not,
   the results won't accurate, but they'll still mean something useful in some sense.

   If hM_K1g is NULL, we aren't including the term corresponding to the initial states.
   If chainsPrior is NULL, we are assuming the default value of 1/K everywhere.

   The "1" in hM_K1g is meant to represent the first state, so it's actually the zeroth column.
*/
double chain_negKL_from_twoslice_marginals(double * Ag, double * eta_KKTg, double * hM_K1g, double * chainsPrior, int K, int Tminus1) {

    double entropy = 0;
    int t,j,k;
    double safety;

    if (hM_K1g != NULL) {
        for(k=0;k<K;++k) {
            if (hM_K1g[k] > 0)
                entropy -= hM_K1g[k] * log(hM_K1g[k]);

            if (chainsPrior == NULL)
                entropy += hM_K1g[k] * -log(K);
            else if (chainsPrior[k] > 0)
                entropy += hM_K1g[k] * log(chainsPrior[k]);
        }
    }

    for(t=0;t<Tminus1;++t) {
        for(j=0;j<K;++j) {
            for(k=0;k<K;++k) {
                if (eta_KKTg[j + K*k + K*K*t] > 0)
                    entropy -= eta_KKTg[j + K*k + K*K*t] * log(eta_KKTg[j + K*k + K*K*t]);

                if (Ag[j+K*k] > 0)
                    entropy += eta_KKTg[j + K*k + K*K*t] * log(Ag[j+K*k]);
                }
        }
    }

    /* which isn't really the entropy */
    return entropy;

}


void normalizeRows(double * A, int M, int N, double padding) {

    int m,n;
    double sum;

    for (m=0;m<M;++m) {
        /* find the sum for that row m*/
        sum = 0;
        for (n=0;n<N;++n) {
            A[m+M*n] += padding;
            sum += A[m+M*n];
        }

        /* then divide */
        for (n=0;n<N;++n) {
            if (sum > 0)
                A[m+M*n] /= sum;
            else
                A[m+M*n] = 1/N;
        }
    }
    return;
}


void normalizeColumns(double * A, int M, int N, double padding) {

    int m,n;
    double sum;

    for (n=0;n<N;++n) {
        /* find the sum for that column n*/
        sum = 0;
        for (m=0;m<M;++m) {
            A[m+M*n] += padding;
            sum += A[m+M*n];
        }

        /* then divide */
        for (m=0;m<M;++m) {
            if (sum > 0)
                A[m+M*n] /= sum;
            else
                A[m+M*n] = 1/M;
        }
    }
    return;
}




/* Samples N sequences of T values in {0,...,K-1} from a single HMM with observations
   defined from beliefs on the initial state and on the specific transitions at every
   step (obtained by running forwards-backwards earlier).

   The output goes into the first argument, and we use the last argument to get uniform
   random values in [0,1] from the calling function.
 
   It's very important to note that we're often using the two-slice marginals whose rows do not
   sum to 1. We can't naively sample with a cumsum because of that. We have to renormalize. This
   is why there is an "n" for "normalized" in front of the two-slice marginals argument name.
 
*/

void sampleFromSingleHiddenChain(int * samples_TN, double * smoothedFirstStates_K, double * nxi_KKTm1, int K, int T, int N, double * uniformRandomValues_TN) {

    int t,n,k;

    for (n=0; n<N; ++n) {

        samples_TN[0 + T*n] = sampleVectorForIndex(smoothedFirstStates_K, K, uniformRandomValues_TN[0 + T*n], 1);

        for (t=1;t<T;++t) {
             /* a certain row of that normalize xi at slice (t-1) going to time t */
            samples_TN[t + T*n] = sampleVectorForIndex(nxi_KKTm1 + ( samples_TN[t-1 + T*n ] + K*K*(t-1) ), K, uniformRandomValues_TN[t + T*n], K);
        }
    }

    return;
}

/* Does the thing you'd expect from values
   V = [0.1, 0.4, 0.3, 0.2]
   K = 4
   u = rand(1,1)

   I added the stride argument because I was going to use this function to
   sample along a row of xi_KKTm1 and I needed a stride of K. The default
   value for the stride is 1 and not 0.
*/
int sampleVectorForIndex(double * V, int K, double u, int stride) {

    int k;
    double cumsum = 0;

    for (k=0; k<K; ++k) {
	cumsum += V[k*stride];
	if (cumsum >= u)
	    return k;
    }
    /* This should never happen, but it's not like I'm checking for errors here.
       If this function returns K, it's because V wasn't normalized, because K was
       wrong or because u wasn't in [0,1]
    */
    return K;
}

/* let's write B = T-1 */
void normalize_xi_KKBG(double * xi_KKBG, int K, int B, int G) {
    int t,g;
    
    for (t=0;t<B;++t)
        for (g=0;g<G;++g)
            normalizeRows(xi_KKBG + (K*K*t + K*K*B*g), K, K, 0);
            
    return;
}

/* I'm not sure what is the equivalent in Numpy of the Matlab "find" function
 * that can be applied to multiple lines to find, for example, what are the first
 * columns that have a nonzero value (and do that for every row).
 */

/*void sampleRowMultinomial(int * sampledIndices, double * E, int N, int D, double * uniformRandomValues_N) {
    / We expect E to be a table with N rows, D columns, with the columns normalized. /
    int n, d;

    for (n=0;n<N;++n) {
        
    
    }
    return;
}*/

void findFirstPositiveIndexInEveryRow(int * returnedIndices, double * E, int N, int D) {
    /* E is N-by-D, fortran ordering */
	
	int n,d;
	for(n=0;n<N;++n)
		for(d=0;d<D;++d)
			if(E[n+N*d] > 0) {
				returnedIndices[n] = d;
				break;
			}
	return;
}

double vectorMinimum(double *A, int L) {
    double result = A[0];
    int n;
    for(n=0; n<L; ++n)
        if (A[n] < result)
            result = A[n];
    return result;
}

double vectorMaximum(double *A, int L) {
    double result = A[0];
    int n;
    for(n=0; n<L; ++n)
        if (A[n] < result)
            result = A[n];
    return result;
}

