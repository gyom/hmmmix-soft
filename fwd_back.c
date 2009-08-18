#include <math.h>
#include <stdlib.h>
#include <assert.h>

static double vectorMaximum(double * A, int M);
static double vectorRescaleLogs(double * A, int M);
static double vectorLogsumexp(double * A, int M);
static void vectorAdd(double * destination, double * A, double * B, int L);
static void vectorMultiply(double * destination, double * A, double * B, int L);
static void vectorAdd_withStrides(double * destination, double * A, double * B, int L, int strideA, int strideB);
static void vectorMultiply_withStrides(double * destination, double * A, double * B, int L, int strideA, int strideB);

/* This function writes to the already allocated memory regions pointed to by the pointers OUT_logSmoothedMarginals and OUT_logTwoSliceMarginals. If the argument twoSliceOption is set to 0, the second output argument won't be written to (nor computed).

   We are using the Fortran (and Matlab) format for array indexing, where the first index runs faster than the last index.
   e.g. A(1,1), A(2,1), A(3,1) are consecutive addresses in memory. 
   All the arrays are referenced as (double *) instead of multi-dimensional arrays so we have to keep track of the strides.

   The dimensions of the variables are the following.
      K = number of discrete states
      T = length of the time sequence
      OUT_logSmoothedMarginals : (K,T)
      OUT_logTwoSliceMarginals : (K,K,T)
      logPrior : (K)
      logTransmat : (K,K) or (K,K,T-1)
      logObslik : (K,T)

   It's possible to use a different transition matrix at every step of the sequence by specifying the argument logTransmat as an array of size (K,K,T-1). In that case, we need to set the value of 'transmatStride' to K*K instead of 0. This way, at step t we will use the K-by-K transition matrix found at logTransmat + K*K*(t-1).

   Note that all the arguments are in "log form". It would have been really nice to be able to support values of -Inf in the transition matrices so that we'd have exp(-Inf)=0 but the #defines from math.h don't seem to work properly. I'd rather not support those than rely on something that's not portable. Putting values of -100000 would be a reasonable to do instead of having -Inf.

Finally, no attempt is made here to keep track of the would-be normalizing constants that we'd get if we didn't work with logarithms. When doing the logsumexp trick, we simply work with a maximum function and we don't really normalize.
 */

void fwd_back(double * OUT_logSmoothedMarginals, double * OUT_logTwoSliceMarginals, double * logPrior, double * logTransmat, int transmatStride, double * logObslik, int K, int T, int twoSliceOption) {

    double * alpha;
    double * beta;
    int j,k,t;
    double * workbench;
    double paddingfactor;

    /* We call those "alpha" and "beta", but they're actually the logs of those quantities. */
    alpha = (double *)malloc(sizeof(double)*K*T);
    beta = (double *)malloc(sizeof(double)*K*T);
    workbench = (double *)malloc(sizeof(double)*K);

    /* debugging
    for(k=0;k<K;++k) {
        mexPrintf("logPrior[%d] = %f\n", k, logPrior[k]);
    }
    for(t=0;t<T;++t) {
        for(k=0;k<K;++k) {
             mexPrintf("logObslik[%d+K*%d] = %f\n", k, t, logObslik[k+K*t]);
        }
    }
    for(t=0;t<T-1;++t)
        for(k=0;k<K;++k)
            for(j=0;j<K;++j)
                mexPrintf("logTransmat[%d+K*%d+transmatStride*%d] = %f\n", j, k, t, logTransmat[j + K*k + transmatStride*t]);
    */
   
    /* Forwards. */
    /*mexPrintf("Initializing forwards.\n");*/
    t=0;
    for(k=0;k<K;++k) {
        alpha[k+K*t] = logPrior[k] + logObslik[k+K*t];
    }
    vectorRescaleLogs(alpha + K*t, K);

    /*mexPrintf("Running forwards.");*/
    /* formula (26.17) from my pdf of Kevin's book */
    for(t=1;t<T;++t) {
	for(k=0;k<K;++k) {
	    /* getting column k of the transition matrix */
	    vectorAdd(workbench, logTransmat + K*k + transmatStride*(t-1), alpha + K*(t-1), K);
	    paddingfactor = vectorRescaleLogs(workbench, K);
	    alpha[k+K*t] = logObslik[k+K*t] + paddingfactor + vectorLogsumexp(workbench,K);
	}
	/* in the spirit of normalization, but not quite */
	vectorRescaleLogs(alpha + K*t, K);
    }

    /* Backwards. */
    /*mexPrintf("Initializing backwards.");*/
    t=T-1;
    for(k=0;k<K;++k)
	beta[k+K*t] = 0; /* log(1) = 0 */
    
    /*mexPrintf("Running backwards.");*/
    for(t=T-2;t>=0;--t) {
	for(k=0;k<K;++k) {
	    /* Getting row k of the transition matrix. */
	    /* Double-check that thing with the transition matrices that run one less index. */
	    vectorAdd_withStrides(workbench, logTransmat + k + transmatStride*t, beta + K*(t+1), K, K, 1);
	    vectorAdd(workbench, workbench, logObslik + K*(t+1), K);
	    paddingfactor = vectorRescaleLogs(workbench, K);
	    beta[k+K*t] = paddingfactor + vectorLogsumexp(workbench,K);
	}
    }

    /* Combine the two into the marginals. */
    /*mexPrintf("Combining forwards and backwards.");*/
    for(t=0;t<T;++t) {
	vectorAdd(OUT_logSmoothedMarginals + K*t, alpha + K*t, beta + K*t, K);
	vectorRescaleLogs(OUT_logSmoothedMarginals + K*t, K);
    }

    /* Get the two-slice marginals if we want them. */
    /*mexPrintf("Getting twoslice marginals.");*/
    if (twoSliceOption != 0) {
	for(t=0;t<T-1;++t) {
	    for(k=0;k<K;++k)
		for(j=0;j<K;++j)
		    OUT_logTwoSliceMarginals[j + K*k + K*K*t] = logTransmat[j+K*k+transmatStride*t] + alpha[j+K*t] + logObslik[k+K*(t+1)] + beta[k+K*(t+1)];
	    vectorRescaleLogs(OUT_logTwoSliceMarginals + K*K*t, K*K);
	}
    }

    free(alpha);
    free(beta);
    free(workbench);
    return;
}

static double vectorMaximum(double * A, int M) {
    assert(M>0);
    double result = A[0];
    int m;
    for(m=0;m<M;++m) {
        /*mexPrintf("vectorMaximum : A[%d] = %f\n", m, A[m]);*/
        if (A[m] > result)
            result = A[m];
    }
    return result;
}

/* Update A = A - max(A) and returns that old value max(A). */
static double vectorRescaleLogs(double * A, int M) {
    double maxA;
    int m;
    maxA = vectorMaximum(A, M);
    for(m=0;m<M;++m)
	A[m] -= maxA;
    return maxA;
}

/* You should call vectorRescaleLogs before doing this. */
static double vectorLogsumexp(double * A, int M) {
    double total = 0;
    int m;
    for(m=0;m<M;++m)
	total += exp(A[m]);
    /*mexPrintf("vectorLogsumexp : log(%f)\n", total);*/
    /*assert(total > 0);*/
    return log(total);
}

static void vectorAdd(double * destination, double * A, double * B, int L) {
    int m;
    for(m=0;m<L;++m)
	destination[m] = A[m] + B[m];
}

static void vectorAdd_withStrides(double * destination, double * A, double * B, int L, int strideA, int strideB) {
    int m;
    for(m=0;m<L;++m)
	destination[m] = A[m*strideA] + B[m*strideB];
}

static void vectorMultiply(double * destination, double * A, double * B, int L) {
    int m;
    for(m=0;m<L;++m)
	destination[m] = A[m] * B[m];
}

static void vectorMultiply_withStrides(double * destination, double * A, double * B, int L, int strideA, int strideB) {
    int m;
    for(m=0;m<L;++m)
	destination[m] = A[m*strideA] * B[m*strideB];
}
