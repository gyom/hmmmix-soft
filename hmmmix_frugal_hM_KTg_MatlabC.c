#include <math.h>
#include "mex.h"
#include <string.h>
#include "hmmmix_common.h"
#include <assert.h>

/* 
    This function is a simpler version of hmmmix_optimal_M_KTG that requires
    the values of R_TgK (see my report) and processes only one group.
 
    This function is meant to play really well with constants 'S' that
    control the temperature for the optimization, but I realized that it
    wasn't really desirable to actually use the 'S' inside here.
 
    It's convenient for 'S' because of the way by which we can just multiply
    the log evidences, pseudocounts and all when calling this function.

    Instead of having the last values of xi_KKTg(:,:,T) be junk,
    I decided to cut the legacy support and have that variable be of
    size (K,K,T-1).
 
    IN :
        R_KTg, chainsPrior_Kg, A_KKg, xi_pseudocounts_KK
    OUT :
        hM_KTg, xi_KKBg, updated_chainsPrior_Kg, updated_A_KKg, divergenceContributions
*/
 
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double * R_KTg, * chainsPrior_Kg, * A_KKg, * xi_pseudocounts_KK;
    
    int K, T, B;
    int j, k, t;

    double * hM_KTg, * xi_KKBg, * updated_chainsPrior_Kg, * updated_A_KKg, * divergenceContributions;
    
    /* I will need the logs of chainsPrior_Kg, A_KKg, but the log of A_KKg 
       is the only one for which I'd say computations by memorizing the results.
     */
    double * log_chainsPrior_Kg, * log_A_KKg;
    double * rownormalized_xi_KKBg;
    
    
    double * outputToolPtr;
    int trplarr_ndim = 3;
    int * trplarr_dims = mxMalloc(trplarr_ndim*sizeof(int));
    
    double logeps;
    
    if ( !(3<=nrhs && nrhs<=4) || !(1<=nlhs && nlhs<=5))
        mexErrMsgTxt("hmmmix_frugal_hM_KTg_MatlabC requires 3-4 inputs and 1-5 outputs");

    R_KTg = mxGetPr(prhs[0]);
    chainsPrior_Kg = mxGetPr(prhs[1]);
    A_KKg = mxGetPr(prhs[2]);

    
    K = mxGetM(prhs[0]);
    T = mxGetN(prhs[0]);
    
    /* I decided to define B=T-1 because a needed a single letter name for
     * that quantity T-1. I thought naming xi_KKTg when it actually had T-1
     * elements in the last dimension was more confusing than naming it
     * xi_KKBg. I'm not sure if this is a good call.
     */
    B = T-1;
    assert(mxGetM(prhs[1])*mxGetN(prhs[1]) == K);
    assert(mxGetM(prhs[2]) == K);
    assert(mxGetN(prhs[2]) == K);
    
    /* We're going to be using this value to make sure that the zeros in
     * the transition matrices correspond logs whose values are very small
     * compared the the loglikelihoods that we're dealing with in R_KTg.
     */
    logeps = T*vectorMinimum(R_KTg, K*T);
    if (logeps > 0)
        logeps = -T;
    
    /* Our implementation of fwd_back works with logs so we need to make
     * the conversions somewhere. It's a bit stupid to do it in C instead
     * of doing in Matlab with a wrapper, though.
     */
    
    log_A_KKg = mxMalloc(K*K*sizeof(double));
    for (k=0;k<K*K;++k) {
        if (A_KKg[k] < 0) {
            mexPrintf("Found a value < 0 in a transition matrix.\n");
            for (k=0;k<K*K;++k)
                mexPrintf("A_KKg[%d]=%f.\n", k, A_KKg[k]);
            mexErrMsgTxt("transition matrix fault");
        } else if (A_KKg[k] == 0)
            log_A_KKg[k] = logeps;
        else
            log_A_KKg[k] = log(A_KKg[k]);
    }

    log_chainsPrior_Kg = mxMalloc(K*sizeof(double));
    for (k=0;k<K;++k) {
        if (chainsPrior_Kg[k] < 0) {
            mexPrintf("Found a value < 0 in a the initial state priors.\n");
            for (k=0;k<K*K;++k)
                mexPrintf("chainsPrior_Kg[%d]=%f.\n", k, chainsPrior_Kg[k]);
            mexErrMsgTxt("initial state prior fault");
        } else if (chainsPrior_Kg[k] == 0)
            log_chainsPrior_Kg[k] = logeps;
        else
            log_chainsPrior_Kg[k] = log(chainsPrior_Kg[k]);
    }
    
    /* If no pseudocounts were supplied, we fill up the values with zeros. */
    if (4<=nrhs)
        xi_pseudocounts_KK = mxGetPr(prhs[3]);
    else {
        xi_pseudocounts_KK = mxMalloc(K*K*sizeof(double));
        for (k=0; k<K*K; ++k)
            xi_pseudocounts_KK[k] = 0;
    }
    

    /*  assign the memory for the variables in which we'll construct
        the arrays of values that we want to output
     */
    hM_KTg = mxMalloc(K*T*sizeof(double));
    xi_KKBg = mxMalloc(K*K*B*sizeof(double));
    updated_chainsPrior_Kg = mxMalloc(K*sizeof(double));
    updated_A_KKg = mxMalloc(K*K*sizeof(double));

    /* DEBUGGING
    for(k=0;k<K;++k)
        mexPrintf("log_chainsPrior_Kg[%d]=%f\n", k,log_chainsPrior_Kg[k]);
    for(k=0;k<K*K;++k)
        mexPrintf("log_A_KKg[%d]=%f\n", k,log_A_KKg[k]);
    */

    /* Drowned in all the code, this is the line that makes the call to the
     * fwd_back code that does the work.
     */
    
    fwd_back(hM_KTg, xi_KKBg, log_chainsPrior_Kg, log_A_KKg, 0, R_KTg, K, T, 1);
    /* The C function spits out values in log form so we have to exponentiate
     * and then normalize them. We could have used intermediary variables, but
     * this is C so we're being greedy with memory.
     */
    
    for(k=0;k<K*T;++k)
        hM_KTg[k] = exp(hM_KTg[k]);
    normalizeColumns(hM_KTg, K, T, 0);
    for(k=0;k<K*K*B;++k)
        xi_KKBg[k] = exp(xi_KKBg[k]);
    normalizeColumns(xi_KKBg, K*K, B, 0);
    
    
    
    /* find the MLE for A_KKg using the twoslice marginals in xi_KKBg */
    if (nlhs >= 2) {
        /* the T-1 is because the last values are junk */
        transition_matrix_MLE_from_twoslice_marginals_pseudocounts_normalized(updated_A_KKg, xi_KKBg, K, T-1, xi_pseudocounts_KK);
    }

    /* It's a bit pointless to have another variable for that, but I'm not
     *  completely sure right now if the "pi" vector for HMM corresponds to
     *  the first filtered states or it's not that at all.
     */
    for (k=0;k<K;++k)
        updated_chainsPrior_Kg[k] = hM_KTg[k];
    
    /* Now comes the divergence contributions, the complicated part.
     * I'll return 3 values. The first term without the "S" factor in front,
     * the hM log(hM) entropy and then the result of
     *      S*(first term) - (second term)
     */
    if (nlhs >= 5) {
        
        /* initialization */
        
        divergenceContributions = mxMalloc(2*sizeof(double));
        
        rownormalized_xi_KKBg = mxMalloc(K*K*B*sizeof(double));
        memcpy(rownormalized_xi_KKBg, xi_KKBg, K*K*B*sizeof(double));
        for (t=0; t<T-1; ++t)
            normalizeRows(rownormalized_xi_KKBg + K*K*t, K, K, 1e-16);
        
        divergenceContributions[0] = 0;
        divergenceContributions[1] = 0;

        
        /* on R_KTg */
        for (t=0;t<T;++t) {
            for (k=0; k<K; ++k) {
                divergenceContributions[0] += hM_KTg[k+K*t] * R_KTg[k+K*t];
            }
        }
        
        /* log f(M|A,pi) */
        for (k=0;k<K;++k) {
            divergenceContributions[0] += hM_KTg[k + K*0] * log_chainsPrior_Kg[k];
            divergenceContributions[1] += hM_KTg[k + K*0] * log(hM_KTg[k + K*0] + 1e-16);
        }
            
        for (t=0;t<T-1;++t) {
            for (k=0; k<K; ++k) {
                for (j=0; j<K; ++j) {
                    /* newer, maybe right ? */
                    divergenceContributions[0] += hM_KTg[k + K*t] * rownormalized_xi_KKBg[k+K*j+K*K*t] * log_A_KKg[k+K*j];
                    divergenceContributions[1] += hM_KTg[k + K*t] * rownormalized_xi_KKBg[k+K*j+K*K*t] * log(rownormalized_xi_KKBg[k+K*j+K*K*t] + 1e-16);
                    
                    /* older, probably wrong
                    divergenceContributions[0] += xi_KKTg[k+K*j+K*K*t] * log_A_KKg[k+K*j];
                    divergenceContributions[1] += xi_KKTg[k+K*j+K*K*t] * log(rownormalized_xi_KKTg[k+K*j+K*K*t] + 1e-16);
                     */
                }
            }
        }

        /* the pseudocounts. Don't forget to call the function 
         * hmmmix_pseudoCounts_logNormalizingConstant(A,pseudoCounts)
         * in Matlab after to get the correct divergence term.
         */
        for (k=0; k<K; ++k) {
            for (j=0; j<K; ++j) {
                divergenceContributions[0] += xi_pseudocounts_KK[k+K*j] * log_A_KKg[k+K*j];
            }
        }
        
    }
    
    
    /* output to Matlab */
    
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(K,T,mxREAL);
        outputToolPtr = mxGetPr(plhs[0]);
        memcpy(outputToolPtr, hM_KTg, K*T*sizeof(double));
    }

    if (nlhs >= 2) {
        trplarr_dims[0] = K;
        trplarr_dims[1] = K;
        trplarr_dims[2] = T-1;
        
        plhs[1] = mxCreateNumericArray(trplarr_ndim, trplarr_dims, mxDOUBLE_CLASS, mxREAL);
        outputToolPtr = mxGetPr(plhs[1]);
        memcpy(outputToolPtr, xi_KKBg, K*K*B*sizeof(double));
    }

    if (nlhs >= 3) {
        plhs[2] = mxCreateDoubleMatrix(K,1,mxREAL);
        outputToolPtr = mxGetPr(plhs[2]);
        memcpy(outputToolPtr, updated_chainsPrior_Kg, K*sizeof(double));
    }

    if (nlhs >= 4) {
        plhs[3] = mxCreateDoubleMatrix(K,K,mxREAL);
        outputToolPtr = mxGetPr(plhs[3]);
        memcpy(outputToolPtr, updated_A_KKg, K*K*sizeof(double));
    }
    
    if (nlhs >= 5) {
        plhs[4] = mxCreateDoubleMatrix(2,1,mxREAL);
        outputToolPtr = mxGetPr(plhs[4]);
        memcpy(outputToolPtr, divergenceContributions, 2*sizeof(double));
    }
    
    /* cleaning up */
    
    mxFree(trplarr_dims);
    mxFree(hM_KTg);  mxFree(xi_KKBg);
    mxFree(updated_A_KKg);  mxFree(updated_chainsPrior_Kg);
    mxFree(log_A_KKg); mxFree(log_chainsPrior_Kg);
    
    /* because we allocate xi_pseudocounts_KK when (4<=nrhs) fails */
    if (4>nrhs)
        mxFree(xi_pseudocounts_KK);

        
    if (nlhs >= 5) {
        mxFree(divergenceContributions);
        mxFree(rownormalized_xi_KKBg);
    }
    
        
    return;
}
