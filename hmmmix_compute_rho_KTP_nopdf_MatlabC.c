#include "mex.h"
#include <string.h>
#include <math.h>
#include "hmmmix_common.h"

/* a C version of
 * function rho_KTP = hmmmix_compute_rho_KTP(hM_KTG, hC, Y_PT, mu_KP, lambda_KP, alpha_K, nu_K)
 * now without useless parameters
 *
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

    double * hM_KTG, *hC, *alpha_K;

    int K, P, T, G;

    int k,p,t,g,j;

    int ind;

    double *S; /* caching some computations*/

    int trplarr_ndim = 3;
    int * trplarr_dims = mxMalloc(trplarr_ndim*sizeof(int));

    double *rho_KTP;

    double * outputToolPtr;

    if (nrhs!= 3 || nlhs != 1)
    	mexErrMsgTxt("hmmmix_optimal_hM_KTG_MatlabC requires 3 inputs and 1 output. If you used 7 elements, you probably included the useless arguments that are no longer wanted.");

    hM_KTG = mxGetPr(prhs[0]);
    hC = mxGetPr(prhs[1]);
    alpha_K = mxGetPr(prhs[2]);

    K = mxGetM(prhs[0]); /* size(hM_KTG,1) */
    G = mxGetM(prhs[1]); /* size(hC,1) */
    P = mxGetN(prhs[1]); /* size(hC,2) */
    T = round(((double)mxGetN(prhs[0]))/(double)G); /* voodoo against roundoff errors */

    rho_KTP = mxMalloc(K*T*P*sizeof(double));
    trplarr_dims[0] = K;
    trplarr_dims[1] = T;
    trplarr_dims[2] = P;

    S = mxMalloc(K*K*sizeof(double));

    for(k=0;k<K;++k) {
        for(j=0;j<K;++j) {
            if (j==k)
                S[j+K*k] = alpha_K[j] / (alpha_K[j] + 2);
            else
                S[j+K*k] = 1 / (alpha_K[j] + 2);
        }
    }

    for (t=0;t<T;++t)
        for (p=0;p<P;++p) {
            for (k=0;k<K;++k) {
                ind = k+K*t + K*T*p;
                rho_KTP[ind] = 0;
                for (g=0;g<G;++g)
                    for (j=0;j<K;++j) {
                        rho_KTP[ind] += hC[g+G*p] * hM_KTG[j + K*t + K*T*g] * S[j+K*k];
                    }
            }
            normalizeColumns(rho_KTP + K*t + K*T*p, K, 1, 0);
    }

    if (nlhs >= 1) {
        plhs[0] = mxCreateNumericArray(trplarr_ndim, trplarr_dims, mxDOUBLE_CLASS, mxREAL);
        outputToolPtr = mxGetPr(plhs[0]);
        memcpy(outputToolPtr, rho_KTP, K*T*P*sizeof(double));
    }

    mxFree(S); mxFree(rho_KTP);
    mxFree(trplarr_dims);
    return;
}





