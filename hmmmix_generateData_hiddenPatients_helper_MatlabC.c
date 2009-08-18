#include <math.h>
#include "mex.h"
#include <string.h>

/*
 *    (indexM_1Tg, alphaK, rand(1,T)) -> (Z_KTp, indexZ_1Tp)
 *
 *  A = normalize(ones(3) + eye(3));
 *  [E1,E2] = hmmmix_generateData_hiddenChains_helper_MatlabC(A, [1], rand(1,100));
 *  alphaK = [4,4,4];
 *  [F1, F2] = hmmmix_generateData_hiddenPatients_helper_MatlabC(E2, alphaK, rand(1,100));
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* All the indices passed to this function have to be from 0. Same goes
     * for the array of indices returned from this function.
     * It's just a matter of convention.
     */

    double * indexM_1Tg, * randomHelper_T, *alphaK;
    double * Z_KTp, *indexZ_1Tp;
    
    double alpha;
    
    double * outputToolPtr;
    
    int K,T,k,t;
    
    indexM_1Tg = mxGetPr(prhs[0]);
    alphaK = mxGetPr(prhs[1]);
    randomHelper_T = mxGetPr(prhs[2]);

    K = mxGetM(prhs[1])*mxGetN(prhs[1]);
    T = mxGetM(prhs[2])*mxGetN(prhs[2]);

    Z_KTp = mxMalloc(K*T*sizeof(double));
    indexZ_1Tp = mxMalloc(T*sizeof(double));
    
    for (t=0;t<T;++t) {
        for (k=0;k<K;++k)
            Z_KTp[k+K*t] = 0;
        
        alpha = alphaK[(int)round(indexM_1Tg[t])];
        
        /* check to see if we don't need any offset from the state of the chain */
        if (randomHelper_T[t] <=  alpha / (alpha + (K-1))) {
            indexZ_1Tp[t] = indexM_1Tg[t];
            Z_KTp[ (int)round(indexZ_1Tp[t]) +K*t] = 1;
        } else {
            /* Otherwise we're going for another offset.
             * This is the "hey I can't do math" approach, but K is so small
             * that I'm fine doing this to make the code clearer.
             */
            for (k=1;k<K;++k) {
                /* if we're there, or if we exceed the last index (shouldn't happen unless something is wrong */
                if ((randomHelper_T[t] <= (alpha + k) / (alpha + K-1)) || (k==K-1)) {
                    indexZ_1Tp[t] = indexM_1Tg[t] + k;
                    /* avoid modulo nightmares that could crash the thing if we
                       got an out of bound index at the end
                     */
                    if (indexZ_1Tp[t] >= K)
                        indexZ_1Tp[t] -= K;
                    /* no need for a +k here because we're using indexZ_1Tp[t] */
                    Z_KTp[ (int)round(indexZ_1Tp[t]) + K*t] = 1;
                    break;
                }
            }
        }
    }
    
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(K,T,mxREAL);
    	outputToolPtr = mxGetPr(plhs[0]);
        memcpy(outputToolPtr, Z_KTp, K*T*sizeof(double));
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(1,T,mxREAL);
    	outputToolPtr = mxGetPr(plhs[1]);
        memcpy(outputToolPtr, indexZ_1Tp, T*sizeof(double));
    }    
    
    mxFree(Z_KTp); mxFree(indexZ_1Tp);
    
    return;
}

