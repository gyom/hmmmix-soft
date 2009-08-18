#include <math.h>
#include "mex.h"
#include <string.h>

/*
 *    (A, ind_initstate, rand(1,T)) -> (hM_KTg, indexM_1Tg)
 *
 *  A = normalize(ones(3) + eye(3));
 *  [E1,E2] = hmmmix_generateData_hiddenChains_helper_MatlabC(A, [1], rand(1,100));
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* All the indices passed to this function have to be from 0. Same goes
     * for the array of indices returned from this function. Essentially, this
     * just means that we'd like to have that "ind_initstate" value shifted
     * by -1 calling this function from Matlab.
     * It's just a matter of convention.
     */
    double * A_KKg, *randomHelper_T;
    
    double *ind_initstate;

    double * hM_KTg, * indexM_1Tg;
    double accum;
    
    double * outputToolPtr;
    
    int K,T,k,t;
    
    A_KKg = mxGetPr(prhs[0]);
    ind_initstate = mxGetPr(prhs[1]);
    randomHelper_T = mxGetPr(prhs[2]);

    K = mxGetM(prhs[0]);
    T = mxGetN(prhs[2]);

    hM_KTg = mxMalloc(K*T*sizeof(double));
    indexM_1Tg = mxMalloc(T*sizeof(double));
    
    /* setting the initial state */
    t=0;
    for (k=0;k<K;++k)
        hM_KTg[k + K*t] = 0;
    
    hM_KTg[(int)round(ind_initstate[0]) + K*t] = 1;
    indexM_1Tg[t] = (int)round(ind_initstate[0]);

    /* let's get rid of the task of putting zeros everywhere */
    for (t=1;t<T;++t)
        for (k=0;k<K;++k)
            hM_KTg[k+K*t] = 0;
    
    /* iterate over the chain */
    for (t=1;t<T;++t) {
        /* find the transition state */
        accum = 0;
        
        for (k=0;k<K;++k) {
            /* cumsumming over the row hM_Ktg[t-1]
               This could be done before calling the function. I know.
             */
            accum += A_KKg[(int)round(indexM_1Tg[t-1]) + K*k];
            /* see when the probabilities exceed the cumsum, or we reach the last value */
            if ((accum >= randomHelper_T[t]) | (k==K-1)) {
                indexM_1Tg[t] = k;
                hM_KTg[k+K*t] = 1;
                break;
            }
        }
    }
    
    
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(K,T,mxREAL);
    	outputToolPtr = mxGetPr(plhs[0]);
        memcpy(outputToolPtr, hM_KTg, K*T*sizeof(double));
        /*
        for (t=0;t<T;++t)
            for (k=0;k<K;++k)
                outputToolPtr[k+K*t] = hM_KTg[k+K*t];
         */
    }

    if (nlhs >= 2) {
        plhs[1] = mxCreateDoubleMatrix(1,T,mxREAL);
    	outputToolPtr = mxGetPr(plhs[1]);
        memcpy(outputToolPtr, indexM_1Tg, T*sizeof(double));
        /*
        for (t=0;t<T;++t)
            for (k=0;k<K;++k)
                outputToolPtr[k+K*t] = indexM_1Tg[k+K*t];
         */
    }    
    
    mxFree(hM_KTg); mxFree(indexM_1Tg);
    
    return;
}

