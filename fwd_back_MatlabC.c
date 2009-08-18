#include "mex.h"
#include <string.h>

void fwd_back(double * OUT_logSmoothedMarginals, double * OUT_logTwoSliceMarginals, double * logPrior, double * logTransmat, int transmatStride, double * logObslik, int K, int T, int twoSliceOption);

/*
 *	This function is meant as a wrapper for the above C function.
 *  We are working with logarithms of everything.
 *  NaN values or -Inf are not acceptable.
 *
 *  [logSmoothedSequence, logTwoSliceMarginals] = fwd_back_MatlabC(logPrior, logObslik, logTransitionMatrix)
 *	
 *  where   logPrior is of size K
 *          logObslik is of size [K,T]
 *          logTransitionMatrix is of size [K,K ] or [K,K,T-1]
 *
 *  This implementation has been checked against the previous implementation of fwd_back,
 *  which was itself checked against Kevin's Matlab implementation. The version with
 *  non-stationary matrices hasn't been checked for correctness against a basic implementation,
 *  but the results given are the same when a constant matrix is used, and they are what we
 *  should expect them to be when we provide transition matrices that force transitions at
 *  certain places.
 */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double * logPrior;
    double * logObslik;
    double * logTransitionMatrix;
    
    double * logSmoothedSequence;
    double * logTwoSliceMarginals;
    
    int K, T;
    int nonStationaryTransitionMatrices;
    int transmatStride;
    int computeTwoSliceMarginals;

    /* artificial variables for outputting stuff back to Matlab */
    double * outputToolPtr;
    int trplarr_ndim = 3;
    int * trplarr_dims = mxMalloc(trplarr_ndim*sizeof(int));
    
    if ( !(3==nrhs) || !(1<=nlhs && nlhs<=2))
        mexErrMsgTxt("hmmmix_frugal_hM_KTg_MatlabC requires 3 inputs and 1-2 outputs");
    
    if (2==nlhs)
        computeTwoSliceMarginals = 1;
    else
        computeTwoSliceMarginals = 0;
    
    logPrior = mxGetPr(prhs[0]);
    logObslik = mxGetPr(prhs[1]);
    logTransitionMatrix = mxGetPr(prhs[2]);
    
    K = mxGetM(prhs[0])*mxGetN(prhs[0]);
    T = mxGetN(prhs[1]);
    if(K != mxGetM(prhs[1]))
        mexErrMsgTxt("Wrong dimensions for the second argument : logObslik\n");
    if(K != mxGetM(prhs[2]))
        mexErrMsgTxt("Wrong dimensions for the third argument : logTransitionMatrix\n");
    
    if( K != mxGetN(prhs[2]) ) {
    	if( K*(T-1) != mxGetN(prhs[2]) )
            mexErrMsgTxt("Wrong dimensions for the third argument : logTransitionMatrix\n");
        nonStationaryTransitionMatrices = 1;
        transmatStride = K*K;
    } else {
        nonStationaryTransitionMatrices = 0;
        transmatStride = 0;
    }
    
    logSmoothedSequence = mxMalloc(K*T*sizeof(double));
    /* allocate that memory only if it's going to be used */
    if (computeTwoSliceMarginals == 1)
        logTwoSliceMarginals = mxMalloc(K*K*T*sizeof(double));
    else
        logTwoSliceMarginals = (void *)0;
    
    
    fwd_back(logSmoothedSequence, logTwoSliceMarginals, logPrior, logTransitionMatrix, transmatStride, logObslik, K, T, computeTwoSliceMarginals);
    
    /*
    mexPrintf("logSmoothedSequence[0]=%f,logSmoothedSequence[1]=%f\n",logSmoothedSequence[0],logSmoothedSequence[1]);
    */
    
    if (nlhs >= 1) {
        plhs[0] = mxCreateDoubleMatrix(K,T,mxREAL);
        outputToolPtr = mxGetPr(plhs[0]);
        memcpy(outputToolPtr, logSmoothedSequence, K*T*sizeof(double));
    }
    
	if (nlhs >= 2) {
        trplarr_dims[0] = K;
        trplarr_dims[1] = K;
        trplarr_dims[2] = T-1;
        
        plhs[1] = mxCreateNumericArray(trplarr_ndim, trplarr_dims, mxDOUBLE_CLASS, mxREAL);
        outputToolPtr = mxGetPr(plhs[1]);
        memcpy(outputToolPtr, logTwoSliceMarginals, K*K*(T-1)*sizeof(double));
    }
    mxFree(logSmoothedSequence);
    /* free that variable only if we allocated it */
    if (computeTwoSliceMarginals == 1)
    	mxFree(logTwoSliceMarginals);
    mxFree(trplarr_dims);
    return;
}
