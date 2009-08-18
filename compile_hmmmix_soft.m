% Written by Guillaume Alain, summer of 2009.
%
% I usually compile everything by going into the appropriate directory and
% then running all these commands in succession. There are a few commands
% in the last portion that compile functions that are not really necessary
% unless we want to use the hack that overrides some of Sohrab Shah's
% functions from "CNAhmmer".
%
% No harm is done if we don't add the optional_kmedoids_wrapper directory
% to the Matlab path.

% This is the forwards-backwards algorithm, working with logarithms for all
% the quantities. The "fwd_back_MatlabC" file provides a way to call
% fwd_back directly from Matlab if we want. Otherwise, the fwd_back code is
% linked later on with the hmmmix_frugal_hM_KTg_MatlabC code.
mex -c fwd_back.c
mex -c fwd_back_MatlabC.c
mex fwd_back_MatlabC.o fwd_back.o

% Provides a few interesting functions for many of the other C files that
% want to evaluate things like the log-likelihood of the observations.
mex -c hmmmix_common.c

% This is the "frugal" version of the code that computes the updates for
% the hM quantities. The "non-frugal" version computed loglikelihoods by
% itself but that turned out to be a bad way to factor the code. I also
% realized after some time that it was better to separate out the groups
% and not process them all at once. The value of T is so large that we
% don't get penalized for that in terms of speed.
mex -c hmmmix_frugal_hM_KTg_MatlabC.c
mex hmmmix_frugal_hM_KTg_MatlabC.o hmmmix_common.o fwd_back.o

% These two functions are used to generate artificial data faster. There
% were a few problems with non-vectorizable code so these were written.
mex hmmmix_generateData_hiddenChains_helper_MatlabC.c
mex hmmmix_generateData_hiddenPatients_helper_MatlabC.c

mex viterbi_path_MatlabC.c

mex -c hmmmix_compute_rho_KTP_nopdf_MatlabC.c
mex hmmmix_compute_rho_KTP_nopdf_MatlabC.o hmmmix_common.o

fprintf('Now we''re compiling functions that are not really necessary unless we want to use the code for kmedoids.\n')

% I honestly don't remember where repmatC could have been required. I don't
% think it's in my code, but I've had this problem come up one day so I
% decided to include the file in here. 
cd optional_kmedoids_wrapper
mex -c mexutil.c
mex -c repmatC.c
mex repmatC.o mexutil.o

% Functions used for the kmedoids hack with Sohrab Shah's code. These are
% used to speed up the functions that couldn't be vectorized. Since I've
% used many versions of the fwd_back.c code, I had to include here the
% "old" version that behaved exactly like Sohrab's matlab version of the
% 'fwd_back' function. Same for viterbi_path_SSC.
mex fwd_back.c
mex viterbi_path_SSC.c


