fprintf('Unit test for hmmmix_frugal_hM_KTg_MatlabC.\n');

% We already compared hmmmix_frugal_hM_KTg_MatlabC with the previous
% version written in C to make sure that the values not tested here are the
% same.

K = 3;
T = 9;

R_KTg = rand(K,T);
chainsPrior_Kg = normalize(rand(K,1));
A_KKg = normalize(rand(K,K),2);
xi_pseudocounts_KK = ones(K,K);


[hM_KTg, xi_KKBg, updated_chainsPrior_Kg, updated_A_KKg, divergenceContributions] = hmmmix_frugal_hM_KTg_MatlabC(R_KTg, chainsPrior_Kg, A_KKg, xi_pseudocounts_KK);

tolerance = 1e-8;

[logSmoothedSequence, logTwoSliceMarginals] = fwd_back_MatlabC(log(chainsPrior_Kg), R_KTg, log(A_KKg));
E = normalize(exp(logSmoothedSequence),1);
F = exp(logTwoSliceMarginals);
for t=1:T-1
    F(:,:,t) = normalize(F(:,:,t));
end
assert( max(abs(E(:) - hM_KTg(:)))  <  tolerance);
assert( max(abs(F(:) - xi_KKBg(:)))  <  tolerance);

fprintf('\tTest passed.\n');