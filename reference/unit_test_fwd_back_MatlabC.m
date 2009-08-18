fprintf('Unit test for fwd_back_MatlabC compared to Kevin Murphy''s fwd_back implementation.\n');

K = 3;
T = 10;

init_state_distrib = rand(K,1);
transmat = normalize(rand(K,K),2);
obslik = 50*rand(K,T);

[gamma,alpha,beta,loglik] = fwd_back(init_state_distrib, transmat, obslik);

[logSmoothedSequence, logTwoSliceMarginals] = fwd_back_MatlabC(log(init_state_distrib+eps), log(obslik+eps), log(transmat+eps));
E = normalize(exp(logSmoothedSequence),1);

% There's some error to the order of 1e-7. It's not the +eps that I
% introduce. I don't know what it is and it doesn't really matter.
tolerance = 1e-4;

assert(max(abs(E(:) - gamma(:))) < tolerance);

fprintf('Outputs match those of the reference fwd_back.m\n');

%% Now we're forcing the transitions to be from {state 1 in timestep 4} to
%  {state 2 in timestep 5} using the non-stationary transition matrices. I
%  think the random evidence will be pretty much ignored.

K = 3;
T = 10;

init_state_distrib = rand(K,1);
transmat = cat(3, repmat([1,0,0;1,0,0;1,0,0], [1,1,4]), repmat([0,1,0;0,1,0;0,1,0], [1,1,5]));
obslik = 50*rand(K,T);

[logSmoothedSequence, logTwoSliceMarginals] = fwd_back_MatlabC(log(init_state_distrib+eps), log(obslik+eps), log(transmat+eps));
E = normalize(exp(logSmoothedSequence),1);

assert( max(abs( E(:,5) - [1;0;0])) < tolerance );
assert( max(abs( E(:,6) - [0;1;0])) < tolerance );

fprintf('\tTest passed.\n');