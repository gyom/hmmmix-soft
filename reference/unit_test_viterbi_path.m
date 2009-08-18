fprintf('Starting unit test for viterbi_path_MatlabC.\n');

% There is a potential problem here because we're just testing against
% Sohrab Shah's implementation in Matlab. However, his Matlab code is much
% more easy to check against Kevin Murphy's formulas in his book.

K = 3;
T = 10;

obslik = rand(K,T);
transmat = normalize(rand(K,K),2);
prior = normalize(rand(K,1));

[path loglik seg] = viterbi_path_SS(prior, transmat, obslik);

[path2 loglik2 seg2] = viterbi_path_MatlabC(log(prior+eps), log(transmat+eps), log(obslik+eps));
path3 = viterbi_path_MatlabC(log(prior+eps), log(transmat+eps), log(obslik+eps));

% Unfortunately, we don't have a base implementation of Viterbi that we
% trust with non-stationary matrices.
[path4 loglik4 seg4] = viterbi_path_MatlabC(log(prior+eps), repmat(log(transmat+eps), [1,1,T-1]), log(obslik+eps));

% whatever you think is fine
tolerance = 1e-8;

assert( max(abs(path(:) - path2(:))) < tolerance );
assert( max(abs(loglik(:) - loglik2(:))) < tolerance );
assert( max(max(abs(seg(:,1:3) - seg2(:,1:3)))) < tolerance );

assert( max(abs(path(:) - path3(:))) < tolerance );

assert( max(abs(path(:) - path4(:))) < tolerance );
assert( max(abs(loglik(:) - loglik4(:))) < tolerance );
assert( max(max(abs(seg(:,1:3) - seg4(:,1:3)))) < tolerance );

fprintf('\tUnit test for viterbi_path_MatlabC passed.\n');