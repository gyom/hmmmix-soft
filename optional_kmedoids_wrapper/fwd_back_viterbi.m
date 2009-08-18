function [path,gamma,eta,loglik] = fwd_back_viterbi(init_state_distrib,transmat,obslik)

    % My way of circumventing the calls to Sohrab's slow implementation of
    % fwd_back_viterbi. My version of fwd_back wasn't returning the gamma
    % values but, after inspecting Sohrab's code, that values is just
    % returned by 'hmmMix' as first output argument and I know for sure that
    % I'm not using that output. I can then safely put junk in 'gamma' for
    % this hack.
    
    [path loglik junk] = viterbi_path_SSC(init_state_distrib, transmat, obslik);
    
    eps = 10e-200;
    [logSmoothedSequence, logTwoSliceMarginals] = fwd_back_MatlabC(log(init_state_distrib+eps), log(obslik+eps), log(transmat+eps));
    %E = normalize(exp(logSmoothedSequence),1);

    gamma = zeros(size(logSmoothedSequence));

    E = exp(logTwoSliceMarginals);
    
    [K,K0,B] = size(E);
    T = B+1;
    assert(K==K0);
    
    F = reshape(E, [K*K,B]);
    sliceSums = sum(F,1);

    eta0 = E ./ repmat(reshape(sliceSums, [1,1,B]), [K,K,1]);
    eta = cat(3, eta0, zeros(K,K,1));


% check the correctness with this :
%     K = 3;
%     T = 10;
% 
%     obslik = rand(K,T);
%     transmat = normalize(rand(K,K),2);
%     init_state_distrib = normalize(rand(K,1));
% 
%     [path,gamma,eta,loglik] = fwd_back_viterbi(init_state_distrib,transmat,obslik);
%       % CHANGE PATH IN-BETWEEN TO CALL THE OTHER VERSION
%     [path0,gamma0,eta0,loglik0] = fwd_back_viterbi(init_state_distrib,transmat,obslik);
    
end