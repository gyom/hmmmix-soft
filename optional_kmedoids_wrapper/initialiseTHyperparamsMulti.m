function [mu_0,lambda_0,eta,m,gamma,S] = initialiseTHyperparamsMulti(Y,f,K, dataSt, initValuesSt)

    % Guillaume's hack of Sohrab's HmmMix code.

    [numS,N] = size(Y);
    mu_0 = zeros(numS,K);
    lambda_0 = zeros(numS,K);
    eta = zeros(numS,K);
    m = zeros(numS,K);
    gamma = zeros(numS,K);
    S = zeros(numS,K);
    
    
    for s = 1:numS
        [mu_0(s,:),lambda_0(s,:),eta(s,:),m(s,:),gamma(s,:),S(s,:)] = ...
            initialiseTHyperparams(Y(s,:),f,K, dataSt, initValuesSt);
    end
    
    