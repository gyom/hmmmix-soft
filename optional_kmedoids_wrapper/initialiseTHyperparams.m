function [mu_0,lambda_0,eta,m,gamma,S] = initialiseTHyperparams(y,f,K, dataSt, initValuesSt)

    % Guillaume's hack of Sohrab's HmmMix code.

    mu_0 = mean(initValuesSt.m_KP,2)';
    lambda_0 = mean(initValuesSt.lambda_KP,2)';
    
    eta = mean(initValuesSt.eta_KP,2)';
    m = mean(initValuesSt.m_KP,2)';
    gamma = mean(initValuesSt.gamma_KP,2)';
    S = mean(initValuesSt.S_KP,2)';

    %[size(mu_0); size(lambda_0); size(eta); size(m); size(gamma); size(S)]
    
end