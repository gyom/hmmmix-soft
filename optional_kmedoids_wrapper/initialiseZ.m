function [mus lambdas eta m gamma S Z pZ] = initialiseZ(Y,chrIndex,f,K,kappa,nu,maxiter, dataSt, initValuesSt)
% function [mus lambdas eta m gamma S Z] = initialiseZ(Y,chrIndex,f,K,kappa,nu,maxiter)
% Initialised the hidden states Z for each corresponding observed logratio
% Y
% Y(i,j) is the log ratio for patient i, probe j
% Z(i,j) is the hidden state in 1:K for patient i, probe j
% 
% Z and the model parameters mus and lambdas are estimated by a mixture of
% Student-t distribution HMM
% 
% Also returned are the hyperparamters eta, m, gamma, S

% Guillaume's hack of Sohrab's HmmMix code.

[numS, N] = size(Y);
Z = zeros(numS,N);
pZ = zeros(numS,K,N);
% hyperparameters
eta = zeros(numS,K);
m = zeros(numS,K);
gamma = zeros(numS,K);
S = zeros(numS,K);

% model parameters
mus = zeros(numS,K);
lambdas = zeros(numS,K);
for i=1:numS
  disp(['Initialising patient ', int2str(i), ' of ', int2str(numS)]);  
  y = Y(i,:);%-nanmedian(Y(i,:));
  [mu_0,lambda_0,eta(i,:),m(i,:),gamma(i,:),S(i,:)] = initialiseTHyperparams(y,f,K, dataSt, initValuesSt);
  
  [musAll lambdasAll pi loglik Z(i,:) pZ(i,:,:)] = ...
        THMMAcghViterbi(y,chrIndex,mu_0,lambda_0,nu,eta(i,:),m(i,:),gamma(i,:),S(i,:),kappa,maxiter);
  mus(i,:) = musAll(:,end);
  lambdas(i,:) = lambdasAll(:,end); 
end



