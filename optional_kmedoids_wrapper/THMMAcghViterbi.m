function [mus lambdas pi loglik Z rho segs] = ...
    THMMAcghViterbi(Y,chrIndex,mu_0,lambda_0,nu,eta,m,gamma,S,kappa,maxiter)
    % EM algorithm for Student-t HMM MAP estimate
    %
    % function [mu lambda pi loglik] = THMMAcgh(Y,mu_0,lambda_0,nu,maxiter)
    %
    % Y          - data
    % mu_0       - initial value for Student-t means (1-by-K)
    % lambda_0   - initial values for Student-t precision (1-by-K)
    % maxiter    - maximum number of EM iterations
    % nu         - fixed values for Student-t dof (1-by-K)
    % eta        - prior precision on mu
    % m          - prior mean on mu
    % gamma      - prior shape on lambda (Gamma(shape,scale))
    % S          - prior scale on lambda (Gamma(shape,scale))
    % N          - size of data
    % K          - number of discrete states

K = length(mu_0);
N = length(Y);
py = zeros(K,N);                 % local evidence
mus = zeros(K,maxiter);          % state means
lambdas = zeros(K,maxiter);       % state variances
pi = kappa;               % initial state distribution   
converged = 0;                   % flag for convergence
Nk = zeros(1,K);                 % number of Zs in each state
Z = zeros(N,1);
loglik = zeros(1,maxiter);
% SET UP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up the chromosome indicies
% make cell array of chromosome indicies
chrs = unique(chrIndex);
chrsI = cell(1,length(chrs));

piZ = cell(1,length(chrs)); % need 1 initial dist per chromosome
% initialise the chromosome index and the init state distributions 
for c = 1:length(chrs)
    chrsI{c} = find(chrIndex == chrs(c));
    piZ{c} = ones(1,K)/K;
end


% INITIALISATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i = 1;
mus(:,1) = mu_0;
lambdas(:,1) = lambda_0;
% re-calculate the likelihood
for k=1:K
    py(k,:)=tdistPDF(Y,mus(k,i),lambdas(k,i),nu(k));
end
% initialise transition matrix to the prior:
strength=100;
e = 0.999;
A = zeros(K);
for j = 1:K
    A(j,:) = (1-e)/(K-1);
    A(j,j) = e;
end
A_prior = A;%[0.999 0.00005 0.00005; (1e-10)/2 1-1e-10 (1e-10)/2; 0.00005 0.00005 0.999];
dirPrior = A*strength;
pi_prior = kappa;
loglik(i) = -Inf;

% Expectation Maximization
rho = zeros(K,N);    % marginal p(Z_t|Y)
xi  = zeros(K,K,N);  % marginal p(Z_t,Z_{t-1}|Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (~converged && (i < maxiter))
    disp(['THMMAcghViterbi: EM iteration:', int2str(i), ' loglik: ',num2str(loglik(i))]);
    i = i+1;
    % E-step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for c = 1:length(chrsI)
        piZ{c} = pi;
        I = chrsI{c};
       % [Z(I) ll] = viterbi_path_SS(piZ{c}, A, py(:,I));
        [rho(:,I),alpha,beta,xi(:,:,I),ll] = fwd_back(piZ{c}, A, py(:,I));       
        loglik(i) = loglik(i)+ll;
    end
    Zcounts = reshape(sum(xi,3),[K,K]);
%    disp(Zcounts);
    % M-step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % update the noise hyperparams
    mu_i = reshape(mus(:,i-1),[1,K]);
    lambda_i = reshape(lambdas(:,i-1),[1,K]);
    [mus(:,i),lambdas(:,i), pi, u] = ...
        estimateTNoiseParamsMap(Y,mu_i,lambda_i,nu,rho,eta,m,gamma,S,kappa);

    % re-calculate the likelihood
    for k=1:K
        py(k,:)=tdistPDF(Y,mus(k,i),lambdas(k,i),nu(k));
    end
 
    priorA = 0;
    % Update transition matrix A
    for k = 1:K
        A(k,:) = Zcounts(k,:)+dirPrior(k,:);
        A(k,:) = normalise(A(k,:));
        priorA = priorA + log(dirichletpdf(A_prior(k,:),A(k,:)));
    end

    % compute log-likelihood and check convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k = 1:K
        % calculate the p(mu,sigma) = the prior of the normal-Gamma evaluated
        % at the estimated quantities.  See Bishop : pattern recognition and
        % machine learning, pg 101 for how to calculate this
       priorMu(k)      = log(normpdf(mus(k,i),mu_0(k),1));
  
       %priorLambda(k)   = log(gampdf(lambdas(k,i),kappa(k),S(k)));
    end
%    disp(sum(priorMu));
%    disp(sum(priorSigma));
%    disp(priorA);
    loglik(i) = loglik(i) + priorA + sum(priorMu);% + sum(priorSigma);
    if approxeq(loglik(i),loglik(i-1),1e-1) || (loglik(i) < loglik(i-1))
        converged = 1;
    end
end
if converged
    i = i-1;
end
segs = cell(1,length(chrs));
% do the final viterbi pass
for c = 1:length(chrsI)
    piZ{c} = pi;
    I = chrsI{c};
    [Z(I) ll segs{c}] = viterbi_path_SSC(piZ{c}, A, py(:,I));
end
disp('Transition matrix:');
disp(A);
disp(ll);
disp(i);

mus = mus(:,1:i);
lambdas = lambdas(:,1:i);
loglik = loglik(1:i);
