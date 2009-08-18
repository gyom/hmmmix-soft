function [pM M pZ G piG dirPrior loglikEnd logpostEnd mus lambdas Am] = hmmMix(Y,chrIndex,Z,pZ,G_init, ...
    M_init,alpha_prior,G_prior,kappa,nu,mus,lambdas,eta,m,gamma,S,beta, ...
    fixedZ,markovM,maxiter, dataSt, initValuesSt)
%     function [G theta M] = clusteringDMM(Z,G,
% clustering model using iid Dirichlet mixture priors
% Z(p,t)        = discrete state in {1..K} of patient p, probe t (P-by-N)
% G_init        = initial cluster membership vector 1-by-P
% M_init(g,t)   = initial dirichlet mixture component for group g, probe t
% alpha_prior(c,:,t) = non-stationary dirichlet prior for component c of
%                      the mixture for probe t
% M_prior(g,k)  = dirichlet prior on the dirichlet mixture weights for group g
% G_prior(g)    = dirichlet prior on the group membership
%
%This model assumes that Z's are fixed

% Guillaume's hack of Sohrab's HmmMix code.

[P,K,N] = size(pZ);
[numG,N] = size(M_init);
[K,C] = size(alpha_prior);

piG = zeros(size(G_prior));
% initialise the stationary distributions to their priors
epsilon = 0.75;
pG = ones(numG,P)*((1-epsilon)/(numG-1));
pM = ones(numG,C,N)*((1-epsilon)/(C-1));
M = M_init;
for p=1:P
    pG(G_init(p),p)=epsilon;
end
for g = 1:numG
    piG(g) = sum(G_init==g);
%     for n=1:N
%         pM(g,M_init(g,n),n)=epsilon;
%     end
end
piG = normalize(piG);

A=[0.9 0.05 0.05; 0.05 0.9 0.05; 0.05 0.05 0.9];
Am = zeros(numG,C,C);
for g=1:numG
    Am(g,:,:)=A;
end

%piM = [0.1 0.8 0.1];
piM = ones(1,K)/K; % gyom
piM = reshape(repmat(piM,[numG,1]),[numG,K]);
piM = reshape(repmat(piM,[1,22]),[numG,K,22]);
%piM_prior = [1 8 1];
%Am_prior = [100 1 1; 1 100 1; 1 1 100];
piM_prior = [K K K]; % gyom
Am_prior = initValuesSt.pseudoCounts; %gyom

chrs = unique(chrIndex);
localEv = [];
if ~fixedZ
    localEv = tdistPDFMult(Y,mus,lambdas,nu);
else
    localEv = calls2Probs(Z,K);
end
dirPrior = preComputeDirPrior(alpha_prior);

% START EM
%
%
iter = 1;
converged = 0;
maxll = -Inf;
% calculate local evidence for raw data
logpost = -inf(1,maxiter*length(beta));
loglik = inf(1,maxiter*length(beta));
iindex=2;
cll = 0;
G = G_init;
% if isempty(T)
%     plotPatientsAndProfilesDMM(pZ,G,pM,chrIndex);
% else
%     plotPatientsAndProfilesWithTruth(pZ,G,pM,chrIndex,T);
% end

%plotPatientsAndProfilesWithTruth(pZ,G,pM,chrIndex,T);
for B=1:length(beta)
    iter=1;
    while iter<=maxiter && ~converged
        disp(['HMM-MIX iteration ', int2str(iindex-1),' log posterior: ', num2str(logpost(iindex-1))]);
        logPrior = 0;
        noisePrior = 0;
        aPrior = 0;
        % E-step
        %pG = updateGDmm(pM,pZ,piG,alpha_prior,beta(B));

        % update the clusters
        pG = updateGHmmMix(localEv,dirPrior,M,piG);
        [jnk G] = max(pG);
        piG = normalise(sum(pG,2)+G_prior);

        % E-step
        b = hmmMixEmission(localEv,G,dirPrior,numG);
        if markovM
            [pM,M,Am,ll] = updateMMarkovDmm(b,chrIndex,piM,Am,Am_prior);
        else
            ll = 0;
            pM = b;
            [jnk M] = max(pM,[],2); 
        end
        loglik(iindex) = ll;
        aPrior = 0;
        for g = 1:numG
            for s=1:C
                aPrior = aPrior + ...
                    log(dirichletpdf(reshape(Am(g,s,:),[C,1]),Am_prior(s,:)'));
            end
        end 
        
        if ~fixedZ
            pZ = updateZHmmMix(G,M,dirPrior,localEv);
            [jnk,zi] = max(pZ,[],2); 
            [P,N]
            Z = reshape(zi,[P,N]);
            rhoZ = calls2Probs(Z,K);
            % reestimate noise parameters
            for p = 1:P
                rho = reshape(rhoZ(p,:,:),[K,N]);
                [mus(p,:), lambdas(p,:)] = estimateTNoiseParamsMap(Y(p,:),mus(p,:), ...
                    lambdas(p,:),nu,rho,eta(p,:),m(p,:),gamma(p,:),S(p,:),kappa);
            end
            % compute log prior
            stddev = ones(size(eta))./(sqrt(eta));
            noisePrior = sum(sum(log(gampdf(lambdas,gamma,ones(size(S))./S)))) + ...
                 sum(sum(log(normpdf(mus,m,stddev))));
            
           
            localEv = tdistPDFMult(Y,mus,lambdas,nu);

            %reestimate background distribution
            alpha = updateBackground(pZ,G,M);
            alpha_prior(:,2) = alpha;
            dirPrior = preComputeDirPrior(alpha_prior);
        else
            %reestimate background distribution
            alpha = updateBackground(localEv,G,M);
            alpha_prior(:,2) = alpha;
            dirPrior = preComputeDirPrior(alpha_prior);
        end
        
        gPrior =  log(dirichletpdf(piG,G_prior));
        postProb = loglik(iindex)+ aPrior + noisePrior + gPrior;
%         disp(['APrior: ', num2str(aPrior)]);
%         disp(['Noise Prior: ',num2str(noisePrior)]);
%         disp(['GPrior: ',num2str(gPrior)]);
        logpost(iindex) = postProb;
        
        % check for convergence
        if (postProb<=logpost(iindex-1))
            %converged=1;
        end
        iindex=iindex+1;
        iter=iter+1;

        %  figure; imagesc(pG);
%         if isempty(T)
%             plotPatientsAndProfilesDMM(pZ,G,pM,chrIndex);
%         else
%             plotPatientsAndProfilesWithTruth(pZ,G,pM,chrIndex,T);
%         end
    end
end
loglikEnd = ll;
logpostEnd = postProb;
% figure; plot(loglik);hold on; plot(logpost,'-r');hold off;
% legend({'log-lik','log-posterior'});


end