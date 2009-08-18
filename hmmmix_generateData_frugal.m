function dataSt = hmmmix_generateData_frugal(G, P, T, K)

    % G : number of groups
    % P : number of patients
    % T : length of chains
    % K : number of components (centers)
    
    % the data is single-dimensional

    
    % This will be a hard assignment, and I'm using 1:G encoding instead of
    % having G binary values for every patient.
    pi_multinomial_parameters = ones(1,G)/G;
    [patientsAssignedCohort_1K, patientsAssignedCohort_bin] = sample_multinomial(pi_multinomial_parameters, P);
    patientsAssignedCohort = toRow(patientsAssignedCohort_1K);
    
    transitionMatrix_stabilityFactor = T/4;%100; % arbitrary
    transitionMatrix_prior = ones(K,K) + transitionMatrix_stabilityFactor*eye(K);

    A_KKG = zeros(K,K,G);
    for g=1:G
        for k=1:K
            A_KKG(k,:,g) = sample_dirichlet(transitionMatrix_prior(k,:),1);
        end
    end
    
    % simulating the chains M(1:K,1:T) put inside a cell array 1:G
    pi_M_multinomial_parameters = ones(1,K)/K; %normalize([1,10,1]);
    
    alpha_stretch = 10; % arbitrary
    alpha_K = alpha_stretch*ones(1,K);

    %% hidden chains
    hM_KTG = nan(K,T,G);
    indexM_1TG = nan(1,T,G);
    for g=1:G
        [init_index, init_junk] = sample_multinomial(pi_M_multinomial_parameters,1);
        [hM_KTg, indexM_1Tg] = hmmmix_generateData_hiddenChains_helper_MatlabC(squeeze(A_KKG(:,:,g)), init_index-1, rand(1,T));
        hM_KTG(:,:,g) = hM_KTg;
        indexM_1TG(1,:,g) = indexM_1Tg;
    end
    
    
    % Sample the parameters mu, lambda, m, eta for the patients.
    %
    % There is a good deal of confusion around the role of S_K(1:K) because
    % of the two ways to define the gamma pdf. Matlab decided to go with
    % the exp(t/b) and I decided to go with the exp(b*t) to be more
    % compatible with Sohrab's existing code. As a result of this, we have
    % to call the Matlab's gamrnd function with the second argument
    % inverted.
    
    %m_K = sort(trnd(10, K, 1)) + linspace(-1,1,K)';
    eta_K = 4+10*rand(K,1);
    gamma_K = 4+10*rand(K,1); % no idea
    S_K = 4+10*rand(K,1); %0.1 * rand(K,1); % % no idea
    nu_KP = 4+repmat(10*rand(K,1), [1,P]); % no idea

    m_K = (5*abs(trnd(10, K, 1))+3).*[-1;0;1]; %sort(trnd(10, K, 1)) + linspace(-1,1,K)';
    %eta_K = 5000*ones(K,1);%4+10*rand(K,1);
    %gamma_K = ones(K,1); %4+10*rand(K,1); % no idea
    %S_K = 0.0001*ones(K,1); %1000*rand(K,1); %4+10*rand(K,1); % no idea
    %nu_KP = 4+repmat(10*rand(K,1), [1,P]); % no idea

    
    %tau_hidden_param = zeros(K,P);    
    lambda_KP = zeros(K,P);
    mu_KP = zeros(K,P);
    for p=1:P
        for k=1:K
            % 1/S_K has to be the second parameter because it's the one that
            % finds its way as S_K in the later formulas that estimate
            % the parameters. Matlab use exp(-x/b) and I work with
            % exp(-b*x). This is why I put the inverse 1/S_K here.
            lambda_KP(k,p) = gamrnd(gamma_K(k), 1/S_K(k));
            
            % Be careful about this eta_K here. Sohrab's paper would have
            % this line be sqrt(eta_K(k)/lambda_KP(k,p)) but Archambeau's
            % thesis and other sources would have it be
            % sqrt(1/(eta_K(k)*lambda_KP(k,p))), like here.
            mu_KP(k,p) = m_K(k) + sqrt(1/eta_K(k)/lambda_KP(k,p))*randn(1,1);
            
        end
    end
    
    
    %% the chains for the patients
    Z_KTP = nan(K,T,P);
    indexZ_1TP = nan(1,T,P);
    Y_PT = nan(P,T);
    for p=1:P
        g = patientsAssignedCohort_1K(p);
        [Z_KTp, indexZ_1Tp] = hmmmix_generateData_hiddenPatients_helper_MatlabC(indexM_1TG(:,:,g), alpha_K, rand(1,T));
        Z_KTP(:,:,p) = Z_KTp;
        indexZ_1TP(1,:,p) = indexZ_1Tp;
        
        Y_PT(p,:) = hmmmix_generateData_Y_pT_from_Z_KTp(Z_KTp, lambda_KP(:,p), mu_KP(:,p), nu_KP(:,p));
        
    end
    
    dataSt.Y_PT = Y_PT;
    % we could return also the values for the Z and M, but right now I'd
    % like to thing about other things
    
    % more consistent notation
    dataSt.m_K = m_K;
    dataSt.eta_K = eta_K;
    dataSt.gamma_K = gamma_K;
    dataSt.S_K = S_K;
    dataSt.nu_KP = nu_KP;
    dataSt.lambda_KP = lambda_KP;
    dataSt.mu_KP = mu_KP;
    dataSt.alpha_K = alpha_K;
    
    
    dataSt.hM_KTG = hM_KTG;
    dataSt.indexM_1TG = indexM_1TG + 1; % because those start at 0
    dataSt.indexZ_1TP = indexZ_1TP + 1;
    
    % one prior for all the chains
    dataSt.chainsPrior = pi_M_multinomial_parameters;
    dataSt.A_KKG = A_KKG;
    
    % to have this in the same format as the variational "C" component that
    % controls the partial cluster assignments
    dataSt.trueC = patientsAssignedCohort_bin';
    
end