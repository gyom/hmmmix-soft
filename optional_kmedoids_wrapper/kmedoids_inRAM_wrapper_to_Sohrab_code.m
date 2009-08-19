function kmedoidsSt = kmedoids_inRAM_wrapper_to_Sohrab_code(dataSt, initValuesSt)
    % This function is some sort of hack to use Sohrab's code to perform
    % the k-medoids initialization. One of the problems that I had was that
    % I wanted to override some of the hyperparameters and that meant that
    % I had to write different versions of his functions.
    % 
    % Ideally, people using my code for hmmmix-soft would have their own
    % way of initializing the patients (setting initValues.hC_GP), but
    % since I wrote this to use kmedoids, it can be used.
    %
    % This function isn't ready to handle large datasets, unlike
    % 'runHmmsoft_onHD', but I might be possible to rewrite it for that
    % purpose. Since I'm using k-medoids more or less like a black box (to
    % some degree), I'm not touching this. I'd recommend a completely
    % random initialization when using large datasets that require the use
    % of 'runHmmsoft_onHD' instead of 'runHmmsoft_inRAM'.
    %
    % This function returns a structure whose fields contain variables of
    % potential interest that were produced during the k-medoids
    % intialization. The most important field should be kmedoidsSt.G_init
    % where kmedoidsSt.G_init(p) contains the index of the group to which
    % we assign patient p. Hard values for the hidden chains are found in
    % kmedoidsSt.M_init, but they are not used by the hmmsoft algorithm.
    %
    % To use these functions, the directory 'optional_kmedoids_wrapper' has
    % to be included in the path so that the functions in it will be found
    % before those coded by Sohrab Shah.

    [K,P] = size(initValuesSt.nu_KP);
    %T = size(dataSt.Y_PT,2);
    
    if ~hasField(dataSt, 'G')
        numG = size(initValuesSt.hC_GP,1);
    else
        numG = dataSt.G;
    end

    %nu=[3,3,3];
    nu=toRow(mean(initValuesSt.nu_KP,2));
    kappa=[0.1 0.8 0.1]*1000;
    maxiter = 100;

        % impute the Z. Taken from the 'convertToWecca' file.
        % The problem is that the values for lambda, for example, is buried
        % somewhere deep in there and is set to something around 1000. Then the
        % mus are set to somewhere around [-15, 0, 15] and this is completely
        % independant of our actual data. One wonders about the quality of the
        % imputed values for Z after that.
        [mus lambdas eta m gamma S Z pZ] = ...
            initialiseZ(dataSt.Y_PT,dataSt.chromosomeIndices,2,K,kappa,nu,maxiter, dataSt, initValuesSt);

        % now we're into runHmmMixFromInit

        % k-medoids to get a good initialization
        R=30;
        dMetric='Hamming';
        maxiter = 30; % instead of R=100, maxiter=1000 which takes an hour
        [numP,N] = size(Z);
        [E, pM] = calculateEntropy(Z,K);
        En = (E-mean(E))/std(E);
        W = ones(size(En))./(1+exp(-En/0.25));
        [G_init medoids] = multRestartKmedoidsWeighted(Z,W,R,numG,dMetric,maxiter);
    
        % because it breaks something if not all the groups are used, so I
        % remap them all to use integers as little as possible without
        % changing the clustering
            unique_values_of_G = unique(G_init);
            new_G_init = zeros(P,1);
            for p=1:P
                new_G_init(p) = find(unique_values_of_G == G_init(p));
            end
            G_init = new_G_init;
                
        f=2; % as Sohrab asked

        [mus,lambdas,eta,m,gamma,S] = initialiseTHyperparamsMulti(dataSt.Y_PT,f,K, dataSt, initValuesSt);
        pZ = calls2Probs(Z,K);
        Ethresh = 0.75;
        M_init = initialiseMWithEntropy(Z,G_init,K,Ethresh);
        [numP,N] = size(Z);
        
    kmedoidsSt.mu_KP = mus';
    kmedoidsSt.lambda_KP = lambdas';
    kmedoidsSt.eta_KP = eta';
    kmedoidsSt.m_KP = m';
    kmedoidsSt.gamma_KP = gamma';
    kmedoidsSt.S_KP = S';
    
    kmedoidsSt.M_init = M_init;
    % M_init(g,t) = the value of {1,...,K} imputed initially for patient p
    %               at time step t
    kmedoidsSt.G_init = G_init';
    % G_init(p) = the group g in {1,...,G} to which we initially assign
    %               patient p
    
    kmedoidsSt.hC_GP = zeros(numG,P);
    for p=1:P
        kmedoidsSt.hC_GP(G_init(p), p) = 1;
    end
    
end