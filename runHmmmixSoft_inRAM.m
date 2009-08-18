function reportSt = runHmmmixSoft_inRAM(dataSt, initValuesSt, algorithmParamsSt)
    % This function performs the hmmmix-soft algorithm without
    % consideration for memory requirement. The parameters consist of three
    % structures with certain fields explained below. Some default values
    % are provided by this function for unspecified fields. The output
    % 'reportSt' is also a structure with many fields.
    % 
    %
    %%% fields for dataSt : %%%
    %   Y_PT                 (mandatory, the observations, the data, all as one sequence)
    %   chromosomeIndices    (default = zeros(1,T))
    %   G                    (default = size(initValues.hC_GP,1),
    %                           allows for initValues.hC_GP to be left unspecified)  
    %   nu_KP                (mandatory, part of the emission model. Refer to thesis.)
    %    
    %%% fields for initValues. Refer to thesis for the roles of the fields : %%%
    %   initValues.hC_GP      (mandatory, unless we have dataSt.G)
    %   m_KP                  (mandatory, part of the emission model)
    %   eta_KP                (mandatory, part of the emission model)
    %   gamma_KP              (mandatory, part of the emission model)
    %   S_KP                  (mandatory, part of the emission model)
    %   mu_KP                 (default = initValuesSt.mu_KP, part of the emission model)
    %   lambda_KP             (default = 100*rand(K,P), part of the emission model)
    %   pseudoCounts   (default = T*eye(K) + ones(K,K), prior for the transition matrices)
    %   alpha_KK              (mandatory, part of transitions P(Z|M) in hmmmix.
    %                          Arguments with K elements are turned into diagonal matrices.
    %                          With hindsight, these should have been
    %                          normalized from the start.)
    %    
    %
    %%% fields for algorithmParamsSt : %%%
    %
    %   algorithmParamsSt.include_final_hard_assignment (default=false)
    %
    %   number_loops_to_stabilize_with_fixed_assignments
    %   verbose
    %   inner_maximal_number_of_iterations (default = 10)
    %   outer_maximal_number_of_iterations (default = 20)
    %   assignments_log_scaling (default = 1/sqrt(T))
    %   inner_tolerance (default = 0.1)
    %       This controls the convergence of the hidden chains distribution
    %       hM_KTG by determining when we should stop the cycle that
    %       updates the distributions and the parameters.
    %   include_final_hard_assignment (default = false)
    %       If this is set to 'true', after convergence we perform one Viterbi
    %       pass to get hard values for the hidden chains. We then get hard patient
    %       assignments in the same way.
    %
    %   thetaUpdates_forbidGoingBack   (default = true)
    %       The updates to the mu,lambda parameters are done in an
    %       iterative fashion to minimize a lower bound on the quantity
    %       that we really wish to minimize. With this option set, we are
    %       computing the loglikelihoods at every step to evaluate if we
    %       should reject updates that don't improve the loglik. This
    %       procedure becomes the bottleneck of the algorithm if we set
    %       this value to 'true'.
    %   thetaUpdates_whichLogLikelihoodToUse   (default = 'original')
    %       This argument is only used if thetaUpdates_forbidGoingBack is
    %       set to 'true'. Refer to 'hmmmix_mu_lambda_MLE_nopdf_loop.m' for
    %       more explanations. The two acceptable values are 'original' and
    %       'second'.
    %
    % The values of dataSt.chromosomeIndices are expected to be integers
    % indices starting from 1 that indicate the chromosomes to which the
    % values of dataSt.Y_PT belong. The observations dataSt.Y_PT themselves
    % contain all the concatenated sequences for the chromosomes. We assume
    % that, although the information for multiple chromosome can be mixed up,
    % we at least have that, within a given chromosome, the right order is determined
    % by the order in which they occur in dataSt.Y_PT.
    %
    %
    % Main output fields for hmmmix-soft :
    %   reportSt.hM_KTG
    %   reportSt.hC_GP
    %   reportSt.mu_KP
    %   reportSt.lambda_KP
    %   reportSt.transitionMatrices
    %   reportSt.initStates
    %   reportSt.loglikelihood_theta
    %
    % Output fields for monitoring (mostly for debugging) :
    %   reportSt.time
    %   reportSt.globalEnt
    %   reportSt.avLocalEnt
    %
    % Output fields for hard assignments,
    % (computed only if algorithmParamsSt.include_final_hard_assignment=true) :
    %    reportSt.loglikelihood_chains
    %    reportSt.loglikelihood_individual_chains
    %    reportSt.viterbiPathsForAllGroups
    %    reportSt.hC_GP_convertedToHard
    %    reportSt.hM_KTG_convertedToHard
    %    reportSt.M
    %    reportSt.patientsAssignmentIndices
    %    reportSt.loglikelihood_for_hard_projection
    %    reportSt.loglikelihood_for_hard_projection_individual_patients
        
    % To be completely honest, I think there could be a way to add some
    % kind of regularization to the initial priors for the hidden chains.
    % I'm uneasy about the fact that sometimes the initial priors can achive a
    % value of zero, but I don't want to hack something just like that.
    
    %% Read the input arguments. 
    
    Y_PT = dataSt.Y_PT;
    [P,T] = size(Y_PT);
    
    chromosomeIndices = getFieldOrDefault(dataSt, 'chromosomeIndices', ones(1,T));
    
    if ~hasField(dataSt, 'G')
        hC_GP = initValuesSt.hC_GP;
        G = size(hC_GP,1);
    else
        G = dataSt.G;
        hC_GP = getFieldOrDefault(initValuesSt, 'hC_GP', normalize(rand(G,P),1));
    end
    
    nu_KP = dataSt.nu_KP;
    
    m_KP = initValuesSt.m_KP;
    eta_KP = initValuesSt.eta_KP;
    gamma_KP = initValuesSt.gamma_KP;
    S_KP = initValuesSt.S_KP;
   
    mu_KP = getFieldOrDefault(initValuesSt, 'mu_KP', m_KP);
    K = size(mu_KP,1);
    
    lambda_KP = getFieldOrDefault(initValuesSt, 'lambda_KP', 100*rand(K,P));    
    pseudoCounts = getFieldOrDefault(initValuesSt, 'pseudoCounts', T*eye(K) + ones(K,K));
    
    alpha_KK = initValuesSt.alpha_KK;
    if numel(alpha_KK) ~= K*K
        % if the alpha_KK value isn't of size (K,K), it's probably because
        % we specified only the diagonal values and the rest should be ones
        alpha_KK = diag(alpha_KK(:)-1) + ones(K,K);
    end
    
    % Just a sanity check to make sure that there hasn't been a basic
    % misunderstanding.
    assert(all([K,P] == size(nu_KP)));
    assert(all([K,P] == size(mu_KP)));
    assert(all([K,P] == size(lambda_KP)));
    assert(all([K,P] == size(S_KP)));
    assert(all([K,P] == size(m_KP)));
    assert(all([K,P] == size(gamma_KP)));
    assert(all([K,P] == size(eta_KP)));
    assert(all([K,K] == size(alpha_KK)));
    assert(all([G,P] == size(hC_GP)));
    assert(all([P,T] == size(Y_PT)));
    
    %% Set up the transition matrices and initial states variables based on
    %  the number of chromosomes that the data represents.
    
    uniqChromInd = unique(chromosomeIndices);
    transitionMatrices = {};
    initStates = {};
    if K~=3
        fprintf('The only smart values for the initial states are when K=3 and we want to encourage the middle state.\n');
    end
    for g=1:G
        for c = uniqChromInd
            transitionMatrices{g,c} = normalize(pseudoCounts,2);
            % Specific to the case of K=3 with Sohrab's experiment.
            if K==3
                initStates{g,c} = toCol(normalize([1,100,1]));
            else
                initStates{g,c} = ones(K,1)/K;
            end
        end
    end
    
    % Read the parameters for the control of the algorithm. These are not
    % the parameters of the hmmmix model, but the knobs that control the
    % execution of the algorithm, more or less (like the number of
    % iterations to do or whether or not we should check for the
    % loglikelihood increases while optimizing with respect to the
    % parameters of the model (mu_KP, lambda_KP)). This includes the constant
    % 'tau' from my thesis.
    
    number_loops_to_stabilize_with_fixed_assignments = getFieldOrDefault(algorithmParamsSt, 'number_loops_to_stabilize_with_fixed_assignments', 5);
    verbose = getFieldOrDefault(algorithmParamsSt, 'verbose', false);
    inner_maximal_number_of_iterations = getFieldOrDefault(algorithmParamsSt, 'inner_maximal_number_of_iterations', 10);
    outer_maximal_number_of_iterations = getFieldOrDefault(algorithmParamsSt, 'outer_maximal_number_of_iterations', 20);
    assignments_log_scaling = getFieldOrDefault(algorithmParamsSt, 'assignments_log_scaling', 1/sqrt(T));
    inner_tolerance = getFieldOrDefault(algorithmParamsSt, 'inner_tolerance', 0.1);
    
    thetaUpdates_forbidGoingBack = getFieldOrDefault(algorithmParamsSt, 'thetaUpdates_forbidGoingBack', true);
    thetaUpdates_whichLogLikelihoodToUse = getFieldOrDefault(algorithmParamsSt, 'thetaUpdates_whichLogLikelihoodToUse', 'original');
    
    include_final_hard_assignment = getFieldOrDefault(algorithmParamsSt, 'include_final_hard_assignment', false);
    
    %% The main loop of the algorithm. It runs until we set the variable
    %  'done' to 'false' because the counter exceeds the maximal number of
    %  iterations (or until some other condition has been broken, but there
    %  isn't any other at the moment).
    
    done = false;
    niter = 1;
    % Sometimes we want to keep track of the time taken for the execution,
    % or of the entropy the patient assignments have.
    time_History = []; globalEnt_History = []; avLocalEnt_History = [];
    while ~done
        tic
        if verbose, fprintf('Starting an iteration of the main loop.\n'); end
       
        YgivenZ_KTP = hmmmix_YgivenZ_KTP(Y_PT, nu_KP, mu_KP, lambda_KP);
        YgivenM_KTP = hmmmix_YgivenM_KTP(YgivenZ_KTP, alpha_KK);

        R_KTG = hmmmix_R_KTG(YgivenM_KTP, hC_GP);
        hM_KTG = nan(K,T,G);
        for g=1:G
            [smoothedSequencesM_KT, updatedMLE_transitionMatrices, updatedMLE_initStates] = ...
                hmmmix_multiChromosome_singleGroup_smoothing(R_KTG(:,:,g), chromosomeIndices, transitionMatrices(g,:), initStates(g,:), pseudoCounts, inner_tolerance, inner_maximal_number_of_iterations, verbose);
            hM_KTG(:,:,g) = smoothedSequencesM_KT;
            for c = uniqChromInd
                transitionMatrices{g,c} = updatedMLE_transitionMatrices{c};
                initStates{g,c} = updatedMLE_initStates{c};
            end
        end
        
        hC_GP = hmmmix_hC_GP_withFreeEnergy(YgivenM_KTP, hM_KTG, assignments_log_scaling);
        [globalEnt,avLocalEnt] = twoEntropiesOfAssignments(hC_GP);
        if verbose, fprintf('The entropies : '), disp([globalEnt,avLocalEnt]); end
        globalEnt_History = [globalEnt_History, globalEnt];
        avLocalEnt_History = [avLocalEnt_History, avLocalEnt];

        % I should add here some options about evaluating the
        % log-likelihood because we might be more interested in tracking
        % the log-likelihood with hard imputed values than the lower bound
        % (or the lower bound on the lower bound, which is even less
        % direct).
        if number_loops_to_stabilize_with_fixed_assignments < niter
            [mu_KP, lambda_KP] = hmmmix_mu_lambda_MLE_nopdf_loop(hM_KTG, hC_GP, Y_PT, mu_KP, lambda_KP, alpha_KK, nu_KP, m_KP, gamma_KP, eta_KP, S_KP, thetaUpdates_forbidGoingBack, false, thetaUpdates_whichLogLikelihoodToUse);
        end
        
        time_for_this_pass = toc;
        time_History = [time_History,time_for_this_pass];
        if verbose, fprintf('\tIt took %f seconds for this pass in the main loop.\n', time_for_this_pass); end
        niter = niter + 1;
        if niter > outer_maximal_number_of_iterations
            done = true;
        end
    end

    reportSt.hM_KTG = hM_KTG;
    reportSt.hC_GP = hC_GP;
    reportSt.mu_KP = mu_KP;
    reportSt.lambda_KP = lambda_KP;
    reportSt.time = time_History;
    reportSt.globalEnt = globalEnt_History;
    reportSt.avLocalEnt = avLocalEnt_History;
    
    
    % Do one last step to get hard values for patients, just for fun and
    % then maybe impute hard values for their corresponding Z also. We are
    % throwing away here the good values that we have for hC_GP and the
    % hidden chains to replace them with hard values.
    
    if verbose, fprintf('Average time for iterations in main loop is %f seconds.\n', mean(time_History)); end
    
    if include_final_hard_assignment
        viterbiPathsForAllGroups = nan(G,T);
        for g=1:G
            viterbiPathIndices_concatenatedChromosomes = hmmmix_multiChromosome_singleGroup_viterbi(R_KTG(:,:,g), chromosomeIndices, transitionMatrices, initStates);
            viterbiPathsForAllGroups(g,:) = viterbiPathIndices_concatenatedChromosomes;
        end
        
        for g=1:G
            hM_KTG(:,:,g) = convertSequence_Integer_to_1K(viterbiPathsForAllGroups(g,:),K);
        end
        
        hC_GP = hmmmix_hC_GP_withFreeEnergy(YgivenM_KTP, hM_KTG, assignments_log_scaling);
        for p=1:P
            [junk, ind] = max(hC_GP(:,p));
            hC_GP(:,p) = 0;
            hC_GP(ind,p) = 1;
        end
        
        patientsAssignmentIndices = convertSequence_1K_to_Integer(hC_GP);
        
        % Compute the original log-likelihood.
        loglikelihood_for_hard_projection = 0;
        loglikelihood_for_hard_projection_individual_patients = zeros(1,P);
        for p=1:P
            E = log(YgivenM_KTP(:,:,p)) .* hM_KTG(:,:,patientsAssignmentIndices(p));
            loglikelihood_for_hard_projection = loglikelihood_for_hard_projection + sum(E(:));
            loglikelihood_for_hard_projection_individual_patients(p) = sum(E(:));
        end
        
        reportSt.viterbiPathsForAllGroups = viterbiPathsForAllGroups;
        reportSt.hC_GP_convertedToHard = hC_GP;
        reportSt.hM_KTG_convertedToHard = hM_KTG;
        
        % This is somewhat useless since I'm just converting back the data
        % for viterbiPathsForAllGroups.
        %   reportSt.M = zeros(G,T);
        %   for g=1:G
        %       reportSt.M(g,:) = convertSequence_1K_to_Integer(reportSt.hM_KTG_convertedToHard(:,:,g));
        %   end
        
        reportSt.patientsAssignmentIndices = patientsAssignmentIndices;
        reportSt.loglikelihood_for_hard_projection = loglikelihood_for_hard_projection;
        reportSt.loglikelihood_for_hard_projection_individual_patients = loglikelihood_for_hard_projection_individual_patients;


        % This thing isn't really fast, but it's done way after the
        % algorithm finished. You might want to turn this off manually if
        % you're not interested in the actual loglikelihood values
        % returned.        
        loglikelihood_individual_chains = zeros(1,G);
        for c = uniqChromInd
            ind = find(dataSt.chromosomeIndices == c);
            for g=1:G
               logT = log(transitionMatrices{g,c});
               %Pi = initStates{g,c};
               %loglikelihood_individual_chains = loglikelihood_individual_chains + log(Pi(ind(1)));
               for t=ind(2:end)
                   loglikelihood_individual_chains(g) = loglikelihood_individual_chains(g) + logT(reportSt.viterbiPathsForAllGroups(g,t-1),reportSt.viterbiPathsForAllGroups(g,t));
               end
            end
        end
        loglikelihood_chains = sum(loglikelihood_individual_chains);
        reportSt.loglikelihood_chains = loglikelihood_chains;
        reportSt.loglikelihood_individual_chains = loglikelihood_individual_chains;

    end

    reportSt.transitionMatrices = transitionMatrices;
    reportSt.initStates = initStates;

    reportSt.loglikelihood_theta = hmmmixsoft_loglik_theta(reportSt.mu_KP, reportSt.lambda_KP, initValuesSt.m_KP, initValuesSt.eta_KP, initValuesSt.gamma_KP, initValuesSt.S_KP);

    

    
end




