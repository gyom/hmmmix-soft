
trueG=10;
G = 10;
P = 10;
T = 100;
K = 3;

dataSt = hmmmix_generateData_frugal(trueG, P, T, K);

dataSt.chromosomeIndices = ones(1,T);

dataSt.m_KP = repmat(dataSt.m_K, [1,P]);
dataSt.eta_KP = repmat(dataSt.eta_K, [1,P]);
dataSt.gamma_KP = repmat(dataSt.gamma_K, [1,P]);
dataSt.S_KP = repmat(dataSt.S_K, [1,P]);
dataSt.alpha_KK = diag(dataSt.alpha_K-1)+ones(K,K);

% prepare hmmmix-soft
    initValuesSt.nu_KP = dataSt.nu_KP;
    initValuesSt.m_KP = dataSt.m_KP;
    initValuesSt.eta_KP = dataSt.eta_KP;
    initValuesSt.gamma_KP = dataSt.gamma_KP;
    initValuesSt.S_KP = dataSt.S_KP;
    initValuesSt.alpha_KK = dataSt.alpha_KK;
    % A bit cheating would be to select
    %   initValuesSt.lambda_KP = mean(dataSt.lambda_KP)*(1 - mean(randn(10,1)))
    % but both methods would benefit from this ...
    initValuesSt.lambda_KP = 10*rand(K,P);
    %initValuesSt.lambda_KP = mean(dataSt.lambda_KP)*(1 - mean(randn(10,1)))
    initValuesSt.pseudoCounts = T*eye(K) + ones(K,K);
    
    initValuesSt.hC_GP = normalize(rand(G,P),1);

    algorithmParamsSt.number_loops_to_stabilize_with_fixed_assignments = 1;
    algorithmParamsSt.verbose = false;
    algorithmParamsSt.assignments_log_scaling = 1/sqrt(T);
    algorithmParamsSt.include_final_hard_assignment = true;

reportSt_inRAM = runHmmmixSoft_inRAM(dataSt, initValuesSt, algorithmParamsSt);

% use dataSt.Y_PT instead of passing the data by files
algorithmsParamsSt.scratchPath = '/tmp';
reportSt_onHD = runHmmmixSoft_onHD(dataSt, initValuesSt, algorithmParamsSt);

% Now let's make sure that both versions returned the same thing. This is
% because the algorithm is deterministic. It makes a nice way to check the
% equivalence. Don't run the 'runHmmmixSoft_inRAM' with a large T, though.

assert(max(abs(reportSt_inRAM.hC_GP(:) - reportSt_onHD.hC_GP(:))) < 10e-5);
assert(max(abs(reportSt_onHD.avLocalEnt(:) - reportSt_inRAM.avLocalEnt(:)))  < 10e-5);

viterbiPaths = nan(G,T);
for g=1:G
    load(reportSt_onHD.viterbiPaths_GT_cachedFilenames{g}, 'viterbiPaths_gT');
    viterbiPaths(g,:) = viterbiPaths_gT;
end

assert(max(abs(reportSt_inRAM.viterbiPathsForAllGroups(:) - viterbiPaths(:)))  < 10e-5);


