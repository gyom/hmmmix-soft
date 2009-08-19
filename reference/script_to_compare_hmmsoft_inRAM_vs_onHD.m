%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This file is a test to make sure that the two versions of the hmmmix-soft
% algorithm (in memory vs on the hard drive) give the same results.
%
% At the same time, it shows how the code should be used.
%
% Three structures are constructed : dataSt, initValuesSt and algorithmParamsSt.
% Some data from the hmmmix-soft model is also generated at first, so that we can
% use it to see if the model learned is sensible.
%
% For a description about the fields of the structure obtained from the
% core runHmmmixSoft_inRAM function, refer to the source of that function.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% the true number of groups in the data generated
trueG=10;
% the number of groups that we'll use to model that data
G = 10;
% P is the number of patients, all throughout the code
P = 100;
% T is the number of time steps. For this example we have only one
% chromosome, so we set 'dataSt.chromosomeIndices' to be just ones.
T = 1000;
% we have three states for the hidden discrete values
K = 3;

% generate some data from the hmmmix model to use
dataSt = hmmmix_generateData_frugal(trueG, P, T, K);

% if the first 127 time steps represented the data first chromosome, and
% the next 51 time steps for the second chromosome, this would look like
% dataSt.chromosomeIndices = [repmat(1, [1,127]), repmat(2, [1,51])];
dataSt.chromosomeIndices = ones(1,T);

% refer to thesis for the meaning of the quantities
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
    initValuesSt.pseudoCounts = T*eye(K) + ones(K,K);
    initValuesSt.hC_GP = normalize(rand(G,P),1);

    % To use code from Sohrab's k-medoids, you need to add the
    % 'optional_kmedoids_wrapper' directory to your Matlab path.
    % You also need Sohrab's CNAhmmer code, but it shouldn't be found
    % before the functions from 'optional_kmedoids_wrapper' in the path.
    % Uncomment the two following lines to initialize the patient assignments with
    % k-medoids :
    %       kmedoidsSt = kmedoids_inRAM_wrapper_to_Sohrab_code(dataSt, initValuesSt);
    %       initValuesSt.hC_GP = kmedoidsSt.hC_GP;
    % It would be possible to use also the updated estimates for
    % 'kmedoidsSt.S_KP' and others, but it's not really relevant here.
    
    % Controls the number of times that we iterate with the initial patient
    % assignments set. Set to 1 if you don't initialize hC_GP to something
    % meaningful, and ~5 if you do.
    algorithmParamsSt.number_loops_to_stabilize_with_fixed_assignments = 3;
    algorithmParamsSt.verbose = false;
    algorithmParamsSt.assignments_log_scaling = 1/sqrt(T); % heuristic
    algorithmParamsSt.include_final_hard_assignment = true;

fprintf('Running hmmmix-soft in memory\n');
reportSt_inRAM = runHmmmixSoft_inRAM(dataSt, initValuesSt, algorithmParamsSt);

% use dataSt.Y_PT instead of passing the data by files
algorithmsParamsSt.scratchPath = '/tmp';
fprintf('Using the directory %s as scratch space to perform hmmmix-soft on the hard drive.\nYou might want to delete all the generated .mat files afterwards.\nIf you run hmmmix-soft on the hard drive using the same scratch directory all the time,\nthese files will be overwritten and nothing will accumulate.\n', algorithmsParamsSt.scratchPath);
reportSt_onHD = runHmmmixSoft_onHD(dataSt, initValuesSt, algorithmParamsSt);

% Now let's make sure that both versions returned the same thing. This is
% because the algorithm is deterministic. It makes a nice way to check the
% equivalence. Don't run the 'runHmmmixSoft_inRAM' with a large T, though.

assert(max(abs(reportSt_inRAM.hC_GP(:) - reportSt_onHD.hC_GP(:))) < 0.04);
assert(max(abs(reportSt_onHD.avLocalEnt(:) - reportSt_inRAM.avLocalEnt(:)))  < 0.01);

% I'm checking here the values for the imputed chains because it's one of
% the relevant quantities.
viterbiPaths = nan(G,T);
for g=1:G
    load(reportSt_onHD.viterbiPaths_GT_cachedFilenames{g}, 'viterbiPaths_gT');
    viterbiPaths(g,:) = viterbiPaths_gT;
end

assert(max(abs(reportSt_inRAM.viterbiPathsForAllGroups(:) - viterbiPaths(:)))  < 10e-5);


