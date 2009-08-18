function [mu_KP, lambda_KP] = hmmmix_mu_lambda_MLE_nopdf_loop_onHD(hM_KTG_cachedFilenames, hC_GP, Y_PT_cachedFilenames, mu_KP, lambda_KP, alpha_KK, nu_KP, m_KP, gamma_KP, eta_KP, S_KP, rho_KTP_cachedFilenames, YgivenZ_KTP_cachedFilenames, forbidGoingBack, verbose, whichLogLikelihoodToUse)

    if nargin < 14
        % The default value is 'true', but it's a costly thing to do.
        forbidGoingBack = true;
    end

    if nargin < 16
        whichLogLikelihoodToUse = 'original'; % 'original' or 'second'
        % Determines if we'll use 
        %    sum(sum(log(sum(rho_KTP .* YgivenZ_KTP,1))))
        % or
        %    sum(sum(sum(rho_KTP .* log(YgivenZ_KTP))))
        % as loglikelihood value that we want to monitor.
        % With 'original' we use the first one.
    end    


    [G,P] = size(hC_GP);
    K = size(mu_KP,1);
    assert(all([K,P] == size(mu_KP)));
    assert(all([K,P] == size(lambda_KP)));
    assert(all([K,P] == size(m_KP)));
    assert(all([K,P] == size(S_KP)));
    assert(all([K,P] == size(gamma_KP)));
    assert(all([K,P] == size(eta_KP)));
    
    debug = false;

    % computes the rho_KTP and update the files
    hmmmix_compute_rho_KTP_nopdf_onHD(hM_KTG_cachedFilenames, hC_GP, alpha_KK, rho_KTP_cachedFilenames);
    
    % backup
    original_mu_KP = mu_KP;
    original_lambda_KP = lambda_KP;
    
    % this computes  original_loglik
    % and updates the values for YgivenZ_KTp
    original_loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_KP, eta_KP, S_KP, gamma_KP);
    for p=1:P
        load(Y_PT_cachedFilenames{p}); % loads Y_pT
        YgivenZ_KTp = hmmmix_YgivenZ_KTP(Y_pT, nu_KP(:,p), mu_KP(:,p), lambda_KP(:,p));
        save(YgivenZ_KTP_cachedFilenames{p}, 'YgivenZ_KTp');
        load(rho_KTP_cachedFilenames{p}); % loads rho_KTp
        switch whichLogLikelihoodToUse
            case 'original'
                original_loglik = original_loglik + sum(sum(log(sum(rho_KTp .* YgivenZ_KTp,1))));
            case 'second'
                original_loglik = original_loglik + sum(rho_KTp(:) .* log(YgivenZ_KTp(:)));
        end
        %original_loglik = original_loglik + sum(rho_KTp(:) .* log(YgivenZ_KTp(:)));
        % I need that value of T later on so I might as well get it here.
        T = size(YgivenZ_KTp,2);
        clear YgivenZ_KTp;
        clear Y_pT;
    end
    
    % This updates the slices  mu_KP(:,p), lambda_KP(:,p) for p=1:P.
    % It can be done independantly for all the patients without problems.
    for p=1:P
        previous_u_KTp = zeros(K,T,1);
        paramsDone = false;
        safetyCounter = 1;
        load(rho_KTP_cachedFilenames{p}); % loads rho_KTp
        load(Y_PT_cachedFilenames{p}); % loads Y_pT
        
        while ~paramsDone
            u_KTp = hmmmix_compute_u_KTP(rho_KTp, Y_pT, mu_KP(:,p), lambda_KP(:,p), nu_KP(:,p));
            
            [mu_Kp, lambda_Kp] = hmmmix_mu_lambda_MLE(rho_KTp, Y_pT, u_KTp, m_KP(:,p), nu_KP(:,p), gamma_KP(:,p), eta_KP(:,p), S_KP(:,p));
            mu_KP(:,p) = mu_Kp;
            lambda_KP(:,p) = lambda_Kp;
            
            if (max(abs(previous_u_KTp(:) - u_KTp(:))) < 0.01) || safetyCounter > 50
                paramsDone = true;
                %fprintf('optimizing lambda and mu took %d iterations\n', safetyCounter);
            end
            safetyCounter = safetyCounter + 1;
            previous_u_KTp = u_KTp;
        end
    end
    
    if forbidGoingBack
        
        % We want to compute the current loglik to see if it went up or not,
        % but it's so costly that we should only do it if we're interested in
        % using these values.
        loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_KP, eta_KP, S_KP, gamma_KP);
        for p=1:P
            load(Y_PT_cachedFilenames{p}); % loads Y_pT
            YgivenZ_KTp = hmmmix_YgivenZ_KTP(Y_pT, nu_KP(:,p), mu_KP(:,p), lambda_KP(:,p));
            save(YgivenZ_KTP_cachedFilenames{p}, 'YgivenZ_KTp');
            load(rho_KTP_cachedFilenames{p}); % loads rho_KTp
            switch whichLogLikelihoodToUse
                case 'original'
                    loglik = loglik + sum(sum(log(sum(rho_KTp .* YgivenZ_KTp,1))));
                case 'second'
                    loglik = loglik + sum(rho_KTp(:) .* log(YgivenZ_KTp(:)));
            end
            %loglik = loglik + sum(rho_KTp(:) .* log(YgivenZ_KTp(:)));
            clear YgivenZ_KTp;
            clear Y_pT;
        end
        
        % if this check fails, we want to revert to the original parameters
        if loglik < original_loglik

            mu_KP = original_mu_KP;
            lambda_KP = original_lambda_KP;
            if verbose
                fprintf('F : The loglikelihood of the parameters went from %f to %f so we forget about this round.\n', original_loglik, loglik);
            end
        else
            if verbose
                fprintf('S : The loglikelihood of the parameters went from %f to %f so we keep the suggestions.\n', original_loglik, loglik);
            end
        end
    end
    
end