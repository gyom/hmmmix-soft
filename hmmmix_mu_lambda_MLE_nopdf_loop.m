function [mu_KP, lambda_KP] = hmmmix_mu_lambda_MLE_nopdf_loop(hM_KTG, hC_GP, Y_PT, mu_KP, lambda_KP, alpha_K, nu_KP, m_K, gamma_K, eta_K, S_K, forbidGoingBack, verbose, whichLogLikelihoodToUse)

    [K,T,G] = size(hM_KTG);
    P = size(hC_GP,2);

    % As mentioned in my thesis, we are using a lower bound on the original
    % lower bound. Nothing guarantees that the updates here will be
    % beneficial in terms of maximizing the first (tightest) lower bound,
    % so I added the option here of rejecting updates when the values for
    % the loglikelihood of the observations is going down.
    if nargin < 12
        % The default value is 'true', but it's a costly thing to do.
        forbidGoingBack = true;
    end
    
    if nargin < 13
        verbose = false;
    end
    
    if nargin < 14
        whichLogLikelihoodToUse = 'original'; % 'original' or 'second'
        % Determines if we'll use 
        %    sum(sum(log(sum(rho_KTP .* YgivenZ_KTP,1))))
        % or
        %    sum(sum(sum(rho_KTP .* log(YgivenZ_KTP))))
        % as loglikelihood value that we want to monitor.
        % With 'original' we use the first one.
    end
    
    debug = false;


    % I had this alpha_KK in the main code for the FL data, but now I'm
    % hitting this function 'hmmmix_compute_rho_KTP_nopdf_MatlabC' that I
    % won't be rewriting to allow alpha_KK to be square. I also realized
    % that I was picking alpha_KK = [20,1,1; 1,10,1; 1,1,20]; anyways, and
    % I can reproduce that value with alpha_K = [20,10,20]. So there it is.
    % This function 'hmmmix_mu_lambda_MLE_nopdf_loop' now takes alpha_KK
    % also because I'm downgrading it.
    if all([K,K] == size(alpha_K))
        alpha_K = diag(alpha_K);
        % Make sure someone hasn't had the bad (incompatible) idea of
        % normalizing alpha already. Those alpha matrices are always
        % renormalized on the spot by some function.
        assert(min(alpha_K(:)) >= 1)
    end
    assert(length(alpha_K) == K);
    
    rho_KTP = hmmmix_compute_rho_KTP_nopdf_MatlabC(hM_KTG, hC_GP, alpha_K);    
    
    % Backup the original quantities so that, if 'forbidGoingBack' is true
    % and the loglikelihood went down, we can revert to them.
    if forbidGoingBack || debug
        original_mu_KP = mu_KP;
        original_lambda_KP = lambda_KP;
        YgivenZ_KTP = hmmmix_YgivenZ_KTP(Y_PT, nu_KP, mu_KP, lambda_KP);
        switch whichLogLikelihoodToUse
            case 'original'
                original_loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K) + sum(sum(log(sum(rho_KTP .* YgivenZ_KTP,1))));
            case 'second'
                original_loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K) + sum(rho_KTP(:) .* log(YgivenZ_KTP(:)));
        end
    end
    
    previous_u_KTP = zeros(K,T,P);
    paramsDone = false;
    safetyCounter = 1;

    % not used
    if debug
        loglikHistory = [original_loglik];
    end
        
    while ~paramsDone
        
        u_KTP = hmmmix_compute_u_KTP(rho_KTP, Y_PT, mu_KP, lambda_KP, nu_KP);
        [mu_KP, lambda_KP] = hmmmix_mu_lambda_MLE(rho_KTP, Y_PT, u_KTP, m_K, nu_KP, gamma_K, eta_K, S_K);

        if (max(abs(previous_u_KTP(:) - u_KTP(:))) < 0.01) || safetyCounter > 50
            paramsDone = true;
            %fprintf('optimizing lambda and mu took %d iterations\n', safetyCounter);
        end
        safetyCounter = safetyCounter + 1;
        previous_u_KTP = u_KTP;

        % Computing the loglikelihood is expensive, so don't do it unless you
        % want to track the values.
        if debug
            YgivenZ_KTP = hmmmix_YgivenZ_KTP(Y_PT, nu_KP, mu_KP, lambda_KP);
            switch whichLogLikelihoodToUse
                case 'original'
                    loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K) + sum(sum(log(sum(rho_KTP .* YgivenZ_KTP,1))));
                case 'second'
                    loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K) + sum(rho_KTP(:) .* log(YgivenZ_KTP(:)));
            end
            loglikHistory = [loglikHistory, loglik];
        end        
        
        % This would be ok if we weren't looking for a fixed point, I
        % think, but now I'm not sure it's reasonable to expect that every
        % time we set u or mu,lambda we should expect an increase in
        % loglik.
        %if previous_loglik > loglik
        %    mu_KP = previous_mu_KP;
        %    lambda_KP = previous_lambda_KP;
        %    loglik = previous_loglik;
        %    paramsDone = true;
        %end
        
    end
    
    if debug
        loglikHistory
    end   
    
    % Computing the loglikelihood is expensive, so don't do it unless you
    % have to.
    if debug || forbidGoingBack
        YgivenZ_KTP = hmmmix_YgivenZ_KTP(Y_PT, nu_KP, mu_KP, lambda_KP);
        switch whichLogLikelihoodToUse
            case 'original'
                loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K) + sum(sum(log(sum(rho_KTP .* YgivenZ_KTP,1))));
            case 'second'
                loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K) + sum(rho_KTP(:) .* log(YgivenZ_KTP(:)));
        end
        %loglik = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K,
        %eta_K, S_K, gamma_K) + sum(rho_KTP(:) .* YgivenZ_KTP(:));
    end
    
    % if this check fails, we want to revert to the original parameters
    if forbidGoingBack
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