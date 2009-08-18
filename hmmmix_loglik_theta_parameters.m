function LL = hmmmix_loglik_theta_parameters(mu_KP, lambda_KP, m_K, eta_K, S_K, gamma_K)

    % This function may look a bit confusing, but the essential idea is
    % that we are creating these variables of size (K,P) from either
    % arguments of size K or size (K,P).
    % eg : m_K is supposed to have K elements so we replicate these
    %       values P times and then use the result in a vectorized formula.
    %       If m_K starts out with (K,P) elements already, we don't do the
    %       replication.

    [K,P] = size(mu_KP);
    assert(all([K,P] == size(lambda_KP)));

    % after all this, all the _K variables will be of size _KP with copied
    % values
    if numel(m_K) == K
        m_K = repmat(toCol(m_K), [1,P]);
    end
    assert(all([K,P] == size(m_K)));

    if numel(eta_K) == K
        eta_K = repmat(toCol(eta_K), [1,P]);
    end
    assert(all([K,P] == size(eta_K)));
    
    if numel(S_K) == K
        S_K = repmat(toCol(S_K), [1,P]);
    end
    assert(all([K,P] == size(S_K)));
    
    if numel(gamma_K) == K
        gamma_K = repmat(toCol(gamma_K), [1,P]);
    end
    assert(all([K,P] == size(gamma_K)));
    
    % up to a constant
    %E = 1/2*log(eta_K .* lambda_KP) - 1/2*(eta_K .* lambda_KP).*(mu_KP - m_K).^2 + (gamma_K/2 - 1).*log(lambda_KP) - 1/2./S_K .* lambda_KP;
    %LL = sum(E(:));
    
    % Exact value. It may seem wasteful to compute those gammaln(...)
    % calls, but it's happening only K*P times and that's nothing compared
    % to the K*P*T costs of other operations.
  
    Ewish = gamma_K .* log(S_K) - gammaln(gamma_K) + (gamma_K - 1) .* log(lambda_KP) - 0.5 * S_K .* lambda_KP;
    Enorm = -0.5 * log(pi) + 0.5 * log(eta_K .* lambda_KP) - 0.5 * eta_K .* lambda_KP .* (mu_KP - m_K).^2;
    
    LL = sum(Ewish(:) + Enorm(:));

    % For some reason, I had this : 
    %   Ewish = -gamma_K / 2 .* log(2*S_K) - gammaln(gamma_K/2) + ...
    %         (gamma_K - 2)/2 .* log(lambda_KP) - 0.5 ./ S_K .* lambda_KP;
    % before I decided to go ahead with the conversion and change the 1/S_K
    % for S_K. I'm not sure why there were all these divisions by 2 in
    % there. Maybe I copied some code over from the formulas with
    %   log Gamma(u | nu/2, nu/2).
    
end