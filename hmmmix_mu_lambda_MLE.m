function [mu_KP, lambda_KP] = hmmmix_mu_lambda_MLE(rho_KTP, Y_PT, u_KTP, m_K, nu_K, gamma_K, eta_K, S_K)

    [K,T,P] = size(rho_KTP);

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
    
    if numel(nu_K) == K
        nu_K = repmat(toCol(nu_K), [1,P]);
    end
    assert(all([K,P] == size(nu_K)));
    
    m_KP = m_K;
    eta_KP = eta_K;
    S_KP = S_K;
    gamma_KP = gamma_K;
    nu_KP = nu_K;
    
    mu_KP = zeros(K,P);
    lambda_KP = zeros(K,P);

    % I came to a realization while writing code like this. It's not worth
    % complicating the code to vectorize over all dimensions at the same
    % time. If T is large, then it's totally acceptable to iterate slowly
    % over K,P because it won't affect the execution by any large factor.
    for k=1:K
        for p=1:P
            % for mu_KP
            numer = sum(rho_KTP(k,:,p) .* u_KTP(k,:,p) .* Y_PT(p,:)) + eta_KP(k,p)*m_KP(k,p);
            denom = sum(rho_KTP(k,:,p) .* u_KTP(k,:,p)) + eta_KP(k,p);
            mu_KP(k,p) = numer / denom;
       
            % for lambda_KP
            numer1 = sum( rho_KTP(k,:,p) .* u_KTP(k,:,p) .* ( Y_PT(p,:) - mu_KP(k,p) ).^2);
            denom = sum( rho_KTP(k,:,p) ) + gamma_KP(k,p) - 1;
            
            numer2 = nu_KP(k,p)*(mu_KP(k,p) - m_KP(k,p))^2 + S_KP(k,p);

            % we want (numer1/denom + numer2/denom)^-1
            
            % from paper
            lambda_KP(k,p) = 1 / (numer1/denom + numer2/denom);
            
            % test
            % lambda_KP(k,p) = (numer1/denom + numer2/denom);
            
            % test 2. Looks like it's working relatively well, but fails
            % sometimes by a small margin.
            %lambda_KP(k,p) = sqrt(denom/(numer1 + numer2));
            
        end
    end
            
                
end