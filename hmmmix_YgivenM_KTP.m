function YgivenM_KTP = hmmmix_YgivenM_KTP(YgivenZ_KTP, alpha)

    % 'alpha' can have K elements or K^2. If it has K components, we assume
    % that the full corresponding matrix is of the form
    %
    % ones(K,K) + diag(alpha-1)
    
    [K,T,P] = size(YgivenZ_KTP);
    
    if numel(alpha) == K
        alpha = ones(K,K) + diag(alpha-1);
    end

    assert(all([K,K] == size(alpha)));
    
    % Now we assume that alpha represents some sort of transition matrix
    % where the odds of going from M=j to Z=k are given by
    %   alpha(j,k) / sum(alpha(j,:),2)
    
    alpha = normalize(alpha,2);
    
    YgivenM_KTP = reshape(alpha' * reshape(YgivenZ_KTP, [K,P*T]), [K,T,P]);

end