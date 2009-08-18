function Y_pT = hmmmix_generateData_Y_pT_from_Z_KTp(Z_KTp, lambda_Kp, mu_Kp, nu_K)

    % We're using Z_pT here and not the version with indices instead of
    % 1-of-K notation. The last operation wouldn't make sense otherwise.

    K = size(lambda_Kp,1);
    assert(K == size(mu_Kp,1));
    assert(K == size(nu_K,1));

    Z_KTp = squeeze(Z_KTp);
    T = size(Z_KTp,2);
    assert(K == size(Z_KTp,1));
    
    V = repmat(toCol(nu_K), [1,T]);

    KY = trnd(V) ./ repmat(sqrt(lambda_Kp), [1,T]) + repmat(mu_Kp, [1,T]);
    
    % this picks one value from each column
    Y_pT = sum(KY .* Z_KTp);

end