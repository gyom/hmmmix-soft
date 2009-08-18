function YgivenZ_KTP = hmmmix_YgivenZ_KTP(Y_PT, nu_KP, mu_KP, lambda_KP)

    [P,T] = size(Y_PT);
    K = size(nu_KP,1);
    
    assert(all([K,P]==size(nu_KP)));
    assert(all([K,P]==size(mu_KP)));
    assert(all([K,P]==size(lambda_KP)));
    
    % For us to be able to use vectorized calls, we need to set up a few
    % matrices to be fed to the final call to tpdf.
    V = permute(repmat(nu_KP, [1,1,T]), [1,3,2]);
    B_Y = permute(repmat(Y_PT, [1,1,K]), [3,2,1]);
    B_mu = permute(repmat(mu_KP, [1,1,T]), [1,3,2]);
    B_sqrtlambda = permute(repmat(sqrt(lambda_KP), [1,1,T]), [1,3,2]);
    
    X = (B_Y - B_mu) .* B_sqrtlambda;
    YgivenZ_KTP = tpdf(X,V) .* B_sqrtlambda;

    if P == 1
        assert(all([K,T] == size(YgivenZ_KTP)));
    else
        assert(all([K,T,P] == size(YgivenZ_KTP)));
    end
    
end