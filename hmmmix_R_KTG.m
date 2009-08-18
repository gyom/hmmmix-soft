function R_KTG = hmmmix_R_KTG(YgivenM_KTP, hC_GP)

    [K,T,P] = size(YgivenM_KTP);
    G = size(hC_GP,1);
    assert(P == size(hC_GP,2));

    R_KTG = reshape( log(reshape(YgivenM_KTP, [K*T,P])) * hC_GP', [K,T,G]);

end