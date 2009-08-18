function hmmmix_compute_rho_KTP_nopdf_onHD(hM_KTG_cachedFilenames, hC_GP, alpha_KK, rho_KTP_cachedFilenames)

    K = size(alpha_KK,1);
    assert(K == size(alpha_KK,2));

    S = normalize(alpha_KK,2);

    load(hM_KTG_cachedFilenames{1});
    T = size(hM_KTg,2);
    clear hM_KTg;
    
    [G,P] = size(hC_GP);
    
    for p=1:P
        rho_KTp = zeros(K,T,1);
        for g=1:G
            load(hM_KTG_cachedFilenames{g}); % loads hM_KTg
            accum = zeros(K,T,1);
            for k=1:K
                for j=1:K
                    accum(k,:,1) = accum(k,:,1) + S(j,k)*hM_KTg(j,:,1);
                end
            end
            rho_KTp = rho_KTp + hC_GP(g,p) * accum;
        end
        save(rho_KTP_cachedFilenames{p}, 'rho_KTp');
    end

end