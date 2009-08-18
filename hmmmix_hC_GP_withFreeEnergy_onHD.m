function hC_GP = hmmmix_hC_GP_withFreeEnergy_onHD(logYgivenM_KTP_cachedFilenames, hM_KTG_cachedFilenames, S)

    P = length(logYgivenM_KTP_cachedFilenames);
    G = length(hM_KTG_cachedFilenames);
    
    E = zeros(G,P);
    
    for g=1:G
        load(hM_KTG_cachedFilenames{g}, 'hM_KTg'); % loads hM_KTg
        for p=1:P
            load(logYgivenM_KTP_cachedFilenames{p}, 'logYgivenM_KTp'); %loads logYgivenM_KTp
            %size(hM_KTg)
            %size(logYgivenM_KTp)
            E(g,p) = sum(hM_KTg(:) .* logYgivenM_KTp(:));
            clear logYgivenM_KTp;
        end
        clear hM_KTg;
    end
    
    % Because I think I'm getting values that are too small and mess up the
    % normalize call.
    E = E - repmat(max(E), [G,1]);
    
    hC_GP = normalize(exp(S*E),1);


end