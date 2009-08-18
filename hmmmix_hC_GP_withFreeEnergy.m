function hC_GP = hmmmix_hC_GP_withFreeEnergy(YgivenM_KTP, hM_KTG, S)

    [K,T,G] = size(hM_KTG);
    P = size(YgivenM_KTP,3);


    E = reshape(hM_KTG, [K*T,G])' * reshape(log(YgivenM_KTP), [K*T,P]);
    
    % Because I think I'm getting values that are too small and mess up the
    % normalize call.
    E = E - repmat(max(E), [G,1]);
    
    hC_GP = normalize(exp(S*E),1);

    % Just a small test to make sure that the reshape magic is correct, or
    % to make it clearer.
    if false
        E_test = zeros(G,P);
        for g=1:G
            for p=1:P
                E_test(g,p) = sum(sum(hM_KTG(:,:,g) .* log(YgivenM_KTP(:,:,p))));
            end
        end
        fprintf('Testing compact formula for hC_GP.');
        
        assert(   max(abs(E(:) - E_test(:))) < 1e-8    );
        fprintf('Passed.\n');
    end
end