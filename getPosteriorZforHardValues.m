function [posteriorZ_KTP, MLE_Z_KTP, indices_MLE_Z_TP] = getPosteriorZforHardValues(hM_KTG_hard, hC_GP_hard, alpha_KK, Y_PT, nu_KP, mu_KP, lambda_KP)
    
    [K,T,G] = size(hM_KTG_hard);
    [G,P] = size(hC_GP_hard);
    
    
    %YgivenM_KTP = hmmmix_YgivenM_KTP(YgivenZ_KTP, alpha_KK);
    
    YgivenZ_KTP = hmmmix_YgivenZ_KTP(Y_PT, nu_KP, mu_KP, lambda_KP);
    rho_KTP = hmmmix_compute_rho_KTP_nopdf_MatlabC(hM_KTG_hard, hC_GP_hard, diag(alpha_KK));
    
    posteriorZ_KTP = normalize(YgivenZ_KTP .* rho_KTP, 1);
    
    MLE_Z_KTP = zeros(K,T,P);
    indices_MLE_Z_TP = zeros(T,P);
    
    % This is slow because we're in Matlab, but right now I don't care
    % about this and it's just a post-processing step run just once.
    for t=1:T
        for p=1:P
            [junk, ind] = max(posteriorZ_KTP(:,t,p));
            MLE_Z_KTP(ind,t,p) = 1;
            indices_MLE_Z_TP(t,p) = ind;
        end
    end
    
end




