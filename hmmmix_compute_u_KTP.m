function u_KTP = hmmmix_compute_u_KTP(rho_KTP, Y_PT, mu_KP, lambda_KP, nu_KP)

    [K,T,P] = size(rho_KTP);

    u_KTP = zeros(K,T,P);

%     for t=1:T
%         for p=1:P
%             
%             % my pooled version
%             %tmp = (2*mean(nu_KP(:,p)) - 1)/( sum(((Y_PT(p,t) - mu_KP(:,p)).^2).*lambda_KP(:,p).*squeeze(rho_KTP(:,t,p))) + mean(nu_KP(:,p)));
%             % version 2, non-regularized
%             %tmp = (mean(nu_KP(:,p)) - 1)/( sum(((Y_PT(p,t) - mu_KP(:,p)).^2).*lambda_KP(:,p).*squeeze(rho_KTP(:,t,p))));
%             for k=1:K            
%                 %paper and Archambeau
%                 %u_KTP(k,t,p) = (1+nu_KP(k,p))/( ((Y_PT(p,t) - mu_KP(k,p))^2)*lambda_KP(k,p) + nu_KP(k,p));
% 
%                 % mine
%                 u_KTP(k,t,p) = (-1+nu_KP(k,p))/( ((Y_PT(p,t) - mu_KP(k,p))^2)*lambda_KP(k,p) + nu_KP(k,p));
%                 
%                 % my pooled version
%                 %u_KTP(k,t,p) = tmp;
%                 
%                 %test
%                 %u_KTP(k,t,p) = (1+nu_K(k))/( (lambda_KP(k,p)*(Y_PT(p,t) - mu_KP(k,p)))^2 + nu_K(k));
%             end
%         end
%     end

    for p=1:P
        for k=1:K
            u_KTP(k,:,p) = (-1+nu_KP(k,p)) ./ ( ((Y_PT(p,:) - mu_KP(k,p)).^2)*lambda_KP(k,p) + nu_KP(k,p));
        end
    end

end