function LL = hmmmixsoft_loglik_theta(mu_KP, lambda_KP, m_KP, eta_KP, gamma_KP, S_KP)

    logNormal = -0.5*log(2*pi)+0.5*log(lambda_KP .* eta_KP)-0.5*(lambda_KP .* eta_KP) .* (mu_KP - m_KP).^2;
    logGamma = -gammaln(gamma_KP) + gamma_KP.*log(S_KP) + (gamma_KP-1).*log(lambda_KP) - S_KP .* lambda_KP;
    
    LL = sum(logNormal(:) + logGamma(:));

end