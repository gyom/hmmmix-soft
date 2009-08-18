function hmmmix_R_KTG_onHD(logYgivenM_KTP_cachedFilenames, hC_GP, R_KTG_cachedFilenames, K, T)

    %[K,T,P] = size(YgivenM_KTP);
    [G,P] = size(hC_GP);
    assert(P == length(logYgivenM_KTP_cachedFilenames));

    if nargin < 5
        % some any file just to get the value of K,T
        load(logYgivenM_KTP_cachedFilenames{1}); % loads 'logYgivenM_KTp'
        K = size(logYgivenM_KTp,1);
        T = size(logYgivenM_KTp,2);
    end
        
    % This is MADNESS !! Having G*P load/writes just for that ... that's
    % what you get for not thinking about memory. I'm not even trying to
    % predict how much of that is necessary and how much is just wasteful.
    for g=1:G
        % delete the previous file
        %delete(R_KTG_cachedFilenames{g});
        R_KTg = zeros(K,T,1);
        for p=1:P
            load(logYgivenM_KTP_cachedFilenames{p}); % loads 'logYgivenM_KTp'
            R_KTg = R_KTg + logYgivenM_KTp * hC_GP(g,p);
        end
        save(R_KTG_cachedFilenames{g}, 'R_KTg');
    end

end