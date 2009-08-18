function E = convertSequence_Integer_to_1K(A,K)

    % A = [1,2,1,1,3]
    % K = 3;
    %
    % gives E = [1,0,1,1,0;
    %            0,1,0,0,0;
    %            0,0,0,0,1]
    
    if nargin < 2
        K = max(A);
    end

    T = numel(A);
    E = zeros(K,T);

    offsets = (0:(T-1))*K;
    
    E(A + offsets) = 1;
    
end