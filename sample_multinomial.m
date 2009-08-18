function [D_1K, D_bin] = sample_multinomial(theta,N)

    % not vectorial because I don't need it now
    p = cumsum(toRow(normalize(theta)));
    
    D_1K = zeros(N,1);
    D_bin = zeros(N, length(theta));
    
    for n=1:N
        [junk, ind] = find( rand(1,1) <= p, 1, 'first');
        D_1K(n,1) = ind;
        D_bin(n,ind) = 1;
    end

end