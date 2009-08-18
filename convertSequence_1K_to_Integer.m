function A = convertSequence_1K_to_Integer(E)

    % It's slow, but I don't really care about this for now.

    % E = [1,0,1,1,0;
    %      0,1,0,0,0;
    %      0,0,0,0,1]
    % gives
    %    A = [1,2,1,1,3]
    
    [K,T] = size(E);
    
    A = zeros(1,T);

    % make sure there is only one value of 1 per column
    assert(sum(E) == 1);
    
    for t=1:T
        A(t) = find(E(:,t));
    end
    
end