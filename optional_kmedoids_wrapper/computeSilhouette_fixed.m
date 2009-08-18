function s = computeSilhouette_fixed(d, clusters, numG)

    if nargin < 3
        numG = max(clusters(:));
    end

    N = size(d,1);
    assert(N==size(d,2));
    
    b = zeros(1,numG);
    s = zeros(1,N);
    
    for i=1:N
        c = clusters(i);
        cI = find(clusters==c);

        % compute the mean distance of a to other members of the cluster
        a=0;
        if length(cI)>1
            a = mean(d(i,setdiff(cI,i)));
        end
        
        kI = setdiff(1:numG,c);
        b=inf(size(kI));
        for k=1:length(kI)   
            b(k)=mean(d(i,find(clusters==kI(k))));
        end
        minb=min(b);
        s(i) = (minb-a)/max(a,minb);
    end
    
end