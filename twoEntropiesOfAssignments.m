function [globalEnt,avLocalEnt] = twoEntropiesOfAssignments(hC_GP)
    [G,P] = size(hC_GP);

    globalEnt = -entropy(mean(hC_GP,2));
    avLocalEnt = -entropy(hC_GP)/P;
end

function e = entropy(A)
    e = sum(A(:).*log(A(:)+eps));
end