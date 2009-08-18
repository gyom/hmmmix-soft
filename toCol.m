function v = toCol(v0)
    if(size(v0,1) == 1)
        v = v0';
    else
        v = v0;
    end
end
        