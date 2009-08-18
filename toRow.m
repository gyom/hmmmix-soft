function v = toRow(v0)
    if(size(v0,2) == 1)
        v = v0';
    else
        v = v0;
    end
end
        