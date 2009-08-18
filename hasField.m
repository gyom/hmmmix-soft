function v = hasField(struct, fieldname)

    % a arguably not-that-useful function that checks if a structure has a
    % field of name 'fieldname'

    % I forgot about the existence of strcmp when I wrote this function.
    
    v = 0;
    thefields = fields(struct);
    
    for r=1:length(thefields)
        % skip it if they don't have the same length
        if any(size(thefields{r}) ~= size(fieldname))
            continue;
        end
        
        v = v | all(thefields{r}==fieldname);
    end

end