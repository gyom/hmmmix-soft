function value = getFieldOrDefault(inSt, fieldname, defaultvalue)

    % If the input structure "inSt" doesn't have a field with name
    % "fieldname", then return "defaultvalue". Otherwise, return the value
    % of that field.
    %
    % Meant to be used as :
    %       initialAssignments = getFieldOrDefault(dataSt, 'init', ones(M,N));
    
    if ~hasField(inSt, fieldname)
        value = defaultvalue;
    else
        value = inSt.(fieldname);
    end

end