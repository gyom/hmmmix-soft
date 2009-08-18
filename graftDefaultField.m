function outSt = graftDefaultField(inSt, fieldname, defaultvalue)

    % If the input structure "inSt" doesn't have a field with name
    % "fieldname", then add it (setting it to "defaultvalue") and return
    % the new structure.
    %
    % Meant to be used as :
    %       dataSt = graftDefaultField(dataSt, 'iterations', 20);

    outSt = inSt;
    
    if ~hasField(outSt, fieldname)
        outSt.(fieldname) = defaultvalue; % = setfield(outSt, fieldname, defaultvalue);
    end

end