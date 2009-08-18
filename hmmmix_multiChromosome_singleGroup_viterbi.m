function viterbiPathIndices_concatenatedChromosomes = hmmmix_multiChromosome_singleGroup_viterbi(logobslik, chromosomeIndices, transitionMatrices, initStates)

    % This function is a viterbi version of the
    % hmmmix_multiChromosome_singleGroup_smoothing function and it's just
    % mostly a remix of the code. If you ever want to correct a mistake in
    % the code, correct it in the previous function and then carry the
    % changes in this one. You are not supposed to be maintaining these two
    % functions separately.

    [K,T] = size(logobslik);
    assert(length(chromosomeIndices) == T)
    uniqChrInd = unique(chromosomeIndices);
    
    for c = uniqChrInd
        assert(all([K,K] == size(transitionMatrices{c})));
        initStates{c} = toCol(initStates{c});
        assert(all([K,1] == size(initStates{c})));
    end
    
    viterbiPathIndices_concatenatedChromosomes = nan(1,T);
    
    for c = uniqChrInd
        
        % if we have only one value for that chromosome, don't do anything
        currentChromsome_length = sum(chromosomeIndices==c);
        if currentChromsome_length == 1
            updatedMLE_transitionMatrices{c} = transitionMatrices{c};
            updatedMLE_initStates{c} = initStates{c};
            continue;
        end
        
        currentChromosome_logobslik = logobslik(:,chromosomeIndices==c);
        currentChromosome_initStates = initStates{c};
        currentChromosome_transitionMatrix = transitionMatrices{c};
        
        [currentChromosome_viterbi, junk2, junk3] = viterbi_path_MatlabC(log(currentChromosome_initStates), log(currentChromosome_transitionMatrix+eps), currentChromosome_logobslik);
        
        %I'm following Sohrab's standard of putting all the
        % chromosomes in the same array, tagging them with another
        % variable. I'm assuming that they are ordered, though, in that
        % array, and that we are interested in having the smoothed
        % marginals be in the same locations.
        
        viterbiPathIndices_concatenatedChromosomes(1,chromosomeIndices==c) = currentChromosome_viterbi;
    end
    
    
    
    
end