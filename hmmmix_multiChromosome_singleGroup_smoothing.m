function [smoothedSequencesM_KT, updatedMLE_transitionMatrices, updatedMLE_initStates] = hmmmix_multiChromosome_singleGroup_smoothing(logobslik, chromosomeIndices, transitionMatrices, initStates, pseudoCounts, tolerance, maxIter, verbose)

    % This function is intended to deal with only one group, but that group
    % can have more than one chromosome in it. There won't be any
    % parallelization at this level because, if there is any to be had, it
    % will be at the level of doing this for all groups g=1:G.

    % logobslik is K-by-T
    % chromosomeIndices is 1-by-T
    % transitionMatrices is a cell array with C entries of K-by-K
    % initStates is a cell array with C entries of 1-by-K
    % pseudoCounts is either
    %   a cell array of C entries of K-by-K matrices or
    %   a K-by-K matrix
    
    % There an implementation decision to be made concerning the
    % potentially unused chromosomes. Should I allocate the space for 10
    % transition matrices if we have chromosomes 7,8,10 only ? The other
    % alternative is to have an index, but this is a recipe for headaches.
    % I decided to use cell arrays as arguments, but it's not obvious if
    % it's the best decision.
    
    if nargin < 8
        verbose = false;
    end
    
    [K,T] = size(logobslik);
    assert(length(chromosomeIndices) == T)
    uniqChrInd = unique(chromosomeIndices);
    
    for c = uniqChrInd
        assert(all([K,K] == size(transitionMatrices{c})));
        initStates{c} = toCol(initStates{c});
        assert(all([K,1] == size(initStates{c})));
    end
    
    smoothedSequencesM_KT = nan(K,T);
    
    for c = uniqChrInd
        
        % if we have only one value for that chromosome, don't do anything
        currentChromsome_length = sum(chromosomeIndices==c);
        if currentChromsome_length == 1
            updatedMLE_transitionMatrices{c} = transitionMatrices{c};
            updatedMLE_initStates{c} = initStates{c};
            continue;
        end
        
        % set up the initial values for that chromosomes
        currentChromosome_logobslik = logobslik(:,chromosomeIndices==c);
        
        currentChromosome_initStates = initStates{c};
        currentChromosome_transitionMatrix = transitionMatrices{c};
        
        if verbose
            fprintf('Chromosome %d. Initialization. \n', c);
            disp(currentChromosome_initStates);
            disp(currentChromosome_transitionMatrix)
        end
 
        if iscell(pseudoCounts)
            currentChromosome_pseudoCounts = pseudoCounts{c};
        else 
            currentChromosome_pseudoCounts = pseudoCounts;
        end
        all([K,K] == size(currentChromosome_pseudoCounts));
        
        % Now we want to iterate until we reach an acceptably stable
        % solution. This is evaluated in terms of the marginal
        % distributions have stabilized.
        
        niter = 0;
        old_currentChromosome_smoothedMarginals = zeros(K,currentChromsome_length);
        while(niter < maxIter)
        
            %currentChromosome_logobslik, currentChromosome_obslik, currentChromosome_initStates, currentChromosome_transitionMatrix, currentChromosome_pseudoCounts
            
            [new_currentChromosome_smoothedMarginals, junk2, currentChromosome_initStates, currentChromosome_transitionMatrix, junk5] = ...
                hmmmix_frugal_hM_KTg_MatlabC(currentChromosome_logobslik, currentChromosome_initStates, currentChromosome_transitionMatrix, currentChromosome_pseudoCounts);    
            
            % I'm not sure if I should stick to testing the stability of
            % the transition matrices since this test looks a bit
            % expensive. It'd be a shame if we ended up padding the
            % iterations by a factor of about (K+1)/K just because of this long test. 
            if max(abs(new_currentChromosome_smoothedMarginals(:) - old_currentChromosome_smoothedMarginals(:))) < tolerance
                break;
            end
            old_currentChromosome_smoothedMarginals = new_currentChromosome_smoothedMarginals;
            niter = niter + 1;
            
            if verbose
                fprintf('Chromosome %d. Going for one more iteration. \n', c);
                disp(currentChromosome_initStates);
                disp(currentChromosome_transitionMatrix)
            end
        end

        if niter == maxIter
            fprintf('We reached the maximum number of iterations %d for chromosome %c. This might be a sign that there is a problem.\n', maxIter, c);
        end
        
        % Now that we've iterated, we write the values in the output
        % arguments. I'm following Sohrab's standard of putting all the
        % chromosomes in the same array, tagging them with another
        % variable. I'm assuming that they are ordered, though, in that
        % array, and that we are interested in having the smoothed
        % marginals be in the same locations.
        
        smoothedSequencesM_KT(:,chromosomeIndices==c) = new_currentChromosome_smoothedMarginals;
        updatedMLE_transitionMatrices{c} = currentChromosome_transitionMatrix;
        updatedMLE_initStates{c} = currentChromosome_initStates;
    end
    
    
end