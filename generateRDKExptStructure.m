function expt = generateRDKExptStructure(preset_name)
    
   
    if strcmp(preset_name, 'Expt1')
        % Set number of blocks. Block types are as follows:
        %   1  :  all prosaccade, 50% chance primer in target loc, 50% no primer
        %   2  :  all antisacccade, 50% chance primer in target loc, 50% no primer
        
        expt.numTrials = 120;
       
        % 
        
        expt.presVFsup = zeros(1, expt.numTrials);  % logical vector indicating whether to present an RDK in superior field
        expt.presVFinf = zeros(1, expt.numTrials);  % logical vector indicating whether to present in inferior field
        right_trials = binaryrandomize(expt.block(1).numTrials);        % select a random half of trials to cue right
        expt.block(1).trialOrder(right_trials, 1) = 1;          
        expt.block(1).trialOrder(:, 2) = 0;             % all prosaccade trials 
        primer_trials = binaryrandomize(expt.block(1).numTrials);   % select random half of trials for primer presentation
        expt.block(1).trialOrder(primer_trials, 3) = 1;             % present primer on 50% of trials

        % Block 2
        expt.block(2).numTrials = 60;   
        expt.block(2).trialOrder = zeros(expt.block(2).numTrials, 3);  
        right_trials = binaryrandomize(expt.block(2).numTrials);        % select a random half of trials to cue right
        expt.block(2).trialOrder(right_trials, 1) = 1; 
        expt.block(2).trialOrder(:, 2) = 1;             % all antisaccade trials
        primer_trials = binaryrandomize(expt.block(2).numTrials);
        expt.block(2).trialOrder(primer_trials, 3) = 1;             % present primer on 50% of trials
    end


end

function idx = binaryrandomize(numTrials)      
    % outputs a random half of indices from the vector 1:numTrials
    shuffle = randperm(numTrials);
    idx = shuffle(1:floor(length(shuffle)/2));
end

function idx = tertiaryrandomize(numTrials)
    % outputs a vector of length numTrials where 1/3 of indices are 0, 1, or 2, evenly distributed
    onethird = floor(numTrials/3);
    r = mod(numTrials, 3);
    extra0 = 0;
    extra1 = 0;
    extra2 = 0;

    while r > 0
        a = randi(3);
        switch a 
        case 1 
            if extra0 ~= 1
                extra0 = 1;
                r = r-1;
            end
        case 2 
            if extra1 ~= 1
                extra1 = 1;
                r = r-1;
            end
        case 3
            if extra2 ~= 1
                extra2 = 1;
                r = r-1;
            end
        end
    end

    numzeros = onethird + extra0;
    numones = onethird + extra1;
    numtwos = onethird + extra2;

    shuffle = randperm(numTrials);
    idx(shuffle(1:numzeros)) = 0;
    idx(shuffle(numzeros+1:numzeros+numones)) = 1;
    idx(shuffle(numzeros+numones+1:end)) = 2;
end
