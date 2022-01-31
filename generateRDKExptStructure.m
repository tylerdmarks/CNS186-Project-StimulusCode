function expt = generateRDKExptStructure(preset_name)
    
   
    if strcmp(preset_name, 'Sequence report')
        direction_pool = [0, 90, 180, 360];
        lengths = [4 6 8];
        
        expt.numTrials = 120;
        expt.field = zeros(1, expt.numTrials);      % inferior VF = 0, superior VF = 1
        sup_trials = binaryrandomize(expt.numTrials);        % select a random half of trials to present upper
        expt.field(sup_trials) = 1;
        expt.sequence = cell(1, expt.numTrials);       % contains sequence for every trial
        sequence_lengths = zeros(1, expt.numTrials);    % lengths of each of the sequences (random thirds)
        random_lengths = tertiaryrandomize(expt.numTrials);
        sequence_lengths(random_lengths == 1) = lengths(2);
        sequence_lengths(random_lengths == 2) = lengths(3);
        sequence_lengths(sequence_lengths == 0) = lengths(1);
        for ss = 1:expt.numTrials
            curr_length = sequence_lengths(ss);
            random_directions = datasample(direction_pool, curr_length);
            expt.sequence{ss} = random_directions;
        end
        
    elseif strcmp(preset_name, 'Orientation report')
        % two direction pools: left, right & up, down
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
