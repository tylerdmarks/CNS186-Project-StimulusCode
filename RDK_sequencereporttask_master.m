%% Setup
clear
KbName('UnifyKeyNames');        % fixes an error calling KbName('ESCAPE'), might be a Windows thing

screenID = 0;
% File name, change for each new participant       'subjectnumber_nameoftask_date'
ID = 'S001_sequencetask_220120'; 
fr = Screen('NominalFrameRate', screenID);          % get monitor's framerate

% set up save directory and save file for experiment data
fullpath = ('C:\Experiment Data\CNS186 project');            % operational folder for this computer
data_filename = sprintf('%s.mat', ID);

% set up timings and parameters
dot_properties.size = 10;
dot_properties.speed = 80;
dot_properties.ndots = 200;
dot_properties.lifetime = 45;
dot_properties.coherence(1) = 20;       % coherence for inferior visual field (first index)
dot_properties.coherence(2) = 20;       % coherence for superior visual field (second index)
dot_properties.color = [0 0 0];
aperture.size = [1000 500];             % size of aperture (ellipse, [x, y] dimensions)
RDKduration = 1;                        % duration of each element in RDK sequence (seconds)
ISI = 0.3;                              % duration of interval in between RDK sequence elements

% Fixation parameters
fix.size = 16; % pixels, cross
fix.width = 3; % pixels, line width
fix.color = [250 250 250]; % usual color of fixation (gray)
fix.duration = 1;       % duration of initial fixation

% display properties
display.dist = 3;   % cm        distance from screen?
display.width = 14.1;   % cm    width of screen?
display.resolution = displaySize;
display.framerate = fr;

% Generate experiment information
expt = generateRDKExptStructure('Expt1');

%% Initialize

% Shift focus to command window so typing doesn't go to code editor if something
% happens to break
commandwindow
HideCursor;

% initialize Psychtoolbox
% can't run PTB on laptop without skipping sync tests
Screen('Preference', 'SkipSyncTests', 1);  
AssertOpenGL;
[ptr,winrect] = Screen('OpenWindow', screenID, 0);
Screen('BlendFunction', ptr, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

displaySize = winrect(3:4);    % get screen size
winCtr = displaySize/2;    % get center position of the screen

% set fix and aperture positions
fix.pos = env.winCtr;          % set fixation position to center of monitor
X = [-fix.size fix.size 0 0];
Y = [0 0 -fix.size fix.size];
fix.coords = [X; Y];

aperture.pos(1) = [displaySize(1)/2 displaySize(2)*0.25];     % superior field RDK aperture (first index)
aperture.pos(2) = [displaySize(1)/2 displaySize(2)*0.75];     % inferior field RDK aperture (second index)

% create data storage for responses
response.time = zeros(1, expt.numTrials);       % store time it takes to complete response
response.sequence = cell(1, expt.numTrials);    % store response sequences

%% Task preparation
while true
    % Draw the fixation cross
    Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
    % Draw title message
    DrawFormattedText(ptr, 'Look at the fixation cross, then press the space bar.', winCtr(1), winCtr(2) - displaySize(2)/4, fix.color);
    % Render to screen
    Screen('Flip', ptr);
    
    [~, ~, keycode, ~] = KbCheck(-1);
    if keycode(KbName('Space'))
        break
    end
end

%% TASK LOOP

for tt = 1:expt.numTrials
    
    % every 10 trials give a break
    if mod(tt, 10) == 1 && tt ~=1
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
        % Draw title message
        DrawFormattedText(ptr, 'Look at the fixation cross, then press the space bar.', winCtr(1), winCtr(2) - displaySize(2)/4, fix.color);
        % Render to screen
        Screen('Flip', ptr);

        [~, ~, keycode, ~] = KbCheck(-1);
        if keycode(KbName('Space'))
            break
        end
    end
    
    % Set state info
    field = expt.field(tt);   % which visual field (inferior = 0, superior = 1)
    sequence = expt.sequence{tt}; % sequence of directions (orientations)

    %% INITIAL FIXATION
    preFixtime = GetSecs;
    while GetSecs - preFixTime < fix.duration
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
        % Render to screen
        Screen('Flip', ptr);
    end
    
    %% SEQUENCE PRESENTATION
    coherence = dot_properties.coherence(field+1);
    for ss = 1:length(sequence)
        direction = sequence(ss);
        dot_properties = initializeDotProperties(dot_properties, direction, coherence); 
        
        preElementTime = GecSecs;

        while GetSecs-preElementTime < RDKduration
            % present the RDK
        end
        
        preISITime = GetSecs;
        while GetSecs - preISITime < ISI
            % present ISI blank
        end
        
    end
            
    
    
    
    %% RESPONSE PERIOD
    % log duration of response period and response sequence

    
    
end

% Save all workspace variables
try cd('data');
catch
    mkdir('data');
    cd('data');
end
save(data_filename);
  
    
% RDK subfunctions
function dot_properties = initializeDotProperties(dot_properties, direction, coherence)

    direction_pool = 0:45:315;      % possible random directions
    
    dot_properties.random_directions = datasample(direction_pool, dot_properties.ndots);

    dot_properties.coherence_direction = direction;
    dot_properties.coherence = coherence;

    is_already_direction = (dot_properties.random_directions == dot_properties.coherence_direction);
    coherent_dot_idx = datasample(1:dot_properties.ndots,round(dot_properties.coherence * dot_properties.ndots), 'Replace', false);
    is_coherent = false(1, dot_properties.ndots);
    is_coherent(coherent_dot_idx) = true; 

    dot_properties.direction = dot_properties.random_directions;
    dot_properties.direction(is_coherent & ~is_already_direction) = dot_properties.coherence_direction;
end

function pix = angle2pix(display_properties, ang)
    %Calculate pixel size
    pixSize = display_properties.width/display_properties.resolution(1);   %cm/pix

    sz = 2*display_properties.dist*tan(pi*ang/(2*180));  %cm

    pix = round(sz/pixSize);   %pix

end


