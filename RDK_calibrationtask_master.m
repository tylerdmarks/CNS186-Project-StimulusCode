%% Setup
clear
KbName('UnifyKeyNames');        % fixes an error calling KbName('ESCAPE'), might be a Windows thing

screenID = 0;
testing = 0;
% File name, change for each new participant       'subjectnumber_nameoftask_date'
ID = 'S001test_RDKcalibrationtask_220203_2'; 
% Get monitor's refresh rate
fr = Screen('NominalFrameRate', screenID);          
% set background color for drawing
background = [150 150 150];

% set up save directory and save file for experiment data
fullpath = ('C:\Experiment Data\CNS186 project');            % operational folder for this computer
data_filename = sprintf('%s.mat', ID);

% set up timings and parameters
dot_properties.size = 32;
dot_properties.speed = 150;
dot_properties.ndots = 200;
dot_properties.lifetime = 45;
start_coherence = 0.95;
coherence_increment = 0.05;
dot_properties.color = [0 0 0];
aperture.size = [2400 900];             % size of aperture (ellipse, [x, y] dimensions)
RDKduration = 0.6;                        % duration of each element in RDK sequence (seconds)
ISI = 0.1;                              % duration of interval in between RDK sequence elements

% Fixation parameters
fix.size = 24; % pixels, cross
fix.width = 4; % pixels, line width
fix.color = [0 0 0]; % color of fixation
fix.duration = 0.8;       % duration of initial fixation

% Generate experiment information
expt = generateRDKExptStructure('Calibration');

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

% display properties struct
display.dist = 3;   % cm        distance from screen?
display.width = 14.1;   % cm    width of screen?
display.resolution = displaySize;
display.framerate = fr;

% set fix and aperture positions
fix.pos = winCtr;          % set fixation position to center of monitor
X = [-fix.size fix.size 0 0];
Y = [0 0 -fix.size fix.size];
fix.coords = [X; Y];

% these valueus are distance from center of screen (origin)
aperture.pos(1, :) = [0 -displaySize(2)*0.25];     % superior field RDK aperture (first index)
aperture.pos(2, :) = [0 displaySize(2)*0.25];     % inferior field RDK aperture (second index)

% create data storage for responses
response.inf_CohTraj = [];       % store coherence trajectory for inferior VF
response.sup_CohTraj = [];    % store coherence trajectory for superior VF
% Set text size
Screen('TextSize', ptr, 64);
%% Task preparation
while true
    % Gray background
    Screen('FillRect', ptr, background)
    
    % Draw the fixation cross
    Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
    % Draw title message
    DrawFormattedText(ptr, 'Look at the fixation cross, then press the space bar.', 'center', winCtr(2) - displaySize(2)/4, fix.color);
    % Render to screen
    Screen('Flip', ptr);
    
    [~, ~, keycode, ~] = KbCheck(-1);
    if keycode(KbName('Space'))
        break
    end
end

%% TASK LOOP
reversal_threshold = 6;
threshold_reached = false;  
num_reversal_inf = 0;           % counts number of reversals in inferior visual field
num_reversal_sup = 0;           % counts number of reversals in superior visual field
trial_counter = 1;
coherence(1) = start_coherence;     % start both fields at starting coherence (0 = inf, 1 = sup)
coherence(2) = start_coherence;
coherence_log = cell(2, 1);     % initialize records of coherence trace for both hemifields
inf_threshold_reached = false;      % track these separately, when threshold is reached then stop incrementing for the hemifield
sup_threshold_reached = false;
while ~threshold_reached    
    % every 10 trials give a break
    if mod(trial_counter, 10) == 1 && trial_counter ~=1
        % Gray background
        Screen('FillRect', ptr, background)
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
        % Draw title message
        DrawFormattedText(ptr, 'Look at the fixation cross, then press the space bar.', 'center', winCtr(2) - displaySize(2)/4, fix.color);
        % Render to screen
        Screen('Flip', ptr);

        [~, ~, keycode, ~] = KbCheck(-1);
        if keycode(KbName('Space'))
            break
        end
    end
    
    % Set state info
    field = expt.field(trial_counter);   % which visual field (inferior = 0, superior = 1)
    direction = expt.sequence{trial_counter}; % direction for this trial

    %% INITIAL FIXATION
    preFixTime = GetSecs;
    while GetSecs - preFixTime < fix.duration
        % Gray background
        Screen('FillRect', ptr, background)
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
        % Render to screen
        Screen('Flip', ptr);
    end
    
    %% RDK PRESENTATION
    curr_coherence = coherence(field+1);
    aperture_center = aperture.pos(field+1, :);

    dot_properties = initializeDotProperties(dot_properties, direction, curr_coherence); 
        
    dot_properties.x = (rand(1,dot_properties.ndots)-.5)*aperture.size(1) + aperture_center(1);
    dot_properties.y = (rand(1,dot_properties.ndots)-.5)*aperture.size(2) + aperture_center(2);

    pixpos.x = dot_properties.x;
    pixpos.y = dot_properties.y;

    pixpos.x = pixpos.x + winCtr(1);
    pixpos.y = pixpos.y + winCtr(2);

    % dot movement
    dx = dot_properties.speed*sin(dot_properties.direction*pi/180)/fr;
    dy = -dot_properties.speed*cos(dot_properties.direction*pi/180)/fr;

    dx = angle2pix(display,dx);
    dy = angle2pix(display,dy);
    l = aperture_center(1)-aperture.size(1)/2;
    r = aperture_center(1)+aperture.size(1)/2;
    b = aperture_center(2)-aperture.size(2)/2;
    t = aperture_center(2)+aperture.size(2)/2;

    dot_properties.life = ceil(rand(1,dot_properties.ndots)*dot_properties.lifetime);
    preElementTime = GetSecs;
    while GetSecs-preElementTime < RDKduration
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);

        %convert from degrees to screen pixels
        pixpos.x = dot_properties.x+ display.resolution(1)/2;
        pixpos.y = dot_properties.y+ display.resolution(2)/2;

        %Use the equation of an ellipse to determine which dots fall inside.
        goodDots = (dot_properties.x-aperture_center(1)).^2/(aperture.size(1)/2)^2 + ...  
            (dot_properties.y-aperture_center(2)).^2/(aperture.size(2)/2)^2 < 1;

        % if you're using an ellipse, add pixpos.x(goodDots) and for y too
        Screen('DrawDots', ptr, [pixpos.x(goodDots); pixpos.y(goodDots)],  dot_properties.size, dot_properties.color, [0, 0], 3);
        %update the dot position
        dot_properties.x = dot_properties.x + dx;
        dot_properties.y = dot_properties.y + dy;

        %move the dots that are outside the aperture back one aperture
        %width.
        dot_properties.x(dot_properties.x<l) = dot_properties.x(dot_properties.x<l) + aperture.size(1);
        dot_properties.x(dot_properties.x>r) = dot_properties.x(dot_properties.x>r) - aperture.size(1);
        dot_properties.y(dot_properties.y<b) = dot_properties.y(dot_properties.y<b) + aperture.size(2);
        dot_properties.y(dot_properties.y>t) = dot_properties.y(dot_properties.y>t) - aperture.size(2);

        %increment the 'life' of each dot
        dot_properties.life = dot_properties.life+1;

        %find the 'dead' dots
        deadDots = mod(dot_properties.life,dot_properties.lifetime)==0;

        %replace the positions of the dead dots to a random location
        dot_properties.x(deadDots) = (rand(1,sum(deadDots))-.5)*aperture.size(1) + aperture_center(1);
        dot_properties.y(deadDots) = (rand(1,sum(deadDots))-.5)*aperture.size(2) + aperture_center(2);

        Screen('Flip', ptr);
    end
        
    % ISI
    preISITime = GetSecs;
    while GetSecs - preISITime < ISI
        % Gray background
        Screen('FillRect', ptr, background)
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);

        % Render to screen
        Screen('Flip', ptr);
    end
        
    [~, ~, keycode, ~] = KbCheck(-1);
    if keycode(KbName('ESCAPE')) || keycode(KbName('q'))
        % Clear screen
        sca

        % Save all workspace variables
        save(data_filename)
        return
    end

    %% RESPONSE PERIOD
   
    % Gray background
    Screen('FillRect', ptr, background)
    % Instructions
    DrawFormattedText(ptr, 'Which direction was the motion?', 'center', winCtr(2) - displaySize(2)/4, fix.color);

    % Draw response options (arrows at set locations)
    drawArrow(ptr, 'up', [winCtr(1) winCtr(2)-displaySize(2)/30])
    drawArrow(ptr, 'left', [winCtr(1)-displaySize(1)/50 winCtr(2)])
    drawArrow(ptr, 'down', [winCtr(1) winCtr(2)+displaySize(2)/30])
    drawArrow(ptr, 'right', [winCtr(1)+displaySize(1)/50 winCtr(2)])

    % render everything
    Screen('Flip', ptr);

    % check response
    [~, keycode, ~] = KbStrokeWait(-1);
    if keycode(KbName('LeftArrow'))
        curr_response = 270;
    elseif keycode(KbName('RightArrow'))
        curr_response = 90;
    elseif keycode(KbName('UpArrow'))
        curr_response = 0;
    elseif keycode(KbName('DownArrow'))
        curr_response = 180;
    elseif keycode(KbName('ESCAPE'))
        sca
        % Save all workspace variables
        save(data_filename)
        return
    end  
    
    if field == 0 && ~inf_threshold_reached             %only increment if threshold hasn't been reached yet
        if curr_response == direction           % if response is correct
            coherence(field+1) = coherence(field+1)-coherence_increment;    % decrease coherence for the hemifield
        else        % if response is incorrect
            coherence(field+1) = coherence(field+1)+coherence_increment;    % increase coherence for the hemifield
        end
        coherence_log{field+1} = [coherence_log{field+1} coherence(field+1)];       % log the coherence value
    elseif field == 1 && ~sup_threshold_reached
        if curr_response == direction           % if response is correct
            coherence(field+1) = coherence(field+1)-coherence_increment;    % decrease coherence for the hemifield
        else        % if response is incorrect
            coherence(field+1) = coherence(field+1)+coherence_increment;    % increase coherence for the hemifield
        end
        coherence_log{field+1} = [coherence_log{field+1} coherence(field+1)];       % log the coherence value
    end
    
    % correct for dipping below coherence of 0.05 or above coherence of 1
    if coherence(field+1) < 0.05
        coherence(field+1) = 0.05;
        coherence_log{field+1}(end) = 0.05;
    elseif coherence(field+1) > 1
        coherence(field+1) = 1;
        coherence_log{field+1}(end) = 1;
    end
       
    % draw subject's response
    switch curr_response
        case 0
            drawArrow(ptr, 'up', [winCtr(1) winCtr(2)-displaySize(2)/30])
        case 90
            drawArrow(ptr, 'right', [winCtr(1)+displaySize(1)/50 winCtr(2)])
        case 180
            drawArrow(ptr, 'down', [winCtr(1) winCtr(2)+displaySize(2)/30])
        case 270
            drawArrow(ptr, 'left', [winCtr(1)-displaySize(1)/50 winCtr(2)])
    end
    DrawFormattedText(ptr, 'Which direction was the motion?', 'center', winCtr(2) - displaySize(2)/4, fix.color);
    Screen('Flip', ptr);
    pause(0.5)
    
    % determine if reversal
    % if last value equals 2 values ago and does not equal 1 value ago
    if length(coherence_log{field+1,:}) > 3
        if coherence_log{field+1}(end) == coherence_log{field+1}(end-2) && coherence_log{field+1}(end) ~= coherence_log{field+1}(end-1)
            if field == 0  
                num_reversal_inf = num_reversal_inf + 1;
            elseif field == 1 
                num_reversal_sup = num_reversal_sup + 1;
            end
        else            % written as: must be consecutive reversals, get rid of this else to just make it total reversals
            if field == 0
                num_reversal_inf = 0;
            elseif field == 1
                num_reversal_sup = 0;
            end
        end
    end
    
    % if reversal threshold has been reached for the hemifield, indicate as such
    if num_reversal_inf == reversal_threshold
        inf_threshold_reached = true;
    end
    
    if num_reversal_sup == reversal_threshold
        sup_threshold_reached = true;
    end
    
    % if threshold for both hemifields has been reached, terminate the loop
    if inf_threshold_reached && sup_threshold_reached
        threshold_reached = true;
    end
    trial_counter = trial_counter+1;
    
end
% clear screen

sca

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

    is_already_direction = (dot_properties.random_directions == dot_properties.coherence_direction);
    coherent_dot_idx = datasample(1:dot_properties.ndots,round(coherence * dot_properties.ndots), 'Replace', false);
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

% Other subfunctions
function drawArrow(ptr, direction, headloc)
    % Draws an arrow using FillPoly and DrawLines
    % direction = direction the arrowheadloc is facing
    % location = [x,y] location of arrowheadloc
    
    % create a triangle
    arrowwidth  = 60;           % width of arrow headloc
    linelength = 50;
    switch direction
        case 'up'
        tripoints = [ headloc-[arrowwidth/2,0]         
                   headloc+[arrowwidth/2,0]         
                   headloc-[0,arrowwidth] ];
        linecoords = [headloc(1) headloc(2)+linelength headloc(1) headloc(2)];
        case 'down'
        tripoints = [ headloc-[arrowwidth/2,0]       
                   headloc+[arrowwidth/2,0]        
                   headloc+[0,arrowwidth] ];
        linecoords = [headloc(1) headloc(2)-linelength headloc(1) headloc(2)];
        case 'left'
        tripoints = [ headloc-[0,arrowwidth/2]       
                   headloc+[0,arrowwidth/2]        
                   headloc-[arrowwidth,0] ];    
        linecoords = [headloc(1) headloc(2) headloc(1)+linelength headloc(2)];
        case 'right' 
        tripoints = [ headloc-[0,arrowwidth/2]         
                   headloc+[0,arrowwidth/2]         
                   headloc+[arrowwidth,0] ]; 
        linecoords = [headloc(1)-linelength headloc(2) headloc(1) headloc(2)];
    end
    Screen('FillPoly', ptr, [0 0 0], tripoints);
    Screen('DrawLine', ptr, [0 0 0], linecoords(1), linecoords(2), linecoords(3), linecoords(4), 10);
        
end


