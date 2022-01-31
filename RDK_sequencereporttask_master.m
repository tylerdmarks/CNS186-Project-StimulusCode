%% Setup
clear
KbName('UnifyKeyNames');        % fixes an error calling KbName('ESCAPE'), might be a Windows thing

screenID = 0;
testing = 1;
% File name, change for each new participant       'subjectnumber_nameoftask_date'
ID = 'S001_RDKsequencetask_220120'; 
% Get monitor's refresh rate
fr = Screen('NominalFrameRate', screenID);          
% set background color for drawing
background = [150 150 150];

% set up save directory and save file for experiment data
fullpath = ('C:\Experiment Data\CNS186 project');            % operational folder for this computer
data_filename = sprintf('%s.mat', ID);

% set up timings and parameters
dot_properties.size = 20;
dot_properties.speed = 80;
dot_properties.ndots = 200;
dot_properties.lifetime = 45;
dot_properties.coherence(1) = 0.90;       % coherence for inferior visual field (first index)
dot_properties.coherence(2) = 0.90;       % coherence for superior visual field (second index)
dot_properties.color = [0 0 0];
aperture.size = [1600 800];             % size of aperture (ellipse, [x, y] dimensions)
RDKduration = 1;                        % duration of each element in RDK sequence (seconds)
ISI = 1;                              % duration of interval in between RDK sequence elements

% Fixation parameters
fix.size = 24; % pixels, cross
fix.width = 4; % pixels, line width
fix.color = [0 0 0]; % color of fixation
fix.duration = 1;       % duration of initial fixation

% Generate experiment information
expt = generateRDKExptStructure('Sequence report');

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
response.time = zeros(1, expt.numTrials);       % store time it takes to complete response
response.sequence = cell(1, expt.numTrials);    % store response sequences

%% Task preparation
while true
    % Gray background
    Screen('FillRect', ptr, background)
    
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
        % Gray background
        Screen('FillRect', ptr, background)
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
        % Draw title message
        DrawFormattedText(ptr, 'Look at the fixation cross, then press the space bar.', winCtr(1)-displaySize(1)/4, winCtr(2) - displaySize(2)/4, fix.color);
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
    preFixTime = GetSecs;
    while GetSecs - preFixTime < fix.duration
        % Gray background
        Screen('FillRect', ptr, background)
        % Draw the fixation cross
        Screen('DrawLines', ptr, fix.coords, fix.width, fix.color, fix.pos, 2);
        % Render to screen
        Screen('Flip', ptr);
    end
    
    %% SEQUENCE PRESENTATION
    coherence = dot_properties.coherence(field+1);
    aperture_center = aperture.pos(field+1, :);
    for ss = 1:length(sequence)
        direction = sequence(ss);
        dot_properties = initializeDotProperties(dot_properties, direction, coherence); 
        
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
            
            drawArrow(ptr, 'left')
            drawArrow(ptr, 'right') 

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
            
            drawArrow(ptr, 'left')
            drawArrow(ptr, 'right')
            
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
        
    end

    %% RESPONSE PERIOD
    % log duration of response period and response sequence
    if ~testing
    response_complete = false;
    curr_response = [];
    feedback_startpoint_x = [];     % starting x location for feedback display, based on sequence length
    feedback_startpoint_y = [];     % starting y location for feedback
    offset = [];            % element distancing in x direction
    preResptime = GetSecs;
    while ~response_complete
        % if response vector is equal to length of sequence, response is complete
        if length(curr_response) == length(sequence)
            response_complete = true;
        end
        
        % Gray background
        Screen('FillRect', ptr, background)
        % Instructions
        DrawFormattedText(ptr, 'Repeat the sequence with the arrow keys.', winCtr(1)-displaySize(1)/4, winCtr(2) - displaySize(2)/4, fix.color);
        % Draw response options (arrows at set locations)
        
        % Draw 
        for rr = 1:length(response)
            % get x location (add set distance per element)
            
            % draw element
        end
        
        % render everything
        Screen('Flip', ptr);
        
        % check response
        [~, ~, keycode, ~] = KbCheck(-1);
        if keycode(KbName('LeftArrow'))
            curr_response = [curr_response 270];
        elseif keycode(KbName('RightArrow'))
            curr_response = [curr_response 90];
        elseif keycode(KbName('UpArrow'))
            curr_response = [curr_response 0];
        elseif keycode(KbName('DownArrow'))
            curr_response = [curr_response 180];
        elseif keycode(KbName('ESCAPE'))
            sca
            % Save all workspace variables
            save(data_filename)
            return
        end
        
        if response_complete
            pause(1)
            curr_resptime = GetSecs-preResptime;
        end
    end
    response.time(tt) = curr_resptime;
    response.sequence{tt} = curr_response;
    end
    
   
    
    
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
function drawArrow(ptr, direction, location)
    % Draws an arrow using FillPoly and DrawLines
    % direction = direction the arrowhead is facing
    % location = [x,y] location of arrowhead
    
    % create a triangle
    head   = [ 600, 600 ]; % coordinates of head
    width  = 40;           % width of arrow head
    switch direction
        case 'up'
        points = [ head-[width/2,0]         % left corner
                   head+[width/2,0]         % right corner
                   head-[0,width] ];      % vertex
        case 'down'
        points = [ head-[width/2,0]         % left corner
                   head+[width/2,0]         % right corner
                   head+[0,width] ];      % vertex
               
        case 'left'
        points = [ head-[0,width/2]         % left corner
                   head+[0,width/2]         % right corner
                   head-[width,0] ];      % vertex
        case 'right' 
        points = [ head-[0,width/2]         % left corner
                   head+[0,width/2]         % right corner
                   head+[width,0] ];      % vertex
    end
    Screen('FillPoly', ptr, [0 0 0], points);
        
end


