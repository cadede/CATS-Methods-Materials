%% Identify simple lunge cues and foraging strategies
% David Cade, James Fahlbusch and Ross Nichols
% version 2.0.0
% Written in Matlab 2014a
% Tested in 2014a and Matlab 2020b
% Goldbogen Lab and Friedlaender Lab
% Stanford University & University of California, Santa Cruz


% Necessary Scripts and Packages
% 	CATSMatlabTools package
% Matlab Statistical Package
% Subfunctions: savedetection.m


% Import Files:
% prh.mat
% lunge.mat
% strategy.mat

% Export Files and Variables:
% lunges.mat
% LungeDN - Lunge Date Number
% LungeI - Lunge Index
% LungeDepth - Lunge Depth (m)
% Lunge C - Lunge Confidence Scores
% 1 - Maybe
% 2 - Likely
% 3 - High Confidence
% creator - taken from step 1, enter manually in code
% primary_cue - taken from step 1, enter manually in code
% prh_fs - sampling rate of prh file
% starttime - Time of first datapoint in prh file
% created_on - datetime creation of lunge.mat file
% progressIndex - index value of detection progression
% notes - taken from step 1, enter manually in code
% strategy.mat
% StrategyS - Index of Start Markers for Bubble Netting
% StrategyE - Index of End Markers for Bubble Netting
% StrategyC - Confidence score Index for StrategyS and StrategyE
% 1 - High Confidence
% 2 - Low COnfidence
% creator - taken from step 1, enter manually in code
% primary_cue - taken from step 1, enter manually in code
% prh_fs - sampling rate of prh file
% starttime - Time of first datapoint in prh file
% creates_on - datetime creation of lunge.mat file
% progressIndex - index value of detection progression
% notes - taken from step 1, enter manually in code



% Using the Detection Tool - Lunges and Strategy
% In addition to lunge detection, this tool was originally made to denote periods of interest related
% to foraging strategy, namely, bubble net feeding in humpbacks.
% However, this tool allows for manual marking of any specific behavior
% and can be used for more purposes than bubble net feeding detection.
% Lunges may also be detected, as single point markers simultenous to
% strategy markers and is compatible with SimpleLungeDetection.m in
% export and import files.

% TO USE:
% 1.) Fill our data fields (primary_cue,creator,notes), then run Section 1 and select the PRH file you wish to process. If
% you wish to continue where you left on either a lunges.mat or
% strategy.mat file, then those files must be contained within the
% same folder as your PRH for them to load correctly.

% 2.) Once the data is loaded, proceed to Section 2 and run the Plot and Process Section
% Use Mouse and keyboard shortcuts listed on the Left hand side of
% the screen to navigate and place marks. Marks will be plotted
% automatically

% 3.) For Strategy Marking, you can perform 3 main functions:
%.) Q - High Confidence Marks
% Pressing Q will place the High Confidence START  Marker (green),
% indicating the start of the behavior. You may then press
% any key to place the END marker (red).
%.) E - Low Confidence Marks
% Pressing E will place the Low Confidence START  Marker (cyan),
% indicating the start of the behavior. You may then press
% any key to place the END marker (magenta).
%.) R - Delete Set of Marks - Pressing R will delete the nearest
%set of strategy marks, both START and END.

% 4.) For Lunge Marking
%.) Left Click - High Confidence Lunge (Score of 3)
% Indicates a high confidence lunge
%.) L - Likely a Lunge (Score of 3)
% Indicates a moderate confidence score lunge
%.) M - Maybe a Lunge (Score of 3)
%  Indicates a low confidence score lunge
%.) Right Click - Delete Nearest Lunge Marker

% System Commands
%.) S - Save. This function saves a lunge.mat and strategy.mat
%file in the folder containing the loaded prh. This overwrites
%any existing file of the same name.
%.) Return - Go Forward.
%.) B - Go Backward.
%.) Numpad 1-9 - Change Zoom (x10)

%% 1. Load Data
clear; % clears the workspace
%Start where you left off?
atLast = true; %this will look for a variable called progressIndex
M = 10; % number of minutes

% Variables that will be saved in the Lunge and Strategy file
notes = '';
creator = '';
primary_cue = 'speedJJ';


try
    drive = 'Gold Field'; % name of drive where files are located. (can be drive name or drive letter)
    folder = 'CATS/tag_data'; % folder in the drive where the cal files are located (and where you want to look for files) %'Users\Dave\Documents\Programs\MATLAB\Tagging\CATS cal';%
    % make finding the files easier
    a = getdrives;
    for i = 1:length(a)
        [~,vol]=system(['vol ' a{i}(1) ':']);
        if strfind(vol,drive); vol = a{i}(1); break; end
    end
catch
end

cf = pwd; try cd([vol ':/' folder]); catch; end
[filename,fileloc]=uigetfile('*.*', 'the PRH file to analyze'); %look below to save time to make truncated files okay
cd(fileloc);
disp('Loading Data, will take some time');
load([fileloc filename(1:end-3) 'mat']);
whaleName = INFO.whaleName;
ii = strfind(filename,' ');
clear LungeI LungeDN LungeC LungeDepth StI StC StDN StIE StEDN % clears lunge data in case a previous file was still loaded
lungename = [filename(1:ii-1) 'lunges.mat'];
strategyname = [filename(1:ii-1) 'strategy.mat'];
load([fileloc filename]);
try load([fileloc lungename]); catch; end
try load([fileloc strategyname]); catch; end
try
    speedFN = speed.FN;
    speedFN(isnan(speedFN)) = min(speedFN); speedFN = runmean(speedFN,round(fs/2)); speedFN(isnan(speed.FN)) = nan;
catch; speedFN = nan(size(p));
end
if ~exist('head', 'var')
    head = nan(size(p));
end
speedJJ = speed.JJ;
speedJJ(isnan(speedJJ)) = min(speedJJ); speedJJ = runmean(speedJJ,round(fs/2)); speedJJ(isnan(speed.JJ)) = nan;
J = njerk(Aw,fs); J(end+1) = J(end);
if exist('LungeDN','var')
    L = LungeDN;
    LI = LungeI;
    if exist('LungeC', 'var')
        LC = LungeC;
    else
        LC = nan(length(LI));
    end
elseif exist('time','var');
    L = time;
    for ii = 1:length(L)
        [~,LI(ii)] = min(abs(DN-L(ii)));
    end
    if size(LI,2)>1;
        LI = LI';
    end
    if exist('LungeC', 'var')
        LC = LungeC;
    else
        LC = nan(length(LI));
    end
else
    L = nan(0,0);
    LI = nan(0,0);
    LC = nan(0,0);
end

% Convert Loaded Variables to temp variables
if exist('StrategyS','var')
    StI = StrategyS; StDN = DN(StrategyS);end
if exist('StrategyE','var')
    StEI = StrategyE; StEDN = DN(StrategyE);end
if exist('StrategyC','var')
    StC = StrategyC; StEC =  StrategyC; end
if exist('StrategySC','var')
    StC = StrategySC; end
if exist('StrategyEC','var')
    StEC = StrategyEC; end

% Create Foraging Strategy Temp Variables
if ~exist('StI','var') % High Confidence Variables
    StI = nan(0,0);
    StC = nan(0,0);
    StDN = nan(0,0);
    StEI = nan(0,0);
    StEC = nan(0,0);
    StEDN = nan(0,0);
end

if ~exist('StrategyI','var');
    StrategyI = [];end

disp(['TagTurnedOn: ' datestr(DN(1),'mm/dd/yy HH:MM:SS.fff')]);
disp(['TagOnAnimal: ' datestr(DN(find(tagon,1)),'mm/dd/yy HH:MM:SS.fff')]);
disp(['EndData: ' datestr(DN(end),'mm/dd/yy HH:MM:SS.fff')]);
disp(['TagOffAnimal: ' datestr(DN(find(tagon,1,'last')),'mm/dd/yy HH:MM:SS.fff')]);

disp('Section 1 finished: Data Loaded');

%% 2. Plot and Process

% Check to see if we should start at beginning or at a saved index (i)
if ~exist('progressIndex','var') || ~atLast
    i = find(tagon,1); %Sets i as the progress index and automatically sets to Tagon if no index alrady exists
else
    i = progressIndex; % If it exists, set i to progressindex
    % PrograssIndex is only the progress at the time of saving the
    % file. i is the progress throughout the script
end

for iii = 1 % A single loop where iii = 1?
    instructions = sprintf('Controls:\nLeftClick: High Confidence\nL: Likely a Lunge\nM: Maybe a Lunge\nRightClick: Delete Point\n1-9: Change Zoom(x10)\nB: Move Back\nEnter: Move Forward\nS: Save\nQ: High Confidence Bubble Net \nW: Low Confidence Bubble Net\nR: Delete Bubble Net Marker');
    while i<find(tagon,1,'last') % Main while loop that tracks progres (i). Will continue as long as i is smaller than the last tagon index
        figure(101); clf % creates figure and clears the current window
        annotation('textbox', [0, 0.5, 0, 0], 'string', instructions,'FitBoxToText','on') % adds annotation box
        
        e = min(find(p(i+M*60*fs:end)<10,1,'first')+i+(M+1)*60*fs-1,length(p));
        % progress index (i) + minute window (M - 10min) x 60 x fs (sampling rate - 10Hz)
        % p(minutes *60 -> Seconds * sampling rate -> 6000 (10minutes
        % of index) + Progress Index
        % p(ten minutes past progress index: end of pressure)
        % find(first pressure index less than 10 meters from the
        % progress index + 10 minutes
        
        % e variable is the end of the window most likely, spits out an
        % index 11 minutes after i.
        if isempty(e)||isnan(e); e = length(p); end
        I = max(i-60*fs,1):e; % Not sure why there is a max here. I is an index 1 minute before i : e
        tagonI = false(size(p)); tagonI(I) = true; % creates a tagon index using I
        tagonI = tagon&tagonI; % where both are true.
        s1 = subplot(3,1,1); %adds subplot to figure
        [ax1,~,h2] = plotyy(DN(I),p(I),DN(I),J(I)); %plots Jerk and pressure in top subplot
        set(ax1(1),'ydir','rev','nextplot','add','ylim',[-5 max(p(tagonI))]); %Sets axis details
        ylabel('Jerk','parent',ax1(2));
        ylabel('Depth','parent',ax1(1));
        set(ax1(2),'ycolor','m','ylim',[0 1.2*max(J(tagonI))]);
        set(h2,'color','m');
        set(ax1,'xlim',[DN(I(1)) DN(I(end))]);
        set(ax1,'xticklabel',datestr(get(gca,'xtick'),'mm/dd HH:MM:SS'));
        title(filename(1:end-11));
        s2 = subplot(3,1,2); % Add second subplot
        uistack(ax1(1)); %restacks axis in subplot?
        set(ax1(1), 'Color', 'none');
        set(ax1(2), 'Color', 'w')
        [ax2,h1,h2] = plotyy(DN(I),pitch(I)*180/pi,DN(I),roll(I)*180/pi); set(ax2(2),'nextplot','add','ycolor','k','ylim',[-180 180]);
        % Plot pitch, roll ^
        ylabel('Roll and Head','parent',ax2(2));
        plot(ax2(2),DN(I),head(I)*180/pi,'b.','markersize',4);
        % Plot heading ^
        set(ax2(1),'ycolor','g','nextplot','add','ylim',[-90 90]);
        ylabel('pitch','parent',ax2(1));
        set(h1,'color','g'); set(h2,'color','r','linestyle','-','markersize',4);
        set(ax2,'xlim',[DN(I(1)) DN(I(end))]);
        set(ax2,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        s3 = subplot(3,1,3); % add third subplot
        ax3 = plot(DN(I),speedJJ(I),'b',DN(I),speedFN(I),'g');
        % plot speed JJ vs FN? Jiggle versus Flownoise?
        set(s3,'nextplot','add');
        marks = nan(1,3); % creates a marks variable with 1x3 NaN
        stratmarks = nan(1,4);
        if ~isempty(L) % L is Lunge time? This just checks if any lunges have been detected and plots marks with associated confidence scores
            %change color based on confidence
            colors = 'rbk';
            for c=1:3
                II = find(LC==c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                if ~isempty(II)
                    marks(1) = plot(ax1(1),L(II),p(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                    marks(2) = plot(ax2(1),L(II),pitch(LI(II))*180/pi,[colors(c) 's'],'markerfacecolor',colors(c));
                    marks(3) = plot(s3,L(II),speedJJ(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                end
            end
        end
        
        if ~isempty(StI); % Strategy Start Marks
            startcolors = 'gc';
            for c=1:2
                JJ = find(StC==c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                if ~isempty(JJ)
                    stratmarks(1) = plot(ax1(1),StDN(JJ),p(StI(JJ)),['' 's'],'markerfacecolor',startcolors(c));
                    stratmarks(2) = plot(ax2(2),StDN(JJ),head(StI(JJ))*180/pi,['k' 's'],'markerfacecolor',startcolors(c));
                end
            end
        end
        if ~isempty(StEI); % Strategy End Marks
            endcolors = 'rm';
            for c=1:2
                KK = find(StEC==c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                if ~isempty(KK)
                    stratmarks(3) = plot(ax1(1),StEDN(KK),p(StEI(KK)),['k' 's'],'markerfacecolor',endcolors(c));
                    stratmarks(4) = plot(ax2(2),StEDN(KK),head(StEI(KK))*180/pi,['k' 's'],'markerfacecolor',endcolors(c));
                end
            end
        end
        
        % Set graphical settings to plot 3 (Speed plot)
        set(s3,'xlim',[DN(I(1)) DN(I(end))]);
        set(s3,'ylim',[0 1.1*max(speedJJ(tagonI))],'xlim',[DN(I(1)) DN(I(end))]);
        set(s3,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        ylabel('Speed');
        
        button = 1; % set button to 1
        redraw = false; % set redraw to false
        close2lunge = 15; % A value determining the lunge point adjustment
        if strcmp(whaleName(1:2),'bb'); close2lunge = 5; end % Change the adjustment to 5 if it's a minke
        
        while ~isempty(button)
            redraw = false; % set redraw to false
            if ishandle(ax2) % if ax2 is working? begin the gninput command
                [x,~,button] = ginput(1); % graphical selection tool command that returns the x value and the button (that you pressed) and its value
                if isempty(button); continue; end % if button is empty, continue
                switch button % Selectional script base on the button input
                    case 1
                        [~,xI] = min(abs(DN-x));
                        [~,mI] = max(speedJJ(xI-5*fs:xI+5*fs)); %find the max within 5 seconds
                        mI = mI +xI-5*fs - 1;
                        if any(abs(LI-xI)<close2lunge*fs); %if it's close, change one that exists
                            [~,delI] = min(abs(L-x));
                            L(delI) = []; LI(delI) = []; LC(delI) = [];
                            LC = [LC;3];
                            [L,II] = sort([L; x]);
                            LI = sort([LI;xI]);
                            LC = LC(II);
                        else
                            LC = [LC; 3];
                            [L,II] = sort([L;DN(mI)]); LI = sort([LI;mI]); LC = LC(II);
                        end
                    case 108 %l selected - likely lunge
                        [~,xI] = min(abs(DN-x));
                        [~,mI] = max(speedJJ(xI-5*fs:xI+5*fs)); %find the max within 5 seconds
                        mI = mI +xI-5*fs - 1;
                        if any(abs(LI-xI)<close2lunge*fs); %if it's close, change one that exists
                            [~,delI] = min(abs(L-x));
                            L(delI) = []; LI(delI) = []; LC(delI) = [];
                            LC = [LC;2];
                            [L,II] = sort([L; x]);
                            LI = sort([LI;xI]);
                            LC = LC(II);
                        else
                            LC = [LC; 2];
                            [L,II] = sort([L;DN(mI)]); LI = sort([LI;mI]); LC = LC(II);
                        end
                    case 109 %m selected - maybe lunge
                        [~,xI] = min(abs(DN-x));
                        [~,mI] = max(speedJJ(xI-5*fs:xI+5*fs)); %find the max within 5 seconds
                        mI = mI +xI-5*fs - 1;
                        if any(abs(LI-xI)<close2lunge*fs); %if it's close, change one that exists
                            [~,delI] = min(abs(L-x));
                            L(delI) = []; LI(delI) = []; LC(delI) = [];
                            LC = [LC;1];
                            [L,II] = sort([L; x]);
                            LI = sort([LI;xI]);
                            LC = LC(II);
                        else
                            LC = [LC;1];
                            [L,II] = sort([L;DN(mI)]); LI = sort([LI;mI]); LC = LC(II);
                        end
                        
                    case 113 %q to mark High Confidence Bubble netting
                        [~, StrategyI] =  min(abs(DN-x)); % This takes the index of x and saves it to Strategy I
                        StrategyItemp = StrategyI; % For Quality Control on end point selection
                        StDN = sort([StDN;DN(StrategyI)]);
                        [StI,JJ] = sort([StI;StrategyI]);
                        StC = [StC;1]; % Confidence Level 1 (HC)
                        StC = StC(JJ); % Match sort to StDn and StI
                        
                        % Graph Start Point before selecting end point
                        if ~isempty(StI) % Strategy Start Marks
                            startcolors = 'gc';
                            for c = 1:2
                                JJ = find(StC == c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                                if ~isempty(JJ)
                                    stratmarks(1) = plot(ax1(1),StDN(JJ),p(StI(JJ)),['k' 's'],'markerfacecolor',startcolors(c));
                                    stratmarks(2) = plot(ax2(2),StDN(JJ),head(StI(JJ))*180/pi,['k' 's'],'markerfacecolor',startcolors(c));
                                end
                            end
                        end
                        % Select End Point with Quality Controls
                        qc = false;
                        while qc == false
                            [x,~,button] = ginput(1); % New selection
                            [~, StrategyI] =  min(abs(DN-x)); % This takes the index of x and saves it to Strategy I
                            if StrategyItemp >= StrategyI
                                warning('End Point Cannot be before or on Start Point, choose your end point again');
                                continue
                            else
                                qc = true;
                            end
                        end
                        
                        [StEDN,KK] = sort([StEDN;DN(StrategyI)]);
                        StEI = sort([StEI;StrategyI]);
                        StEC = [StEC;1];
                        StEC = StEC(KK);
                        
                    case 119 %w to mark Low Confidence Bubble Netting
                        [~, StrategyI] =  min(abs(DN-x)); % This takes the index of x and saves it to Strategy I
                        StrategyItemp = StrategyI; % For Quality Control on end point selection
                        StDN = sort([StDN;DN(StrategyI)]);
                        [StI,JJ] = sort([StI;StrategyI]);
                        StC = [StC;2]; % Confidence Score 2
                        StC = StC(JJ); % Sort StC to match StI and StDN
                        
                        % Graph Start Point before selecting End Point
                        if ~isempty(StI) % Strategy Start Marks
                            startcolors = 'gc';
                            for c=1:2
                                JJ = find(StC==c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                                if ~isempty(JJ)
                                    stratmarks(1) = plot(ax1(1),StDN(JJ),p(StI(JJ)),['k' 's'],'markerfacecolor',startcolors(c));
                                    stratmarks(2) = plot(ax2(2),StDN(JJ),head(StI(JJ))*180/pi,['k' 's'],'markerfacecolor',startcolors(c));
                                end
                            end
                        end
                        % Select End Point with Quality Controls
                        qc = false;
                        while qc == false
                            [x,~,button] = ginput(1); % New selection
                            [~, StrategyI] =  min(abs(DN-x)); % This takes the index of x and saves it to Strategy I
                            if StrategyItemp >= StrategyI
                                warning('End Point Cannot be before or on Start Point, choose your end point again');
                                continue
                            else
                                qc = true;
                            end
                        end
                        [StEDN, KK] = sort([StEDN;DN(StrategyI)]);
                        StEI = sort([StEI;StrategyI]);
                        StEC = [StEC;2];
                        StEC = StEC(KK);
                        
                    case 3 % delete an x
                        [~,delI] = min(abs(L-x)); % edited to choose the x where L is closest
                        L(delI) = []; LI(delI) = []; LC(delI) = [];
                        redraw = true; button = [];
                        
                    case 114 %r delete Strategy marker
                        [SS,delI] = min(abs(StDN-x));% Find Closest Start Marker
                        [EE,delEI] = min(abs(StEDN-x)); % Find Closest End Marker
                        
                        if length(StI) >0 & length(StEI) >0 % Check if points exist
                            if SS<EE
                                StI(delI) = []; StDN(delI) = []; StC(delI) = [];
                                StEI(delI) = []; StEDN(delI) = []; StEC(delI) = [];
                            else if EE<SS
                                    StI(delEI) = []; StDN(delEI) = []; StC(delEI) = [];
                                    StEI(delEI) = []; StEDN(delEI) = []; StEC(delEI) = [];
                                end
                            end
                        else if length(StI) == 0 || length(StEI) == 0 % Added this to be able to remove marker if only 1 marker exists, clears all variables if so.
                                StI = []; StDN = []; StEI = []; StEDN = []; StC = []; StEC = [];
                            end
                        end
                        
                        SS = []; EE = []; % Clear min values
                        redraw = true; button = [];
                        
                    case 98 %if b, go backwards
                        i = max(find(tagon,1),i-M*60*fs);
                        redraw = true; button = [];
                    case num2cell(49:57) %if you press a number, change drawing to 10*that number
                        M = 10*(button-48);
                        redraw = true; button = [];
                    case 115 %s selected - save progress
                        savedetection(DN,fs,L,LI,p,LC,i,creator,primary_cue,notes,fileloc,filename,StI,StEI,StC,StEC);
                end
                
                if ~isempty(L) % This is the main plotter for the selection points
                    try delete(marks); catch; end
                    %change color base on confidence
                    colors = 'rbk'; % define colors
                    for c=1:3
                        II = find(LC==c); % find where LC equals confidence 1:3
                        if ~isempty(II)
                            marks(1) = plot(ax1(1),L(II),p(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                            marks(2) = plot(ax2(1),L(II),pitch(LI(II))*180/pi,[colors(c) 's'],'markerfacecolor',colors(c));
                            marks(3) = plot(s3,L(II),speedJJ(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                        end
                    end
                end
                
                if ~isempty(StI); % Strategy Start Marks
                    try delete(stratmarks(1:2)); catch; end
                    startcolors = 'gc';
                    for c=1:2
                        JJ = find(StC==c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                        if ~isempty(JJ)
                            stratmarks(1) = plot(ax1(1),StDN(JJ),p(StI(JJ)),['k' 's'],'markerfacecolor',startcolors(c));
                            stratmarks(2) = plot(ax2(2),StDN(JJ),head(StI(JJ))*180/pi,['k' 's'],'markerfacecolor',startcolors(c));
                        end
                    end
                end
                
                if ~isempty(StEI); % Strategy End Marks
                    try delete(stratmarks(3:4)); catch; end
                    endcolors = 'rm';
                    for c=1:2
                        KK = find(StEC==c); % Looks for confidence scores and plots the Lunge Index with the confidence score colors
                        if ~isempty(KK)
                            stratmarks(3) = plot(ax1(1),StEDN(KK),p(StEI(KK)),['k' 's'],'markerfacecolor',endcolors(c));
                            stratmarks(4) = plot(ax2(2),StEDN(KK),head(StEI(KK))*180/pi,['k' 's'],'markerfacecolor',endcolors(c));
                        end
                    end
                end
            end
        end  % while the button value is not empty
        
        if redraw
            continue;
        else
            i = e;
        end
    end % End of First While Loop
    
    
end

% Save Strategy and Lunge File
savedetection(DN,fs,L,LI,p,LC,i,creator,primary_cue,notes,fileloc,filename,StI,StEI,StC,StEC);


