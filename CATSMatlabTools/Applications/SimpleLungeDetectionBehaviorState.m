% % Identify simple lunge cues (just at mouth opening)
% David Cade, James Fahlbusch, and Jake Linsky
% version 2.0
% Goldbogen Lab
% Stanford University
% 
% Input - 10hz PRH mat file
% Outputs - Section 2: Lunge .mat file, behavior state .mat file
%           Section 3: Error check variable and command window report for behavior states.
%           Returns indices and description of suspected errors.
%           Section 4: Behavior state .xlsx file (same format as .mat file)
%           Section 5: Behavior state .xlsx file with durations (if all
%           behaviors have start and finish)
%           
% 
% 
%% 1. Load Data
clear; % clears the workspace
%Start where you left off?
atLast = true; %this will look for a variable called progressIndex
M = 10; % number of minutes
% Variables that will be saved in the Lunge file
notes = '';
creator = 'DEC';
primary_cue = 'speedJJ';

%Set behaviroal state names
states = cellstr(['state1';'state2';'state3';'state4']);

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
    clear LungeI LungeDN LungeC LungeDepth % clears lunge data in case a previous file was still loaded
    behaviorname = [filename(1:ii-1) 'BehaviorState.mat'];
    lungename = [filename(1:ii-1) 'lunges.mat'];
    load([fileloc filename]);
    try load([fileloc behaviorname]); catch; end
    try load([fileloc lungename]); catch; end
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
%  Behavior variables
      BehaviorI = [];
      BehaviorState = [];
      BehaviorText = [];
      Bcolors = ['k'];
   if exist ('BehaviorIndex','var')
       behI = BehaviorIndex;
       behT = BehaviorTime;
       behS = Behavior;
       behSS = StartStop;
   else
       behI = [];
       behT = [];
       behS = [];
       behSS = [];
   end
       

    disp(['TagTurnedOn: ' datestr(DN(1),'mm/dd/yy HH:MM:SS.fff')]);
    disp(['TagOnAnimal: ' datestr(DN(find(tagon,1)),'mm/dd/yy HH:MM:SS.fff')]);
    disp(['EndData: ' datestr(DN(end),'mm/dd/yy HH:MM:SS.fff')]);
    disp(['TagOffAnimal: ' datestr(DN(find(tagon,1,'last')),'mm/dd/yy HH:MM:SS.fff')]);

disp('Section 1 finished: Data Loaded');
    
%% 2. Plot and Process

%Lunge data denoted in all plots as squares.
%Behavior data plotted over dive profile
%O indicates beginning of behavioral state, X indicates termination of
%behavioral state
%If current state (denoted at bottom of window) is active, beginning a new
%state will terminate the previous state. "G" can be used to toggle states
%if this is not desired.
%Y will terminate the current state without placing a new state. For later
%export sections this should be used to finish the PRH file. 
%"H" will set the progressIndex to whatever either the lunge or
%behaviorstate that is furthest along in the file, i.e jump to point of
%furthest progress.

    % Check to see if we should start at beginning or at a saved index (i)
    if ~exist('progressIndex','var') || ~atLast
        i = find(tagon,1);
    else
        i = progressIndex;
    end
    %
    for iii = 1
    instructions = sprintf('Controls:\nLeftClick: High Confidence\nL: Likely a Lunge\nM: Maybe a Lunge\nRightClick: Delete Lunge\nQ: State1\nW: State2\nE: State3\nR: State4\nT: Delete State\nY: End State\nG:Toggle State\n1-9: Change Zoom(x10)\nB: Move Back\nEnter: Move Forward\nH: Go To Last\nS: Save');
    while i<find(tagon,1,'last')
        figure(101); clf
        annotation('textbox', [0, 0.5, 0, 0], 'string', instructions,'FitBoxToText','on')
        
        if ~isempty(BehaviorState)
        annotation('textbox',[.5, 0.1, 0, 0],'string',BehaviorText,'FitBoxToText','on','tag','ba','color',Bcolors(BehaviorState));
        end

        e = min(find(p(i+M*60*fs:end)<10,1,'first')+i+(M+1)*60*fs-1,length(p));
        if isempty(e)||isnan(e); e = length(p); end
        I = max(i-60*fs,1):e;
        tagonI = false(size(p)); tagonI(I) = true;
        tagonI = tagon&tagonI;
        s1 = subplot(3,1,1);
        [ax1,~,h2] = plotyy(DN(I),p(I),DN(I),J(I));
        set(ax1(1),'ydir','rev','nextplot','add','ylim',[-5 max(p(tagonI))]);
        ylabel('Jerk','parent',ax1(2));
        ylabel('Depth','parent',ax1(1));
        set(ax1(2),'ycolor','m','ylim',[0 1.2*max(J(tagonI))]);
        set(h2,'color','m');
        set(ax1,'xlim',[DN(I(1)) DN(I(end))]);
        set(ax1,'xticklabel',datestr(get(gca,'xtick'),'mm/dd HH:MM:SS'));
        title(filename(1:end-11));
        s2 = subplot(3,1,2);
        uistack(ax1(1));
        set(ax1(1), 'Color', 'none');
        set(ax1(2), 'Color', 'w')
        [ax2,h1,h2] = plotyy(DN(I),pitch(I)*180/pi,DN(I),roll(I)*180/pi); set(ax2(2),'nextplot','add','ycolor','k','ylim',[-180 180]);
        ylabel('Roll and Head','parent',ax2(2));
        plot(ax2(2),DN(I),head(I)*180/pi,'b.','markersize',4);
        set(ax2(1),'ycolor','g','nextplot','add','ylim',[-90 90]);
        ylabel('pitch','parent',ax2(1));
        set(h1,'color','g'); set(h2,'color','r','linestyle','-','markersize',4);
        set(ax2,'xlim',[DN(I(1)) DN(I(end))]);
        set(ax2,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        s3 = subplot(3,1,3);
        ax3 = plot(DN(I),speedJJ(I),'b',DN(I),speedFN(I),'g');
        set(s3,'nextplot','add');
        marks = nan(1,3);
        if ~isempty(L)
            %change color based on confidence
            colors = 'rbk';
            for c=1:3
                II = find(LC==c);
                if ~isempty(II)
                    marks(1) = plot(ax1(1),L(II),p(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                    marks(2) = plot(ax2(1),L(II),pitch(LI(II))*180/pi,[colors(c) 's'],'markerfacecolor',colors(c));
                    marks(3) = plot(s3,L(II),speedJJ(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                end
            end 
        end

%plot colors and o start x end for behavior states
        if ~isempty(behI)
        redraw = true;
                 try delete(Bmarks); catch; end;
                  Bcolors = 'rgkc';
                  Bshapes = 'ox';
                  for ss = 1:2;
                          SS = find (behSS==ss);
            %change cc end value to # of behavioral states used
                     for cc = 1:4
                      QQ = find (behS==cc);
                          CS = intersect(QQ,SS);
                              if ~isempty(CS)
                                   Bmarks = plot(ax1 (1),behT(CS),p(behI(CS)),[Bcolors(cc) Bshapes(ss)])
                              end
                     end
                  end
            end

        set(s3,'xlim',[DN(I(1)) DN(I(end))]);
        set(s3,'ylim',[0 1.1*max(speedJJ(tagonI))],'xlim',[DN(I(1)) DN(I(end))]);
        set(s3,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
        ylabel('Speed');

        button = 1;
        redraw = false;
        close2lunge = 15;
        if strcmp(whaleName(1:2),'bb'); close2lunge = 5; end
        
        while ~isempty(button)
            redraw = false;
            if ishandle(ax2)
                [x,~,button] = ginput(1);
                if isempty(button); continue; end
                switch button
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

                    case 103 %g toggle BehaviorState
                        if isempty (BehaviorState)
                            BehaviorState = 1;
                        else if BehaviorState < 4
                                BehaviorState = BehaviorState + 1;
                            else
                                BehaviorState = [];
                            end
                        end
                            
                    case 113 %q mark state1
                        [~, BehaviorI] = min(abs(DN-x));
                        if ~isempty (BehaviorState)
                            behS = [behS;BehaviorState];
                            behSS = [behSS;2];
                            [behT,II] = sort([behT;DN(BehaviorI-1)]); behI = sort([behI;BehaviorI-1]); behS = behS(II); behSS = behSS(II);              
                        end
                        BehaviorState = 1;
                        behS = [behS;BehaviorState];
                            behSS = [behSS;1];
                            [behT,II] = sort([behT;DN(BehaviorI)]); behI = sort([behI;BehaviorI]); behS = behS(II); behSS = behSS(II); 
                             
                    case 119 %w mark state2
                        [~, BehaviorI] = min(abs(DN-x));
                         if ~isempty (BehaviorState)
                      behS = [behS;BehaviorState];
                            behSS = [behSS;2];
                            [behT,II] = sort([behT;DN(BehaviorI-1)]); behI = sort([behI;BehaviorI-1]); behS = behS(II); behSS = behSS(II);              
                        end
                        BehaviorState = 2;
                        behS = [behS;BehaviorState];
                            behSS = [behSS;1];
                            [behT,II] = sort([behT;DN(BehaviorI)]); behI = sort([behI;BehaviorI]); behS = behS(II); behSS = behSS(II); 
                             
                    case 101 %e mark state3
                        [~, BehaviorI] = min(abs(DN-x))
                         if ~isempty (BehaviorState)
                            behS = [behS;BehaviorState];
                            behSS = [behSS;2];
                            [behT,II] = sort([behT;DN(BehaviorI-1)]); behI = sort([behI;BehaviorI-1]); behS = behS(II); behSS = behSS(II);              
                        end
                        BehaviorState = 3;
                        behS = [behS;BehaviorState];
                            behSS = [behSS;1];
                            [behT,II] = sort([behT;DN(BehaviorI)]); behI = sort([behI;BehaviorI]); behS = behS(II); behSS = behSS(II); 
                             
                        
                    case 114 %r mark state 4
                        [~, BehaviorI] = min(abs(DN-x))
                         if ~isempty (BehaviorState)
                            behS = [behS;BehaviorState];
                            behSS = [behSS;2];
                            [behT,II] = sort([behT;DN(BehaviorI-1)]); behI = sort([behI;BehaviorI-1]); behS = behS(II); behSS = behSS(II);              
                        end
                        BehaviorState = 4;
                        behS = [behS;BehaviorState];
                            behSS = [behSS;1];
                            [behT,II] = sort([behT;DN(BehaviorI)]); behI = sort([behI;BehaviorI]); behS = behS(II); behSS = behSS(II); 
                             
                        
                    case 116 %t delete state
                        [~, delB] = min(abs(DN-x))
                        [~, btI] = min(abs(behI-delB))
                        behI(btI) = [];
                        behT(btI) = [];
                        behS(btI) = [];
                        behSS(btI) = [];
                        %clear delB btI
                        redraw = true; button = [];
                        
                    case 121 %y place end state marker (for final open state or corrections)
                        [~, BehaviorI] =  min(abs(DN-x));
                        if ~isempty (BehaviorState)
                            behS = [behS;BehaviorState];
                            behSS = [behSS;2];
                            [behT,II] = sort([behT;DN(BehaviorI)]); behI = sort([behI;BehaviorI]); behS = behS(II); behSS = behSS(II);              
                        end
                        
                    case 104 % set progressIndex to last point
                        i = max([behI;LI]);
                        redraw = true; button = [];

                    case 3 % delete a lunge
                        [~,delI] = min(abs(L-x));
                        L(delI) = []; LI(delI) = []; LC(delI) = [];
                          redraw = true; button = [];           

                    case 98 %if b, go backwards
                        i = max(find(tagon,1),i-M*60*fs);
                        redraw = true; button = [];
                    case num2cell(49:57) %if you press a number, change drawing to 10*that number
                        M = 10*(button-48);
                        redraw = true; button = [];
                    case 115 %s selected - save progress
                        % set the created on date vector
                        d1 = datevec(now());
                        created_on = [d1(2) d1(3) d1(1)];
                        clearvars d1;
                        % store temp variables to lunge file 
                        starttime = DN(1);
                        prh_fs = fs;
                        LungeDN = L;
                        depth = p(LI);
                        time = L;
                        LungeI = LI;
                        LungeC = LC;
                        LungeDepth = depth;
                        progressIndex = i;
                        save([fileloc getWhaleID(filename) 'lunges.mat'],'LungeDN','LungeI','LungeDepth','LungeC','creator','primary_cue','prh_fs','starttime','created_on', 'progressIndex', 'notes');

                        %store behavior variables to BehaviorState file
                        BehaviorIndex = behI;
                        BehaviorTime = behT;
                        Behavior = behS;
                        StartStop = behSS;
                        save([fileloc getWhaleID(filename) 'BehaviorState.mat'],'starttime','prh_fs','progressIndex','creator','created_on','notes','BehaviorIndex','BehaviorTime','Behavior','StartStop');
                     
                end
                if ~isempty(L)
                    try delete(marks); catch; end
                    %change color base on confidence
                    colors = 'rbk';
                    for c=1:3
                        II = find(LC==c);
                        if ~isempty(II)
                            marks(1) = plot(ax1(1),L(II),p(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                            marks(2) = plot(ax2(1),L(II),pitch(LI(II))*180/pi,[colors(c) 's'],'markerfacecolor',colors(c));
                            marks(3) = plot(s3,L(II),speedJJ(LI(II)),[colors(c) 's'],'markerfacecolor',colors(c));
                        end
                    end 
                end

%plot behavior points
            if ~isempty(behI)
                 try delete(Bmarks); catch; end;
                  Bcolors = 'rgkc';
                  Bshapes = 'ox';
                  for ss = 1:2;
                          SS = find (behSS==ss);
            %change cc end value to # of behavioral states used
                     for cc = 1:4
                      QQ = find (behS==cc);
                          CS = intersect(QQ,SS);
                              if ~isempty(CS)
                                   Bmarks = plot(ax1 (1),behT(CS),p(behI(CS)),[Bcolors(cc) Bshapes(ss)])
                              end
                     end
                  end
            end
            %create textbox for auditor to keep track of state
            BehaviorText = states(BehaviorState);
            if exist('ba')==1; delete(findall(gcf,'tag','ba')); end;
            if ~isempty(BehaviorState);
            ba = annotation('textbox',[.5, 0.1, 0, 0],'string',BehaviorText,'FitBoxToText','on','tag','ba','color',Bcolors(BehaviorState));
            end

            end
        end
        if redraw
            continue;
        else
            i = e;
        end
    end
    % set the created on date vector
    d1 = datevec(now());
    created_on = [d1(2) d1(3) d1(1)];
    clearvars d1;
    % store temp variables to lunge file 
    starttime = DN(1);
    prh_fs = fs;
    LungeDN = L;
    depth = p(LI);
    time = L;
    LungeI = LI;
    LungeC = LC;
    LungeDepth = depth;
    progressIndex = i;
    save([fileloc getWhaleID(filename) 'lunges.mat'],'LungeDN','LungeI','LungeDepth','LungeC','creator','primary_cue','prh_fs','starttime','created_on', 'progressIndex', 'notes');

                        %store behavior variables to BehaviorState file
                        BehaviorIndex = behI;
                        BehaviorTime = behT;
                        Behavior = behS;
                        StartStop = behSS;
                        save([fileloc getWhaleID(filename) 'BehaviorState.mat'],'starttime','prh_fs','progressIndex','creator','created_on','notes','BehaviorIndex','BehaviorTime','Behavior','StartStop','states');

    end
    
%% 3. Clean up variables (behavior audit only)
% Check for errors in cell. Returns CheckResult variable with indices of
% suspected errors along with a description. Use progressIndex = I where
% I is an idex returned in column one of the error check and run section 2
% again to view in plot. 
clear CheckResult;
for i = 1:length(behI)
    if behSS(i) == 1
        if i<length(behI) && behSS(i+1) == 1
            if ~exist('CheckResult','var')
           CheckResult = {behI(i),'previous behavior ongoing'};
            else
                CheckResult(end+1,1) = {behI(i)}; CheckResult(end,2)= {'previous behavior ongoing'}
            end
    elseif i == length(behI)
             if ~exist('CheckResult','var')
                CheckResult = {behI(i),'ongoing final behavior'};
            else
                CheckResult(end+1,1) = {behI(i)}; CheckResult(end,2)= {'ongoing final behavior'}
             end
        end
    end

    if behSS(i) == 2
        if i<length(behI) behSS(i+1) == 2
             if ~exist('CheckResult','var')
           CheckResult = {behI(i),'double stop'};
            else
                CheckResult(end+1,1) = {behI(i)}; CheckResult(end,2)= {'double stop'}
             end
        end
        if i<length(behI) && behSS(i-1) == 1 && behS(i-1)~=behS(i)
             if ~exist('CheckResult','var')
           CheckResult = {behI(i),'end behavior does not match open behavior'};
            else
                CheckResult(end+1,1) = {behI(i)}; CheckResult(end,2)= {'end behavior does not match open behavior'}
             end
        end
    end
        
end

%% 4. Export behavior state excel file
% Exports excel file in same format as .mat file
SStext ={'start','stop'};
for i = 1:length(behS)
    sText(i) = states(behS(i));
    ssText(i)= SStext(behSS(i));
end
    sText = sText(:); ssText = ssText(:);
    depName(1:length(behS)) = {convertCharsToStrings(whaleName)}; depName = depName(:);
    Headings = {'Deployment','Time','Index','Behavior','Start/Stop'};
    Table = table(depName,datestr(BehaviorTime),BehaviorIndex,sText,ssText,'VariableNames',Headings);
    writetable(Table,[whaleName 'BehaviorStates.xlsx']);  
%% 5. Calculate durations and export excel file
% exports excel file with durations for each state
clear durTable newTable
Headings = {'Deployment','State','Start time','End time','Duration'};
for i = 1:length(states)
    clear newTable depName statename
      I  = find(behS == i);
      Istart = find(behSS(I) == 1);
      Istop = find(behSS(I) == 2);
      depName(1:length(Istart)) = {convertCharsToStrings(whaleName)}; depName = depName(:);
      statename(1:length(Istart)) = states(i);statename = statename(:);
      if ~exist('durTable')
         durTable = table(depName,statename,datestr(behT(Istart)),datestr(behT(Istop)),between(datetime(behT(Istart),'ConvertFrom','datenum'),datetime(behT(Istop),'ConvertFrom','datenum'),'time'),'VariableNames',Headings);
      else
          newTable = table(depName,statename,datestr(behT(Istart)),datestr(behT(Istop)),between(datetime(behT(Istart),'ConvertFrom','datenum'),datetime(behT(Istop),'ConvertFrom','datenum'),'time'),'VariableNames',Headings);
          durTable = [durTable;newTable];
      end
end
    writetable(durTable,[whaleName 'BehaviorDurations.xlsx']);  
%% for just flownoise
i = 1;
fs = 10;
   if exist('LungeDN','var')
        L = LungeDN;
        LI = LungeI;
    elseif exist('time','var');
        L = time;
        for ii = 1:length(L)
            [~,LI(ii)] = min(abs(DN-L(ii)));
        end
        if size(LI,2)>1;
            LI = LI';
        end
    else
        L = nan(0,0);
        LI = nan(0,0);
    end
while i<length(DN)
    figure(1); clf
    e = min(i+(M+1)*60*fs-1,length(DN));
    %         if isempty(e)||isnan(e); e = length(p); end
    I = max(i-60*fs,1):e;
    %         tagonI = false(size(p)); tagonI(I) = true;
    %         tagonI = tagon&tagonI;
    %         s1 = subplot(3,1,1);
    %         [ax1,~,h2] = plotyy(DN(I),p(I),DN(I),J(I));
    %         set(ax1(1),'ydir','rev','nextplot','add','ylim',[-5 max(p(tagonI))]);
    %         ylabel('Jerk','parent',ax1(2));
    %         ylabel('Depth','parent',ax1(1));
    %         set(ax1(2),'ycolor','m','ylim',[0 1.2*max(J(tagonI))]);
    %         set(h2,'color','m');
    %         set(ax1,'xlim',[DN(I(1)) DN(I(end))]);
    %         set(ax1,'xticklabel',datestr(get(gca,'xtick'),'mm/dd HH:MM:SS'));
    %         s2 = subplot(3,1,2);
    %         uistack(ax1(1));
    %         set(ax1(1), 'Color', 'none');
    %         set(ax1(2), 'Color', 'w')
    %         [ax2,h1,h2] = plotyy(DN(I),pitch(I)*180/pi,DN(I),roll(I)*180/pi); set(ax2(2),'nextplot','add','ycolor','k','ylim',[-180 180]);
    %         ylabel('Roll and Head','parent',ax2(2));
    %         plot(ax2(2),DN(I),head(I)*180/pi,'b.','markersize',4);
    %         set(ax2(1),'ycolor','g','nextplot','add','ylim',[-90 90]);
    %         ylabel('pitch','parent',ax2(1));
    %         set(h1,'color','g'); set(h2,'color','r','linestyle','.','markersize',4);
    %          set(ax2,'xlim',[DN(I(1)) DN(I(end))]);
    %         set(ax2,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    %         s3 = subplot(3,1,3);
    ax3 = plot(DN(I),flownoise(I));%,'b',DN(I),speedFN(I),'g');
    set(gca,'nextplot','add');
    marks = nan(1,3);
    if ~isempty(L)
        %             marks(1) = plot(ax1(1),L,p(LI),'rs','markerfacecolor','r');
        %             marks(2) = plot(ax2(1),L,pitch(LI)*180/pi,'rs','markerfacecolor','r');
        marks(3) = plot(gca,L,flownoise(LI),'rs','markerfacecolor','r');
    end
    set(gca,'xlim',[DN(I(1)) DN(I(end))]);
    set(gca,'ylim',[1.1*min(flownoise(I)) 0.8*max(flownoise(I))],'xlim',[DN(I(1)) DN(I(end))]);
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    ylabel('Speed');
    button = 1;
%     title(filename(1:end-12));
    redraw = false;
    while ~isempty(button)
        try [x,~,button] = ginput(1); catch; end %JAF added TryCatch
        if isempty(button); continue; end
        switch button
            case 1
                [~,xI] = min(abs(DN-x));
                [~,mI] = max(flownoise(xI-5*fs:xI+5*fs)); %find the max within 5 seconds
                mI = mI +xI-5*fs - 1;
                if any(abs(LI-xI)<15*fs); %if it's close, change one that exists
                    [~,delI] = min(abs(L-x));
                    L(delI) = []; LI(delI) = [];
                    L=sort([L; x]); LI = sort([LI;xI]);
                else
                    L = sort([L;DN(mI)]); LI = sort([LI;mI]);
                end
            case 3 % delete an x
                [~,delI] = min(abs(L-x));
                L(delI) = []; LI(delI) = [];
            case 98 %if b, go backwards
                i = max(find(tagon,1),i-M*60*fs);
                redraw = true; button = [];
            case num2cell(49:57) %if you press a number, change drawing to 10*that number
                M = 10*(button-48);
                redraw = true; button = [];
        end
        if ~isempty(L)
            try delete(marks); catch; end
            %                 marks(1) = plot(ax1(1),L,p(LI),'rs','markerfacecolor','r');
            %                 marks(2) = plot(ax2(1),L,pitch(LI)*180/pi,'rs','markerfacecolor','r');
            marks(3) = plot(gca,L,speedJJ(LI),'rs','markerfacecolor','r');
        end
    end
    if redraw
        continue;
    else
        i = e;
    end
end
LungeDN = L;
depth = p(LI);
time = L;
LungeI = LI;
oi2 = strfind(fileloc,'\');
socalname = [fileloc(oi2(end-1)+1:oi2(end-1)+13) '-' filename(1:2)];
% save(['B:\Dropbox\Shared\AcouDarts\' socalname 'lunges'],'depth','time'); %
% save([fileloc filename(1:end-12) 'lunges.mat'],'LungeDN','LungeI');

