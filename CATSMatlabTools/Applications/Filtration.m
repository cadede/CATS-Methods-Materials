 % %%% Performance on a Dive Scale %%%
% %%% For this we need:
% %%% 1)individual
% %%% 2) total number of lunges
% 3a) hours of day w/tag on
% 3b) number of lunges during day
% 4a) hours of night w/ tag on
% 4b) number of lunges during night
% 5a) hours of dawn/dusk w/ tag on
% 5b) number of lunges during dawn/dusk
% 6)  Filter time for each lunge
% 6b) mean, sd filter time for each whale
%stroke_glide, signal processing, runcirc_mean, Sun_Angle needed for work 


%% How "click" through the plotted figure generated in Section 4
% The lunge will already be identified and appear as a black square
% Your first click should be at mouth opening:the peak in speed that occurs
% during the lunge (the lunge marker is usually very close to this)
% Your second click should be at mouth closing: the sharpest point of
% deceleration following the lunge 
% Your third click should encompass the 'purge time': you click just before
% the gyroscope picks up again, indicating fluke movement
% Then you use keyboard key 'g' or 'b' to recognized this lunge as a 'good'
% or 'bad' lunge. Important when multiple people are generating lunge
% detection files. 
% To move onto the next lunge, you hit SPACEBAR



%% 1) Select an individual whale, open prh
clear;

cf = pwd; 
[filename,fileloc]=uigetfile('*prh.mat', 'the PRH file to analyze'); %look below to save time to make truncated files okay
cd(fileloc);
disp('Loading Data, will take some time'); 
load([fileloc filename(1:end-3) 'mat']); 
    ii = strfind(filename,' ');
    lungename = [filename(1:ii-1) 'lunges.mat'];
    load([fileloc filename]);
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
    disp(['TagTurnedOn: ' datestr(DN(1),'mm/dd/yy HH:MM:SS.fff')]);
    disp(['TagOnAnimal: ' datestr(DN(find(tagon,1)),'mm/dd/yy HH:MM:SS.fff')]);
    disp(['EndData: ' datestr(DN(end),'mm/dd/yy HH:MM:SS.fff')]);
    disp(['TagOffAnimal: ' datestr(DN(find(tagon,1,'last')),'mm/dd/yy HH:MM:SS.fff')]);

    
disp('Section 1 finished: Data Loaded');
%% 2) Total Number of Lunges (if lunges already found), as well as other elements
% Dependent on lunges.mat file

%Upload lunges.mat file
[lungefile, fileloc] = uigetfile('*lunges.mat','multiselect','on');
load([fileloc lungefile]);

%Prepare some variables of interest 
numberoflunges = length(LungeI);
depth = p(LungeI); length(depth);
timeoflunge = datestr(DN(LungeI), 'HH:MM:SS');
avgdepth = mean(depth);
sdfordepth = std(depth);
Max = max(depth);
Min = min(depth);
t1 = find(tagon,1);
t2 = find(tagon, 1, 'last');

% Calculcate Solar Elevations for Lunges using first GPS point and Datetime
% for each lunge
[dielLunge el]  = DielLunge_v2(LungeDN,GPS,INFO);


%Two tables are made: 
% 1) A metadata table. One line per individual that you upload
% 2) A lunge table. Each lunge is it's own line. The index, the depth, the
% time, and the time of day are recorded for each lunge.

FiltrationMetadata = table({INFO.whaleName}, t1, t2);

LungeTable = table(LungeI, timeoflunge, depth, dielLunge, el);

%Save the tables using something like: 
%xlswrite('C:\Users\Shirel\Documents\Shirel\Antarctica\FiltrationMetadata.xls', table2cell(FiltrationMetadata));
%xlswrite('C:\Users\Shirel\Documents\Goldbogen Lab\Thesis\Chapter 2- Filtration\Lunges\LungeTable.xls',table2cell(LungeTable)); 

%Actually want both saved directly to the individual they came from

disp('Section 2 finished: Lunge Table and Filtration Metadata xls made');



%% 3) Smooth speed
%

speedJJ = speed.JJ;
speedFN = speed.FN;
speedJJ(isnan(speedJJ)) = min(speedJJ); speedJJ = runmean(speedJJ,ceil(fs/2)); speedJJ(isnan(speed.JJ)) = nan; %Smooths 
speedFN(isnan(speedFN)) = min(speedFN); speedFN = runmean(speedFN,ceil(fs/2)); speedFN(isnan(speed.FN)) = nan; %Smooths
SmoothGw=runcirc_mean(Gw(:,2),4); 
P = runcirc_mean(SmoothGw,4*fs);
[glds,stks] = stroke_glide([SmoothGw-P roll head],fs,1,4*pi/180,7);
%% 4) Figure set up
figure(); clf;
J = njerk(Aw,fs); 
cd(fileloc);
purgeloc = dir('*PurgeTable.mat');
try load(purgeloc.name); catch; end
%this loads the saved version from the last empty place, does not overwrite

if exist('purgeProgess', 'var')
    i = purgeProgress
else     
    i = 1;
end 

while i<=length(LungeI)
    
    %subplot 1: Plots depth with the lunge marked on it 
    sp1 = subplot(12,1,1:3); % sSbplot divides figure into grid
    e2 = LungeI(i); % The pre-detected lunge will show up as an X
    I = e2-fs*120:e2+fs*120; % Display 2 minutes before and 2 minutes after lunge, can be changed by modifying the '120' (it's in seconds)
    plot(DN(I),p(I),'b',DN(round(stks(:,1)*fs)),p(round(stks(:,1)*fs)),'rx'); hold on; % Adds little red Xs to represent the fluke strokes 
    for ii = 1:length(glds(:,1)) %plots green line for the glide long the depth axis
        II = round(glds(ii,1)*fs):round(fs*glds(ii,2));
        plot(DN(II),p(II),'g');
    end
    set(gca,'xlim',DN([I(1) I(end)])); % Limits plot view to the instruction given in the I line (line 185)
    set(gca,'ydir','rev','xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    
    %subplot 2: Plots speed and jerk, with two different Y axes
    sp2 = subplot(12,1,4:6);
    [axs,hs1, hs2] = plotyy(DN(I),speed.SP(I),DN(I),J(I)); % speed (hs1), jerk (hs2)
    set(axs,'nextplot','add');
    set(hs1,'color','k');
    plot(axs(1),DN(I),speedJJ(I),'g');
    plot(axs(1), DN(e2), speedJJ(e2),'rs', 'markerfacecolor', 'r');
    set(axs(1),'ylim',[0 5.6]);
    set(hs2,'color','m');
    set(axs,'xlim',DN([I(1) I(end)]));
    set(axs,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    

  %subplot 3 % Plots the gyroscope
    sp3 = subplot(12,1,7:9);
%     plot(DN(I),pitch(I)-P(I));
    plot(DN(I),Gw(I,2),DN(I),P(I),'k');
    %set(axx,'xlim',DN([I(1) I(end)]));
    xlim(DN([I(1) I(end)]));
   
    
   % subplot 4: Plot pitch (green), roll (red), and heading (blue)
    sp4 = subplot(12,1,10:12);
    [ax,h1,h2] = plotyy(DN(I),pitch(I)*180/pi,DN(I),roll(I)*180/pi);
    set(ax,'nextplot','add'); plot(ax(2),DN(I),head(I)*180/pi,'b.');
    set(h1,'color','g','linestyle','-');
    set(h2,'color','r','linestyle','-');
    set(ax,'xlim',DN([I(1) I(end)]));
    set(ax,'xticklabel',datestr(get(gca,'xtick'),'HH:MM:SS'));
    plot(DN(I),P(I)*180/pi,'k');
    %linkaxes([sp1 sp2 sp3 sp4],'x'); 
   
    % The section below determines how many times you can click and what
    % keyboard key you click to represent the filtration variables 
    % The filtration variables are purge2, purge2, and whether the lunge is
    % deep (keyboard key d) or shallow (keyboard key s). In theory, these could be changed to whatever
    % is of interest to the user. 
    % In your table, the lungecat column will have 100, representing d, and 115 representing s 

    clickthislunge = true;
    tic;
    while clickthislunge;
        [x,y,button] = ginput(4); %the four clicks are opening, closing, purge1, and bad or good if you disagree with the lunge detection 
        xi = nan(3, 1);
        for ii = 1:3 %only need 3 bc d/s is not based on speed or timing
            [~,xi(ii)] = min(abs(DN-x(ii))); %fills in xi with values that you've clicked
        end
        looking4purge= plot(axs(1),DN(xi),speedJJ(xi),'ks','markerfacecolor','m');
        if button(end) == 103; % 103 is g, for good. is this a good lunge? 
            lungecat = 103;
            D=text(DN(e2+25),speedJJ(e2), 'G', 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'parent', sp2);
        elseif button(end) == 098;
            lungecat = 098;
            D=text(DN(e2+25),speedJJ(e2),'B', 'FontSize', 12, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'parent', sp2);
            xi = nan(3, 1); % makes xi filled with NaNs. Used for when you disagree with the lunge, b is bad
        else lungecat = nan;
            D=text(DN(e2+25),speedJJ(e2),'Click F. You Messed Up.'); %allows you to re start that lunge 
        end

        [~,~,button] = ginput(1);
        if button == 102;
            clickthislunge = true;
            delete (looking4purge);
            delete (D); % SPACEBAR MOVES YOU FORWARD
        elseif button == 032;
            clickthislunge = false;
        else clickthislunge = true;
            text(DN(e2+25),speedJJ(e2),'You''ve really done it now. Scream in frustration');
            delete (looking4purge);
            delete (D); %reload your purgetable 
        end
    end
    LungeTable.MO(i, 1) = xi(1); %index of mouth open
    LungeTable.MC(i, 1) = xi(2); % index of mouth closed 
    LungeTable.purge1(i, 1) = (xi(3)-xi(2))/fs; %calculating purge1 which is bsed on mouth close to beginning of gyro signal
    LungeTable.lungecat(i, 1) = lungecat; % adds the lunge category to the table  
    purgeProgress = i;
    save([fileloc getWhaleID(filename) 'PurgeTable.mat'], 'LungeTable', 'purgeProgress')
    i = i+1;
end 