function [speedJJ,speedTable,speedstats] = SpeedFromRMS3(RMS1,lab1,RMS2,lab2,fs,p,pitch,roll,DN,speedper,tagslip,tagon,binSize,filterSize,minDepth0,minPitch0,minSpeed,minTime)

global label1 label2 numRMS minPitch minDepth maxDepth maxRR% only used in the context of this function
label1 = lab1; label2 = lab2; minPitch = minPitch0; minDepth = minDepth0;
% Based on and available with Cade et al. 2018, Determining forward speed from
% accelerometer jiggle in aquatic environments, Journal of Experimental
% Biology, http://jeb.biologists.org/lookup/doi/10.1242/jeb.170449

if nargin <18 || isempty(minTime); minTime = 2/fs; end
if nargin <17 || isempty(minSpeed); minSpeed = 0; end
if nargin <16 || isempty(minPitch); minPitch = 45; end
if nargin <15 || isempty(minDepth); minDepth = 5; end
if nargin <14 || isempty(filterSize); filterSize = 1; end
if nargin <13 || isempty(binSize); binSize = 1/fs; end
if nargin <12 || isempty(tagon); tagon = true(size(p)); end
if nargin <11 || isempty(speedper); speedper = length(p); end
if nargin <10 || isempty(DN); DN = (0:length(p)-1)'/24/60/60/fs; end
if nargin < 9; help SpeedFromRMS; error('Must have 9 inputs, RMS2 can be empty'); end
if length(DN) == 1; DN = (DN:1/fs/24/60/60:DN+((length(p)-1)/fs)/24/60/60)'; end %if DN is just a starttime

warning('off','stats:nlinfit:IllConditionedJacobian');
warning('off','stats:nlinfit:ModelConstantWRTParam');
warning('off','MATLAB:legend:IgnoringExtraEntries');
warning('off','stats:nlinfit:IterationLimitExceeded');
warning('off','MATLAB:rankDeficientMatrix');

if size(RMS1,1)<size(RMS1,2); RMS1 = RMS1'; end
if isempty(RMS2); numRMS = 1; else numRMS = 2; end
% if size(RMS1,2) >3; error('This function takes a maximum of 3 axes of RMS data from which to calculate speed'); end
if numRMS == 2 && size(RMS2,2)>1; error('For comparison, RMS2 must have only 1 axis'); end

if size(RMS1,2)>1; % allows for multivariate modeling
    model = true;
    J = RMS1;
    RMS1 = 20*log10(nanmean(10.^(RMS1(:,1:end)/20),2)) ; % mean across all axes, used as a default value and then will be replaced later using the multivariate model
    func = '';
    for i = 2:size(J,2)
        func = [func '+a(' num2str(i+2) ')*X(:,' num2str(i) ')']; %creates another variable for each column
    end
    func = [func '))'];
    for i = 1:size(J,2)
        func = [func '*(1-(a(' num2str(i+2) ')<0|a(' num2str(i+2) ')>1))'];
    end
    func = [func ')*(sum(a(' num2str(3) ':' num2str(i+2) '))>.9&sum(a(' num2str(3) ':' num2str(i+2) '))<1.1)'];
    func = ['fun = @(a,X)(a(1)*exp(a(2)*(a(3)*X(:,1)' func ';'];
    eval(func);
    OJ = J;
else model = false; fun = []; J = RMS1; OJ = J;
end
if numRMS == 1; RMS2 = RMS1; end; % use the same RMS1, just don't graph it

maxRR = 15; maxRR2 = maxRR; % the thresholds2 are the newly adjusted thresholds
%error checking
if length(RMS1)~=length(p) || length(p) ~=length(pitch) || length(RMS1) ~= length(RMS2); error('Vectors are not the same length'); end
if 1-.05<binSize*fs && binSize*fs < 1+.05; binData = false; else binData = true; end
if binData && abs(round(binSize*fs)-binSize*fs) > .01; error('binSize (in seconds) * sample rate must be an integer'); end

%Ovar = Original variable, before downsampling into a bin size
OJigRMS = RMS1; Ospeedper = speedper; Otagon = tagon; ODN = DN; Op = p; 
OFNRMS = RMS2;

%depth deviation
dd = runmean(p,round(filterSize/2*fs)); %smooth depth to the same filterSize size as your RMS data
dd = [diff(dd); 0];
spitch = runcirc_mean(pitch,round(filterSize*2*fs)); %smoothed pitch
spitch = [circ_mean([spitch(1:end-1) spitch(2:end)],[],2); nan]; % average each point with the one after it (since each depth deviation is happening between two pitches)
OspeedSP = -dd./sin(spitch)*fs; OspeedSP(abs(spitch*180/pi)<minPitch) = nan; OspeedSP(p<minDepth) = nan;
OspeedJJ = nan(size(RMS1)); % speed from Jiggle
clear dd;

if binData % bin data (so that you are looking at chunks of data instead of each data point on its own)
    bin = round(fs*binSize);
    X = buffer(RMS1(:,1),bin,0,'nodelay'); X(X == 0) = nan; Y = 20*log10(nanmean(10.^(X/20))); RMS1 = Y';
    X = buffer(RMS2(:,1),bin,0,'nodelay'); X(X == 0) = nan; RMS2 = 20*log10(nanmean(10.^(X/20)))';
    X = buffer(pitch,bin,0,'nodelay'); pitch = circ_mean(X)'; spitch = pitch;
    dd = runmean(p,round(filterSize/2*fs)); %smooth depth to the same filterSize size as your RMS data
    dd = [diff(dd); 0];
    X = buffer(dd,bin,0,'nodelay'); dd = sum(X)';
    X = buffer(roll,bin,0,'nodelay'); roll= circ_mean(X)';
    X = buffer(p,bin,0,'nodelay'); p = nanmean(X)';
    X = buffer(DN,bin,0,'nodelay'); X(X==0) = nan; DN = nanmean(X)'; DNf = DN-floor(DN); DNf = round(DNf*24*60*60*fs)/24/60/60/fs; DN = floor(DN)+DNf;
    speedper = ceil(speedper/bin); %speedPerEnd = ceil(speedPerEnd/bin);
    X = buffer(tagon,bin,0,'nodelay'); tagon = logical(min(X)');
    fs = 1/binSize;
    if model
        Y = nan(length(Y),0);
        for i = 1:size(J,2)
            X = buffer(J(:,i),bin,0,'nodelay'); X(X == 0) = nan; Y = [Y 20*log10(nanmean(10.^(X/20)))'];
        end
        J = Y;
    end
end

% from here on out, everything is at the reduced, binned, sizes
speedSP = -dd./sin(spitch)*fs;
% get rid of excluded values
speedSP(abs(spitch*180/pi)<minPitch) = nan;
speedSP(p<minDepth) = nan;
speedSP(speedSP<minSpeed) = nan;

speedJJ = nan(size(RMS1)); % speed from Jiggle (assuming that is RMS1
speedFN = nan(size(RMS2)); % speed from flow noise (assuming that is RMS2)
Istot = []; Istot2 = [];

% only use sections with a number of consecutive good points
S1 = speedper;
todel = [];
for i = 1:length(speedper(:,1))
    I = round(speedper(i,1)):speedper(i,2); % period of time between tag slips
    I = I(tagon(I));
    
    goodspeedI = find(~isnan(speedSP(I)) & sign(spitch(I)) == -sign(dd(I))) +I(1)-1; %the sign part gives you some gaps but ensures the whale is heading in the direction of the dive/ascent
    [s, e] = consec(goodspeedI);
    Ie = (s(2:end)-e(1:end-1))<=3; % if the consecutive sections are missing 3 or fewer samples, don't worry about that gap
    s([false Ie]) = []; e([Ie false]) = [];
    Ie = e-s>minTime*fs; % find sections longer than minTime seconds where depth is greater than minDepth and pitch is greater than minPitch degrees
    s = s(Ie); e = e(Ie);
    if isempty(s); speedFN(I) = nan;
        Is = 1; else Is = []; end % if there aren't any places with enough data to calibrate, out of luck
    for j = 1:length(s); Is = [Is s(j):e(j)]; end
    
    IsJ = Is(~isnan(RMS1(Is,1))&~isnan(speedSP(Is)));
    IsJ2 = Is(~isnan(RMS2(Is,1))&~isnan(speedSP(Is)));
    if length(IsJ)<20; % if you have less than 20 points, don't try to calibrate, just skip this section
        todel = [todel i];
        if i >1; display(['Speed section ' num2str(i) ' does not have enough points (at least 20) to calibrate, combining with prior section']);
        else display(['Speed section ' num2str(i) ' does not have enough points (at least 20)  to calibrate, combining with subsequent section']);
        end
        
    end
    Istot = [Istot IsJ];
    Istot2 = [Istot2 IsJ2];
end

speedper(todel,:) = []; speedper = [[S1(1);speedper(2:end,1)] [speedper(2:end,1); S1(end,2)]];

oi = false(size(p)); oi(Istot) = true; Istot = oi; %convert Istot to an index
oi = false(size(p)); oi(Istot2) = true; Istot2 = oi; %conve

maxDepth = nan;
% adj thresholds for all speed periods to limit data to the best sections
% for calibration
threshes = {'minDepth2' 'maxDepth2' 'minPitch2' 'maxRR2'};
[newIstot,newIstot2,minDepth2,maxDepth2,minPitch2,maxRR2] = adjthresh(Istot,Istot2,speedper,minPitch,minDepth,maxDepth,maxRR,RMS1,RMS2,speedSP,roll,pitch,p,fs);

[speedper,todel] = checksmallsections (speedper,newIstot,S1);
for iii = 1:length(threshes); eval([threshes{iii} '(todel) = [];']); end

restart = true;
while restart
    Ospeedper = round(speedper*bin);
    fJ = cell(size(speedper,1),numRMS); gofJ = fJ;
    c = cell(size(speedper,1),1);
    m = c;
    
    I = [newIstot&tagon newIstot2&tagon];
    if model
        [JigRMSA,OJigRMSA,speedJJ,OspeedJJ,mA] = makemultimodel(I,1:length(tagon),fun,J,RMS1,nan(length(OJ),1),speedSP,speedJJ,OspeedJJ,OJ,1:length(OJ));
    else
        JigRMSA = RMS1; OJigRMSA = OJigRMS;
    end
    
    fitopts = fitoptions('exp1'); fitopts.Robust = 'on';
    
    % make plots using the section's data
    for i = 1:length(speedper(:,1))
        [fJt,gofJt,JigRMSA,mt,ct,speedper,Ospeedper,newIstot,newIstot2,tagon,OJigRMS,RMS1,RMS2,OJigRMSA,OFNRMS,speedSP,speedJJ,OspeedJJ,OJ]...
            = plotsectiondata(i,fitopts,bin,DN,ODN,p,tagslip,model,fJ,gofJ,speedper,Ospeedper,newIstot,newIstot2,tagon,fun,J,RMS1,RMS2,JigRMSA,OJigRMSA,OFNRMS,speedSP,speedJJ,OspeedJJ,OJ,m,c);
        if numRMS == 1; fJt = fJt(:,1); gofJt = gofJt(:,1); end
        fJ(i,:) = fJt(i,:); gofJ(i,:) = gofJt(i,:);
        m(i) = mt(i); c(i) = ct(i);
    end
    
    [fJtot, gofJtot] = fit(RMS1(newIstot&tagon&~isnan(RMS1(:,1)),1),speedSP(newIstot&tagon&~isnan(RMS1(:,1))),'exp1',fitopts);
    [fJtot2, gofJtot2] = fit(RMS2(newIstot2&tagon&~isnan(RMS2)),speedSP(newIstot2&tagon&~isnan(RMS2)),'exp1',fitopts);
    if model
        [fJtot, gofJtot] = fit(JigRMSA(newIstot&tagon&~isnan(JigRMSA)),speedSP(newIstot&tagon&~isnan(JigRMSA)),'exp1',fitopts);
    else
        JigRMSA = RMS1(:,1);
    end
    
    for i = 1:length(speedper(:,1))
        plotalldata(i,fJtot,fJtot2,RMS1,RMS2,newIstot,newIstot2,gofJtot,gofJtot2,JigRMSA,speedSP);
    end
    
    speed = table(OspeedSP,OspeedJJ,nan(size(OspeedJJ)),zeros(size(Op)),cell(size(Op)),...
        'VariableNames',{'SP',label1,'r2','section','sectionUsed'});
    % speed.type(~isnan(OspeedSP)) = {'SP'}; speed.type(isnan(OspeedSP)) = {'NA'};
    r2 = table(nan(size(speedper(:,1))),nan(size(speedper(:,1))),cell(size(speedper(:,2))),'VariableNames',{[label1 'r2'],[label2 'r2'],'sectionUsed'});
    i = 0;
    restart = false;
    while i < length(speedper(:,1)) && ~restart
        i = i+1;
        figure(300+i); oi = get(gcf,'children');
        title({'See instructions';  'GOAL: match blue line to green dots in 2nd plot'},'parent',oi(end),'fontweight','bold','fontsize',12);
        instructions = sprintf('Controls:\nEnter: accept values\nLeft(right) click: zoom in (out)\nt: readjust this section''s thresholds (and use this section)\ns:add new speed period break at closest tag slip\nn:merge current period with next period\np: merge current period with prior period\nw: apply the cal from the whole deployment to this section\nor: press a section number to see that section''s cal applied here');
        try delete(a); catch; end
        a = annotation('textbox', [0.005, .995, 0, 0], 'string', instructions,'FitBoxToText','on');
        tpos = get(a,'position');
        for ii = length(oi)-1:length(oi);
            pos = get(oi(ii),'position'); if pos(1)<tpos(3)+.005; pos(1) = tpos(2)+.005; set(oi(ii),'position',pos); end
        end
        axes(oi(end))
        [x,~,button] = ginput(1);
        
        ci = i; %current i
        sectI = round(speedper(i,1)):speedper(i,2); % period of time between tag slips
        sectI = sectI(tagon(sectI));
        tempIs = [newIstot newIstot2];
        tempIs([1:speedper(i,1)-1,speedper(i,2)+1:end],:) = false;
        
        
        %Ic = sectI; %Ic = I(camondec(I));
        OsectI = round(Ospeedper(i,1)):Ospeedper(i,2); % period of time between tag slips
        OsectI = OsectI(Otagon(OsectI)); OIc = OsectI; %Ic = I(camondec(I));
        I = intersect(sectI,find(newIstot==1));
        I2 = intersect(sectI,find(newIstot2==1));
        fJ{length(speedper(:,1))+1,1} = fJtot; gofJ{length(speedper(:,1))+1,1} = gofJtot;
        fJ{length(speedper(:,1))+1,2} = fJtot2; gofJ{length(speedper(:,1))+1,2} = gofJtot2;
        sp1 = oi(end);  sp2 = oi(end-2); if strcmpi(get(sp2,'type'),'legend'); sp2 = oi(end-1); end
        xs = get(sp1,'xlim'); xsc = xs; %current xs;
              
        if ci>length(speedper(:,1))
            oi = cellstr(repmat('All',size(OIc')));
        else
            oi = cellstr(char(num2str(ci*ones(size(OIc')))));
        end
        speed.sectionUsed(OIc) = oi;
        for ii = 1:numRMS
            speed.r2(OIc,ii) = round(gofJ{ci,ii}.rsquare*100)/100;
        end
        speed.section(OIc) = i;
        
        r2.([label1 'r2'])(i) = round(gofJ{ci,1}.rsquare*100)/100;
        if numRMS > 1; r2.([label2 'r2'])(i) = round(gofJ{ci,2}.rsquare*100)/100; end
        if ci>length(speedper(:,1)); r2.sectionUsed(i) = {'All'}; else  r2.sectionUsed(i) = {num2str(ci)}; end
        restart = false; % restart reapplies all values, including new speedperiods
        while ~isempty(button)
            if button>48 && button <49+length(speedper(:,1))
                ci = button-48;
            elseif button == 119 %w
                ci = length(speedper(:,1))+1;
            elseif button == 1
                set([sp1 sp2], 'xlim',[x-5/24/60 x+5/24/60]); % show the 10 minutes around the zoom section
                set([sp1 sp2],'xticklabel',datestr(get(sp2,'xtick'),'HH:MM:SS'));
                xsc = [x-5/24/50 x+5/24/60];
                [x,~,button] = ginput(1);
                continue;
            elseif button == 3
                set([sp1 sp2], 'xlim',xs); xsc = xs;
                set([sp1 sp2],'xticklabel',datestr(get(sp2,'xtick'),'HH:MM:SS'));
                [x,~,button] = ginput(1);
                continue;
            elseif button == 116 %t
                delete(a);
                a = annotation('textbox', [0.005, .995, 0, 0], 'string', 'adjust thresholds below, then press enter to see the new thresholds applied','FitBoxToText','on');
                % want to use the original restrictions (Istot) that then get reduced within the script, and only want to use this section's data, thus:
%                 (Istot|tempIs(:,1))&(1:length(Istot))'>=speedper(i,1)&(1:length(Istot))'<speedper(i,2)                
                [tempIs(:,1),tempIs(:,2),minDepth2(i),maxDepth2(i),minPitch2(i),maxRR2(i)] = ...
                    adjthresh((Istot|tempIs(:,1))&(1:length(Istot))'>=speedper(i,1)&(1:length(Istot))'<speedper(i,2),(Istot2|tempIs(:,2))&(1:length(Istot2))'>=speedper(i,1)&(1:length(Istot2))'<speedper(i,2),speedper,minPitch2(i),minDepth2(i),maxDepth2(i),maxRR2(i),RMS1,RMS2,speedSP,roll,pitch,p,fs,i);
                newIstot(sectI) = tempIs(sectI,1); newIstot2(sectI) = tempIs(sectI,2);
                [fJt,gofJt,JigRMSA,mt,ct,speedper,Ospeedper,newIstot,newIstot2,tagon,OJigRMS,RMS1,RMS2,OJigRMSA,OFNRMS,speedSP,speedJJ,OspeedJJ,OJ]...
                    = plotsectiondata(i,fitopts,bin,DN,ODN,p,tagslip,model,fJ,gofJ,speedper,Ospeedper,newIstot,newIstot2,tagon,fun,J,RMS1,RMS2,JigRMSA,OJigRMSA,OFNRMS,speedSP,speedJJ,OspeedJJ,OJ,m,c);
                fJ(i,:) = fJt(i,:); gofJ(i,:) = gofJt(i,:);
                m(i) = mt(i); c(i) = ct(i);
                [fJtot, gofJtot] = fit(RMS1(newIstot&tagon&~isnan(RMS1(:,1)),1),speedSP(newIstot&tagon&~isnan(RMS1(:,1))),'exp1',fitopts);
                [fJtot2, gofJtot2] = fit(RMS2(newIstot2&tagon&~isnan(RMS2)),speedSP(newIstot2&tagon&~isnan(RMS2)),'exp1',fitopts);
                if model
                    [fJtot, gofJtot] = fit(JigRMSA(newIstot&tagon&~isnan(JigRMSA)),speedSP(newIstot&tagon&~isnan(JigRMSA)),'exp1',fitopts);
                else
                    JigRMSA = RMS1(:,1);
                end
                plotalldata(i,fJtot,fJtot2,RMS1,RMS2,newIstot,newIstot2,gofJtot,gofJtot2,JigRMSA,speedSP);
                i = i -1; break;
            elseif button == 115 %s
                [~,ii] = min(abs(ODN-x));
                [~,ii] = min(abs(tagslip(2:end-1,1)-ii));
                ii = ii+1;
                si = find(Ospeedper(:,1)<tagslip(ii,1),1,'last');
                Ospeedper = [Ospeedper(1:si-1,:);[Ospeedper(si,1) tagslip(ii,1)]; [tagslip(ii,1) Ospeedper(si,2)]; Ospeedper(si+1:end,:)];
                speedper = ceil(Ospeedper/bin);
                for iii = 1:length(threshes)
                    eval([threshes{iii} ' = [' threshes{iii} '(1:si); ' threshes{iii} '(si); ' threshes{iii} '(si+1:end)];']);
                end
                restart = true;
                break;
            elseif button == 110 % n
                Ospeedper = [Ospeedper(1:i-1,:); [Ospeedper(i,1) Ospeedper(i+1,2)]; Ospeedper(i+2:end,:)];
                speedper = ceil(Ospeedper/bin);
                restart = true;
                for iii = 1:length(threshes)
                    eval([threshes{iii} '(i) = [];']);
                end
                break;
            elseif button == 112 % p
                Ospeedper = [Ospeedper(1:i-2,:); [Ospeedper(i-1,1) Ospeedper(i,2)]; Ospeedper(i+1:end,:)];
                speedper = ceil(Ospeedper/bin);
                restart = true;
                for iii = 1:length(threshes)
                    eval([threshes{iii} '(i) = [];']);
                end
                break;
            else
                [x,~,button] = ginput(1);
                continue;
            end
            if ci>length(speedper(:,1)); sect = 'All'; else sect = num2str(ci); end
            try delete(sp2); catch; end
            sp2 = subplot(412);
            OspeedTJ = fJ{ci,1}.a*exp(fJ{ci,1}.b*OJigRMS(OIc,1));
            if model&&ci>length(speedper(:,1)); OspeedTJ = fJ{ci,1}.a*exp(fJ{ci,1}.b*OJigRMSA(OIc,1)); end
            OspeedTJ2 = fJ{ci,2}.a*exp(fJ{ci,2}.b*OFNRMS(OIc,1));
            hold on; plot(DN(sectI),speedSP(sectI),'r.','markersize',8);
            ylabel('speed (m/s)');
            plot(DN(I),speedSP(I),'g.','markersize',6);
            plot(ODN(OsectI), OspeedTJ);
            if numRMS ==2;   plot(ODN(OsectI), OspeedTJ2,'m--');   end
            legend('OCDR (Default restrictions)','OCDR (Updated restrictions)','Speed from RMS','Speed from RMS2','location','westoutside');
            set(sp2,'ylim',[0 5.5]);
            set(sp2,'xlim',xsc);
            pos2 = get(gca,'position'); pos2(3) = pos(3); pos2(4) = 1.2*pos2(4); set(gca,'position',pos2);
            oi = datestr(get(sp2,'xtick'),'HH:MM:SS'); set([sp1 sp2],'xticklabel',oi);
            linkaxes([sp1,sp2],'x')
            oi = 'Exponential';
            title([oi ' model from section # ' sect '.'],'verticalalignment','top');
            OspeedJJ(OIc,1) = OspeedTJ;
            OspeedJJ(OIc,2) = OspeedTJ2;
            [x,~,button] = ginput(1);
            if ci>length(speedper(:,1))
                oi = cellstr(repmat('All',size(OIc')));
            else
                oi = cellstr(char(num2str(ci*ones(size(OIc')))));
            end
            speed.sectionUsed(OIc) = oi;
            speed.r2(OIc,1) = round(gofJ{ci,1}.rsquare*100)/100;
            speed.r2(OIc,2) = round(gofJ{ci,2}.rsquare*100)/100;
            speed.section(OIc) = i;
            if model&&ci>length(speedper(:,1))&&isempty(button);
                OJigRMS(OIc,1) = OJigRMSA(OIc);
            end
            r2.([label1 'r2'])(i) = round(gofJ{ci,1}.rsquare*100)/100;
            r2.([label2 'r2'])(i) = round(gofJ{ci,2}.rsquare*100)/100;
            if ci>length(speedper(:,1)); r2.sectionUsed(i) = {'All'}; else  r2.sectionUsed(i) = {num2str(ci)}; end
        end
    end
end

speed.(label1) = OspeedJJ;
% everything back to normal length;
speedJJ = OspeedJJ; p = Op; speedper = Ospeedper;  JigRMS = OJigRMS; tagon = Otagon; DN = ODN; if numRMS == 2; JigRMS = [JigRMS OFNRMS]; end
% speedSP(~newIstot) = nan;
R2 = speed.r2;
thresholds = struct('minDepth',minDepth2,'maxDepth',maxDepth2,'minPitch',minPitch2,'maxRollRate',maxRR2,'minTime',minTime,'filterSize',filterSize,'minSpeed',minSpeed);
section = speed.section;
sectionUsed = speed.sectionUsed;
if numRMS== 1; speedJJ(:,2) = nan(size(speedJJ(:,1))); R2(:,2) = nan(size(R2(:,1))); end


P95 = nan(length(p),2);
P68 = P95;
C95 = P95;
P952 = P95; P682 = P68; C952 = C95;

%
% error('note for Dave- this was set at i = 2, why would that be?  Check with next set of data.');

for i = 1:length(speedper(:,1))
    sectI = round(speedper(i,1)):speedper(i,2); % period of time between tag slips
    sectI = sectI(tagon(sectI));
    P95(sectI,:) = predint(fJ{i,1},JigRMS(sectI,1),0.95,'observation','off');
    try    P952(sectI,:) = predint(fJ{i,2},JigRMS(sectI,2),0.95,'observation','off');
    catch; P952(sectI,:) = nan(length(sectI),2); end
    P68(sectI,:) = predint(fJ{i,1},JigRMS(sectI,1),0.68,'observation','off');
    try P682(sectI,:) = predint(fJ{i,2},JigRMS(sectI,2),0.68,'observation','off');
    catch; P682(sectI,:) = nan(length(sectI),2); end
    C = confint(fJ{i,1});
    try C2 = confint(fJ{i,2});
    catch; C2 = nan(2);
    end
    C95(sectI,:) = [C(1,1)*exp(C(1,2)*JigRMS(sectI,1)) C(2,1)*exp(C(2,2)*JigRMS(sectI,1))];
    try     C952(sectI,:) = [C2(1,1)*exp(C2(1,2)*JigRMS(sectI,2)) C2(2,1)*exp(C2(2,2)*JigRMS(sectI,2))];
    catch; C952(sectI,:) = nan(length(sectI),2);
    end
end
if numRMS == 2
    speedTable = table(OspeedSP, speedJJ(:,1), JigRMS(:,1),P68, P95, C95, R2(:,1), speedJJ(:,2),P682, P952, C952, R2(:,2),section,sectionUsed,...
        'VariableNames',{'SP',label1,'JRMS',[label1 'P68'],[label1 'P95'],[label1 '95'],[label1 'r2'],label2,[label2 'P68'],[label2 'P95'],[label2 '95'],[label2 'r2'],'section','sectionUsed'});
else
    speedTable = table(OspeedSP, speedJJ(:,1), JigRMS(:,1), P68, P95, C95, R2(:,1),section,sectionUsed, ...
        'VariableNames',{'SP',label1,[label1 'RMS'],[label1 'P68'],[label1 'P95'],[label1 '95'],[label1 'r2'],'section','sectionUsed'});
end
models = fJ;
modelsFit = gofJ;
if model; multiModels = m; else multiModels = []; end

sectionsendindex = speedper(:,2); sectionsendindex(end) = find(tagon,1,'last'); % ensure that speedper is applied all the way to the end and there is an exact match for apply cals later.
speedstats = struct();
speedstats.Models = models(:,1); 
speedstats.ModelFits = modelsFit(:,1); 
speedstats.r2used = r2;
speedstats.sections_end_index = sectionsendindex;
speedstats.info = 'Each row is the model for a tag orientation section.  The last row in Models is the regression on all data.  The r2 used is the r2 for each section (sometimes the whole data section was used if there was not enough data to get a good regression for an individual section';
if numRMS >1
    speedstats.(label2).Models = models(:,2);
    speedstats.(label2).ModelFits = modelsFit(:,2);
end

for ii = 1:length(multiModels);
    speedstats.multiModels{ii} = multiModels{ii}.Coefficients;
end
speedstats.Thresh = thresholds;

end

function [newIstot,newIstot2,minDepth2,maxDepth2,minPitch2,maxRR2] = adjthresh(Istot,Istot2,speedper,minPitch1,minDepth1,maxDepth1,maxRR1,RMS1,RMS2,speedSP,roll,pitch,p,fs,sectnum)
global label2 numRMS minPitch minDepth maxDepth maxRR
if ~exist('sectnum','var'); horiz = false; else horiz = true; end %horiz means plot horizontally just a single sections worth of data

button = 9;
if ~horiz; Is = Istot; sect = 'All'; else
    Istot(1:speedper(sectnum,1)-1) = false;
    Istot(speedper(sectnum,2)+1:end) = false;
    Is = Istot;
    Istot2(1:speedper(sectnum,1)-1) = false;
    Istot2(speedper(sectnum,2)+1:end) = false;
    %     = false(size(p)); Is(Is1{sectnum}) = true; Is2o = false(size(p)); Is2o(Is2{sectnum}) = true;
    sect = label2;
end

Is1 = Is; Is2 = Is; Is3 = Is; Is4 = Is;
rollrate = abs(runcirc_mean(wrapToPi([diff(roll); 0]),round(fs/2))*180/pi)*fs;
Is1(rollrate>maxRR1)=false;
Is2(abs(pitch*180/pi)<minPitch1) = false;
if(isnan(maxDepth1)); maxDepth1 = ceil(max(p)); end
Is3(p<minDepth1 | p>maxDepth1) = false;
% only uses first RMS column to find bad points
% allow for thresholds to be different for each section
if ~horiz
    minDepth2 = minDepth1*ones(size(speedper,1),1); maxDepth2 = maxDepth1*ones(size(speedper,1),1); minPitch2 = minPitch1*ones(size(speedper,1),1); maxRR2 = maxRR1*ones(size(speedper,1),1);
else minDepth2 = minDepth1; maxDepth2 = maxDepth1; minPitch2 = minPitch1; maxRR2 = maxRR1;
end
while ~isempty(button)
    sp = nan(4,1);
    if ~horiz; figure(1); clf; else subplot(2,1,2); end
    for i = 1:4
        sp(i) = subplot(2,2*(horiz+1),i+4*horiz);
        Is = Is1 & Is2 & Is3 & Is4 &Istot;
        if horiz&&i == 4; Is = Is1 & Is2 & Is3 & Is4 &Istot2; end
        %         X = RMS1(Is,1);
        %         Y = speedSP(Is);
        %         Z = zeros(length(X),1);
        switch i
            case 1
                C = rollrate(Is);
                C = -C;
                T = 'Roll Rate (deg/s)';
                D = ['Max: ' num2str(maxRR1)];
            case 2
                C = abs(pitch(Is)*180/pi);
                T = '|pitch| (deg)';
              try  D = ['Min: ' num2str(minPitch2(sectnum))];
              catch; D = ['Min: ' num2str(minPitch2(1))];
              end
            case 3
                C = p(Is);
                T = 'Depth (m)';
               try  D = ['Min: ' num2str(minDepth2(sectnum)) ', Max: ' num2str(maxDepth2(sectnum))];
               catch;  D = ['Min: ' num2str(minDepth2(1)) ', Max: ' num2str(maxDepth2(1))];
               end
            case 4
                if ~horiz
                    C = zeros(size(p));
                    for ii = 1:length(speedper(:,1))
                        C(speedper(ii,1):speedper(ii,2)) = ii;
                    end
                    C = C(Is);
                    T = 'Calibration Period';
                    D = ['Section: ' sect];
                else
                    %                     Is = Istot2;
                    C = ones(size(p)); C = C(Is);
                    T = label2;
                    D = 'All points';
                end
        end
        [C,I] = sort(C,'ascend');
        X = RMS1(Is,1);
        if i == 4 && horiz; X = RMS2(Is,1); if numRMS == 1; X = nan(sum(Is),1); end; end
        Y = speedSP(Is);
        Z = zeros(length(X),1);
        X = X(I); Y = Y(I); Z = Z(I);
        h = scatter(X,Y,50,C,'marker','.');
        ylabel('Speed (m/s)');
        if i == 1; xlabel('amplitude (dB)'); end
        title([T ', ' D]);
        xs = get(gca,'xlim'); ys = get(gca,'ylim');
        text(xs(1),ys(2),num2str(i),'horizontalalignment','left','verticalalignment','top','fontsize',20);
        if i ==4
            xl = xlabel({'Left click on the Color Bar to cut off bluer values (lower in the bar). Or right click Depth to cut off red values (higher in the bar)'; 'Click a plot to return it to the original default value.';  'Press enter to continue to next section.'},'horizontalalignment','right');
        end
        if ~horiz || i < 4; c = colorbar;
            try if i == 1; set(c,'yticklabel',num2str(-str2num(get(c,'yticklabel')))); end
            catch;  if i == 1; set(c,'yticklabel',num2str(-cellfun(@str2num,get(c,'yticklabel')))); end
            end
        end
    end
    again = true;
    vers = version('-release'); if strcmp(vers(end),'a'); vers = str2num(vers(1:4)); else vers = str2num(vers(1:4))+0.1; end
    axs = get(gcf,'children');

    
    while again
        again = false;
        if vers>2014 %accounts for recent versions of matlab that treat colorbars differently than older versions.
            if horiz; AA = 6:-2:2; else AA = 7:-2:1; end
            for iii = AA
                axes('position',get(axs(iii),'position'),'ylim',get(axs(iii),'Limits'),'color','none','yticklabel',[]);
            end
        end
        [x,y,button] = ginput(1);
        axs1 = get(gcf,'children');
        axs = axs1;
        
        if vers>2014 && ~horiz
            axs([5 7 9 11]) = axs(1:4); axs(1:4) = []; 
        elseif vers>2014
            axs([5 7 9]) = axs(1:3); axs(1:3) = []; 
        end
        ii = find(axs == gca);
        if horiz && ii ~=1; ii = ii+1; end
        secti = 1:length(minDepth2);
        switch ii
            case 1
                if ~horiz
                    Is4([1:speedper(round(y),1)-1 speedper(round(y),2)+1:length(Is)]) = false;
                    sect = num2str(round(y));
                    secti = round(y);
                else again = true;
                    secti = 1:length(minDepth2);
                end
            case 2
                Is4 = true(size(Is));
                sect = 'All';
            case 3
                if button == 1
                    Is3(p<round(y)) = false;
                    minDepth2(secti) = round(y);
                elseif button == 3
                    Is3(p>round(y)) = false;
                    maxDepth2(secti) = round(y);
                end
            case 4
                if button == 1
                    minDepth2(secti) = minDepth;
                    Is3(p>=minDepth) = true;
                elseif button == 3;
                    maxDepth2(secti) = maxDepth;
                    Is3(p>=minDepth) = true;
                end
            case 5
                Is2(abs(pitch*180/pi)<round(y)) = false;
                minPitch2(secti) = round(y);
            case 6
                minPitch2(secti) = minPitch;
                Is2(abs(pitch*180/pi)>=minPitch) = true;
            case 7
                maxRR2(secti) = abs(round(y));
                Is1(rollrate>abs(round(y)))=false;
            case 8
                maxRR2(secti) = maxRR;
                Is1(rollrate<=maxRR) = true;
            otherwise
                again = true;
        end
        if isempty(button); again = false; end
    end
    
end
newIstot = Is1 & Is2 & Is3 & Is4 &Istot;
newIstot2 = Is1 & Is2 & Is3 & Is4 &Istot2;
end

function [speedper,toosmall] = checksmallsections (speedper,newIstot,S1)
toosmall = [];

if length(speedper(:,1))>1
    for i = 1:length(speedper(:,1))
        sectI = newIstot;
        sectI([1:speedper(i,1)-1 speedper(i,2)+1:end])= false;
        if sum(sectI)<20; toosmall = [toosmall; i];
        end
    end
end
speedper(toosmall,:) = []; speedper = [[S1(1);speedper(2:end,1)] [speedper(2:end,1); S1(end,2)]];
% Is1(toosmall) = []; Is2(toosmall) = [];
if ~isempty(toosmall); disp(['Sections ' num2str(toosmall') ' now too small, also joining with adjacent sections']); end
end

function [RMS1,OJigRMS,speedJJ,OspeedJJ,mA,cA] = makemultimodel(I,sectI,fun,J,RMS1,OJigRMS,speedSP,speedJJ,OspeedJJ,OJ,OsectI)
IsJ = I(:,1)&~isnan(RMS1(:,1))&~isnan(speedSP);
opts = statset('nlinfit');
opts.RobustWgtFun = 'bisquare';
try
    mdl = fitnlm(J(IsJ,:),speedSP(IsJ),fun,[1.3 .02 repmat(1/3,1,size(J,2))],'Options',opts);
catch err
    disp('Error: likely not enough data points in a section of data to make a model.  Try increasing the size of the tag orientation periods');
    throw( err);
end
mA = mdl;
cA = mdl.Coefficients.Estimate';
speedJJ(sectI,1) = fun(cA,J(sectI,:));% should sectI be IsJ (that's how it was in the first model usage)
OspeedJJ(OsectI,1) = fun(cA,OJ(OsectI,:));
RMS1(sectI,1) = log(fun([1 1 cA(3:end)],J(sectI,:)));
OJigRMS(OsectI,1) = log(fun([1 1 cA(3:end)],OJ(OsectI,:)));
if any(isinf(RMS1(sectI,1))); RMS1(sectI,1) = cA(3)*J(sectI,1) + cA(4)*J(sectI,2) + cA(5)*J(sectI,3); end
if any(isinf(OJigRMS(sectI,1))); OJigRMS(sectI,1) = cA(3)*OJ(sectI,1) + cA(4)*OJ(sectI,2) + cA(5)*OJ(sectI,3); end
end

function [fJ,gofJ,JigRMSA,m,c,speedper,Ospeedper,newIstot,newIstot2,tagon,OJigRMS,RMS1,RMS2,OJigRMSA,OFNRMS,speedSP,speedJJ,OspeedJJ,OJ]...
    = plotsectiondata(i,fitopts,bin,DN,ODN,p,tagslip,model,fJ,gofJ,speedper,Ospeedper,newIstot,newIstot2,tagon,fun,J,RMS1,RMS2,JigRMSA,OJigRMSA,OFNRMS,speedSP,speedJJ,OspeedJJ,OJ,m,c)
global numRMS label1 label2
% for i = 1:length(speedper(:,1))
I = [newIstot&tagon newIstot2&tagon]; % the points to compare in RMS1 and RMS2
I([1:speedper(i,1)-1 speedper(i,2)+1:end],:) = false;
sectI = speedper(i,1):speedper(i,2);
OsectI = Ospeedper(i,1):Ospeedper(i,2);
if model
    [RMS1,OJigRMS,speedJJ,OspeedJJ,m{i},c{i}] = makemultimodel(I,sectI,fun,J,RMS1,OJigRMSA,speedSP,speedJJ,OspeedJJ,OJ,OsectI);
else
    OJigRMS = OJigRMSA;
end

for ii = 1:2
    switch ii; case 1; RMS = RMS1; case 2; RMS = RMS2; end
    IsJ = I(:,ii)&~isnan(RMS)&~isnan(speedSP);
    try
        [fJ{i,ii}, gofJ{i,ii}] = fit(RMS(IsJ),speedSP(IsJ),'exp1',fitopts);
    catch err
        if ii == 2 && length(sectI) == sum(isnan(RMS(sectI))) % if RMS2 just has no values (i.e. the hydrophone is off for flow noise), nan out all resulting speed but continue
            disp(['Section ' num2str(i) ' of RMS2 has no points, skipping this section for this variable.']);
            [fJ{i,ii}, gofJ{i,ii}] = fit(zeros(50,1),zeros(50,1),'exp1',fitopts);
        else
            disp(['Error: likely not enough data points section ' num2str(i) ' of RMS' num2str(ii) ' data to make a model.  Try increasing the size of the tag orientation periods']);
            throw(err);
        end
    end
    speedJJ(sectI,ii) = fJ{i,ii}.a*exp(fJ{i,ii}.b*RMS(sectI));
    if ii == 1
        OspeedJJ(OsectI,ii) = fJ{i,ii}.a*exp(fJ{i,ii}.b*OJigRMS(OsectI));
    else
        OspeedJJ(OsectI,ii) = fJ{i,ii}.a*exp(fJ{i,ii}.b*OFNRMS(OsectI));
    end
end
IsJ2 = IsJ;
IsJ = I(:,1)&~isnan(RMS1(:,1))&~isnan(speedSP);%&s
fig = figure(300+i); clf;
set(fig,'windowStyle','docked');
sp(1) = subplot(411);
h1 = plot(DN(sectI),p(sectI));
ylabel('Depth (m)');
ax = gca;
set(ax,'xlim',[DN(sectI(1)) DN(sectI(end))],'nextplot','add');
set(ax(1),'ydir','rev','ylim',[0 max(p(sectI))]);
% pos = get(gca,'position');
oi = datestr(get(ax(1),'xtick'),'HH:MM:SS'); set(ax(1),'xticklabel',oi);
plot(ax(1),[DN(round(tagslip(:,1)/bin)) DN(round(tagslip(:,1)/bin))]', repmat([-10 1000],length(tagslip(:,1)),1)','k','linewidth',2);
TEX = text(DN(round(tagslip(:,1)/bin)),repmat(max(p(sectI)),length(tagslip(:,1)),1),'Tag Slip','parent',ax(1),'color','k','rotation',90,'verticalalignment','top');
%     legend(h1,'Depth','location','eastoutside');

sp(2) = subplot(412);
hold on; plot(DN(sectI),speedSP(sectI),'r.','markersize',8);
ylabel('speed (m/s)');
plot(DN(I(:,1)),speedSP(I(:,1)),'g.','markersize',6);
%     plot(DN(sectI),speedSP(sectI),'g.','markersize',6);
cols = {'b-';'m--'};
for ii = 1:numRMS
    plot(ODN(OsectI), OspeedJJ(OsectI,ii),cols{ii});
end
legend('OCDR (Default restrictions)','OCDR (Updated restrictions)',['Speed from ' label1 'RMS'],['Speed from ' label2],'location','westoutside');
set(sp(2),'ylim',[0 5.5]);
set(sp(2),'xlim',DN([sectI(1) sectI(end)]));
pos2 = get(gca,'position');  pos2(4) = 1.2*pos2(4); set(gca,'position',pos2);% pos2(3) = pos(3);
oi = datestr(get(gca,'xtick'),'HH:MM:SS'); set(gca,'xticklabel',oi);
linkaxes([sp(2),ax],'x')
pos1 = get(sp(1),'position'); pos1([1 3]) = pos2([1 3]); pos1(4) = 0.8*pos1(4); set(sp(1),'position',pos1);

sp(5) = subplot(4,2,[5 7]); legend;%
hold off;
plot(RMS1(IsJ,1),speedSP(IsJ),'.');
hold on
xs = get(sp(5),'xlim'); xs = xs(1):diff(xs)/1000:xs(2); ys = fJ{i,1}.a*exp(fJ{i,1}.b*xs); %ys = fJ2{i}.p1*xs.^2 + fJ2{i}.p2*xs + fJ2{i}.p3;
plot(xs,ys,'m','linewidth',2);
B = predint(fJ{i,1},xs,0.68,'observation','off');
xt = get(sp(5),'xlim');  yt = get(sp(5),'ylim');
plot(xs,B(:,1),'m--'); plot(-1000,1000,'g.'); plot(-1000,1000,'color',[0 100 0]/255,'linewidth',4);
plot(xs,B(:,2),'m--');
set(sp(5),'xlim',xt,'ylim',yt);
legend('Sect. data','Exp','68% pred. int','All data','Regr on all','location','Northwest');
xt = xt(1)+diff(xt)*5/6;
yt = [yt(1)+diff(yt)*.05 yt(1)+diff(yt)*.12];
text(xt,yt(1),['R^{2} = ' num2str(round(100*gofJ{i,1}.rsquare)/100)],'color','m','fontweight','bold');
ylabel('OCDR ({\Delta}d/sin(pitch))'); xlabel([label1 ' RMS']);

if numRMS == 2
    sp(4) = subplot(4,2,[6 8]); legend;%
    hold off;
    plot(RMS2(IsJ2),speedSP(IsJ2),'.'); hold on
    xs = get(sp(4),'xlim'); xs = xs(1):diff(xs)/1000:xs(2); ys = fJ{i,2}.a*exp(fJ{i,2}.b*xs); %ys = fJ2{i}.p1*xs.^2 + fJ2{i}.p2*xs + fJ2{i}.p3;
    plot(xs,ys,'m','linewidth',2);
    xt = get(sp(4),'xlim');  yt = get(sp(4),'ylim');
    try B = predint(fJ{i,2},xs,0.68,'observation','off');
        plot(xs,B(:,1),'m--'); plot(-1000,1000,'g.'); plot(-1000,1000,'color',[0 100 0]/255,'linewidth',4);
        plot(xs,B(:,2),'m--');
    catch; disp(['Too few observations for a prediction interval in section ' num2str(i) ' for JRMS2']);
    end
    
    set(sp(4),'xlim',xt,'ylim',yt);
    legend('Sect. data','Exp','68% pred. int','All data','Regr on all','location','Northwest');
    xt = xt(1)+diff(xt)*5/6;
    yt = [yt(1)+diff(yt)*.05 yt(1)+diff(yt)*.12];
    text(xt,yt(1),['R^{2} = ' num2str(round(100*gofJ{i,2}.rsquare)/100)],'color','m','fontweight','bold');
    ylabel('OCDR ({\Delta}d/sin(pitch))'); xlabel([label2 ' RMS']);
end
% end
end

function plotalldata(i,fJtot,fJtot2,RMS1,RMS2,newIstot,newIstot2,gofJtot,gofJtot2,JigRMSA,speedSP)
% for i = 1:length(speedper(:,1))
global numRMS
figure(300+i); oi = get(gcf,'children');
sp(5) = oi(2);
if numRMS == 2; sp(5) = oi(4); sp(4) = oi(2);
    xs = get(sp(4),'xlim'); xs = xs(1):diff(xs)/1000:xs(2); ys = fJtot2.a*exp(fJtot2.b*xs);
    plot(sp(4),RMS2(newIstot2),speedSP(newIstot2),'g.'); plot(sp(4),xs,ys,'linewidth',4,'color',[0 100 0]/255);
    ch = get(sp(4),'children'); ch = [ch(3:end); ch(1:2)];
    set(sp(4),'children',ch);
    yt = get(sp(4),'ylim'); xt = get(sp(4),'xlim');
    xt = xt(1)+diff(xt)*5/6;
    yt = [yt(1)+diff(yt)*.05 yt(1)+diff(yt)*.12];
    text(xt,yt(2),['R^{2} = ' num2str(round(100*gofJtot2.rsquare)/100)],'color',[0 100 0]/255,'fontweight','bold','parent',sp(4));
end
xs = get(sp(5),'xlim'); xs = xs(1):diff(xs)/1000:xs(2); ys = fJtot.a*exp(fJtot.b*xs);
plot(sp(5),JigRMSA(newIstot,1),speedSP(newIstot),'g.'); plot(sp(5),xs,ys,'linewidth',4,'color',[0 100 0]/255);
ch = get(sp(5),'children'); ch = [ch(3:end); ch(1:2)];
set(sp(5),'children',ch);
yt = get(sp(5),'ylim'); xt = get(sp(5),'xlim');
xt = xt(1)+diff(xt)*5/6;
yt = [yt(1)+diff(yt)*.05 yt(1)+diff(yt)*.12];
text(xt,yt(2),['R^{2} = ' num2str(round(100*gofJtot.rsquare)/100)],'color',[0 100 0]/255,'fontweight','bold','parent',sp(5));
end

%public subfunctions
% circ_mean.m
% (c) Phillipp Berens, 2009
% https://www.mathworks.com/matlabcentral/fileexchange/10676-circular-statistics-toolbox--directional-statistics-

function [mu ul ll] = circ_mean(alpha, w, dim)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_mean(alpha, [], dim)
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit 
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al. 
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 3
  dim = 1;
end

if nargin < 2 || isempty(w)
  % if no specific weighting has been specified
  % assume no binning has taken place
	w = ones(size(alpha));
else
  if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1) 
    error('Input dimensions do not match');
  end 
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
  t = circ_confmean(alpha,0.05,w,[],dim);
  ul = mu + t;
  ll = mu - t;
end
end

%
% consec
% (c) Dave Mellinger, 1998
% https://www.mathworks.com/matlabcentral/fileexchange/70-osprey/content/osprey/utils/consec.m
function [y,z] = consec(x)
%CONSEC		Run-length encode consecutive integers as (start,end) runs.
%
% y = consec(x)
%   Given a vector of numbers x that includes sequences of consecutive
%   ascending integers, encode it as (start,end) pairs, where each pair
%   stands for one sequence of consecutive numbers.  Each pair forms a
%   column of the return value y -- that is, y is a 2-by-N array.
%
%   Example:
%          consec([4 5 6 17 18 2 11 12 13]) ==> [4 17 2 11]
%                                               [6 18 2 13]
%
% [starts,stops] = consec(x)
%   As above, but return two 1-by-N vectors instead of one 2-by-N array.

n = length(x);
if (n)
  s = find(x(1:n-1) ~= x(2:n)-1);
  if (nCols(x) == 1)
    y = x([1 s'+1])';			% x is a column vector
    z = x([s' n])';
  else
    y = x([1 s+1]);			% x is a row vector
    z = x([s n]);
  end
else
  y = [];
  z = [];
end

if (nargout < 2)
  y = [y; z];
end
end

%
% nCols
% (c) Dave Mellinger, 1998
% https://www.mathworks.com/matlabcentral/fileexchange/70-osprey
function n = nCols(array)
%   nCols(array) returns the number of columns in the array.  This is the
%   size of the second dimension, i.e., it is shorthand for size(array,2).  
%   For n-dimensional arrays, this is different from the 'n' value given
%   by "[m,n] = size(array)", as 'size' with two output arguments rolls
%   together all the dimensions after the first one.
%
% See also nRows, size.

n = size(array,2);
end

%
% runcirc_mean.m
% (c) David Cade, 2016
function y = runcirc_mean (x,m) %input in radians, not fast

if sum(size(x)>1)<2
X = buffer(x,2*m,2*m-1,x(1)*ones(2*m-1,1));
X = X(:,m:end);
X2 = rot90(buffer(flipud(x(end-(2*m-1):end)),2*m,2*m-1,x(end)*ones(2*m-1,1)),2);
X = [X X2(:,2:m)];
y = circ_mean(X);
if size(x,1)>size(x,2); y = y'; end
else
    y = nan(size(x));
    for i = 1:size(x,2)
        x2 = x(:,i);
        X = buffer(x2,2*m,2*m-1,x2(1)*ones(2*m-1,1));
        X = X(:,m:end);
        X2 = rot90(buffer(flipud(x2(end-(2*m-1):end)),2*m,2*m-1,x2(end)*ones(2*m-1,1)),2);
        X = [X X2(:,2:m)];
        y(:,i) = circ_mean(X)';
    end
end
end

%
% runmean.m
% (c) Jos van der Geest 2006
% http://www.mathworks.com/matlabcentral/fileexchange/10113-runmean
function Y = runmean(X, m, dim, modestr) ;
% RUNMEAN - Very fast running mean (aka moving average) filter
%   For vectors, Y = RUNMEAN(X,M) computes a running mean (also known as
%   moving average) on the elements of the vector X. It uses a window of
%   2*M+1 datapoints. M an positive integer defining (half) the size of the
%   window. In pseudo code: 
%     Y(i) = sum(X(j)) / (2*M+1), for j = (i-M):(i+M), and i=1:length(X) 
%
%   For matrices, Y = RUNMEAN(X,M) or RUNMEAN(X,M,[]) operates on the first   
%   non-singleton dimension of X. RUNMEAN(X,M,DIM) computes the running
%   mean along the dimension DIM.
%
%   If the total window size (2*M+1) is larger than the size in dimension
%   DIM, the overall average along dimension DIM is computed.
%
%   As always with filtering, the values of Y can be inaccurate at the
%   edges. RUNMEAN(..., MODESTR) determines how the edges are treated. MODESTR can be
%   one of the following strings:
%     'edge'    : X is padded with first and last values along dimension
%                 DIM (default)
%     'zero'    : X is padded with zeros
%     'mean'    : X is padded with the mean along dimension DIM 
%
%   X should not contains NaNs, yielding an all NaN result. NaNs can be
%   replaced by using, e.g., "inpaint_nans" created by John D'Errico.
%
%   Examples
%     runmean([1:5],1) 
%       % ->  1.33  2  3  4 4.67
%     runmean([1:5],1,'mean') 
%       % ->  2 2 3 4 4
%     runmean([2:2:10],1,1) % dimension 1 is larger than 2*(M=1)+1 ...
%       % -> 2 4 6 8 10
%     runmean(ones(10,7),3,2,'zero') ; % along columns, using mode 'zero'
%     runmean(repmat([1 2 4 8 NaN 5 6],5,1),2,2) ; 
%       % -> all NaN result
%     A = rand(10,10) ; A(2,7) = NaN ;
%     runmean(A,3,2) ; 
%       % -> column 7 is all NaN
%     runmean(1:2:10,100) % mean
%       % -> 5 5 5 5 5
%
%   This is an incredibly fast implementation of a running mean, since
%   execution time does not depend on the size of the window.
%
%   See also MEAN, FILTER

% for Matlab R13
% version 3.0 (sep 2006)
% Jos van der Geest
% email: jos@jasen.nl

% History:
%   1.0 (2003) created, after a snippet from Peter Acklam (?)
%   1.1 (feb 2006) made suitable for the File Exchange (extended help and
%       documentation)
%   1.2 (feb 2006) added a warning when the window size is too big
%   1.3 (feb 2006) improved help section
%   2.0 (sep 2006) working across a dimension of a matrix. 
%   3.0 (sep 2006) several treatments of the edges. 

% Acknowledgements: (sep 2006) Thanks to Markus Hahn for the idea of
% working in multi-dimensions and the way to treat edges. 

error(nargchk(2,4,nargin)) ;

if ~isnumeric(m) || (numel(m) ~= 1) || (m < 0) || fix(m) ~= m,
    error('The window size (M) should be a positive integer') ;
end

if nargin == 2,
    dim = [] ;
    modestr = 'edge' ;
elseif nargin==3,
    if ischar(dim),
        % no dimension given
        modestr = dim ;
        dim = [] ;
    else 
        modestr = 'edge' ;
    end
end

modestr = lower(modestr) ;

% check mode specifier
if ~ismember(modestr,{'edge','zero','mean'}),
    error('Unknown mode') ;
end

szX = size(X) ;
if isempty(dim),
    dim = min(find(szX>1)) ;
end

if m == 0 || dim > ndims(X),
    % easy
    Y = X ;
else
    mm = 2*m+1 ;
    if mm >= szX(dim),
        % if the window is larger than X, average all
        sz2 = ones(size(szX)) ; 
        sz2(dim) = szX(dim) ;
        Y = repmat(mean(X,dim),sz2) ;
    else
        % here starts the real stuff
        % shift dimensions so that the desired dimensions comes first
        [X, nshifts] = shiftdim(X, dim-1); 
        szX = size(X) ;
        % make the rest of the dimensions columns, so we have a 2D matrix
        % (suggested of Markus Hahn)
        X = reshape(X,szX(1),[]) ; 
        % select how to pad the matrix
        switch (modestr),
            case 'edge'
                % pad with first and last elements
                Xfirst = repmat(X(1,:),m,1) ;
                Xlast = repmat(X(end,:),m,1) ;
            case 'zero'
                % pad with zeros
                Xfirst = zeros(m,1) ;
                Xlast= zeros(m,1) ;
            case 'mean',
                % pad with the average
                Xfirst = repmat(mean(X,1),m,1) ;
                Xlast = Xfirst ;
        end        
        % pad the array
        Y = [zeros(1,size(X,2)) ; Xfirst ; X ; Xlast] ;       
        % the cumsum trick (by Peter Acklam ?)
        Y = cumsum(Y,1) ;
        Y = (Y(mm+1:end,:)-Y(1:end-mm,:)) ./ mm ;
        
        % reshape into original size
        Y = reshape(Y,szX)   ;
        % and re-shift the dimensions
        Y = shiftdim(Y,ndims(Y)-nshifts) ;
    end
end

% =====================
%  CODE OF VERSION 1.3 
% =====================

% function Y = runmean(X,m) ;
% % RUNMEAN - Very fast running mean filter for vectors
% %   Y = RUNMEAN(X,M) computes a running mean on vector X using a window of
% %   2*M+1 datapoints. X is a vector, and M an positive integer defining
% %   (half) the size of the window. In pseudo code:
% %     Y(i) = sum(X(j)) / (2*M+1), for j = (i-M):(i+M), and i=1:length(X)
% %
% %   If the total window size (2M+1) is larger than the length of the vector, the overall
% %   average is returned.
% %
% %   Example:
% %     runmean(1:10,1) % ->
% %     [1.3333 2 3 4 5 6 7 8 9 9.6667]
% %
% %   This is an incredibly fast implementation of a running average, since
% %   execution time does not depend on the size of the window.
% %
% %   X should not contains NaNs (a NaN will result in a all NaN result)
% %   At both ends the values of Y can be inaccurate, as the first and last
% %   values of X are used multiple times. 
% %
% %   See also MEAN
% 
% % for Matlab R13
% % version 1.3 (feb 2006)
% % Jos van der Geest
% % email: jos@jasen.nl
% 
% % History:
% % 1.0 (2003) created, after a snippet from Peter Acklam (?)
% % 1.1 (feb 2006) made suitable for the File Exchange (extended help and
% % documentation)
% % 1.2 (feb 2006) added a warning when the window size is too big
% % 1.3 (feb 2006) improved help section
% 
% error(nargchk(2,2,nargin)) ;
% 
% sz = size(X) ;
% 
% if numel(sz) ~= 2 || (min(sz) ~= 1),
%     error('X should be a vector') ;
% end
% 
% if any(isnan(X)),
%     error('NaNs cannot be dealt with') ;
% end
% 
% if ~isnumeric(m) || (numel(m) ~= 1) || (m < 0) || fix(m) ~= m,
%     error('The window size (M) should be a positive integer') ;
% elseif m == 0,
%     Y = X ;
%     return ;
% end
% 
% mm = 2*m+1 ;
% 
% if mm >= prod(sz),
%     % if the window is larger than X, average all
%     warning('Window size is larger than the length of the vector.')
%     Y = repmat(mean(X),sz) ;
% else
%     % the cumsum trick ...
%     Y = [repmat(X(1),m,1) ; X(:) ; repmat(X(end),m,1)] ;
%     Y = [0 ; cumsum(Y)] ;
%     Y = (Y(mm+1:end)-Y(1:end-mm)) / mm ;
%     Y = reshape(Y,sz) ;
% end
end
