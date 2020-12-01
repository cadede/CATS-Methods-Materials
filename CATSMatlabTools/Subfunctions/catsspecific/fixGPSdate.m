function [GDN,fixjump] = fixGPSdate(DN,GDN,DHz,GHz,fixjump)
global fixweird
if nargin<5; fixjump = false; end

if length(DN)~=length(GDN); error('Lengths not the same'); end
if nargin<4||isempty(GHz)
    GHz = diff(GDN);
    GHz(GHz==0) = nan;
    GHz = round(nanmedian(1./(GHz*24*60*60)));
    if isnan(GHz)||GHz<1; warning('GPS sample rate not detected, using 10Hz'); GHz = 10; end
end
if nargin<3||isempty(DHz)
    DHz = diff(DN);
    DHz(DHz==0) = nan;
    DHz = round(nanmedian(1./(DHz*24*60*60)));
end

if length(GDN)~=length(DN); error('GDN & DN need to be the same length'); end

I = find(diff(GDN)~=0);
tdiff = nanmedian(GDN(I)-DN(I));
I = find(abs(GDN-DN-tdiff)>DN(end)-DN(1)); % find where the jump is more than 5 seconds
GDN(I) = nan; %GDN = fixgaps(GDN); 
disp(['Removed ' num2str(length(I)) ' points that had bad date stamps (more off than the length of the deployment and replaced with an interpolation']);


bigjumps = find(abs(diff(GDN))>1);
if ~isempty(bigjumps)
    if fixjump; button = 50;
    else
        figure
        subplot(2,1,1);
        plot(GDN); set(gca,'yticklabel',datestr(get(gca,'ytick'),'yyyy-mmm-dd HH:MM:SS'));
        title ('Original times: press 1 for this one');
    end
    GDN2 = GDN;
    if mod(length(bigjumps),2) == 1; bigjumps(end+1) = length(GDN); end
    for i = 1:2:length(bigjumps)
        oldDN = floor(GDN2(bigjumps(i)));
        newDN = floor(GDN2(bigjumps(i)+1));
        GDN2(bigjumps(i)+1:bigjumps(i+1)) = GDN2(bigjumps(i)+1:bigjumps(i+1)) - (newDN-oldDN);
    end
    if ~fixjump
        subplot(2,1,2);
        plot(GDN2); set(gca,'yticklabel',datestr(get(gca,'ytick'),'yyyy-mmm-dd HH:MM:SS'));
        title ('new times: press 2 for this one');
        [~,~,button] = ginput(1);
    end
    if button == 50; GDN = GDN2; fixjump = true;
    elseif button ~=49
        error('Must choose the original or corrected GPS date numbers');
    else fixjump = false;
    end
else fixjump = false;
end
    
sames = find(diff(GDN) == 0)+1;
GDN(sames) = nan;
II = find(~isnan(GDN));
if any(diff(II)./(diff(GDN(II))*24*60*60)>3*DHz/GHz); % if there are any points that have more than 3 seconds off from where they should be
    III = II(find(diff(II)./(diff(GDN(II))*24*60*60)>3*DHz/GHz)+1);
    numadj = 0;
    while ~isempty(III)
        GDN(III(1)) = nan;
        II = find(~isnan(GDN));
        III = II(find(diff(II)./(diff(GDN(II))*24*60*60)>3*DHz/GHz)+1);
        numadj = numadj+1;
        if numadj>length(GDN); error('?'); end
    end
    if isempty(fixweird)
    fixweird = 1; % use to allow for user input on whether to error out of this input([num2str(numadj) ' points seem to be misplaced.  Okay to erase the time stamps (1 = yes, 2 = no)?']);
    warning([num2str(numadj) ' points seem to be misplaced. Removing them from the data set']);
    if fixweird == 2; error('Misplaced timestamps'); end
    end
end

if isnan(GDN(end)); ii = find(~isnan(GDN),1,'last'); GDN(end) = round((DN(end)-DN(ii))*GHz*24*60*60)/24/GHz/60/60 + GDN(ii); end
GDN = fixgaps(GDN);
GDN = round(GDN*GHz*24*60*60)/GHz/24/60/60;



% sames = [find(diff(GDN) == 0)+1; I]; % find repeated values and also where the nans from the big spikes were, they will be removed
% I = 1:length(GDN); I = I';
% I(sames) = nan;
% ii = I(1);
% for i = 1:length(I)
%     if isnan(I(i))
%         I(i) = ii;
%     else
%         ii  = I(i);
%     end
% end
% DNdiff = floor((DN-DN(I))*24*60*60*GHz)/GHz/24/60/60;
% GDN(sames) = GDN(sames)+DNdiff(sames);

% badsames = sames(sames-[zeros(DHz/GHz,1); sames(1:length(sames
% firstcorrect = [find(~isnan(GDN),1); sames(find(diff(sames)>1)+1)-1];
% timesince = DN
% 
% for i = ceil(DHz/GHz)+1:length(GDN)
%     if GDN(i)<=GDN(i-DHz/GHz)
%         GDN(i) = GDN(i)+round(diff([DN(i) DN(find(GDN<=GDN(i),1,'last'))])*24*60*60*GHz)/24/60/60/GHz;
%     end
% end
%     