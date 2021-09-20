
function [vidframes,numbad,badsect,DROP] = checkbadframes(vidframes,insecs,replsingles,vidnum)
if nargin<2; insecs = false; end
if nargin<3; replsingles = false; end
if size(vidframes,1)>size(vidframes,2); vidframes = vidframes'; end
if ~insecs
    vidstart = vidframes(find(~isnan(vidframes),1));
    vidstartI = find(~isnan(vidframes),1);
    vidframes = (vidframes - vidstart)*24*60*60; % get to units of seconds
end
vidframeso=vidframes;
% first guess at fixes for smoothing later
fdif = mean([prctile(diff(vidframes),45) prctile(diff(vidframes),55)]);
badframes = isnan(vidframes);% come back to frame 1 if it is bad
numbad = sum(badframes);
for i = find(badframes);
    if i~=1; vidframes(i) = vidframes(i-1)+fdif; end
end
bf2 = fliplr(find(isnan(vidframes))); % now go backwards in case the first frames were bad
for i = bf2
    vidframes(i) = vidframes(i+1)-fdif;
end
dV = diff(vidframes);
% meandif = mean(dV(abs(dV)<4*fdif));
meandif = mean(dV);

% find a line that goes closest to the most points
optfun = @(b) sum(abs(meandif*(1:length(vidframes))+b-vidframes)>.5*meandif);
b = -100*meandif:meandif/10:100*meandif; v = nan(size(b));
for i = 1:length(b);
    v(i) = optfun(b(i)); 
end
[~,i] = min(v);
b = b(i);
fun = @(t) meandif*t+b; % best fit function
I = 1:length(vidframes);
meanvid = fun(I);


badsect = 0;
drops = find(diff(vidframes)<0);
DROP = nan(2,length(drops));
i = 1;
while i<=length(drops)
      ldif = abs(vidframes(drops(i))-vidframes(drops(i)+1));
    if replsingles && vidframes(drops(i)+2)>vidframes(drops(i))
        vidframes(drops(i)+1:drops(i)+2) = nan; vidframes = fixgaps(vidframes);
        disp(['replaced a single bad frame (# ' num2str(drops(i)+1) ') in video # ' num2str(vidnum)]);
        continue;
    end
    sectstart = find(diff(vidframes(1:drops(i)))>fdif*1.5,1,'last')+1; %find the jump up that is bigger than one normal frame, move all down to match.
    if isempty(sectstart); sectstart = 1; end
    sectend = find(diff(vidframes(drops(i)+1:end))>fdif*1.5,1,'first')+drops(i);
    if isempty(sectend); sectend = length(vidframes); end
    
     II =  max(drops(i)-4*(drops(i)-sectstart),1):min(drops(i)+2*(sectend-drops(i)),length(vidframes)); meanvidloc = meanvid;
    if sum(vidframes(II)>meanvid(II))>.8*length(II) || sum(vidframes(II)>meanvid(II))<.2*length(II); % if you have more than 80% of the local points are above the line, recalculate line for this section
        optfun = @(b) sum(abs(meandif*(II)+b-vidframes(II))>.5*meandif);
        b = -100*meandif:meandif/10:100*meandif; v = nan(size(b));  for ii = 1:length(b); v(ii) = optfun(b(ii)); end; [~,ii] = min(v); b = b(ii);
        fun = @(t) meandif*t+b; % best fit function
        meanvidloc(II) = fun(II);
    end
    %     [~,dir] = min([abs(drops(i)-sectstart),abs(drops(i)-sectend)]); % assumes that the drop is associated with the shorter of the section of the drop to the edge or to the drop in the opposite direction
    %     if vidframes(drops(i)+1)<meanvid(drops(i)); dir = 2; else dir = 1; end
    [~,dir] = max(abs(vidframes(drops(i):drops(i)+1)-meanvidloc(drops(i):drops(i)+1)));
    %     xs = max(drops(i)-100,1):min(drops(i)+100,length(vidframes));  figure; subplot(211); plot(I,meanvidloc,'--',I,vidframeso,'m',I,vidframes); xlim([min(xs) max(xs)]); ylim([min([vidframeso(xs) vidframes(xs)]) max([vidframeso(xs) vidframes(xs)])]);
    if sectstart == 1 && sectend == length(vidframes);
       dir = 2; % if it's the whole length with just a drop, assume it started right. 
    end
    if dir == 1; vidframes(sectstart:drops(i)) = vidframes(sectstart:drops(i)) - ldif - fdif;
        badsect = badsect+(drops(i)-sectstart)+1;
        if isnan(DROP(2,i)) || DROP(2,i)<drops(i); DROP(2,i) = drops(i); end
        if isnan(DROP(1,i)) || DROP(1,i)>sectstart; DROP(1,i) = sectstart; end
        if sectstart>1&& diff(vidframes(sectstart-1:sectstart))<0; drops(i) = sectstart-1; i = i-1; end
               
    else   vidframes(drops(i)+1:sectend) = vidframes(drops(i)+1:sectend) + ldif+fdif;
        badsect = badsect+sectend-drops(i);
        if isnan(DROP(1,i)) || DROP(1,i)>drops(i); DROP(1,i) = drops(i); end
        if sectend<length(vidframes) && diff(vidframes(sectend:sectend+1))<0; drops(i) = sectend; i = i-1;
        else DROP(2,i) = sectend;
        end
    end
    % use this and above plot line to examine individual changes in the
    % frame data
%     subplot(212); plot(I,meanvidloc,'--',I,vidframeso,'m',I,vidframes); xlim([min(xs) max(xs)]); ylim([min([vidframes(xs) vidframes(xs)]) max([vidframes(xs) vidframes(xs)])]);
%     pause
    i = i+1;
end
    

% 
% 
%     if vidframes(drops(i)) > fun(drops(i)); % if the drop is after a section that is above the mean line
%         sectstart = find(diff(vidframes(1:drops(i)))>ldif*1.5,1,'last')+1; %theory is that any jump shoule be bigger than the drop
%         if isempty(sectstart); sectstart = 1; end
% %         if (sectstart == drops(i) || ldif < 1.5*meandif) && replsingles; 
% %             vidframes(drops(i)+1:drops(i)+2) = nan; vidframes = fixgaps(vidframes); 
% %             disp(['replaced a single bad frame (# ' num2str(sectstart) ') in video # ' num2str(i)]);
% %         else
%             vidframes(sectstart:drops(i)) = vidframes(sectstart:drops(i)) - ldif-fdif;
% %         end
%         badsect = badsect+(drops(i)-sectstart)+1;
%     else
%         sectend = find(diff(vidframes(drops(i)+1:end))>ldif*1.5,1,'first')+drops(i); %theory is that any jump shoule be bigger than the drop
%         if isempty(sectend); sectend = length(vidframes); end
% %         if (sectstart == drops(i) || ldif < 1.5*meandif) && replsingles;
% %             sectend == drops(i)+1 && replsingles; vidframes(sectend) = nan; vidframes = fixgaps(vidframes);
% %             
% %         else
%             vidframes(drops(i)+1:sectend) = vidframes(drops(i)+1:sectend) + ldif+fdif;
% %         end
%         badsect = badsect+sectend-drops(i);
%     end
% end
if ~insecs; vidframes = vidframes/24/60/60+vidstart-vidframes(vidstartI)/24/60/60; end


