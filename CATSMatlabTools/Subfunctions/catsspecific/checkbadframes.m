
function [vidframes,numbad,badsect] = checkbadframes(vidframes,insecs,replsingles,vidnum)
if nargin<2; insecs = false; end
if nargin<3; replsingles = false; end
if ~insecs
    vidstart = vidframes(find(~isnan(vidframes),1));
    vidstartI = find(~isnan(vidframes),1);
    vidframes = (vidframes - vidstart)*24*60*60; % get to units of seconds
end
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
meandif = mean(dV(abs(dV)<4*fdif));

% find a line that goes closest to the most points
optfun = @(b) sum(abs(meandif*(1:length(vidframes))+b-vidframes)>1.1*meandif);
b = -100:meandif:100; v = nan(size(b));
for i = 1:length(b);
    v(i) = optfun(b(i)); 
end
[~,i] = min(v);
b = b(i);
fun = @(t) meandif*t+b; % best fit function
badsect = 0;
drops = find(diff(vidframes)<0);
for i = 1:length(drops)
    ldif = abs(vidframes(drops(i))-vidframes(drops(i)+1)); 
    if vidframes(drops(i)+2)>vidframes(drops(i)) && replsingles
        vidframes(drops(i)+1:drops(i)+2) = nan; vidframes = fixgaps(vidframes); 
        disp(['replaced a single bad frame (# ' num2str(drops(i)+1) ') in video # ' num2str(vidnum)]);
        continue;
    end
    if vidframes(drops(i)) > fun(drops(i)); % if the drop is after a section that is above the mean line
        sectstart = find(diff(vidframes(1:drops(i)))>ldif*1.5,1,'last')+1; %theory is that any jump shoule be bigger than the drop
        if isempty(sectstart); sectstart = 1; end
%         if (sectstart == drops(i) || ldif < 1.5*meandif) && replsingles; 
%             vidframes(drops(i)+1:drops(i)+2) = nan; vidframes = fixgaps(vidframes); 
%             disp(['replaced a single bad frame (# ' num2str(sectstart) ') in video # ' num2str(i)]);
%         else
            vidframes(sectstart:drops(i)) = vidframes(sectstart:drops(i)) - ldif-fdif;
%         end
        badsect = badsect+(drops(i)-sectstart)+1;
    else
        sectend = find(diff(vidframes(drops(i)+1:end))>ldif*1.5,1,'first')+drops(i); %theory is that any jump shoule be bigger than the drop
        if isempty(sectend); sectend = length(vidframes); end
%         if (sectstart == drops(i) || ldif < 1.5*meandif) && replsingles;
%             sectend == drops(i)+1 && replsingles; vidframes(sectend) = nan; vidframes = fixgaps(vidframes);
%             
%         else
            vidframes(drops(i)+1:sectend) = vidframes(drops(i)+1:sectend) + ldif+fdif;
%         end
        badsect = badsect+sectend-drops(i);
    end
end
if ~insecs; vidframes = vidframes/24/60/60+vidstart-vidframes(vidstartI)/24/60/60; end


