function prh2Acq (fl,fn)
% column order:
% Depth Aw_x Aw_y Aw_z pitch(deg) roll(deg) heading(deg) speedJJ speedFN Jerk MSA
%%

if exist('fl','var')&&exist('fn','var'); fileloc = fl; filename = fn; else
    cf = pwd;
    [filename,fileloc]=uigetfile('*.mat', 'select prh files to convert to Acq','multiselect','on'); %look below to save time to make truncated files okay
    cd(cf);
end
if ~iscell(filename); filename = {filename}; end

for n = 1:length(filename)
   clearvars -except fileloc filename n
   load([fileloc filename{n}]);
   % in 2020, I have to load pitch separately, apparently.
   pitch = load([fileloc filename{n}],'pitch');
   pitch = pitch.pitch;
   try S = speed.JJ; S2 = speed.FN; catch; try S = speedFN; S2 = S; disp('No Jiggle Speed Found'); catch; S = nan(size(p)); S2 = S; disp('No Speed Found'); end; end
   S(isnan(S)) = 1;
   S = runmean(S,round(fs/2));
   njerk = (9.81*fs)*sqrt(diff(Aw).^2*ones(3,1)) ; njerk(end+1) = njerk(end);
%    Sa = 
%    ODBA = 9.81/fs*sum(sum(abs(Sa)));
   MSA = sqrt(sum(Aw.^2,2))-1;
   data = [p Aw [pitch roll head]*180/pi S S2 njerk MSA];
   try data = data(find(tagon,1,'first'):find(tagon,1,'last'),:); catch; end
   try Time = datestr(DN(find(tagon,1,'first')),'HHMMSS'); catch; Time = ''; end
   dlmwrite([fileloc filename{n}(1:end-7) 'ACQ_' Time '.txt'],data,'\t');
end