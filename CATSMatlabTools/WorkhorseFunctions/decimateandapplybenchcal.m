function [fs,Mt,At,Gt,DN,Temp,Light,LightIR,TempInternal,tagondec,camondec,audondec,tagslipdec] = decimateandapplybenchcal(data,Depth,CAL,ofs,DN,df,Hzs,tagon,camon,audon,tagslip)

% this function facilitates stopping and starting of the prh process by
% applying in situ cals and decimation factors to the variables in "data"

nout = length(Depth);

try Temp = (data.Temp-CAL.Tconst)*CAL.Tcal; catch; Temp = data.Temp; end
Temp = decimateM(Temp,ofs,Hzs.THz,df,nout,'THz');

try try Light = decimateM(data.Light,ofs,Hzs.lHz,df,nout,'lHz'); %decdc(data.Light,df); 
        LightIR = nan(size(Light)); 
    catch 
        Light = decimateM(data.Light1,ofs,Hzs.lHz,df,nout,'lHz');
        LightIR = decimateM(data.Light2,ofs,Hzs.lHz,df,nout,'lHz');
%         Light = decdc(data.Light1,df); LightIR = decdc(data.Light2,df); 
    end % IR = infrared
catch; Light = nan(size(Temp)); LightIR = nan(size(Light));
end

try 
    TempInternal = decimateM(data.Temp1,ofs,Hzs.T1Hz,df,nout);
catch
    try
    TempInternal = decimateM(data.TempDepthInternal,ofs,Hzs.TDIHz,df,nout);
    warning('No Temp1, used TempDepthInternal');
    catch
        TempInternal = nan(size(Temp));
        warning('No Temp1, TempInternal is nans');
    end
end
% Temp1 = decdc(data.Temp1,df);
% DV = DV(1:df:end,:); DV = DV(1:length(Temp),:);

fs = ofs/df;
names =fieldnames(CAL);
for ii = 1:length(names);
    eval([names{ii} ' = CAL.' names{ii} ';']);
end
names =fieldnames(Hzs);
for ii = 1:length(names);
    eval([names{ii} ' = Hzs.' names{ii} ';']);
end
numrows = nout; %size(Depth,1);
% Mt = (decdc(filterCATS([data.Comp1 data.Comp2 data.Comp3],ceil(fs/8),round(fs),.05),df)-repmat(magconston,numrows,1))*magcalon; %t indicates tag recorded values
% Mtoff = (decdc(filterCATS([data.Comp1 data.Comp2 data.Comp3],ceil(fs/8),round(fs),.05),df)-repmat(magconstoff,numrows,1))*magcaloff;
% if tagnum == 36 %|| tagnum == 6
%     Mt = filterCATS([data.Comp1 data.Comp2 data.Comp3],ceil(ofs/8),round(ofs),.05);
%     Mt = filterMag(Mt,ofs,tagon);
% else
try
    Mt = filterCATS([data.Comp1 data.Comp2 data.Comp3],ceil(ofs/8),round(ofs),.05); 
catch; Mt = nan(size(Depth,1),3);
    warning ('No magnetometer detected, Mt is nans');
end
% end
I = isnan(Mt);
if sum(I(:,1)) ~=size(I,1)
    Mt = edgenans(Mt); MtNoNan= Mt;
    % Mt = decdc(filterCATS(Mt,ceil(fs/8),round(fs),.05),df);
    Mt = decimateM(Mt,ofs,magHz,df,nout,'magHz'); MtO = Mt;
    % Mt = decdc(Mt,df); MtO = Mt;
    % Mt(188200:end,:) = nan;
    if ~exist('magconstoff','var');
        try magconstoff = magconston; magcaloff = magcalon; catch
            try magconstoff = magconst; magconston = magconst; magcalon = magcal; magcaloff = magcal;
            catch; warning('no magnetometer bench cal, to continue press enter'); pause;
                magconston = [0 0 0]; magcalon = diag(ones(3,1));
                magconstoff = magconston; magcaloff = magcalon;
            end
        end
    end
    
    
    I = I(1:df:end,:); if size(I,1)<size(Mt,1); I(end,:) = I(end-1,:); else I = I(1:size(Mt,1),:); end
    Mtoff = Mt;
    Mtoff = (Mtoff-repmat(magconstoff,numrows,1))*magcaloff;
    Mt = (Mt-repmat(magconston,numrows,1))*magcalon;
    Mt(I) = nan; Mtoff(I) = nan;
end
camondec = camon(1:df:end); camondec = camondec(1:size(Mt,1));
camondec = interp2length(camondec,ofs/df,ofs/df,nout);
audondec = audon(1:df:end); audondec = audondec(1:size(Mt,1));
audondec=interp2length(audondec,ofs/df,ofs/df,nout);
% camoff = ~camon;
%camoffdec = camoff(1:df:end); camoffdec = camoffdec(1:size(Mt,1));
% camoffdec = ~camondec;
tagondec = tagon(1:df:end); tagondec = tagondec(1:size(Mt,1));
tagondec = interp2length(tagondec,ofs/df,ofs/df,nout);
Mt(~camondec,:) = Mtoff(~camondec,:);
tagslipdec = ceil(tagslip/df);

% Gt = (decdc(filterCATS([data.Gyr1 data.Gyr2 data.Gyr3],ceil(fs/8),round(fs),.05),df)-repmat(gyconst,numrows,1))*gycal;
% At = (decdc(filterCATS([data.Acc1 data.Acc2 data.Acc3],ceil(fs/8),round(fs),.05),df)-repmat(aconst,numrows,1))*acal;
try
Gt = [data.Gyr1 data.Gyr2 data.Gyr3];
catch
    Gt = nan(size(Depth,1),3);
    warning('no gyros detected, Gt is nans');
     gyconst = [0 0 0]; gycal = diag(ones(3,1));
end
I = isnan(Gt);
if sum(I(:,1)) ~=size(I,1)
    Gt = edgenans(Gt);
    % Gt = decdc(filterCATS(Gt,ceil(fs/8),round(fs),.05),df);
    Gt = decimateM(Gt,ofs,gyrHz,df,nout,'gyrHz');
    % Gt = decdc(Gt,df);
    I = I(1:df:end,:); if size(I,1)<size(Gt,1); I(end,:) = I(end-1,:); else I = I(1:size(Gt,1),:); end
    Gt = (Gt-repmat(gyconst,numrows,1))*gycal;
    Gt(I) = nan;
end
At = [data.Acc1 data.Acc2 data.Acc3];
I = isnan(At);
At = edgenans(At);
% At = decdc(filterCATS(At,ceil(fs/8),round(fs),.05),df);
At = decimateM(At,ofs,accHz,df,nout,'accHz');

I = I(1:df:end,:); if size(I,1)<size(At,1); I(end,:) = I(end-1,:); else I = I(1:size(At,1),:); end
At = (At-repmat(aconst,numrows,1))*acal;
At(I) = nan;



fs = round(1./mean((DN(50:60)-DN(49:59))*24*60*60));
if abs(round(fs)-fs)<.01; fs = round(fs); end
% try TempI = decdc(data.Temp1,df); catch; end;
% if you have a speed sensor, you could use a variation of this script to
% include speed.
% try
%     ds = diff(data.Speed);
%     oi = diff(find(ds~=0));
%     speedfs = fs*df/min(oi(5:end)); %get rid of the first few values just to be safe;
%     Paddles = decdc(data.Speed*speedfs,df); % paddles is then kind of a smoothed speed in impulses/sec (one impulse is 1/4 rotation of the wheel)
% catch err
%     data.Speed = zeros(size(data.Temp));
%     speedfs = fs;
%     Paddles = decdc(data.Speed*speedfs,df);
% end

% endRow = length(Temp)+1;
% time1 = 0:length(Temp)-1;%data(:,1);
% Times = (DN-DN(1))*24*60*60; % time in seconds from start%data(:,2);


figure(3); clf;
s1 = subplot(3,1,1);
Mt_mag = sqrt(sum(Mt.^2,2));
plot(Mt); hold on; plot(Mt_mag,'k');
legend('X','Y','Z','|Mt|'); title('Mt bench cal'); try ylim([min(Mt(:)) max([Mt(:); Mt_mag])]); catch; end
s2 = subplot(3,1,2);
At_mag = sqrt(sum(At.^2,2));
plot(At); hold on; plot(At_mag,'k');
legend('X','Y','Z','|At|'); title('At bench cal'); ylim([min(At(:)) max([At(:); At_mag])]);
s3 = subplot(3,1,3);
% Mt_mag = sqrt(sum(Mt.^2,2));
ax = plotyy(1:length(Depth),Depth,1:length(Depth),Temp); hold on; %plot(Temp,'k');
set(ax(1),'ydir','rev');
legend('Depth','Temperature'); title('Temp bench cal'); %ylim([min(Mt(:)) max([Mt(:); Mt_mag])]);
linkaxes([s1 s2 s3],'x');
