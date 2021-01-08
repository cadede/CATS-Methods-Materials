function [Depth,At,Mt,Gt,Temp,Temp1,Light,LightIR] = applyCal2(data,DN,CAL,camondec,ofs,Hzs,df)
nout = length(DN);

if nargin<5; df = 1; end
% if you don't need to decimate, assume df = 1;

names =fieldnames(CAL);
for ii = 1:length(names);
    eval([names{ii} ' = CAL.' names{ii} ';']);
end

try pressTemp = data.Temp1; catch; pressTemp = data.Temp; end
try Depth = decimateM((data.Pressure-CAL.pconst)*CAL.pcal+polyval([CAL.pc.tcomp,CAL.pc.poly(2)],pressTemp-CAL.pc.tref),ofs,Hzs.pHz,df,nout,'pHz');
catch; Depth = decimateM((data.Pressure-CAL.pconst)*CAL.pcal,ofs,Hzs.pHz,df,nout,'pHz');
    disp('Bench pressure cal applied');
end
try Temp = (data.Temp-CAL.Tconst)*CAL.Tcal; catch; Temp = data.Temp; end
Temp = decimateM(Temp,ofs,Hzs.THz,df,nout,'THz');
try Temp1 = decimateM(data.Temp1,ofs,Hzs.T1Hz,df,nout,'THz'); catch; Temp1 = nan(size(Temp)); end
try try Light = decimateM(data.Light,ofs,Hzs.lHz,df,nout,'lHz'); %decdc(data.Light,df);
        LightIR = nan(size(Light));
    catch
        Light = decimateM(data.Light1,ofs,Hzs.lHz,df,nout,'lHz');
        LightIR = decimateM(data.Light2,ofs,Hzs.lHz,df,nout,'lHz');
        %         Light = decdc(data.Light1,df); LightIR = decdc(data.Light2,df);
    end % IR = infrared
catch; Light = nan(size(Temp)); LightIR = nan(size(Light));
end
At = decimateM([data.Acc1 data.Acc2 data.Acc3],ofs,Hzs.accHz,df,nout,'accHz');
try Mt = decimateM([data.Comp1 data.Comp2 data.Comp3],ofs,Hzs.magHz,df,nout,'magHz');
catch; Mt = nan(size(Depth,1),3);
    warning ('No magnetometer detected, Mt is nans');
end
try Gt = decimateM([data.Gyr1 data.Gyr2 data.Gyr3],ofs,Hzs.gyrHz,df,nout,'gyrHz');
catch
    Gt = nan(size(Mt));
    warning('no gyros detected, Gt is nans');
    gyconst = zeros(1,size(Gt,2)); gycal = diag(ones(size(Gt,2),1));
end
numrows = size(Gt,1);

Gt = (Gt-repmat(gyconst,numrows,1))*gycal;
if ~exist('Acal','var') || isempty(Acal)
    At = (At-repmat(aconst,numrows,1))*acal;
else
    axA = (acal./abs(acal)); axA(isnan(axA)) = 0;
    At = At*axA;
    At = (At*diag(Acal.poly(:,1))+repmat(Acal.poly(:,2)',[size(At,1),1]))*Acal.cross;
end
if ~exist('magconstoff','var'); 
    try magconstoff = magconston; magcaloff = magcalon; catch
        try magconstoff = magconst; magconston = magconst; magcalon = magcal; magcaloff = magcal;
        catch; warning('no magnetometer bench cal, to continue press enter'); pause;
            magconston = [0 0 0]; magcalon = diag(ones(3,1));
            magconstoff = magconston; magcaloff = magcalon;
        end
    end
end


if ~exist('Mcal','var') || isempty(Mcal)
    Mt(camondec,:) = (Mt(camondec,:)-repmat(magconston,sum(camondec),1))*magcalon;
    Mt(~camondec,:) = (Mt(~camondec,:)-repmat(magconstoff,sum(~camondec),1))*magcaloff;
else
    axM = (magcalon./abs(magcalon)); axM(isnan(axM)) = 0;
    Mt = Mt*axM;
    if any(strcmp('on',fieldnames(Mcal)));
        Mt1 = (Mt*diag(Mcal.on.poly(:,1))+repmat(Mcal.on.poly(:,2)',[size(Mt,1),1]))*Mcal.on.cross;
        Mt2 = (Mt*diag(Mcal.off.poly(:,1))+repmat(Mcal.off.poly(:,2)',[size(Mt,1),1]))*Mcal.off.cross;
        Mt(Mcal.camoffNTdec,:) = Mt2(Mcal.camoffNTdec,:);
        Mt(Mcal.camonNTdec,:) = Mt1(Mcal.camonNTdec,:);
        Mt(~(Mcal.camonNTdec|Mcal.camoffNTdec),:) = nan;
        Mt = fixgaps(Mt);
    else
        Mt = (Mt*diag(Mcal.poly(:,1))+repmat(Mcal.poly(:,2)',[size(Mt,1),1]))*Mcal.cross;
    end
end
