function [Mt,Mcalnew] = calM(data,DN,tagondec,camondec,camon,nocam,ofs,magHz,df,CAL,Temp,b,I,resThresh)


% this function performs an in situ calibration on the magnetometer data,
% based on the spherical_cal scripts at animaltags.org

nout = length(DN);
if nargin < 13 || isempty(I)
    I = find(tagondec);
    % else, could limit I (e.g. I = [40000:60000 90000:100000];)
end
if nargin<14
    resThresh = 0.05;
end

names =fieldnames(CAL);
for ii = 1:length(names);
    eval([names{ii} ' = CAL.' names{ii} ';']);
end
t1 = find(tagondec,1); t2 = find(tagondec,1,'last');

axM = (magcalon./abs(magcalon)); axM(isnan(axM)) = 0;
Mt = fixgaps([data.Comp1 data.Comp2 data.Comp3])*axM; %applyMagcalfirst = false;
Mt = decimateM(Mt,ofs,magHz,df,length(DN),'magHz');
Mcalnew00.poly = [1 0; 1 0; 1 0;]; Mcalnew00.cross = diag([1 1 1]);
% if you don't want to start with a previous calibration, comment out the next line
if exist('Mcal3d0','var'); Mcalnew00 = Mcal3d0; else Mcalnew00.poly = [ones(3,1) (-magconstoff*axM)']; Mcalnew00.cross = axM^-1 * magcaloff; end 

Mt = (Mt*diag(Mcalnew00.poly(:,1))+repmat(Mcalnew00.poly(:,2)',size(Mt,1),1))*Mcalnew00.cross;  % apply old cal first or no cal

% Mt = (decdc([data.Comp1 data.Comp2 data.Comp3],df)-repmat(magconstoffnew,numrows,1))*magcaloffnew; applyMagcalfirst = true;

% I = tagondec;
% if there is trouble converging, choose a limited value for I and uncomment the next line
% I = [93000:98000];%[125000:130000 1210000:1230000];%[]% 187500:192500];%600000:620000 1000000:1020000];% 404000:405000 446000:446200 820000:802000 805000:806000];%1000000:1030000];% 260000:270000];

if isnan(b); error('may be problem with GPS- b does not exist'); end
[Mtt, Mcalnew,output] = spherical_calwk(Mt(I,:),b,'cross');
Mtnew = (Mt*diag(Mcalnew.poly(:,1))+repmat(Mcalnew.poly(:,2)',size(Mt,1),1))*Mcalnew.cross;
oi = norm2(Mtnew); disp(std(oi)/mean(oi)); clear Mtt;
if std(oi)/mean(oi)<0.000001;
    disp('Cross method could not converge, trying again with a bound, if no luck, try limiting I');
    [Mtt, Mcalnew,output] = spherical_calwk(Mt(I,:),b,'cross',true);
end
Mtnew = (Mt*diag(Mcalnew.poly(:,1))+repmat(Mcalnew.poly(:,2)',size(Mt,1),1))*Mcalnew.cross;
Mcal2 = Mcalnew00;
D = Mcalnew00.cross*diag(Mcalnew.poly(:,1));
D2 = Mcalnew.poly(:,2)'*D^-1;
Mcal2.poly(:,2) = Mcalnew00.poly(:,2)+D2';
Mcal2.cross = D*Mcalnew.cross;
Mcalnew = Mcal2;
oi = norm2(Mtnew(tagondec,:)); disp(std(oi)/mean(oi)); clear Mtt;
axB = 100/cond(Mtnew(tagondec,:)'*Mtnew(tagondec,:));

disp('You want residual < 5% and axial balance > 20%');    
figure(3); clf;
s1 = subplot(3,1,1);
Mt_mag = sqrt(sum(Mt.^2,2));
plot(Mt); hold on; plot(Mt_mag,'k');
legend('X','Y','Z','|Mt|'); title('Mt bench spherical cal'); ylim([min(Mt(:)) max([Mt(:); Mt_mag])]);

try delete(s2); catch; end; try delete(s3); catch; end
s2 = subplot(3,1,2);
Mtnew_mag = sqrt(sum(Mtnew.^2,2));
plot(Mtnew);  hold on; plot(Mtnew_mag,'k');
 oi2 = double(camondec); oi2(~camondec) = nan;
  plot(s2,oi2,'m','linewidth',2);
legend('X','Y','Z','Mag','Cam on'); title('Mt Calibrated'); ylim([-max(abs(Mt(:))) max([Mt(:); Mt_mag])]);
linkaxes([s1 s2],'xy')
MtO = Mtnew;
McalnewO = Mcalnew;
%
% if std(oi)/mean(oi)>resThresh || isnan(std(oi)/mean(oi))
%     disp('High Residual, trying with a Temperature calibration')
%     MtO = Mtnew;
%     McalnewO = Mcalnew;
%     [Mtt, McalnewT] = spherical_cal_T(Mt(I,:),b,Temp(I),'cross');
%     Mtnew = applycalT(Mt,Temp,McalnewT);
%     oi = norm2(Mtnew); 
%     if std(oi)/mean(oi)<0.000001;
%         disp('Cross method could not converge, trying again without cross method');
%         [Mtt, McalnewT] = spherical_cal_T(Mt(I,:),b,Temp(I));
%     end
%     Mtnew = applycalT(Mt,Temp,McalnewT);
%     McalnewT.McalOrig = Mcalnew00;
%     Mcalnew = McalnewT;
%         
%     figure(3);
%     s3 = subplot(3,1,3);
%     Mtnew_mag = sqrt(sum(Mtnew.^2,2));
%     plot(Mtnew);  hold on; plot(Mtnew_mag,'k');
%     legend('X','Y','Z'); title('Mt Calibrated'); ylim([min(Mt(:)) max([Mt(:); Mt_mag])]);
%     linkaxes([s1 s2 s3],'xy');
%     oi = norm2(Mtnew); disp(std(oi)/mean(oi)); clear Mtt;
% end
%
% I = sort([find(camondec&tagondec)' 400000:500000]);
axB = 100/cond(Mtnew(tagondec,:)'*Mtnew(tagondec,:));
if (std(oi)/mean(oi)>resThresh || isnan(std(oi)/mean(oi)) || axB<10) && ~nocam 
    disp('High Residual or poor balance, trying with different cam on and cam offs')
    MtT = Mtnew;
    try delete(s4); catch; end
    try CAM = decdc(data.Camera,df);  disp('Assuming cam on indicator value is 14'); 
        CAM = interp2length(CAM,ofs/df,ofs/df,nout);
        isCAM = logical(floor(runmean(CAM>=11.5 & CAM <= 24.5,1)));
        TRANS = ~isCAM&CAM>0;
        camonNTdec = isCAM;
        camoffNTdec = CAM<=0;
    catch
        disp('No Camera variable');
        I2 = I;
        I= find(diff(camon) == 1); camonNoTrans = camon; % get rid of 13 seconds of transition while the camera turns on/off, three seconds from the camoff time
        for i = 1:length(I);  camonNoTrans(I(i):round(I(i)+10*ofs))= false; end; if length(camonNoTrans)>length(camon); camonNoTrans(length(camon)+1:end) = []; end
        I= find(diff(camon) == -1); for i = 1:length(I);  camonNoTrans(round(max(I(i)-10*ofs,1)):I(i))= false; end
        camoff = ~camon;
        I= find(diff(camoff) == 1); camoffNoTrans = camoff;
        for i = 1:length(I);  camoffNoTrans(I(i):round(I(i)+3*ofs))= false; end
        I= find(diff(camoff) == -1); for i = 1:length(I);  camoffNoTrans(round(max(I(i)-3*ofs,1)):I(i))= false; end
        camonNTdec = camonNoTrans(1:df:end); camonNTdec = camonNTdec(1:size(Mt,1));
        camoffNTdec = camoffNoTrans(1:df:end); camoffNTdec = camoffNTdec(1:size(Mt,1));
        CAM = camondec; isCAM = camondec;
        TRANS = ~camonNTdec&~camoffNTdec & tagondec;
        I = I2;
    end
    isI = false(size(tagondec)); isI(I) = true;
    disp('CAM ON:')
    [~, Mcalnewon] = spherical_calwk(Mt(tagondec&camonNTdec&isI,:),b,'cross');
    disp('CAM OFF:')
    [~, Mcalnewoff] = spherical_calwk(Mt(tagondec&camoffNTdec&isI,:),b,'cross');
    Mtnewon = (Mt*diag(Mcalnewon.poly(:,1))+repmat(Mcalnewon.poly(:,2)',size(Mt,1),1))*Mcalnewon.cross;
    Mtnewoff = (Mt*diag(Mcalnewoff.poly(:,1))+repmat(Mcalnewoff.poly(:,2)',size(Mt,1),1))*Mcalnewoff.cross;
    oi = norm2(Mtnewon);
    if std(oi)/mean(oi) < 0.000001;
        disp('Cross method for cam on could not converge, trying again with bounded values');
        disp('CAM ON:')
        [~, Mcalnewon] = spherical_calwk(Mt(tagondec&camonNTdec&isI,:),b,'cross',true);
        Mtnewon = (Mt*diag(Mcalnewon.poly(:,1))+repmat(Mcalnewon.poly(:,2)',size(Mt,1),1))*Mcalnewon.cross;
    end
    oi = norm2(Mtnewoff);
    if std(oi)/mean(oi) < 0.000001;
        disp('Cross method for cam off could not converge, trying again with bounded values');
        disp('CAM OFF:')
        [~, Mcalnewoff] = spherical_calwk(Mt(tagondec&camoffNTdec&isI,:),b,'cross',true);
        Mtnewoff = (Mt*diag(Mcalnewoff.poly(:,1))+repmat(Mcalnewoff.poly(:,2)',size(Mt,1),1))*Mcalnewoff.cross;
    end
    Mcal2 = Mcalnew00;
    D = Mcalnew00.cross*diag(Mcalnewon.poly(:,1));
    D2 = Mcalnewon.poly(:,2)'*D^-1;
    Mcal2.poly(:,2) = Mcalnew00.poly(:,2)+D2';
    Mcal2.cross = D*Mcalnewon.cross;
    Mcalnewon = Mcal2;
    
    Mcal2 = Mcalnew00;
    D = Mcalnew00.cross*diag(Mcalnewoff.poly(:,1));
    D2 = Mcalnewoff.poly(:,2)'*D^-1;
    Mcal2.poly(:,2) = Mcalnew00.poly(:,2)+D2';
    Mcal2.cross = D*Mcalnewoff.cross;
    Mcalnewoff = Mcal2;
    
        
    Mtnewoff(camonNTdec,:) = Mtnewon(camonNTdec,:);
    Mtnewoff( ~(camonNTdec | camoffNTdec),:) = nan;
    Mtnew = Mtnewoff;
%     sum(isnan(Mtnew))
%     figure(3);
%     s2 = subplot(9,1,4:5);
%     Mtnew_mag = sqrt(sum(MtO.^2,2));
%     plot(MtO);  hold on; plot(Mtnew_mag,'k');
%     legend('X','Y','Z'); title('Mt Calibrated Original'); ylim([min(Mt(:)) max([Mt(:); Mt_mag])]);
%     
%     s3 = subplot(9,1,6:7);
%     Mtnew_mag = sqrt(sum(MtT.^2,2));
%     plot(MtT);  hold on; plot(Mtnew_mag,'k');
%     oi= double(TRANS); oi(~TRANS) = nan;
%     plot(s3,isCAM*15-45,'g--','linewidth',2);
%     plot(s3,oi*15-45,'r--','linewidth',2);
%     legend('X','Y','Z','Cam ramping'); title('Mt Calibrated Temp'); ylim([min(Mt(:)) max([Mt(:); Mt_mag])]);
    
%     s4 = subplot(9,1,8:9);
    s4 = subplot(3,1,3);
    plot(Mtnew); hold on;
    legend('X','Y','Z'); title('Mt Calibrated Cam on/off different'); ylim([min(Mt(:)) max([Mt(:); Mt_mag])]);
    oi = double(camondec); oi(~camondec) = nan;
    plot(s2,oi,'m','linewidth',2);   plot(s4,oi,'m','linewidth',2); %plot(s3,oi,'m','linewidth',2);
     linkaxes([s1 s2 s4],'xy');
     % fill in gaps with values from Mt (to keep the curves accurate).
%      Mnan = isnan(Mt(:,1));
     sI = find(diff(TRANS) == 1); %find the transitions 
     eI = find(diff(TRANS) == -1);
     if length(sI)~=length(eI); error('TRANS has a funny length'); end
     Mjumps = percentile(abs(diff(Mt)),.97); %the big jumps are spikes from magnetometer, so just even them out.
     for i = 1:length(sI);
         Msmall = Mt(sI(i)-1:eI(i)+1,:);
         Mstart = Mtnew(sI(i)-1,:);
         Mend = Mtnew(eI(i)+1,:);
         j = find(any(abs(diff(Msmall))>repmat(Mjumps,size(Msmall,1)-1,1),2));
         for ii = 1:length(j) % get rid of jumps;
            jd = diff(Msmall(j(ii):j(ii)+1,:));
            Msmall(j(ii)+1:end,:) = Msmall(j(ii)+1:end,:)-repmat(jd,size(Msmall(j(ii)+1:end,:),1),1);
         end
         Msmall = Msmall+repmat(diff([Msmall(1,:); Mstart]),size(Msmall,1),1); % starts Msmall at Mstart;
         Md = cumsum([zeros(1,3); repmat((Mend - Msmall(end,:))/(size(Msmall,1)-1),size(Msmall,1)-1,1)]); %difference for what it should be.
         Msmall = Msmall + Md;
         if sum(abs(Mend-Msmall(end,:))) > 0.01 || sum(abs(Mstart-Msmall(1,:))) > 0.01; error('Error'); end
         Mtnew(sI(i):eI(i),:) = Msmall(2:end-1,:);
         plot(s4,sI(i):eI(i),Mtnew(sI(i):eI(i),:),'--');
            
     end
     Mtnew_mag = sqrt(sum(Mtnew.^2,2));
     hold on; plot(s4,Mtnew_mag,'k');
     oi = norm2(Mtnew); disp(['Total Residual: ' num2str(std(oi)/mean(oi))]);
     
     xlim([t1 t2])
      z1 = zoom(s1); z2 = zoom(s2);  z4 = zoom(s4); %z3 = zoom(s3);
      set([z1 z2 z4],'enable','on','Motion','both');
    Mchoice = input('Which axis do you want to use? (Enter 1-3) ');

     switch Mchoice
         case 1;
         disp('Rerun section 6 above to get Mt');
         case 2;
             disp('Used Whole Data set Calibration');
             Mt = MtO;
             Mcalnew = McalnewO;
%          case 4;
%              disp('Used Temperature Calibration');
%              Mt = MtT;
%              Mcalnew = McalnewT;
         case 3;
             disp('Used different cals for cam off and cam on');
             Mt = Mtnew;
             Mcalnew = struct();
             Mcalnew.on = Mcalnewon;
             Mcalnew.off = Mcalnewoff;
             Mcalnew.camonNTdec = camonNTdec;
             Mcalnew.camoffNTdec = camoffNTdec;
             if any(any(isnan(Mt)))
                 disp(['Fixing ' num2str(sum(sum(isnan(Mt)))/3) ' nans']);
                 Mt = fixgaps(Mt);
             end
     end
else
    Mt = Mtnew;
end
xlim([t1 t2])