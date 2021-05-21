function [At,Acal] = calA(data,DN,tagondec,ofs,accHz,df,CAL,Depth,I)

% this function performs an in situ calibration on the accelerometer data,
% based on the spherical_cal scripts at animaltags.org

if nargin < 9
    I = find(tagondec);
    % else, could limit I (e.g. I = [40000:60000 90000:100000];)
end

names =fieldnames(CAL);
for ii = 1:length(names);
    eval([names{ii} ' = CAL.' names{ii} ';']);
end


II = I; %find(tagondec);
t1 = find(tagondec,1); t2 = find(tagondec,1,'last');
At = [data.Acc1 data.Acc2 data.Acc3]; I = isnan(At); At = edgenans(At); %At = decdc(At,df); 
At = decimateM(At,ofs,accHz,df,length(DN),'accHz');
numrows = size(At,1);
I = I(1:df:end,:); if size(I,1)<size(At,1); I(end,:) = I(end-1,:); else I = I(1:size(At,1),:); end; At = (At-repmat(aconst,[numrows,1]))*acal; At(I) = nan;
At_mag=sqrt(sum(At.^2,2));
axA = (acal./abs(acal)); axA(isnan(axA)) = 0;
Att = [data.Acc1 data.Acc2 data.Acc3]*axA;
Att = decimateM(Att,ofs,accHz,df,length(DN),'accHz');
[~,Acalnew,~] = spherical_calwk(Att(II,:),1,'cross');
Atnew = (Att*diag(Acalnew.poly(:,1))+repmat(Acalnew.poly(:,2)',[size(Att,1),1]))*Acalnew.cross;
oi = norm2(Atnew); disp(std(oi)/mean(oi)); %clear Mtt;
if std(oi)/mean(oi)<0.000001;
    disp('Cross method could not converge, trying again with bounds, if no luck, try limiting II');
    [~,Acalnew,~] = spherical_calwk(Att(II,:),1,'cross',true);
    Atnew = (Att*diag(Acalnew.poly(:,1))+repmat(Acalnew.poly(:,2)',[size(Att,1),1]))*Acalnew.cross;
end
figure(6); clf; sp1 = subplot(4,1,1:2);
plot(t1:t2,At(t1:t2,:),t1:t2,Atnew(t1:t2,:)); ylim([-2 2]);
legend('X_{bench}','Y_{bench}','Z_{bench}','X_{sphere}','Y_{sphere}','Z_{sphere}');
sp2 = subplot(4,1,3);
plot(t1:t2,At_mag(t1:t2),t1:t2,sqrt(sum(Atnew(t1:t2,:).^2,2)),[t1 t2],[1 1]); legend('At','At_sphericalcal (from animaltags.org)','1');
ylim([0.5 1.5]);
disp(['At median: ' num2str(median(At_mag(t1:t2)))]);
disp(['At_sphericalcal median: ' num2str(median(sqrt(sum(Atnew(t1:t2,:).^2,2))))]);
sp3 = subplot(4,1,4);
plot(t1:t2,Depth(t1:t2)); set(gca,'ydir','rev','ylim',[-5 max(Depth)]); ylabel('Depth (m)');
linkaxes([sp1 sp2 sp3],'x');
z1 = zoom(sp1); z2 = zoom(sp2); z3 = zoom(sp3);% z4 = zoom(s4);
set([z1 z2 z3],'enable','on','Motion','both');
Mchoice = input('Bench cal = 1, spherical cal = 2. Which do you want to use? (Enter 1 or 2) ');
switch Mchoice
    case 1
        Acal = [];
    case 2
        At = Atnew;
        Acal = Acalnew;
    otherwise
        Acal = []; disp('No calibration chosen, so bench test was used by default');
end
         
