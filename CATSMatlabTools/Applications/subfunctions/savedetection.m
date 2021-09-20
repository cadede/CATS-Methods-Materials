function [] = savedetection(DN,fs,L,LI,p,LC,i,creator,primary_cue,notes,fileloc,filename,StI,StEI,StC,StEC)

d1 = datevec(now());
created_on = [d1(2) d1(3) d1(1)];
clearvars d1;
% store temp variables to lunge file
starttime = DN(1);
prh_fs = fs;
LungeDN = L;
depth = p(LI);
%time = L;
LungeI = LI;
LungeC = LC;
LungeDepth = depth;
progressIndex = i;
save([fileloc getWhaleID(filename) 'lunges.mat'],'LungeDN','LungeI','LungeDepth','LungeC','creator','primary_cue','prh_fs','starttime','created_on', 'progressIndex', 'notes');

% Strategy Save Information
StrategyS = StI;
StrategyE = StEI;
%clear StrategyC StrategySC StrategyEC % clears old variables
if length(StC) == length(StEC)
    StrategyC = StC;
    save([fileloc getWhaleID(filename) 'strategy.mat'],'StrategyS','StrategyE','StrategyC','creator','prh_fs','starttime','created_on', 'progressIndex', 'notes');
    disp('Strategy file saved')
else
    StrategySC = StC;
    StrategyEC = StEC;
    save([fileloc getWhaleID(filename) 'strategy.mat'],'StrategyS','StrategyE','StrategySC','StrategyEC','creator','prh_fs','starttime','created_on', 'progressIndex', 'notes');
    disp('Strategy file saved, however, strategy marks are uneven (Different Number of Start and End Points). Check for un-deleted strategy marker.');
end


try aa = strfind(fileloc,'CATS\tag_data\');
    save([fileloc(1:aa+13) 'prhlunges\' getWhaleID(filename) 'lunges.mat'],'LungeDN','LungeI','LungeDepth','LungeC','creator','primary_cue','prh_fs','starttime','created_on', 'progressIndex', 'notes');
    disp('Lunge file saved in fileloc and prhlunges folder');
catch; disp('Lunge file saved, but prhlunges folder not found');
end

end