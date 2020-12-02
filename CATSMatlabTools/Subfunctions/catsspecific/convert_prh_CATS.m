function		ncfile = convert_prh_CATS(prhfile,infostruct,savedir)
% Converts CATS PRH file to NC format
% Saves output NC file in savedir, or pwd if not specified
%
% Danuta Wisniewska, updated by James Fahlbusch
% Version 2.21.18
% Goldbogen Lab
% Stanford University
%
% format:
% ncfile = convert_prh_CATS(prhfile,infostruct,savedir);
% ncfile = convert_prh_CATS(prhfile,infostruct);
%
% returns: name of ncfile created
% ncfile = convert_prh_CATS(...); 
%

VARS = {'At','Aw','Mt','Mw','Gt','Gw','p','T','Light','DN','pitch','roll','head', 'GPS', 'INFO.TempInternal', 'Temp1'} ;
vtype = {'acc','acc','mag','mag','gyr','gyr','press','temp','light','datenum','pitch','roll','head', 'pos','temp','temp'} ;
frm = [0 1 0 1 0 1 0 0 0 0 1 1 1 1 0 0] ;					% 0 = tag, 1 = animal
scf = [9.81,9.81,1,1,1,1,1,1,1,1,1,1,1,1,1,1] ;			% conversions to standard units

X=load(prhfile) ;
dbfile = dir(fullfile(prhfile));
prhCreation = datestr(dbfile.date,'dd-mmm-yyyy HH:MM:SS');
% info=csv2struct(csvfile);
info = infostruct;
utcOffset = info.dephist_utc2loc;
clearvars -except info X scf frm VARS vtype prhfile savedir utcOffset prhCreation

fn = fieldnames(X) ;
if ~isfield(X,'fs'),
    fprintf(' PRH file must contain a varable fs with the sampling rate\n');
    return
end
%check for internal temp in INFO and add to fieldnames (fn)
if isfield(X, 'INFO') && isfield(X.INFO, 'TempInternal')
    fn(end+1) = cellstr('INFO.TempInternal');
end

fs = X.fs(1) ;
ncfile = sprintf('%s_prh%d',info.depid,round(fs));
if nargin>2 
    save_nc([savedir ncfile],info) ;
else
    save_nc(ncfile,info) ;
end

for i=1:length(VARS)
    k = strcmpi(fn,VARS{i}) ;
    % 	if isempty(k), continue, end
    if sum(k)==0, continue, end
    
    if sum(k)~=0 && strcmp(fn{k},'DN')
        %convert DN to UTC
        V=sens_struct(X.(fn{k})*scf(i)-(utcOffset/24),fs,info.depid,vtype{i});
        V.description = 'DN is in UTC. The offset to local time is located in info.dephist_utc2loc';
    elseif sum(k)~=0 && strcmp(fn{k},'GPS')
        V=sens_struct(X.GPS(~isnan(X.GPS(:,1)),:),find(~isnan(X.GPS(:,1)))/fs,info.depid,vtype{i});
    elseif sum(k)~=0 && strcmp(fn{k},'pitch')
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},VARS{i});
        V.description = 'pitch (+ up, - down)';
    elseif sum(k)~=0 && strcmp(fn{k},'roll')
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},VARS{i});
        V.description = 'roll (+ right, - left)';
    elseif sum(k)~=0 && strcmp(fn{k},'head')
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},VARS{i});
        V.description = 'heading (+ right, - left, 0 True North)';
    elseif sum(k)~=0 && strcmp(fn{k},'T')
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},'Texternal');    
        V.description = 'External Temperature';
    elseif sum(k)~=0 && strcmp(fn{k},'INFO.TempInternal')
        fs = X.fs(1) ;
        V=sens_struct(X.INFO.TempInternal,fs,info.depid,vtype{i},'Tinternal');
        V.description = 'Internal Temperature';
    elseif sum(k)~=0 && strcmp(fn{k},'Temp1')
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},'Tinternal');
        V.description = 'Internal Temperature';        
    elseif sum(k)~=0 && strcmp(fn{k},'p')
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i});        
    elseif sum(k)==0
        continue
    elseif sum(k)~=0 && (strcmp(fn{k},'Aw') || strcmp(fn{k},'At') || strcmp(fn{k},'Mw') || strcmp(fn{k},'Mt') || strcmp(fn{k},'Gw') || strcmp(fn{k},'Gt'))
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},VARS{i});
        V.axes = 'FRD (NED)';
    else
        fs = X.fs(1) ;
        V=sens_struct(X.(fn{k})*scf(i),fs,info.depid,vtype{i},VARS{i});
    end

    if frm(i) == 1
        V.frame = 'animal' ;
    else
        V.frame = 'tag' ;
    end
    V.creation_date = prhCreation; %sets date to the creation date of the prh file
    V.author = 'Dave Cade, davecade@stanford.edu';
    if nargin>2 
        add_nc(fullfile(savedir,ncfile),V) ;
    else
        add_nc(ncfile,V) ;
    end
end
