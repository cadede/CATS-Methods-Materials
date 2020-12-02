function [ID,tagtype,species,datenumber] = getWhaleID(STR)
%STR can be a string or a cell array of strings
if ~iscell(STR)    
    str = STR;
    str = char(str);
    ID = str(1:min([regexp(str,' '), regexp(str,'prh'),length(str)-3])-1);
    if ~isempty(strfind(ID,'PAD')); tagtype = 'PAD';
    elseif ~isempty(strfind(ID,'_')); tagtype = 'DTAG';
    else tagtype = 'CATS';
    end
    sp = ID(1:2);
    if strcmp(tagtype,'DTAG')
        DN = datenum([2000+str2num(ID(3:4)) 1 str2num(ID(6:8)) 0 0 0]);
    else
        DN = datenum(ID(3:8),'yymmdd');
    end
else
    ID = cell(size(STR));
    tagtype = cell(size(STR));
    sp = cell(size(STR));
    DN = cell(size(STR));
    for i = 1:prod(size(STR))
        if isempty(STR{i}); continue; end
        str = STR{i};
        str = char(str);
        ID{i} = str(1:min([regexp(str,' '), regexp(str,'prh')])-1);
        if ~isempty(strfind(ID{i},'PAD')); tagtype{i} = 'PAD';
        elseif ~isempty(strfind(ID{i},'_')); tagtype{i} = 'DTAG';
        else tagtype{i} = 'CATS';
        end
        sp{i} = ID{i}(1:2);
        if strcmp(tagtype{i},'DTAG')
            DN = datenum([2000+str2num(ID{i}(3:4)) 1 str2num(ID{i}(6:8)) 0 0 0]);
        else
            DN = datenum(ID{i}(3:8),'yymmdd');
        end
    end
end
species = sp;
datenumber = DN;