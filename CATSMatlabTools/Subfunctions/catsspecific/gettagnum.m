function tagnum = gettagnum(fname)

dashes = regexp(fname,'-'); if isempty(dashes); dashes = 0; else dashes = dashes(end); end
digs=find(isstrprop(fname(dashes+1:end),'digit'),2,'first')+dashes; if diff(digs)~=1; digs = digs(1); end
tagnum = str2num(fname(digs));
