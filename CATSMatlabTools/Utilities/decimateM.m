function M = decimateM(M,ofs,sHz,df,nout,name,removenans)
% M = data to decimate
% ofs = sample rate of M
% sHz = rate at which M was sampled (may be lower than ofs since ofs is usually the sample rate of all data, and values are repeated)
% df = decimation factor 
% nout = number of samples in the resulting vector (accounts for slightly
% different sampling when df does not do exactly evenly into the size of M)
% name = name of sample rate to list if there is an error (optional)
% by default, this function remmoves nans (linearly interpolates), then
% puts them back in

% This function accounts for the original sampling rates of the data.  So,
% for instance, if pressure was sampled at 10 Hz, then upsampled to 50 Hz
% during import, downsampling to 10 Hz is just sampling the data.  But if
% the magnetometer was sampled at 50 Hz, then downsampling involves a
% loss-less decimation using dec_dc from animaltags.org

if nargin<7; removenans = true; end

if removenans
    nI = cell(1,size(M,2));
    for i = 1:length(nI); nI{i} = unique(round(find(isnan(M(:,i)))/df)); nI{i}(nI{i}==0) = 1; end
    M = fixgaps(edgenans(inpaint_nans(M)));
    M(isnan(M)) = 0;
    
end
    
if nargin<5; name = 'sampling rate'; end
if sHz<ofs/df; sHz = ofs/df; end; % if the rate at which data was sampled is smaller than the final rate, don't sample smaller than the final rate.
if sHz>ofs; sHz = ofs; end; % if data has already been downsampled (e.g. the accelerometer), further decimate it without sampling it
if abs(round(ofs/sHz)-ofs/sHz)>.01; error([ name ' does not divide evenly into original sample rate']); end
M = M(1:ofs/sHz:end,:);
if abs(round(df/(ofs/sHz)) - df/(ofs/sHz)) > .01; error(['decimation factor does not divide evenly into ' name]); end
M = decdc(M,df/(ofs/sHz));

if removenans
    for i = 1:length(nI);
        nI{i}(nI{i} == size(M,1)) = size(M,1);
        M(nI{i},i) = nan;
    end
end

if abs(size(M,1)-nout)>1; error(['intended length is ' num2str(size(M,1)-1) ' samples different from the decimated length']);
else
    for i = 1:size(M,2)
        X(1:nout,i) = interp2length(M(:,i),ofs/df,ofs/df,nout);
    end
    M = X;
end

    