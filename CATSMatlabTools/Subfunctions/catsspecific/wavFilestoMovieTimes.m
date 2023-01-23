audioloc = [fileloc 'AudioData\'];
D = dir(audioloc); audfiles = {D.name}'; audfiles = audfiles(~cellfun(@isempty,cellfun(@(x) strfind(x,'.wav'),audfiles,'uniformoutput',false)));
if isempty(audfiles); [~,audioloc]=uigetfile('*.wav', 'select wav file in directory with audiofiles for this deployment'); D = dir(audioloc); audfiles = {D.name}'; audfiles = audfiles(~cellfun(@isempty,cellfun(@(x) strfind(x,'.wav'),audfiles,'uniformoutput',false))); end
vidDN = nan(size(audfiles)); vidDurs = vidDN; frameTimes = {}; movies = {}; vidNam = {}; vid4k = false; frameSize = nan(1,2);
disp('Trying to read timestamps from audiofiles with dateformat yymmdd-HHMMSS.wav');
for i = 1:length(audfiles)
    disp(['Reading audio file ' num2str(i)]);
    clear aud;
    try
        [aud.data,aud.rate,aud.bits] = wavread([audioloc audfiles{i}]);
    catch
        audioI = audioinfo([audioloc audfiles{i}]);
        [aud.data,aud.rate] = audioread([audioloc audfiles{i}]);
        aud.bits = audioI.BitsPerSample;
    end
    if size(aud.data,2) == 2 && all(aud.data(:,2) == 0); aud.data = aud.data(:,1); aud.nrChannels = 1;
    elseif size(aud.data,2) == 2 && all(aud.data(:,1) == 0); aud.data = aud.data(:,2); aud.nrChannels = 1;
    else aud.nrChannels = size(aud.data,2);
    end
    aud.totalDuration = size(aud.data,1)/aud.rate;
    totalDuration = aud.totalDuration;
    vidDN(i) = datenum(audfiles{i}(end-16:end-4),'yymmdd-HHMMSS');
    vidDurs(i) = totalDuration;
    lastwarn('');
    try
        save([audioloc audfiles{i}(1:end-4) 'audio.mat'],'aud','totalDuration');
        if ~isempty(lastwarn)
            error(lastwarn);
        end
    catch %v7.3 allows for bigger files, but makes a freaking huge file if used when you don't need it
        save([audioloc audfiles{i}(1:end-4) 'audio.mat'],'aud','totalDuration','-v7.3');
        disp(['Made a version 7.3 file (audio ' movies{n}(1:end-4) ' was a large file)']);
    end
    movies{i} = audfiles{i}; vidNam{i,1} = audfiles{i};
    disp(['Starts at: ' datestr(vidDN(i)) ', and is ' num2str(vidDurs(i)/60) ' minutes long, sampled at ' num2str(aud.rate) ' Hz.']);
end
save([fileloc filename(1:end-4) 'movieTimes.mat'],'vidDurs','frameTimes','movies','vidDN','vidNam','frameSize');