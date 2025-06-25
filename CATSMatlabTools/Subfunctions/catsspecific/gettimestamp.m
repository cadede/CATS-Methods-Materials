function [time,flag] = gettimestamp(frame,useoldstyle,grayThresh)
%%
if nargin==2 && useoldstyle 
    [time,flag] = gettimestampSEC(frame);
    return;
end

if nargin <3
    incgray = false;
else
    incgray = true;
end

% figure(1); clf;
% I = frame(1054:1064,1675:1772,:);
for k = 1:2
    if size(frame,1) == 720 && size(frame,2) == 1280
        I = frame(695:703,[1037:1055 1060:1080 1083:1100 1106:1132],:);
        incgray = false;
    elseif size(frame,1)  == 1080 && size(frame,2) == 1920
        if k == 1
            I = frame(1055:1064,[1675:1695 1699:1718 1723:1743 1746:1772],:);
        elseif k == 2 % new cats tags have box in a different place
            I = frame(1064:1074,[1820:1838 1842:1863 1867:1886 1890:1920],:);
        % else
        %     dbstop if k>2
        %     error('Timestamp could not be read, try plotting frame and extracting indices of timestamp')
        end
        %     incgray = true;
    else error ('Double check that you have timestamps in the bottom right of your screen in a black box, else set vid4k = true to run make4kmovieTimes instead of this script. Alternatively, could be a frame size issue');
    end
    I = 255-I;
    % image(I);

    imagen = I;
    %turns some gray to black for high resolution frames
    if incgray; imagen(imagen<grayThresh) = 0; end
    if size(imagen,3)==3 %RGB image
        imagen=rgb2gray(imagen);
    end
    % Convert to BW
    threshold = graythresh(imagen);
    re =~im2bw(imagen,threshold);
    % Remove all object containing fewer than 30 pixels
    % imagen = bwareaopen(imagen,10);
    %Storage matrix word from image
    word=[ ];
    % re=imagen;
    %Opens text.txt as file for write
    % fid = fopen('text.txt', 'wt');
    % Load templates
    load CATStemplates
    % templates = CATStemplates;
    % global templates
    % Compute the number of letters in template file
    num_letras=size(templates,2);
    imgn=re;
    %Uncomment line below to see lines one by one
    %imshow(fl);pause(0.5)
    %-----------------------------------------------------------------
    % Label and count connected components
    [L Ne] = bwlabel(imgn);
    try
        if Ne~=9
            for i = 1:100
                if threshold*(1+i/50)>1; Ne = 8; break; end
                re =~im2bw(imagen,threshold*(1+i/50)); %if it's not reading enough letters, try increasing the threshold until you have 9
                imgn=re;
                [L Ne] = bwlabel(imgn);
                if Ne == 9; break; end
            end
        end
    catch
        Ne = 8;
    end
    if Ne~=9; time = 0; flag = true; continue; end
    flag = false;
    for n=1:Ne
        [r,c] = find(L==n);
        if sum(length(r),length(c))<9; flag = true; time = 0; continue; end % if there's just a little piece, you're probably not getting the whole number
        % Extract letter
        n1=imgn(min(r):max(r),min(c):max(c));
        % Resize letter (same size of template)
        img_r=imresize(n1,[42 24]);
        %Uncomment line below to see letters one by one
        %imshow(img_r);pause(0.5)
        %-------------------------------------------------------------------
        % Call fcn to convert image to text
        letter=read_number(img_r,num_letras,templates);
        % Letter concatenation
        word=[word letter];
    end
    if flag; continue; end
    %final check to see if number makes sense
    if str2num(word(1:2))>23 || str2num(word(3:4))>59 || str2num(word(5:6)) > 59 || sum(sum(sum(I<10)))<20;
        time = 0; flag = true; continue;
    end
    time = datenum([word(1:2) ':' word(3:4) ':' word(5:6) '.' word(7:end)]);
    time = time-floor(time);
    if ~flag;break; end
end

