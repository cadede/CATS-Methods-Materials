function [time,flag] = gettimestampSEC(frame)
%% %
% figure(1); clf;
I = frame(1054:1064,1730:1796,:);
I = 255-I;

imagen = I;
if size(imagen,3)==3 %RGB image
    imagen=rgb2gray(imagen);
end
% Convert to BW
threshold = graythresh(imagen);
imagen =~im2bw(imagen,threshold);
% Remove all object containing fewer than 30 pixels
imagen = bwareaopen(imagen,10);
%Storage matrix word from image
word=[ ];
re=imagen;
%Opens text.txt as file for write
% fid = fopen('text.txt', 'wt');
% Load templates
load CATStemplates
% templates = CATStemplates;
% global templates
% Compute the number of letters in template file
num_letras=size(templates,2);
fl = re;
imgn=fl;
%Uncomment line below to see lines one by one
%imshow(fl);pause(0.5)
%-----------------------------------------------------------------
% Label and count connected components
[L Ne] = bwlabel(imgn);
if Ne<6; time = 0; flag = true; return; end
flag = false;
for n=1:Ne
    [r,c] = find(L==n);
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
time = datenum([word(1:2) ':' word(3:4) ':' word(5:6)]);
time = time-floor(time);


