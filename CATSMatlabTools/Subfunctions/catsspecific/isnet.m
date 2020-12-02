function I = isnet(head,pitch,fs)
%%
dhead = head*180/pi;
I = pitch<45*pi/180;
I0 = dhead<5&dhead>-5;
I180 = abs(dhead)>175;
I90 = dhead>85 & dhead<95;
In90 = dhead>-95&dhead<-85;
for i = find(I0)'
    I0(max(1,i-fs*45):min(length(pitch),i+fs*45)) = true;
end
for i = find(I180)'
    I180(max(1,i-fs*45):min(length(pitch),i+fs*45)) = true;
end
for i = find(I90)'
    I90(max(1,i-fs*45):min(length(pitch),i+fs*45)) = true;
end
for i = find(In90)'
    In90(max(1,i-fs*45):min(length(pitch),i+fs*45)) = true;
end
I = I&I0&I180&I90&I-90;
% close gaps of 15 seconds or less
Ic = find(~I);
[s,e] = consec(Ic);
gI =find((e-s)<15*fs);
for i = 1:length(gI)
    I(s(gI(i)):e(gI(i))) = true;
end
% get rid of sections less than 30 s long
Ic = find(I);
[s,e] = consec(Ic);
gI =find((e-s)<30*fs);
for i = 1:length(gI)
    I(s(gI(i)):e(gI(i))) = false;
end
% s1 = subplot(211); plot(head,'.');
% s2 = subplot(212); plot(I);%(1:length(I0),I0,1:length(I180),I180*2,1:length(I90),I90*3,1:length(In90),In90*4); 
% linkaxes([s1 s2],'x');
% ylim([0 2])
