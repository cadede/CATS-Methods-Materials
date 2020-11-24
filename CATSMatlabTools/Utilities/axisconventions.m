function [axA, axM, axG] = axisconventions(tagtype)
%%
% David Cade
% version 9.29.15
% Goldbogen Lab
% Stanford University

% Example for platypus:from Niko: (want to change to be right-handed with NorthEastDown orieentation) for
% orientation in all axes.  So original platypus was:
% Compass: (platypus)
% 1 surge with front facing North giving min value (opposite)
% 2 heave with top facing North giving max value  (opposite and switch with three)
% 3 sway with left facing North giving min value (sway is correct (right is max), but switch to be two)
% Gyro:
% 1 roll with left rotating downwards giving min/neg values -- Right (from tag tutorial): The roll angle in degrees. Roll is defined in the animal frame as a rotation around the caudorostralaxis with a positive roll being a clockwise rotation looking rostrally along the axis.
% 2 pitch with front rotating upwards giving min/neg values -- Wrong: Pitch is defined in the animal frame as a rotation around the leftright axis with a positive pitch being a clockwise rotation looking towards the right along the axis
% 3 yaw with clockwise rotation giving min/neg values -- Wrong: Yaw is defined in the animal frame as a rotation around the ventrodorsal axis with a positive yaw being an anti-clockwise rotation looking dorsally along the axis.
% Acc:
% 1 surge with front facing up giving max/pos value -- Right
% 2 sway with left facing down giving min/neg value -- wrong, this needs to be positive (right up = positive)
% 3 heave with top facing up giving max/pos value -- wrong- in NED, z axis is facing down so -1

x = 1; y = 2; z = 3; 
% in here, put what each column actually is (that is, the [y -x z] implies
% that data in the 1 column is actually y, and data in the second column is
% actually negative x
if strcmpi(tagtype,'Wireless') || strcmpi(tagtype,'motus') || strcmpi(tagtype,'acousonde')
    axAo = [x y z];
    axMo = [x y z];
    axGo = [x y z]; % no gyros in acousonde
elseif strcmpi(tagtype,'platypus')
axAo = [x -y -z];
axMo = [-x -z y];
axGo = [x -y -z];
elseif strcmpi(tagtype,'data')
axAo = [-x -y z];
axMo = [-x -y z];
axGo = [-x -y z];
elseif strcmpi(tagtype,'TDR10')
axAo = [-x y z];
axMo = [x y z];
axGo = [x y z];
elseif strcmpi(tagtype,'TDR10_rotate')
axAo = [x -y z];
axMo = [x y z];
axGo = [x y z];
elseif strcmpi(tagtype,'kitten')
axAo = [y -x z];
axMo = [y -x z];
axGo = [y -x z];
elseif strcmpi(tagtype,'minion')
    %no longer have this tag to check conventions... reconstructed from
    %calibrateCATS1
axAo = [x y z];
axMo = [x -z -y];
axGo = [x y z];    
    
elseif strcmpi(tagtype,'froback')
axAo = [-x y -z];
axMo = [y -x z];
axGo = [-x y -z];
    
elseif strcmpi(tagtype,'PAD')
   axAo = [x -y -z];
   axMo = [-y -x -z];
   axGo = [-x y -z];
      
elseif strcmpi(tagtype,'4k')
   axAo = [y x -z];
   axMo = [x -y -z];
   axGo = [-y -x -z];
   
elseif strcmpi(tagtype,'JellyOld')
   axAo = [z y -x]; %axAo = [-x y -z]; originals over here, now adjusted for jelly orientation of horizontal
   axMo = [-z -y x]; %axMo = [x -y z];
   axGo = [z y -x]; %axGo = [-x y -z];
elseif strcmpi(tagtype,'Jelly')
    axAo = [z y x]; %axAo = [-x y -z]; originals over here, now adjusted for jelly orientation of horizontal
    axMo = [y -z -x]; %axMo = [x -y z]; 
    axGo = [z y x]; %axGo = [-x y -z];
elseif strcmpi(tagtype,'OpenTag')
   axAo = [y -x z];
   axMo = [y -x z]; % unsure
   axGo = [y -x z];  % unsure
else error('TagType not Recognized!');
       
end

axA = diag(zeros(3,1));
for j = 1:3
    axA(abs(axAo) == j,j) = sign(axAo(abs(axAo) == j));
end
axM = diag(zeros(3,1));
for j = 1:3
    axM(abs(axMo) == j,j) = sign(axMo(abs(axMo) == j));
end
axG = diag(zeros(3,1));
for j = 1:3
    axG(abs(axGo) == j,j) = sign(axGo(abs(axGo) == j));
end

