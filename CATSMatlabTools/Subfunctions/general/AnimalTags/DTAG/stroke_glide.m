 function      [GL,KK] = stroke_glide(est,fs,n,J,tmax)
    %
%    [GL,KK] = stroke_glide(pry,Ahf,fs,n,J,tmax)
%    Identify glides (GL) and strokes(KK).
%    
%    INPUT:
%       pry = estimated body rotations around the fundamental axis of 
%            rotation [pitch, roll, yaw] which depending on the 
%            locomotion style it can be around the y, x or z axis. Angles  
%            not estimated are set to 0. Angles are in radians.
%       Ahf = high-pass filtered 3-axis acceleration signal, in m/s2.  
%       fs = sensor sampling rate in Hz 
%       n = if Ahf is used it reffers to the axis of acceleration 
%           1 for accelerations in the x axis, longitudinal acceleration.
%           2 for accelerations in the y axis, lateral acceleration.
%           3 for accelerations in the z axis, dorso-ventral acceleration.
%       n = if pry is used it reffers to the fundamental axis around which 
%           body rotations are analysed. 
%           1 for rotations around the y axis, pitch method.
%           2 for rotations around the x axis, roll method.
%           3 for rotations around the z axis, yaw method.
%       k = sample range over which to analyse.
%       J = magnitude threshold for detecting a fluke stroke in radians.
%           If J is not given, fluke strokes will not be located 
%           but the rotations signal (pry) will be computed.If no J is
%           given or J=[], no GL and KK output will be generated.
%       tmax = maximum duration allowable for a fluke stroke in seconds. 
%           A fluke stroke is counted whenever there is a cyclic variation 
%           in the pitch deviation with peak-to-peak magnitude 
%           greater than +/-J and consistent with a fluke stroke duration
%           of less than tmax seconds, e.g., for Mesoplodon choose tmax=4.
%           If no tmax is given or tmax=[], no GL and KK output will be
%           generated. 
%
%
%    OUTPUT: 
%       GL = matrix containing the start time (first column) and end time
%           (2nd column) of any glides (i.e., no zero crossings in tmax or 
%           more seconds).Times are in seconds.
%       KK = matrix of cues to zero crossings in seconds (1st column) and
%           zero-crossing directions (2nd column). +1 means a 
%           positive-going zero-crossing. Times are in seconds.
%
%    NOTE: Be aware that when using devices that combine different sensors 
%       as in here (accelerometer and magnetometer) their coordinate
%       systems should be aligned. If they are not physically aligned, 
%       invert as necessary for all the sensor's axes to be aligned so that
%       when a positive rotation in roll occurs in one sensor it is also 
% `     positive in all the sensors. 

%
%   Lucia Martina Martin Lopez & Mark Johnson  (June 2013)
                     

 
 
if isempty(J) | isempty(tmax),
   GL = [] ; KK=[] ;
    fprintf('Cues for strokes(KK) and glides (GL) are not given as J and tmax are not set');
   return ;
end

 
 % Find cues to each zero-crossing in vector pry(:,n), rotations around the
 % n axis. 
 
 K = findzc(est(:,n),J,tmax*fs/2) ;
 
 % find glides - any interval between zeros crossings greater than tmax
 k = find(K(2:end,1)-K(1:end-1,2)>fs*tmax) ;
 glk = [K(k,1)-1 K(k+1,2)+1] ;
 
 % shorten the glides to only include sections with jerk < J
 glc = round(mean(glk,2)) ;
 
 for k=1:length(glc),
     kk = glc(k):-1:glk(k,1) ;
     
     test=find (isnan(est(kk,n)));
     if ~isempty(test),
         glc(k)=NaN;
         glk(k,1)=NaN;
         glk(k,2)=NaN;
     else
         glk(k,1) = glc(k) - find(abs(est(kk,n))>=J,1)+1 ;
         kk = glc(k):glk(k,2) ;
         glk(k,2) = glc(k) + find(abs(est(kk,n))>=J,1)-1 ;
     end
 end
 % convert sample numbers to times in seconds
 KK = [mean(K(:,1:2),2)/fs K(:,3)] ;
 
 GL = glk/fs ;
 GL = GL(find(GL(:,2)-GL(:,1)>tmax/2),:) ;
