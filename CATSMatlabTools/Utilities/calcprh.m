function [pitch,roll,head] = calcprh(Aw,Mw,dec)
% dec is optional (declination).  If it is not input, heading will be based
% on magnetic north.

% calculates animal pitch, roll and heading from unfiltered animal-frame
% sensor data

 At_mag = sqrt(sum(Aw.^2,2));
    pitch=asin(Aw(:,1)./At_mag);  % in radians
    
% Descriptions below from animal tag tutorial. Steps 1-5 are done in the script tagframe2whaleframe.m.  See Johnson 2011.
% 3) Roll
% "THE Roll estimator leads to high errors at high pitch angles. So roll errors could be 6 times
% larger when |p|=80 degrees than  when it is horizontal"->animal_orientation_tutorial.pdf
% ROLL HERE: is not as robust at high pitch angles either completely up or completely down.
%
% During each sample time:
%2) Measure the accelerometer and magnetometer vectors, At and Mt, in the
%tag frame
%3) Multiply At by W to get Aa
%4) Estimate pitch and roll using:
% p_hat=-asin(Aa(:,1)/abs(At)), and r_hat=atan(Aa(:,2)/Aa(:,3))  %DC note: our axis conventions differ from Johnson 2011, thus the difference in signs below
%5) Multiply Mt by W to get Ma
%6) Compute the rotation matrices P(p_hat) and R(r_hat) from p_hat and
%r_hat.
%7) Multiply Ma by R(r_hat)-1*P(p_hat)-1 to get Mh.  (could use just the
%transpose here instead of the inverse since the matrix is assumed to be
%symetric
%8) Estimate the heading using h_hat=atan(-mh(:,2)/mh(:,1))

roll=atan2(-Aw(:,2),-Aw(:,3)); % atan2 is used to maintain the correct quadrants
head = nan(size(roll)); head2 = head;
Mh = nan(size(Mw));

% DETERMINE HEADING: TRANSFORM ANIMAL FRAME TO NAVIGATION FRAME:
% Estimator of Heading (used to get magnetometer values with respect to the
% navigation frame). These are the rotation matrices for pitch and roll
% respectively which will be used to go from animal to the navigation
% frame.

% Loop to calculate heading (based on pitch and roll of the animal,
% magnetic field, and offset from True North). 2 Methods for Values with
% respect to True North are calculated. Method 1 (Tyack Johnson), and
% Method 2 is based on Honeywell Mx, My quadrant checks.
if exist('dec','var') && ~isempty(dec)
for i=1:length(Mw(:,1))
    % note that these rotation matrices differ from Johnson 2011 due to
    % differences in axis conventions
    % H = [cos(h) -sin(h) 0; sin(h) cos(h) 0; 0 0 1];
    P = [cos(pitch(i)) 0 sin(pitch(i)); 0 1 0; -sin(pitch(i)) 0 cos(pitch(i))];
    R = [1 0 0; 0 cos(roll(i)) -sin(roll(i)); 0 sin(roll(i)) cos(roll(i))];
    
    % Estimate step 7) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Mh= is to siginify that it comes from the horizontal frame. Independent of pitch and roll.
    %Mh=Mw*R_roll^-1*P_pitch^-1;
    %Mh=b.[cos(i)*(cos(yaw).*cos(d)+sin(yaw)*sin(d)),-cos(i)*(sin(yaw).*cos(d)-cos(yaw)*sin(d)),-sin(i)]; %this is right, which makes me think their rotation matrices are wrong because using them as is should not give this, which is what I get from our NED reference frame
    Mh(i,:) = Mw(i,:)*R'*P';
    % Mh = b*[cos(i)cos(h),-cos(i)sin(h),-sin(i)] using the substitution h = yaw-d.
    head(i) = wrapToPi(atan2(-Mh(i,2),Mh(i,1))+dec);% add dec back in to get the true north bearing
    if head(i) >0; head2(i) = head(i)-pi; else head2(i) = head(i)+pi; end
end

end