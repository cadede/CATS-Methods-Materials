function [Aw,Mw,Gw,pitch,roll,head,W,tagprh] = tagframe2whaleframe(At,Mt,Gt,Depth,dec,calperiodI,W)
% inputs: tag calibrated data At, Mt, Gt, Depth, and either a set of
% calibration indices (stable horizontal periods and a dive), or a W
% matrix.  This function assumes calperiod is indexed to the section of data
% included in At and Mt, Depth etc. If calperiod is empty, W is used, else
% W is calculated from the calperiod entries. head is referenced to true
% north
% Depth only needs to be input if using calperiodI.  If using W, it is not needed

% outputs: whale frame data Aw, Mw, Gw, pitch, roll, heading, and the tag
% orientation data W and tagprh

% W^-1 uses: yaw_o, pitch_o, and roll_o initially, and then is updated for
% the x amount of times that it surfaces and begins to dive again.
% W^-1=(Y_yaw_o*P_pitch_o*R_roll_o)^-1 = R_roll_o'*P_pitch_o'*Y_yaw_o' %DC- so yaw_o, pich_o,roll_o are the animal's with respect to the tag orientation


if ~isempty(calperiodI) % if calperiodI is empty just apply W matrix.  If using calperiodI, W will be rewritten or can be excluded
    % for a calperiod, calculated the W matrix
    % calperiod are paired indices at the start and end of every calibration period
    % for i = 1:length(calperiodI)
    asc = false; %ascent switch
    At_mean = nan(1,3); n = 0;
    for j = 1:2:length(calperiodI)-2
        n = n+1;
        At_mean(n,:) = nanmean(At(calperiodI(j):calperiodI(j+1),:),1);
    end
    At_mean = At_mean./repmat(sqrt(sum(At_mean.^2,2)),1,3); % divide by the magnitude (more robust to dynamic acceleration errors)
    At_mean=nanmean(At_mean,1);
    pitch0 = asin(At_mean(1));
    roll0 = atan2(-At_mean(2),-At_mean(3));
    % check but should be fine:
    if abs(At_mean(2)+cos(pitch0)*sin(roll0))>.05 || abs(At_mean(3)+cos(pitch0)*cos(roll0))>.05
        error('roll0 appears to be in the wrong quadrant (or numbers are way off)');
    end
    % H = [cos(h) -sin(h) 0; sin(h) cos(h) 0; 0 0 1];
    P = [cos(pitch0) 0 sin(pitch0); 0 1 0; -sin(pitch0) 0 cos(pitch0)];
    R = [1 0 0; 0 cos(roll0) -sin(roll0); 0 sin(roll0) cos(roll0)];
    Apr = At(calperiodI(end-1):calperiodI(end),:)*R'*P'; % can be ascent or descent
    %     if mode(sign(Apr(:,2))) == -1
    % yaw0 = atan2(Apr(:,2),-Apr(:,1));  %method 1- find the circ_mean of the yaw that zeros each y-value
    yaw02 = atan2(mean(Apr(:,2)),-mean(Apr(:,1))); %find the yaw that zeros the mean of the y-values
    % the math in the packet suggests that y should be negative, but that's backwards.  the x needs to be negative so that yaw is close to 0 as long as pitch is the way it's supposed to be (negative for a dive), but the calculation of atan2 has to have x positive to be small angle
    %     else
    %         yaw0 = atan2(-Apr(:,2),Apr(:,1));
    %         yaw02 = atan2(-mean(Apr(:,2)),mean(Apr(:,1)));
    %     end
    
    %allow for ascents to be selected
    if Depth(calperiodI(end))<Depth(calperiodI(end-1))
        asc = true;
        %     yaw0 = wrapToPi(yaw0+pi);
        yaw02 = wrapToPi(yaw02+pi);
    end
    
%     yaw0 = circ_mean(yaw0);
    % Y = [cos(yaw0) -sin(yaw0) 0; sin(yaw0) cos(yaw0) 0; 0 0 1];
    Y2 = [cos(yaw02) -sin(yaw02) 0; sin(yaw02) cos(yaw02) 0; 0 0 1];
    % Aw0 = Apr*Y';
    Aw02 = Apr*Y2';
    % testW = [sum(Aw0(:,[1 3]).^2,2) Aw0(:,2) sum(Aw02(:,[1 3]).^2,2) Aw02(:,2)]; % the sum of ax^2 and az^2 should be one since they are the cos(pitch) and sin(pitch), and y should be zero if properly adjusted with 0 roll
    % just for my own curiousity, test which yaw method is better for this
    % case
    % [~,whichy] = min(sum((testW(:,[1 3])-1).^2));
    % [~,whichy2] = min(sum((testW(:,[2 4])).^2));
    % if whichy ~= whichy2
    %     disp('Both yaw methods have advantages, but');
    %     ys = sqrt(sum((testW(:,[1 3])-1).^2));
    %     ys2 = sqrt(sum((testW(:,[2 4])).^2));
    %     if ys(whichy)/ys(whichy2) < ys2(whichy2)/ys2(whichy)
    %         whichy2 = whichy; else whichy = whichy2;
    %     end
    % end
    % disp(['yaw method ' num2str(whichy) ' is better']);
    % if whichy == 2; Y = Y2; yaw0 = yaw02; Aw0 = Aw02; end
    Y = Y2; Aw0 = Aw02; yaw0 = yaw02;
    % just a check to ensure it's in the right quadrant (a check on my sanity)
    if ((mean(Aw0(:,1))>0 || mean(Aw0(:,3))>0) && ~asc) || asc && (mean(Aw0(:,1))<0 || mean(Aw0(:,3))>0); error('yaw0 appears to be in the wrong quadrant'); end  %ensures acceleration is negative in x for the dive by choosing the right quadrant for arctan
    W = (Y*P*R)';
    tagprh = [pitch0 roll0 yaw0]; % could use to determine if there was tagslip
end
    %
    % Relationship between animal and navigation frames:
    % Examples of Rotation or Transformation Matrices:
    % Recap:
    % Q= orientation of the animal with respect to the navigation frame.
    % W= orientation of the animal with respect to the tag
    Aw = At*W;
    Mw = Mt*W;
    Gw = Gt*W;
   
    [pitch,roll,head] = calcprh(Aw,Mw,dec);





% 
