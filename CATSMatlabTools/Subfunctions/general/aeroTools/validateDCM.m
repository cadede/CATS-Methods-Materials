function validateDCM(dcm, varargin )
%VALIDATEDCM internal function to check that the input dcm is orthogonal
%and proper.
% The criteria for this check are:
%      - the transpose of the matrix multiplied by the matrix equals 
%        1 +/- tolerance
%      - determinant of matrix == +1
%
% This function is intended for internal use only.

%   Copyright 2017-2018 The MathWorks, Inc.

%Validate inputs
validateattributes(dcm,{'single','double'},{'real','nonempty','size',[3 3 NaN]},'','dcm');
len = size(dcm,3);

tolerance = eps(2);
if nargin < 2
    action = 'none';
elseif nargin < 3
    action = varargin{1};
else
    action = varargin{1};
    tolerance = varargin{2};   
end

validateattributes(tolerance,{'numeric'},{'real','nonempty','scalar', 'nonnegative'},'','Tolerance');

actionoptions = {'None','Warning','Error'};
if isstring(action) || ischar(action)
    action = lower(validatestring(action, actionoptions,'','Action'));
else
    validateattributes(action,{'numeric'},{'real','nonempty','scalar','<=', 3,'>=', 1},'','Action');
    action = lower(actionoptions{action});
end


if ~strcmp(action,'none')
    % Verify that rotation matrices are proper
    isnotproper = arrayfun(@(k) abs(det(dcm(:,:,k)) - 1) > tolerance,1:len);
    if strcmp(action,'error')
        if sum(isnotproper) > 0 && len == 1
            error(message('aero:validateDCM:isNotProper'));
        elseif sum(isnotproper) == 1 && len > 1
            error(message('aero:validateDCM:isNotProperSingle',num2str(find(isnotproper > 0,1))));
        elseif sum(isnotproper) > 1 && len > 1
            error(message('aero:validateDCM:isNotProperMulti',num2str(find(isnotproper > 0))));
        end
    else
        if sum(isnotproper) > 0 && len == 1
            warning(message('aero:validateDCM:isNotProper'));
        elseif sum(isnotproper) == 1 && len > 1
            warning(message('aero:validateDCM:isNotProperSingle',num2str(find(isnotproper > 0,1))));
        elseif sum(isnotproper) > 1 && len > 1
            warning(message('aero:validateDCM:isNotProperMulti',num2str(find(isnotproper > 0))));
        end
    end
    % Verify that rotation matrices are orthogonal
    isnotorthogonal = arrayfun(@(k) any(reshape(abs(transpose(dcm(:,:,k))* ...
        dcm(:,:,k)-eye(3)) > tolerance,[],1)) && ~isnotproper(k), 1:len);
    if strcmp(action,'error')
        if sum(isnotorthogonal) > 0 && len == 1
            error(message('aero:validateDCM:isNotOrthogonal'));
        elseif sum(isnotorthogonal) == 1 && len > 1
            error(message('aero:validateDCM:isNotOrthogonalSingle',num2str(find(isnotorthogonal > 0,1))));
        elseif sum(isnotorthogonal) > 1 && len > 1
            error(message('aero:validateDCM:isNotOrthogonalMulti',num2str(find(isnotorthogonal > 0))));
        end
    else
        if sum(isnotorthogonal) > 0 && len == 1
            warning(message('aero:validateDCM:isNotOrthogonal'));
        elseif sum(isnotorthogonal) == 1 && len > 1
            warning(message('aero:validateDCM:isNotOrthogonalSingle',num2str(find(isnotorthogonal > 0,1))));
        elseif sum(isnotorthogonal) > 1 && len > 1
            warning(message('aero:validateDCM:isNotOrthogonalMulti',num2str(find(isnotorthogonal > 0))));
        end
    end
    

end
end
