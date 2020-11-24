function x = version4
% Return 1 if this is MATLAB version 4.

x = strncmp(version, '4.', 2);
