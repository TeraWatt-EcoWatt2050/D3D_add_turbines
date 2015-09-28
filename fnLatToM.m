function [ dy ] = fnLatToM( dlat,alat )
%FNLATTOM Converts a difference in latitudes to a distance in metres
% taken from http://www-pord.ucsd.edu/~matlab/coord.htm

% dy   = latitude difference in meters
% dlat = latitude difference in degrees
% alat = average latitude between the two fixes
% Reference: American Practical Navigator, Vol II, 1975 Edition, p 5 

rlat = alat * pi/180;
m = 111132.09 * ones(size(rlat)) - ...
    566.05 * cos(2 * rlat) + 1.2 * cos(4 * rlat);
dy = dlat .* m ;

end

