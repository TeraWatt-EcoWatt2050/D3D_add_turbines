function [ dx ] = fnLonToM( dlon, alat )
%FNLONTOM Converts a difference in longitudes to a distance in metres
% taken from http://www-pord.ucsd.edu/~matlab/coord.htm

% dx = lon_to_m(dlon, alat)
% dx   = longitude difference in meters
% dlon = longitude difference in degrees
% alat = average latitude between the two fixes

rlat = alat * pi/180;
p = 111415.13 * cos(rlat) - 94.55 * cos(3 * rlat);
dx = dlon .* p;

end

