function [ TurbineInfo ] = fnReadTurbinesFile( filename, skip )
%FNREADTURBINESFILE Reads a file containing turbine information
%   At the moment this is simply a CSV. Read it and do some error-checking.
%   Filename is the name of the file as a char.
%   skip is the number of rows to skip at the top of the CSV. CSV columns are as follows:
%
%   x coordinate
%   y coordinate
%   z coordinate (hub height, relative to mean sea level. So -30 is 30m
%           below MSL) - This must be present, but is ignored for D3D - it's in here because
%           it's needed for MIKE and the same CSV files are used.
%           The D3D script will take dz and work out z from it.
%   Height above seabed (dz)
%   (optional): Turbine Orientation.

%   TurbineInfo is a structure containing fields x, y, z and dz, (and optionally o) each of which is
%   a nx1 array.

% Copyright Simon Waldman / Heriot-Watt University, 2014-2015

if (nargin < 1)
    error ('Need a filename');
end
if ~isa(filename, 'char')
    error('First input variable isn''t a char');
end
if ~exist(filename,'file')
    error('File not found');
end
if (nargin < 2)
    skip = 0;
end
%FIXME should probably check that skip is a number somehow?

M = csvread(filename, skip); %yes, matlab makes an exception and uses 0-based indexing here for "convenience", so it is "rows to skip" rather than "row number to start at".

%FIXME what checking of the input should I be doing here? Successful read? 

switch size(M,2)    % how many columns?
    case 4
        TurbineInfo = struct('x', M(:,1), 'y', M(:,2), 'dz', M(:,4));
        TurbineInfo.o = [];
    case 5
        TurbineInfo = struct('x', M(:,1), 'y', M(:,2), 'dz', M(:,4), 'o', M(:,5));
    otherwise
        error('CSV file should have four or five columns: X, Y, Z, dZ, (O).');
end

end

