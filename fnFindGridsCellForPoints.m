function [ Mout, Nout ] = fnFindGridsCellForPoints( gX, gY, pX, pY )
%FNFINDGRIDSCELLFORPOINTS Finds which grid cells a set of points are in
%   Inputs:
%   gX and gY are 2D arrays of the coordinates of grid points representing
%   the centres of cells.
%   pX and pY are 1D arrays giving the coordinates of the query points
%   Outputs:
%   Mout and Nout are 1D arrays giving the grid coordinates of the cells that contain
%   the points. They are of the same length as pX and pY and has the points in
%   the same order.

%   NB this assumes a rectilinear grid (in whichever coordinate system is
%   used here), and will give incorrect results on curvilinear or unstructured grids.

% Copyright (C) Simon Waldman / Heriot-Watt University Dec 2014 - Apr 2015

if nargin < 4
    error('Not enough arguments.');
end
if ~all(size(gX)==size(gY)) || ndims(gX) ~= 2
    error('gX and/or gY are the wrong size. Must be 2D, and same size as each other.');
end
if ~all(size(pX) == size(pY))
    error('pX and pY must be 1D and the same length as each other.');
end

% Let's assume that if the cell whose centre is closest to the current
% point is the one that it's in. Shoudl be valid for rectilinear grids.

%Reform the grid into a nx2 matrix
GridCoords = [ reshape(gX, [], 1) reshape(gY, [], 1) ];  %wow, that's horrible syntax.
% We won't actually get rid of the NaNs, just set them to a silly number.
% This is because it makes it easier to retain a relationship between the
% row number in this array and m,n coords
GridCoords(any(isnan(GridCoords),2),:) = flintmax; %set to the largest value possible in double precision. Shoudl be well out of the way of nearest neighbour!
QueryCoords = [pX pY];

k = dsearchn(GridCoords, QueryCoords);  

%k gives the row number in GridCoords that corresponds to the cell that
%each query point is in. To work back to m and n we'll remember that each
%column of GridCoords is effectively a linear-indexed version of the
%original grid matrices. So we can use ind2sub to get the coordinates back.

xDim = size(gX,1);
yDim = size(gX,2);

siz = [xDim yDim];
[qM,qN] = ind2sub(siz, k);
    
Mout = qM;
Nout = qN;

end
