function [ Orientations ] = fnGuessTurbineOrientations( qM, qN, U, V )
%FNGUESSTURBINEORIENTATIONS Takes a list of m,n coordinates and a set of DA
%velocities, and returns the direction that the current is coming from when
%fastest in each cell. Used to guess at orientations for turbines.
%   Orientations: 1xk matrix of headings in radians, measured clockwise
%   from north.
%   qM, qN: Coordinates of k cells for which results are wanted.
%   U, V: 3d matrices giving u- and v- velocities for all points and all
%   time steps to consider. Dimensions: time, m, n.

% Copyright Simon Waldman / Heriot-Watt University, 2014-2015

%FIXME rather than direction at moment of highest speed, should really use
%direction (direction bin?) of highest energy over the full period.
%(or principal component from u & v?)

if nargin < 4
    error('Not enough arguments.');
end
if length(qM) ~= length(qN)
    error('Lengths of qM and qN must match.');
end
if ~all(size(U)==size(V))
    error('Sizes of U and V must match.');
end
%FIXME more checking of inputs - e.g. how many dimensions, interger values
%for m/n, etc.

NumQueries = length(qM);
Orientations = nan(1,NumQueries);   %preallocate.

for c = 1:NumQueries
    m = qM(c);
    n = qN(c);
    Speeds = sqrt( (U(:,m,n)).^2 + (V(:,m,n)).^2 ); %speed at every time step
    [~, I] = max(Speeds);   %find the index (timestep) of max speed
    dir = atan2(U(I,m,n), V(I,m,n));    % direction the flow is going *to*, in radians anticlockwise from north, -pi:pi. (north, not east, as I've done atan2(x,y) rather than (y,x))
    dir2 = -dir + pi; %direction the flow is coming from, in radians clockwise from north
%   u(c) = U(I,m,n); %for debugging. Allows a quiver of the directions.
%   v(c) = V(I,m,n);
    Orientations(c) = dir2;      
end

end


