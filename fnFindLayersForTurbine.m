function [ LayersIntersected ] = fnFindLayersForTurbine( HubElevation, Diameter, SeabedElevation, NumLayers )
%FNFINDLAYERSFORTURBINE Finds which sigma layers of a model a turbine rotor
%intersects.
%   HubElevation is the elevation of the centre of the turbine with respect
%   to mean sea level. Metres. Should be negative.
%   Diameter is obvious. In metres.
%   SeabedElevation: Elevation of seabed wrt MSL. Ie depth, but negative.
%   NumLayers: Number of equal-spaced sigma layers in the model
%
%   Output: LayersIntersected: Vector listing the layer numbers that are
%   intersected (top layer = 1).

% Copyright Simon Waldman / Heriot-Watt University, 2014-2015

if nargin < 4
    error('Not enough arguments.');
end
if Diameter <= 0
    error('Diameter must be positive');
end
if NumLayers < 1 | round(NumLayers) ~= NumLayers
    error('NumLayers must be a positive integer'); %and an integer!
end
if HubElevation + Diameter / 2 > 0
    error('Turbine will stick out of the water!');
end
if HubElevation - Diameter / 2 < SeabedElevation
    %error('Turbine will hit the seabed!');
    LayersIntersected = [];
    return;
end

LayerBoundaries = linspace(0,SeabedElevation,NumLayers+1); 
% So the boundaries of layer n will be boundaries n and n+1. (top & bottom
% respectively)

%FIXME floating point equalities here.
TopLB = max(find(LayerBoundaries >= (HubElevation + Diameter / 2)));
BottomLB = min(find(LayerBoundaries <= (HubElevation - Diameter / 2)));

% So which layer is TopLB the top of, and which is BottomLB the bottom of?
TopL = TopLB;
BottomL = BottomLB - 1;

%output
LayersIntersected = TopL:BottomL;

end

