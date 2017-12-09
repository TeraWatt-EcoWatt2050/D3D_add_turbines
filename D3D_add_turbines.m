% Insert tidal turbine representations into delft3D model. Requires
% Delft3D-MATLAB toolbox.
%Initial bodgy script by Simon Waldman, Jan 2015.

%% Variables

UseCorrection = false;   % whether to correct for free-stream vs cell velocity. This is only an approximation, and is untested in arrays.

TurbineInfoFile = 'all_longlat.csv'; %csv file as described in fnReadTurbinesFile.m
TurbineFileSkip = 1; %lines of header to skip in this file.

BathyFile = 'D:\temp\initial bed level.mat';    %.mat file exported from Quickplot "initial bed level", using an output of the model without turbines
VelocityFile = 'D:\temp\depth averaged velocity.mat'; %.mat file exported from Quickplot "depth averaged velocity", using an output of the model without turbines.
% unless turbine orientations are set in the CSV TurbineInfo file, they
% will be set from this using the moment of peak velocity in each cell. NB
% use a short period around springs, don't try to export entire runs...
% 1000 time steps is already pushing things on a largeish grid, unless you
% have lots of RAM.

OutputFilename = '1000_turbines_uncorr_rightdepths.ppl';

%grid info
GridFile = 'Orkney.grd';    % D3D grid file.
NumLayers = 10; % sigma layers in model

% Turbine specs
Ct = 0.85;
Diameter = 20;  %metres.

% Other
TidalRange = 4;     %metres. The rotor will be kept half of this value beneath MSL.

%derived stuff
RotorArea = (Diameter / 2)^2 * pi;

%% Read turbine info

disp('Reading turbine info from CSV...');
TurbineList = fnReadTurbinesFile(TurbineInfoFile, TurbineFileSkip);

NumTurbines = length(TurbineList.x);

%% Read grid

disp('Reading grid...');
G = wlgrid('read', GridFile);

%% Find the grid cell that each turbine lies in

disp('Finding grid cells for turbines...');
[ TurbineList.M, TurbineList.N ] = fnFindGridCellsForPoints(G.X, G.Y, TurbineList.x, TurbineList.y);
% The fnFindGrdsCellForPoints function is only valid for rectilinear grids
% - but this is rectilinear in spherical coords.

%% Read the bathymetry

   disp('Loading BathyFile and positioning turbines vertically...');
   
   load(BathyFile);
   if ~exist('data', 'var')
       error('Error reading BathyFile.');
   end
   if ~isfield(data, 'Val')
       error('Format or contents of BathyFile are incorrect.');
   end
   if size(data.Val) ~= size(G.X) + 1
       error('Grid in BathyFile is not the same grid as that in GridFile.');
   end
   
   for t = 1:NumTurbines        %loop through turbinelist.
%        m = TurbineList.N(t);
%        n = TurbineList.M(t);  %THIS WAS WRONG! There's no need to swap
%                               m and n.
       
%        CellBathy = data.Val(m,n);
       CellBathy = data.Val(TurbineList.M(t),TurbineList.N(t));
       TurbineList.Bathy(t) = CellBathy;
       
       z = CellBathy + TurbineList.dz(t);   % dz should be distance from the seabed.
       if z > -(Diameter/2 + TidalRange/2); %if the turbine will stick out of the water 
           z = -(Diameter/2 + TidalRange/2);  %push it down so it's just under the water.
       end
       TurbineList.z(t) = z;
   end

   % if any turbine coordinates were where D3D thinks there is land, they
   % will now have NaN for z-values. Hmm.
   
    clear('data'); 


%% If we haven't been provided orientations for the turbines, guess at them from the direction of peak flow
% NB I think this will give incorrect angles if the M and N grid axes do
% not align with east and north.

if isempty(TurbineList.o)
    disp('Loading VelocityFile and guessing at turbine orientations...');
    load(VelocityFile); %this may be big. Havne't tried with anything really big yet.
    if ~exist('data', 'var')
        error('Error reading VelocityFile.');
    end
    if ~all(isfield(data, {'XComp' 'X'}))
        error('Format or contents of VelocityFile are incorrect.');
    end
    if size(data.X) ~= size(G.X) + 1
        error('Grid in VelocityFile is not the same size as that in GridFile.');
    end
    TurbineList.o = fnGuessTurbineOrientations(TurbineList.M, TurbineList.N, data.XComp, data.YComp);    
    clear('data');  %get rid of that rather large dataset that we don't need any more.
else % ie if we have orientations given in the CSV.
    disp('Using turbine orientations from CSV...');
    TurbineList.o = TurbineList.o * 2 * pi / 360;   %convert degrees in the CSV to radians.
end

%% Set turbine-related drag coeff for each cell in each direction

disp('Finding Closs coefficients for cells with turbines...');
% loop through all grid cells, do work on each one
xDim = size(G.X, 1) - 1;    %don't use the last cell in each direction. we want to be able to look at the distance between any grid location and the one after it. And we shouldn't have turbines at the edge anyway.
yDim = size(G.X, 2) - 1;
numCells = xDim * yDim;

Ax = zeros(xDim, yDim);  %preallocate arrays for rotor effective areas. 
Ay = zeros(xDim, yDim);  % these are the total area of rotor in the horizontal cell *over the full vertical height* that shows when viewed in the x and y directions.
Layers = cell(xDim, yDim);  %the layers that turbines in this cell occupy
HasTurbines = false(xDim, yDim);    %logical matrix showing which cells have turbines in
NuX = zeros(xDim,yDim); % see position paper for explanation of nu and how this relates to Closs (FIXME CITATION)
NuY = zeros(xDim,yDim);
Clossx = zeros(xDim,yDim);
Clossy = zeros(xDim,yDim);
tooshallow = 0;
for m = 1:xDim
    for n = 1:yDim
        cellNo = (m-1)*yDim + n;
        %fprintf('Working on cell %i of %i (%.2f%%). ', cellNo, numCells, (cellNo/numCells * 100));
        %find which (if any) turbines are in this cell.
        TinC = find(TurbineList.M==m & TurbineList.N==n)';  % yay for random transpositions necessary to make matlab behave. Without this, the for loop in a few lines doesn't work.
        %fprintf('%i turbines in cell.\n', length(TinC));
        if ~isempty(TinC) % if there are turbine(s) in this cell (this if may be unnecessary - could leave it to a for loop that never executes)
            HasTurbines(m,n) = true;
            for t = TinC
                % we want to calculate the total "effective area" of rotor
                % in each direction in this cell.
                Ax(m,n) = Ax(m,n) + RotorArea * abs(sin(TurbineList.o(t)));   %abs because we're assuming that turbines work as well backwards as forwards.
                Ay(m,n) = Ay(m,n) + RotorArea * abs(cos(TurbineList.o(t)));
            end 
            % other info that we'll need for this cell in a minute. We'll
            % take mean values from the turbines in the cell.
            Bathy = mean(TurbineList.Bathy(TinC)); %depth to use for this cell
            DeltaZ = -Bathy / NumLayers; %thickness of a vertical layer
            %find which layers the turbines in this cell occupy (we have to
            %use a mean z value for all turbines in the cell, as we can't have more than one overlapping porous
            %plate)
            Layers{m,n} = fnFindLayersForTurbine(mean(TurbineList.z(TinC)), Diameter, Bathy, NumLayers);
  
            if isempty(Layers{m,n})
                fprintf('TOO SHALLOW. %i turbines omitted. m: %i, n: %i.\n', length(TinC), m, n);
                HasTurbines(m,n) = false;   %can we think of something better to do here?
                tooshallow = tooshallow + length(TinC); %add to counter.
                break;
            end
        
            %now convert those areas into the Closs coefficients for each cell.
            
            %We need to know DeltaX and DeltaY in *metres* - which varies
            %across the grid because it's in spherical coordinates. 
            %This would be so much easier if
            %we had the Mapping Toolbox!!! Fortunately I found the
            %functions to do it on the internet ;-)
            DeltaLon = abs(G.X(m,n) - G.X(m+1,n));
            DeltaX = fnLonToM(DeltaLon, G.Y(m,n));
            DeltaLat = abs(G.Y(m,n) - G.Y(m,n+1));
            DeltaY = fnLatToM(DeltaLat, G.Y(m,n));
            
            if UseCorrection % do we correct for free-stream vs cell velocity?
                
                %Calculate the value nu (see position paper for details)
                
                NuX(m,n) = (Ct * Ax(m,n) ) / (DeltaY * DeltaZ * length(Layers{m,n}));  %layers taken into account to split total area of rotors into area per vertical cell covered.
                NuY(m,n) = (Ct * Ay(m,n) ) / (DeltaX * DeltaZ * length(Layers{m,n}));
                
                % calculate Closs
                Clossx(m,n) = (2 * NuX(m,n)) / (1 + sqrt( 1 - NuX(m,n) )).^2;
                Clossy(m,n) = (2 * NuY(m,n)) / (1 + sqrt( 1 - NuY(m,n) )).^2;
            else
                Clossx(m,n) = (Ct * Ax(m,n) ) / (2 * DeltaY * DeltaZ * length(Layers{m,n}));
                Clossy(m,n) = (Ct * Ay(m,n) ) / (2 * DeltaX * DeltaZ * length(Layers{m,n}));
            end
        end
    end
end
fprintf('Total of %i turbines omitted because too shallow.\n', tooshallow);

%% Form and write the ppl file

disp('Writing PPL file...');
FID = fopen(OutputFilename, 'w');
siz = [xDim yDim];
mLen = fix(log10(xDim) + 1);    %number of digits 
nLen = fix(log10(yDim) + 1);

%Loop through all the cells with turbines. We can be clever about this -
%find(HasTurbines) should give linear indices of all of them.
for c = (find(HasTurbines))'
    [m, n] = ind2sub(siz, c);   %get m,n coordinates from the linear index
    fprintf(FID, 'U  %*i  %*i  %*i  %*i  %2i  %2i  %7.5e  .\n', ...
        mLen, m, nLen, n, mLen, m, nLen, n, Layers{m, n}(1), Layers{m,n}(end), Clossx(m,n));
    fprintf(FID, 'V  %*i  %*i  %*i  %*i  %2i  %2i  %7.5e  .\n', ...
        mLen, m, nLen, n, mLen, m, nLen, n, Layers{m, n}(1), Layers{m,n}(end), Clossy(m,n));
end

fclose(FID);

disp('Done!');

%% Plot the horizontal locations of the porous plates

siz = [xDim yDim];
f = figure;
hold on;
cm = flip(parula(101));
maxCloss = max([ max(max(Clossx)) max(max(Clossy)) ]);
set(gcf,'DefaultAxesColorOrder',cm);
for c = (find(HasTurbines))'
    [m, n] = ind2sub(siz, c);   %get m,n coordinates from the linear index
    xcolno = round(100 * Clossx(m,n) / maxCloss) + 1;
    line([m m], [n n+1], 'color', cm(xcolno,:));
    ycolno = round(100 * Clossy(m,n) / maxCloss) + 1;
    line([m m+1], [n n], 'color', cm(ycolno,:), 'LineWidth', 1);
end
colormap(cm);
cb = colorbar;
caxis([0 maxCloss]);
ylabel(cb, 'c loss');
