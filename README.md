# D3D_add_turbines
MATLAB script for representing tidal turbines in Delft3D-Flow models by representing them as porous plates.

This script takes a CSV file of turbine locations and produces a .ppl file of porous plate
specifications for Delft3D-Flow.

Further information on this implementation is available in:
Baston S, Waldman S & Side J (2014) “Modelling energy extraction in tidal flows”, Position Paper, 
output of the TeraWatt UKCMER Grand Challenge project. Rev. 3.1 issued 2015. Available at http://www.masts.ac.uk/about/masts-publications/terawatt-publications/

If you use it in a project that leads to a report or publication I would 
appreciate (a) knowing about it; (b) acknowledgement; and, if appropriate, a citation for the following paper
for which it was developed:

Waldman S, Baston S, Nemalidinne R, Chatzirodou A, Venugopal V, Side J, “Implementation of tidal turbines in MIKE 3 and Delft3D models of Pentland Firth & Orkney Waters” (2017)
Ocean & Coastal Management. http://dx.doi.org/10.1016/j.ocecoaman.2017.04.015

Limitations:
- Only valid for rectilinear grids in lat/lon coordinates (that is, rectilinear if unprojected)
- When converting between lat/lon and metres the script makes the approximation that any given cell 
    is rectangular, but allows for the size of that rectangle to be different across the domain. Be careful
    with very large cells, or with very high latitudes.

Improvements to this script are welcome. The latest version may be found at https://github.com/TeraWatt-EcoWatt2050/D3D_add_turbines.

-Simon Waldman.

Address for correspondance or queries: smw13@hw.ac.uk
