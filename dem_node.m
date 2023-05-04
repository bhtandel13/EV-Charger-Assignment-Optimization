%% CE 295 - Team 4 Project
%   Minimizing EV Charger Driving Distance: Demand Node Coordinates
%   Author: Bhumi Tandel
%   Prof. Arnold
%   last updated: 05/04/2023

% dem_node.m

function [dL, cL] = dem_node(p1, d, row, col)
%FUNCTION: 
         % Creates a grid with the specified number of rows and columns
         % Calculates the coordinates of each grid square centroid and the
         % corner points given the coordinates of the top leftmost demand node
         % Coordinates are backcalculated based on the Haversine formula (grid
         % length is given as a great-circle distance)

%INPUTS: 
         % p1 is a two element vector specificying the latitude and
         % longitude of the first demand node in decimal degrees
         % L is the side length of the desired grid in meters
         % row is how many grids to generate in the horizontal direction
         % col is how many grids to generate in the vertical direction

%OUTPUTS:
         % dL is the set of latitudes and longitudes for each demand
         % node centroid in a matrix of size [row*col, 2]
         % cL is analogous dL, and corresponds to the corner points of the 
         % demand node grids
         % There are (row+1)*(col+1) different corner points


%% PARAMETERS
r = 6371000; %radius of the Earth [m]

%% DEMAND NODE CENTROID COORDINATES

%parse out latitude and longitude of top left demand node & convert to radians
lat1 = p1(1); lon1 = p1(2);
lat1 = lat1*pi/180; lon1 = lon1*pi/180;

%calculate the latitude & longitude of the adjacent (right) point in radians
lat2 = d/r+lat1;
lon2 = 2*asin(sin(d/(2*r))/cos(lat1))+lon1;

%calculate the angular rotation needed to achieve desired grid size
d_lat_diff = lat2-lat1;
d_lon_diff = lon2-lon1;

%create vector of demand node centroid latitudes and longitudes
dy = linspace(lat1, lat1-d_lat_diff*(row-1), row)*180/pi;
dx = linspace(lon1, lon1+d_lon_diff*(col-1), col)*180/pi;

%create grid for deamnd node centroids
[dX, dY] = meshgrid(dx, dy);
dX2 = reshape(dX, [row*col 1]);
dY2 = reshape(dY, [row*col 1]);
dL = [dY2, dX2];

%% GRID CORNER COORDINATES

%calculate the latitude & longitude of the top left corner point of first
%demand centroid
lat3 = (d/2)/r+lat1;
lon3 = 2*asin(sin((d/2)/(2*r))/cos(lat1))+lon1;

%calculate the angular rotation needed to move to top left corner
c_lat_diff = lat3-lat1;
c_lon_diff = lon3-lon1;

%calculate starting coordinates of top left corner of the grid
lat0 = lat1 + c_lat_diff;
lon0 = lon1 - c_lon_diff;

%create vector of grid corner point latitudes and longitudes
cy = linspace(lat0, lat0-d_lat_diff*row, row+1)*180/pi;
cx = linspace(lon0, lon0+d_lon_diff*col, col+1)*180/pi;

%create grid for deamnd node corners
[cY, cX] = meshgrid(cy, cx);
cX2 = reshape(cX, [(row+1)*(col+1) 1]);
cY2 = reshape(cY, [(row+1)*(col+1) 1]);
cL = [cY2, cX2];

end