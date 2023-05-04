%% CE 295 - Team 4 Project
%   Minimizing EV Charger Driving Distance: Station Parameters
%   Author: Bhumi Tandel
%   Prof. Arnold
%   last updated: 05/04/2023

% sta_param.m

function [] = sta_param(mat_file, save_file)
%FUNCTION: 
         % Computes and saves all parameters related to charging stations j = 1, 2, ... J 
         % For each node j, parameters quantify # of charging stalls
         % available, distance between charger j and demand node k
         % Simulates and saves radial transmissions grid parameters by
         % randomly assigning maximum apparent power threshold, resistance,
         % reactance, complex current magnitude to eeach charging station
         % node j = 1, 2, ... J (each node in the simulated transmission grid 
         % coincides with a charging station j)
         % *** Note that transmission grid parameters are randomized and will
         % be different each time this script is run
         % *** Must replace Google API key to calculate distances

%INPUTS: 
         % charging station coordinates (latitude & longitude
         % number of stalls at each station
         % demand node parameters saved in output from dem_param.m

%OUTPUTS:
         % Outputs saved & passed to MILP optimization
         % dist: distance between demand node k and station j
         % c_j: number of parking spots at station j
         % S_max: max apparent power consumption at station j [kVAr]
         % r: resistance on line (i,j)
         % x: reactance on line (i,j)
         % I_max: maximum magnitude of complex current on line (i,j)
         % A: adjacency matrix for simulated transmission grid
         % rho: parent nodes for simulated transmission grid
         % vmin: minimum nodal voltage
         % vmax: maximum nodal voltage
         % t: power demand of a level 2 charger

%% PARAMETERS

%initialize URL request parameters
mode='driving'; %want to calculate driving distance
key = 'PROVIDE YOUR OWN API KEY HERE'; %provide unique API key to use Google distance & geocoding APIs

%load demand node parameters to demand node coordinates
load(mat_file)

%% LOAD STATION INFORMATION

%load station data (Lot ID, Lat, Long, # of spots, binary: 1 if has EV spots already)
station = readmatrix('Parking_data - Copy.csv', 'Range', 'A45:E54'); %Richmond

J = size(station, 1); %J = number of potential charging stations
p_node = station(:, 2:3); %station coordinates
c_j = station(:, 4); %extract # of parking stalls for each station

%% CENTROID TO CHARGER DISTANCES

%store travel distances in 'dist' matrix where rows index stations and columns index demand nodes
dist = nan(J, K);
for k = 1:K
    orig_coord = [num2str(dL(k,1)) ',' num2str(dL(k,2))];
    for j = 1:J
        dest_coord = [num2str(p_node(j,1)) ',' num2str(p_node(j,2))];
        url = ['https://maps.googleapis.com/maps/api/distancematrix/json?origins=',orig_coord,'&destinations=',dest_coord,'&mode=',mode,'&language=en-EN&sensor=fa&key=',key];
        str = urlread(url);
        d = extractBetween(str, '"text" : "', ' km"'); %parse through URL results to extract driving distance
        dist(j, k) = str2double(d);
    end
end

dist = reshape(dist, [J*K, 1]); %reshape to make compatible with cvx optimization code

%Google Distance Matrix API URL Building documentation --> https://developers.google.com/maps/documentation/distance-matrix/distance-matrix

%% STATION TRANISSMION GRID ASSIGNMENT

%load PG&E ICA map data
ICA = readmatrix('ICAMap_Richmond.csv', 'Range', 'E2:E1001'); %Extract PG&E load hosting capacity values to sample from

%randomly assign maximum active power at each demand node based on ICA data
P_max = nan(J, 1);
for j = 1:J
    P_max(j) = ICA(randi(1000,1,1)); %active power [kW]
end

%calculate maximum apparent power assuming that maximum reactive power is 10% of the max active power
S_max = sqrt(P_max.^2 + (0.1*P_max).^2); % [kVA]

%convert from kW to p.u. using following bases:
V_base = 4.17; % [kV]
S_base = 1*10^3; % [kW]
P_base = sqrt(1000^2/1.01); % [kW]

S_max = S_max/S_base;

%min and max nodal voltages for each line (i,j)
v_min = 0.95; v_max = 1.05;

%power drawn from a level 2 charger normalized by P_base
t = 13/P_base; % [kW]

%Extract resistance, reactance, and max current values from HW 3

% r_ij: Resistance [p.u.]
r = readmatrix("Grid_Sim.csv", "range", "B6:N18");

% x_ij: Reactance [p.u.]
x = readmatrix("Grid_Sim.csv", "range", "B24:N36");

% I_max_ij: Maximal line current [p.u.]
I_max = readmatrix("Grid_Sim.csv", "range", "B42:N54");

%generate matrix of potential r, x, I_max values to pull from
grid_rxI = [];
for i = 1:size(r, 1)
    for j = 1:size(r, 2)
        if r(i,j) ~= 0
           grid_rxI(end+1, 1:3) = [r(i,j) x(i,j) I_max(i,j)];
        end
    end
end

%generate random adjacency matrix
A = zeros(J, J);
for j = 2:J % index across all columns except the first one (first node has no parent)
    row = randi(J, 1, 1);
    A(row, j) = 1;
end

%regenerate new resistance, reactance, and max current matrices for new grid by randomly assignning 
%resistance, reactance, and max current values to each station
r = zeros(J, J);
x = zeros(J, J);
I_max = zeros(J, J);
for i = 1:J
    for j = 1:J
        if A(i,j) ~= 0
            row = randi(length(grid_rxI), 1, 1);
            r(i,j) = grid_rxI(row, 1);
            x(i,j) = grid_rxI(row, 2);
            I_max(i,j) = grid_rxI(row, 3);
        end
    end
end

%generate rho matrix based on the adjacency matrix
rho = nan(J, 1);
for i = 1:J
    for j = 1:J
        if A(i,j) ~= 0
            rho(j) = i;
        end
    end
end
rho(1) = 1; %convention: assign node 1 to be the parent of node 1 (power flow from node 1 to node 1 = 0)

%% SAVE DATA

save(save_file, 'J', 'c_j', 'dist', 'S_max', 'r', 'x', 'I_max', 'A', 'rho', 'v_min', 'v_max', 't');
fprintf('%s saved\n', save_file);

end
