%% CE 295 - Team 4 Project
%   Minimizing EV Charger Driving Distance: Demand Node Parameters
%   Author: Bhumi Tandel
%   Prof. Arnold
%   last updated: 05/04/2023

% dem_param.m

function [] = dem_param(p1, L, row, col, save_file)
%FUNCTION: 
         % Computes and saves all parameters related to the demand nodes k = 1, 2, ...K
         % For each node k, parameters quantify EV adoption rate, land use zoning metric, 
         % residual demand
         % *** Must replace Google API key to calculate distances
         % *** Sample demand node inputs to be passed to dem_node (Richmond):
         % p1 = [37.9556027777778, -122.4162583333333]; L = 2000; row = 3; col = 5;
         % ^creates 3x5 grid of demand nodes where p1 is the coordinate of
         % the top left demande node (each grid square is 2x2 km)

%INPUTS: 
         % p1: demand node coordinates [decimal degrees]
         % L: grid spacing [meters]
         % row: Number of rows in the demand node grid
         % col: Number of columns in the deamnd node grid
         % save_file: .mat file name to save demand parameters under

%OUTPUTS:
         % Outputs saved in the saved .mat file: 'dem_param/mat'
         % K: total number of demand nodes
         % dL: coodrinates for each demand node (lat, lon) [decimal degrees]
         % char_dist: distance between existing chargers and demand nodes [km]
         % zipC: zip code at demand node k
         % u_k: EV adoption rate at demand node k [%]
         % v_k: customer class percentage at deamnd node k [%]
         % R_k: residual demand at station k [veh/day]
         % m: number of vehicles supported by a level 2 charger per day [veh/day]

%% PARAMETERS

key = 'PROVIDE YOUR OWN API KEY HERE'; %provide unique API key to use Google distance & geocoding APIs

[dL, cL] = dem_node(p1, L, row, col); %generate coordiantes for entire grid
K = row*col; %total number of demand nodes

%% CENTROID ZIP CODES

%create vector of zipcodes associated with each demand node
zipC = nan(K, 1); % dL = demand node centroids
for k = 1:K
    dem_coord = [num2str(dL(k,1)) ',' num2str(dL(k,2))];
    url = ['https://maps.googleapis.com/maps/api/geocode/json?address=', dem_coord, '&key=', key];
    str = urlread(url);
    z = extractBetween(str, ' CA ', ', USA');
    zipC(k, 1) = str2double(z(1,:));
end

%% CENTROID EV ADOPTION RATE

%load zip code data (2021 PHEV+BEV adoption rates by zip code)
CEC_data = readmatrix('EV Percent.csv', 'Range', 'A2:B47');

%generate percentage of EVs on the road for each demand node based on zip code
u_k = nan(size(zipC));
for k = 1:size(u_k, 1)
    ind = find(CEC_data(:, 1) == zipC(k));
    u_k(k) = CEC_data(ind, 2);
end

%% LAND USE CHARACTERIZATION

%load 2022 PG&E customer data
%column 1 = zipcode, column 2 = month, column 4 = customer type, column 6 = # of customers
Q1 = readtable('Q1.csv'); Q2 = readtable('Q2.csv'); Q3 = readtable('Q3.csv'); Q4 = readtable('Q4.csv');

%combine data for all 4 quarters
Q_all = [Q1; Q2; Q3; Q4];

%clean data (remove unnecessary columns and only keep zip code, month, customer type, and # of customers)
Q_all(:, 3) = []; Q_all(:, 4) = []; Q_all(:, 5:6)=[];
Q_all = table2cell(Q_all);

%convert customer type to numerical variable
% 1 = residential, 2 = commerical, 3 = agricultural, 4 = industrial
class = nan(size(Q_all, 1), 1);
i = 1;
while i <= length(class)
    if contains(Q_all(i, 3), 'Residential', 'IgnoreCase', true) == 1
        class(i) = 1;
    elseif contains(Q_all(i, 3), 'Commercial', 'IgnoreCase', true) == 1
        class(i) = 2;
    elseif contains(Q_all(i, 3), 'Agricultural', 'IgnoreCase', true) == 1
        class(i) = 3;
    elseif contains(Q_all(i, 3), 'Industrial', 'IgnoreCase', true) == 1
        class(i) = 4;
    else 
        Q_all(i, :) = []; class(i) = []; %delete any unclassified rows
    end
    i = i+1;
end

%construct Q_mat to be a copy of Q_all with num type variables
Q_mat = cell2mat(Q_all(:,1:2)); Q_mat(:,3) = class; Q_mat(:,4) = cell2mat(Q_all(:,4));

%rearrange data in Q_3D as (month, customer class, zip code)
%zip code index corresponds to the zip code stored in ZCode under the same index value
ZCode = unique(Q_mat(:,1));
Q_3D = nan(12, 4, length(ZCode));
for z = 1:length(ZCode)
    for i = 1:size(Q_mat, 1)
        if ZCode(z)== Q_mat(i, 1)
            Q_3D(Q_mat(i,2), Q_mat(i, 3), z) = Q_mat(i, 4);
        end
    end
end

%calculate the class percent for each month in each zipcode (ignore nans)
%i fclass percent = nan, indicates 0 customers in that zip code or no data
Q_3D(:, 5, :) = nan();
for i = 1:size(Q_3D, 3)
    for m = 1:12
        Q_3D(m, 5, i) = sum(Q_3D(m, 1:2, i), 'omitnan')/sum(Q_3D(m, 1:4, i), 'omitnan')*100;
    end
end

%average the percentages for each zip code across all 12 months
class_per = [];
for i = 1:length(ZCode)
    class_per(end+1, 1) = mean(Q_3D(:,5,i), 'omitnan');
end

%remove any zip codes with no data
cust_data = [ZCode class_per];
cust_data(find(isnan(cust_data(:,2))), :) =[];

%generate percentage of residential/commercial customers for each demand node based on zip code
%all areas are 100% residential/commercial for this example
v_k = nan(size(zipC)); 
for k = 1:size(u_k, 1)
    ind = find(cust_data(:, 1) == zipC(k)); 
    v_k(k) = cust_data(ind, 2); 
end

%% AVERAGE TRAFFIC LEVELS

m = 5; %maximum vehicle charging cycle for level 2 charger
m_DC = 50; %maximum vehicle charging cycle for DC fast charger

%load Alameda County AADT counts & lat, long cooridnates
AADT = readmatrix('Traffic AADT.csv', 'range', 'N79:P96'); %Richmond

avg_AADT = AADT(:, 1); %AADT values for each post mile station from Caltrans
latY = AADT(:, 2); %post mile latitude
lonX = AADT(:, 3); %post mile longitude

%interpolate to estimate AADT at each demand node 
%nearest neighbor interpolation provided most reasonable results (for Richmond)
T_k = griddata(lonX, latY, avg_AADT, dL(:, 2), dL(:, 1), ['natural']);
T_k = reshape(T_k, [row, col]);

for i = row:-1:1 %index from last coordinate bc closest to original post mile points > smaller error from interpolation
    for j = col:-1:1
        if isnan(T_k(i,j)) %if no interpolated value for T_k, estimate value based on AADT of surounding neighbors
            [neigh_X, neigh_Y] = meshgrid(i-1:i+1, j-1:j+1);
            neigh_X=neigh_X'; neigh_Y=neigh_Y'; %grid of all possible surrounding point indices
            neigh = [reshape(neigh_X, [9 1]) reshape(neigh_Y, [9 1])];
            s = 1;
            while s <= size(neigh, 1)
                if neigh(s, 1) > row | neigh(s, 1) < 1 | neigh(s, 2) > col | neigh(s, 2) < 1 %remove indices that are out of bounds
                   neigh(s, :) = []; 
                else 
                    s = s+1;
                end
            end
            neigh_val = [];
            for k = 1:1:size(neigh, 1)
                neigh_val(end+1) = T_k(neigh(k, 1), neigh(k, 2)); %create vector of neighboring AADT values
            end
            T_k(i,j) = mean(neigh_val, 'omitnan'); %replace nan with average of neighboring traffic estimates
        end
    end
end

T_k = reshape(T_k, [K, 1]);
T_k = T_k.*(u_k/100).*(v_k/100);

%load existing charger coordinates
char_coord = readmatrix('Existing Public Chargers.csv', 'range', 'U283:V314');%Richmond
char_num = readmatrix('Existing Public Chargers.csv', 'range', 'O283:P314');%Richmond

%translate number of chargers (& type) to number of vehicles served
veh_serv = nan(length(char_num), 1);
for i = 1:length(veh_serv)
    if isnan(char_num(i, 1)) %DC fast charger
        veh_serv(i) = char_num(i, 2)*m_DC;
    elseif isnan(char_num(i, 2))
        veh_serv(i) = char_num(i, 1)*m; %level 2 charger
    else
        veh_serv(i) = char_num(i, 2)*m_DC + char_num(i, 1)*m; %both DC fast charger & level 2 charger
    end
end

mode='driving'; %want to calculate driving distance

%calculate all distances between a demand node k and each existing charger
char_dist = nan(K, length(veh_serv));
for k = 1:K
    orig_coord = [num2str(dL(k,1)) ',' num2str(dL(k,2))];
    char_dist_k = [];
    for i = 1:length(veh_serv)
        dest_coord = [num2str(char_coord(i, 1)) ',' num2str(char_coord(i, 2))];
        url = ['https://maps.googleapis.com/maps/api/distancematrix/json?origins=',orig_coord,'&destinations=',dest_coord,'&mode=',mode,'&language=en-EN&sensor=fa&key=',key];
        str = urlread(url);
        d = extractBetween(str, '"text" : "', ' km"');
        char_dist_k(end+1) = str2double(d);
    end
    char_dist(k,:) = char_dist_k';
end

T_k = reshape(T_k, [row, col]);
char_dist_assign = char_dist; %duplicate copy of char_dist that can be modified to assign chargers to nodes
S_k = zeros(size(T_k));
while prod(prod(isnan(char_dist_assign))) == 0 & sum(sum(S_k < T_k)) ~= 0
    for i = row:-1:1 %index from last coordinate because it's closest to original post mile points > smaller interpolation error
        for j = col:-1:1
            for k = 1:K
                if S_k(i,j) < T_k(i,j)
                    min_d = min(char_dist_assign(k, :));
                    ind = find(char_dist_assign(k, :) == min_d, 1);
                    if isempty(ind) == 0 & S_k(i,j) + (veh_serv(ind)) < T_k(i,j) %only assign charger if demand won't be exceeded
                        S_k(i,j) = S_k(i,j) + (veh_serv(ind));
                        char_dist_assign(:, ind) = nan; %set distance to nan so charger is not reassigned somewhere else
                    end
                end
            end
        end
    end
end

R_k = T_k-S_k; %residual demand = total demand - supported demand

%round demand values to mulitples of m [veh/day]
%if remainder of (demand/m) is > (m/2) i.e. 50% of m, then demand is rounded up to the next multiple of m
R_k = R_k(:) - rem(R_k(:), m) + m*round(rem(R_k(:), m)/m);

%Update demand to reflect areas of ~0 demand (coordinates in the water)
R_k(3) = m; R_k(4) = m; R_k(6) = m; %Richmond

%% SAVE DATA

save(save_file, 'K', 'dL', 'char_dist', 'zipC', 'u_k', 'v_k', 'R_k', 'm');
fprintf('%s saved\n', save_file);

end
