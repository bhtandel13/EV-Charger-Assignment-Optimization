%% CE 295 - Team 4 Project
%   Minimizing EV Charger Driving Distance: Bisection Search Algorithm
%   Author: Bhumi Tandel
%   Prof. Arnold
%   last updated: 05/04/2023

% opt_eval.m

function [N, Y, W, D] = opt_eval(dem_P, sta_P, a, W_i)
%FUNCTION: 
         % Recursively solves optimization using the bisection method to
         % find the maximum number of chargers that can be installed

%INPUTS: 
         % dem_P: saved demand node parameters .mat variables
         % sta_P: saved station parameters .mat variables
         % a: percentage of stalls that can be allocated to EV's 
         % 0 < a <= 1
         % W_i: initial guess at the maximum number fo chargers that can be installed

%OUTPUTS:
         % N: number of chargers to install at each station j
         % Y: Binary assignment matrix with rows indexing through stations
         % 1,2, ...J and columns indexing through demand nodes 1,2, ...K
         % W: Maximum number of chargers that can be installed for the
         % given percentage of stalls allocated to EVs
         % D: objective function value (minimized average driving distance) [km driven/veh]

%% SEARCH ALGORITHM

b = 0; %solution found indicator
i = 100; %step size (# of chargers)

while b == 0
    [n, y, d] = opt(dem_P, sta_P, a, W_i);
    if sum(isnan(n)) ~= 0
        i = round(i/2);
        W_i = W_i - i;
    else
        W_i = W_i + i;
    end
    if i == 1
       [n, y, d] = opt(dem_P, sta_P, a, W_i);
       [N, Y, D] = opt(dem_P, sta_P, a, W_i-i);
       if sum(isnan(n)) ~= 0 & sum(isnan(N)) == 0
           b = 1;
       end
    end
end
W = W_i-1;

end

