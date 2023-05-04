%% CE 295 - Team 4 Project
%   Minimizing EV Charger Driving Distance: Optimization Solver
%   Author: Bhumi Tandel
%   Prof. Arnold
%   last updated: 05/04/2023

% opt.m

function [N, Y, D] = opt(dem_P, sta_P, a, W)
%FUNCTION: 
         % Optimization formulation to solve for the number of chargers to
         % install and their assignment to different demand nodes

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
         % D: objective function value (minimized average driving distance) [km driven/veh]

%% PARAMETERS
load(dem_P); %load relevant demand node parameters
load(sta_P); %load relevant station parameters

R_k_vec = []; %create residual demand vector that's the same size as y (repeated demand values)
for k = 1:K
    R_k_vec(end+1:end+J, 1) = R_k(k);
end

%% OPTIMIZATION FORMULATION

cvx_clear
cvx_begin
cvx_solver('gurobi');
cvx_solver_settings('NonConvex', 2);
    variables s(J, 1) p(J, 1) q(J, 1) P(J, J) Q(J, J) V(J) L(J, J) %transmission grid variables
    variable y(J*K, 1) binary %1 if station j is serving demand node k, 0 otherwise
    variable n(J, 1) integer %number of chargers
    variable z(J, J)
    minimize((dist./R_k_vec)'*y)
    
    subject to
    %Transmission Grid Constraints
        % Boundary condition for power line flows
        P( 1 , 1 ) == 0;
        Q( 1 , 1 ) == 0;
        
        % Boundary condition for squared line current
        L( 1 , 1 ) == 0;
        
        % Fix node 0 voltage to be 1 "per unit" (p.u.)
        V(1) == 1;
        
        % Force y_jk to be 0 for trivial demand nodes (nodes in water)
        for kk = 1:K
            if R_k(kk) == m
                y(J*(kk-1)+1 : kk*J, 1) == zeros(size(J, 1)); %set zeros for trivial node k and all stations j
            end
            sum(y(J*(kk-1)+1 : kk*J, 1)) <= 3 %limit number of chargers assigned to a demand node to 2
        end
        
        for jj = 1:J  %Loop over each node

            i = rho(jj); %identify parent node

            %EV chargers installed cannot exceed 10% of the number of parking spaces
            n(jj) <= a*c_j(jj)

            %non-negativity constraints
            n(jj) >= 0
            p(jj) >= 0
            q(jj) >= 0

            %no station should be underutilized
            sum(y(jj:J:end).*R_k) >= n(jj)*m %need this to avoid unassigned chargers
            %indexing > single station j serving all demand nodes k

            % Line Power Flows
            P(i,jj) == (t*n(jj)-p(jj))+sum(A(jj,:).*P(jj, :));
            Q(i,jj) == (0.1*t*n(jj)-q(jj))+sum(A(jj,:).*Q(jj, :));
            
            % Nodal voltage
            V(jj) == V(i) + (r(i,jj)^2+x(i,jj)^2)*L(i,jj)-2*(r(i,jj)*P(i,jj)+x(i,jj)*Q(i,jj))
            
            % Squared current magnitude on lines
            L(i,jj) >= quad_over_lin(P(i,jj), V(jj)) + quad_over_lin(Q(i,jj), V(jj))

            % Compute apparent power from active & reactive power
            norm([p(jj), q(jj)]) <= s(jj)

        end

        %minimum number of chargers to be installed
        sum(n) >= W %need this to be here to avoid trivial solution, 0 chargers

        %assignment constraint
        sum(y.*R_k_vec) >= sum(n)*m %need this to avoid assigned demand exceeding charger capacity
       
        % Squared line current limits
        L <= I_max.^2
            
        % Nodal voltage limits
        v_min^2 <= V 
        V <= v_max^2
        
        % Apparent Power Limits
        s <= S_max

cvx_end

N = round(n); % round to integer values if optimization returns floats
Y = reshape(round(y),[J K]); % round to binary values if optimization returns floats
D = cvx_optval; %objective function value

end