function [F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe
%IDX_PLINE   Defines constants for named column indices to pipeline matrix.
%   Example:
%
%   [F_NODE, T_NODE, C_O, K_O, FMAX_O, COST_O, FG_O] = idx_pipe
% 
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    branch(4, BR_STATUS) = 0;              % take branch 4 out of service
%    Ploss = branch(:, PF) + branch(:, PT); % compute real power loss vector
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-9 must be included in input matrix (in case file)
% 
%     F_NODE      = 1;    %% f, from node number
%     T_NODE      = 2;    %% t, to node number
%     FG_O        = 3;    %% gas flow in pipeline
%     K_O         = 4;    %% kewmouth constant
%     DIAM        = 5;    %% diameter [inch]
%     LNG         = 6;    %% length [km]
%     FMAX_O      = 7;    %% maximun gas flow
%     FMIN_O      = 8;    %% maximun gas flow
%     COST_O      = 9;    %% cost of gas transport
%
%   See also DEFINE_CONSTANTS_GAS.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license 

%% define the indices
F_NODE      = 1;    %% f, from node number
T_NODE      = 2;    %% t, to node number
FG_O        = 3;    %% gas flow in pipeline
K_O         = 4;    %% kewmouth constant
DIAM        = 5;    %% diameter [inch]
LNG         = 6;    %% length [km]
FMAX_O      = 7;    %% maximun gas flow
FMIN_O      = 8;    %% maximun gas flow
COST_O      = 9;	%% cost of gas transport
