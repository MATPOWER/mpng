function [DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
            COST_UNP, GD, NGD] = idx_node
% IDX_NODE   Defines constants for named column indices to node matrix.
%   Example:
%
%   [DEM, WELL, NODE_I, NODE_TYPE, P, PMAX, PMIN, OVP, UNP, COST_OVP, ...
%             COST_UNP, GD, NGD, GD_N] = idx_node
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    Pd = bus(4, PD);       % get the real power demand at bus 4
%    bus(:, VMIN) = 0.95;   % set the min voltage magnitude to 0.95 at all buses
%
%   The index, name and meaning of each column of the bus matrix is given
%   below:
%
%   columns 1-13 must be included in input matrix (in case file)
%   NODE_I      = 1;    %% node number (1 to 29997)
%   NODE_TYPE	= 2;    %% node type (1 - demand node, 2 - extraction node)
%   PR          = 3;    %% pressure ( )
%   PRMAX       = 4;    %% Maximun pressure ( )
%   PRMIN       = 5;    %% Minimun pressure ( )
%   OVP         = 6;    %% overpressure
%   UNP         = 7;    %% underpressure
%   COST_OVP	= 8;    %% overpressure cost
%   COST_UNP	= 9;    %% underpressure cost
%   GD          = 10;	%% Total gas demand 
%   NGD         = 11;	%% Number of different type of demands (positive integer)
% 
%   additional constants, used to assign/compare values in the BUS_TYPE column
%    1  DEM     Gas demand node
%    2  WELL    Gas well node
%
%   See also DEFINE_CONSTANTS_GAS.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license 

%% define node types
DEM         = 1;
WELL        = 2;

%% define the indices
NODE_I      = 1;    %% node number (1 to 29997)
NODE_TYPE	= 2;    %% node type (1 - demand node, 2 - extraction node)
PR          = 3;    %% pressure ( )
PRMAX       = 4;    %% Maximun pressure ( )
PRMIN       = 5;    %% Minimun pressure ( )
OVP         = 6;    %% overpressure
UNP         = 7;    %% underpressure
COST_OVP	= 8;    %% overpressure cost
COST_UNP	= 9;    %% underpressure cost
GD          = 10;	%% Total gas demand 
NGD         = 11;	%% Number of different type of demands (positive integer)


