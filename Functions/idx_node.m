function [DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
            COST_UNP, GD, NGD] = idx_node
% IDX_NODE  Defines constants for named column indices to node info matrix.
%   Example:
%
%   [DEM, WELL, NODE_I, NODE_TYPE, P, PMAX, PMIN, OVP, UNP, COST_OVP, ...
%             COST_UNP, GD, NGD, GD_N] = idx_node;
%
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    gd = node.info(4, GD);     % get the total gas demand at node 4
%    node.info(:, PR) = 1200;   % set the max pressure to 1200 at all nodes
%
%   The index, name and meaning of each column of the node info matrix is given
%   below:
%
%   columns 1-11 must be included in input matrix (in case file)
%    1  NODE_I      node number (positive integer)                
%    2  NODE_TYPE   node type (1 = DEM, 2 = WELL)
%    3  PR          pressure (psi)
%    4  PRMAX       maximum pressure (psi)
%    5  PRMIN       minimum pressure (psi)
%    6  OVP         overpressure (psi)
%    7  UNP         underpressure (psi)
%    8  COST_OVP    overpressure cost ($/psi^2)
%    9  COST_UNP    underpressure cost ($/psi^2)
%    10 GD          total gas demand (MMSCFD)
%    11 NGD         number of different type of gas demands (positive integer)
% 
%   additional constants, used to assign/compare values in the NODE_TYPE column
%    1  DEM     Gas demand node
%    2  WELL    Gas well node
%
%   See also DEFINE_CONSTANTS_GAS.

%   MPNG: MATPOWER - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

%% define node types
DEM         = 1;
WELL        = 2;

%% define the indices
NODE_I      = 1;    %% node number 
NODE_TYPE	= 2;    %% node type (1 - demand node, 2 - extraction node)
PR          = 3;    %% pressure (psi)
PRMAX       = 4;    %% maximun pressure (psi)
PRMIN       = 5;    %% minimun pressure (psi)
OVP         = 6;    %% overpressure (psi)
UNP         = 7;    %% underpressure (psi)
COST_OVP	= 8;    %% overpressure cost  ($/psi^2)
COST_UNP	= 9;    %% underpressure cost  ($/psi^2)
GD          = 10;	%% total gas demand (MMSCFD) 
NGD         = 11;	%% number of different type of demands (positive integer)


