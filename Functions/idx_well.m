function [WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well
%IDX_WELL   Defines constants for named column indices to well matrix.
%   Example:
%
%   [WELL_NODE, G, PW, GMAX, GMIN, G_COST, WELL_STATUS] = idx_well;
% 
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    g = well(4, G);     % get the gas injection of well 4
%    well(:, GMIN) = 0;  % set to zero the minimum gas injection limit of
%                           all wells
% 
%   The index, name and meaning of each column of the well matrix is given
%   below:
%
%   columns 1-7 must be included in input matrix (in case file)
%    1  WELL_NODE       node number      
%    2  G               gas injection (MMSCFD)
%    3  PW              known pressure at well (psi)
%    4  GMAX            maximum gas output (MMSCFD) 
%    5  GMIN            minimum gas output (MMSCFD) 
%    6  WELL_STATUS     well status (not currently used)
%    7  COST_G          gas production cost ($/MMSCFD)
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

%% define the indices
WELL_NODE	= 1;    %% node number      
G           = 2;    %% gas injection (MMSCFD)
PW          = 3;    %% known pressure at well (psi)
GMAX        = 4;    %% maximum gas output (MMSCFD) 
GMIN        = 5;    %% minimum gas output (MMSCFD) 
WELL_STATUS	= 6;    %% well status (not currently used)
COST_G      = 7;    %% gas production cost ($/MMSCFD)
