function [WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well
%IDX_WELL   Defines constants for named column indices to well matrix.
%   Example:
%
%   [WELL_NODE, G, PW, GMAX, GMIN, G_COST, WELL_STATUS] = idx_well
% 
%   The index, name and meaning of each column of the well matrix is given
%   below:
%
%   columns 1-21 must be included in input matrix (in case file)
%     WELL_NODE	  = 1;    %% Well node
%     G           = 2;    %% Gas injection( )
%     PW          = 3;    %% Known pressure at well( )
%     GMAX        = 4;    %% Gmax, maximum gas output ( )
%     GMIN        = 5;    %% Gmin, minimum gas output ( )
%     WELL_STATUS = 6;    %% well status
%     COST_G      = 7;    %% gas production cost
% 
%
%   See also DEFINE_CONSTANTS_GAS.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license 

%% define the indices
WELL_NODE	= 1;    %% Well node
G           = 2;    %% Gas injection( )
PW          = 3;    %% Known pressure at well( )
GMAX        = 4;    %% Gmax, maximum gas output ( )
GMIN        = 5;    %% Gmin, minimum gas output ( )
WELL_STATUS	= 6;    %% well status
COST_G      = 7;    %% gas production cost
