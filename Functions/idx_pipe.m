function [F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe
%IDX_PIPE   Defines constants for named column indices to pipeline matrix.
%   Example:
%
%   [F_NODE, T_NODE, C_O, K_O, FMAX_O, COST_O, FG_O] = idx_pipe;
% 
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    pipe(4, COST_O) = 0;   % set to zero the pipe 4 transport cost
%    k = pipe(:, K_O)       % get weymouth constant vector
% 
%   The index, name and meaning of each column of the branch matrix is given
%   below:
%
%   columns 1-9 must be included in input matrix (in case file)
%    1  F_NODE      f, from node number
%    2  T_NODE      t, to node number
%    3  FG_O        gas flow in pipeline (MMSCFD)
%    4  K_O         kewmouth constant 
%    5  DIAM        diameter (inch)
%    6  LNG         length (km)
%    7  FMAX_O      maximum gas flow (MMSCFD)     
%    8  FMIN_O      minimum gas flow (MMSCFD)     
%    9  COST_O      cost of gas transport ($/MMSCFD)
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
F_NODE      = 1;    %% f, from node number
T_NODE      = 2;    %% t, to node number
FG_O        = 3;    %% gas flow in pipeline (MMSCFD)
K_O         = 4;    %% kewmouth constant 
DIAM        = 5;    %% diameter (inch)
LNG         = 6;    %% length (km)
FMAX_O      = 7;    %% maximum gas flow (MMSCFD)     
FMIN_O      = 8;    %% minimum gas flow (MMSCFD)     
COST_O      = 9;	%% cost of gas transport ($/MMSCFD)
