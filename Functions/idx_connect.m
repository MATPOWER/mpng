function [GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect
%IDX_CONNECT   Defines constants for named column indices to connect struct.
%   Example:
%
%   [GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect;
%       
%    energy = connect.energy(:, MAX_ENER);  % get the maximum energy vector cost
%    connect.interc.term(3, EFF) = 8e-3;    % set to 8e-3 the efficiency of
%                                             the third gas-fired generator
% 
%   The index, name, and meaning of the colums in the diffent matrixes
%   of the "connect" struct are given below:
%
%    maximun daily energy for power plants
%     1  GEN_ID         generator index          
%     2  MAX_ENER       maximum energy available (MWh/d)
% 
%    interconection through compressors
%     1  COMP_ID        compressor index          
%     2  BUS_ID         bus index, where the compressor is connected
% 
%    interconection through gas-fired power plants
%     1  GEN_ID         generator index (already named)
%     2  NODE_ID        node index, where the plant is located
%     3  EFF            plan efficiency (MMSCFD/MW)
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

%% Max daily energy for power plants
GEN_ID      = 1;    % Generator ID
MAX_ENER	= 2;    % Max energy available

%% Interconection through compressors
COMP_ID     = 1;    % Compressor ID
BUS_ID      = 2;    % Bus ID

%% Interconection through gas-fired power plants
% GEN_ID      = 1;  % Generator ID (already called)
NODE_ID     = 2;    % Node where plant is located
EFF         = 3;    % Plant efficiency
