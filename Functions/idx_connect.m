function [GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect
%IDX_WELL   Defines constants for named column indices to connect struct.
%   Example:
%
%   [GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect
% 
%   The index, name, and meaning of the colums in the diffent matrixes
%   of the "connect"  struct is given%   below:
%
%       Max daily energy for power plants
%           GEN_ID      = 1;    % Generator ID
%           MAX_ENER	= 2;    % Max energy available
% 
%       Interconection through compressors
%           COMP_ID     = 1;    % Compressor ID
%           BUS_ID      = 2;    % Bus ID
% 
%       Interconection through gas-fired power plants
%                               % Generator ID (already called)
%           NODE_ID     = 2;    % Node where plant is located
%           EFF         = 3;    % Plant efficiency 
%
%   See also DEFINE_CONSTANTS DEFINE_CONSTANTS_GAS.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license 

%% Max daily energy for power plants
GEN_ID      = 1;    % Generator ID
MAX_ENER	= 2;    % Max energy available

%% Interconection through compressors
COMP_ID     = 1;    % Compressor ID
BUS_ID      = 2;    % Bus ID

%% Interconection through gas-fired power plants
                    % Generator ID (already called)
NODE_ID     = 2;    % Node where plant is located
EFF         = 3;    % Plant efficiency
