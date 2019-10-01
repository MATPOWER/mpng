function [COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp
%IDX_COMP   Defines constants for named column indices to compressor matrix.
%   Example:
%
%   [COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
%     GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp
% 
%   The index, name and meaning of each column of the compressor matrix is given
%   below:
%
%   Compressors type:
% 
%   COMP_G      = 1;    %% compressor working with gas
%   COMP_P      = 2;    %% compressor working with power
% 
%   Define the indices:
%
%   F_NODE      = 1;    %% f, from node number
%   T_NODE      = 2;    %% t, to node number
%   TYPE_C      = 3;    %% type of compressor
%   FG_C        = 4;    %% gas flow through compressor  \ 
%   PC_C        = 5;	%% power consumed by compressor  | -> VARIABLES    
%   GC_C        = 6;	%% gas consumed by compressor   /
%   RATIO_C     = 7;    %% compressor ratio
%   B_C         = 8;    %% B_C
%   Z_C         = 9;    %% Z_C
%   ALPHA       = 10;	%% alpha \    
%   BETA        = 11;	%% beta   | --> gas compsumption parameter of gas-fired compresor 
%   GAMMA       = 12;	%% gamma /
%   FMAX_C      = 13;   %% maximun flow trhough compressors
%   COST_C      = 14;   %% cost of compressor, depended of the gas flow
% 
%   See also DEFINE_CONSTANTS_GAS.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license     
%% compressor type

COMP_P      = 1;    %% compressor working with power
COMP_G      = 2;    %% compressor working with gas


%% define the indices
% COMP        = 1;    %%      
F_NODE      = 1;    %% f, from node number
T_NODE      = 2;    %% t, to node number
TYPE_C      = 3;    %% type of compressor
FG_C        = 4;    %% gas flow through compressor  \ 
PC_C        = 5;	%% power consumed by compressor  | -> VARIABLES    
GC_C        = 6;	%% gas consumed by compressor   /
RATIO_C     = 7;    %% compressor ratio
B_C         = 8;    %% B_C
Z_C         = 9;    %% Z_C
ALPHA       = 10;	%% alpha \    
BETA        = 11;	%% beta   | --> gas compsumption parameter of gas-fired compresor 
GAMMA       = 12;	%% gamma /
FMAX_C      = 13;   %% maximun flow trhough compressors
COST_C      = 14;   %% cost of compressor, depended of the gas flow

