function [COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, ... 
    BETA, GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp
%IDX_COMP   Defines constants for named column indices to compressor matrix.
%   Example:
%
%   [COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, ... 
%     BETA, GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
% 
%   Some examples of usage, after defining the constants using the line above,
%   are:
% 
%    comp(3,TYPE_C) = COMP_P;   % set to power compressor type the
%                                 compressor 3
%    cost_c = comp(:,COST_C);   % get the vector compressor cost of all
%                                 compressors
% 
%   The index, name and meaning of each column of the compressor matrix is given
%   below:
% 
%   columns 1-14 must be included in input matrix (in case file)
%    1  F_NODE      f, from node number            
%    2  T_NODE      t, to node number
%    3  TYPE_C      compressor type (1 = COMP_G, 2 = COMP_P)
%    4  FG_C        gas flow through compressor (MMSCFD) 
%    5  PC_C        power consumed by compressor (MW)
%    6  GC_C        gas consumed by compressor at from node (MMSCFD)  
%    7  RATIO_C     compressor ratio
%    8  B_C         compression design parameter
%    9  Z_C         ratio design parameter
%    10 ALPHA       x, gas compsumption parameter of gas-fired compresor  
%    11 BETA        y, gas compsumption parameter of gas-fired compresor  
%    12 GAMMA       z, gas compsumption parameter of gas-fired compresor  
%    13 FMAX_C      maximum gas flow trhough compressor (MMSCFD)
%    14 COST_C      compression cost ($/MMSCFD)
% 
%   additional constants, used to assign/compare values in the NODE_TYPE column
%    1  COMP_G      gas-fired compressor
%    2  COMP_P      power-fired compressor
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
%% compressor type
COMP_P      = 1;    %% compressor working with power
COMP_G      = 2;    %% compressor working with gas

%% define the indices 
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

