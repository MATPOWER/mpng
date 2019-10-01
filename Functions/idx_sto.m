function [STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto
%IDX_STO  Defines constants for named column indices to storage matrix.
% 
%   The index, name and meaning of each column of the storage matrix 
%   are given below:
%
%   STO_NODE	= 1;    %% storage node
%   STO         = 2;    %% storage level -> (after running program)
%   STO_0       = 3;    %% storage initial level 
%   STOMAX      = 4;    %% maximun storage
%   STOMIN      = 5;    %% minimun storage
%   FSTO        = 6;    %% storage outflow diference \      
%   FSTO_OUT	= 7;    %% storage outflow            | -> (after running program)  
%   FSTO_IN     = 8;    %% storage inflow            /
%   FSTOMAX     = 9;    %% maximun storage outflow diference     
%   FSTOMIN     = 10;	%% minimun storage outflow diference
%   S_STATUS	= 11;	%% storage status 
%   COST_STO	= 12;	%% storage cost 
%   COST_OUT	= 13;   %% storage outflow cost
%   COST_IN     = 14;   %% storage inflow cost
% 
%   See also DEFINE_CONSTANTS_GAS.

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales

%   3-clause bsd license 

%% define the indices
STO_NODE	= 1;    %% storage node
STO         = 2;    %% storage level -> (after running program)
STO_0       = 3;    %% storage initial level 
STOMAX      = 4;    %% maximun storage
STOMIN      = 5;    %% minimun storage
FSTO        = 6;    %% storage outflow diference
FSTO_OUT	= 7;    %% storage outflow 
FSTO_IN     = 8;    %% storage inflow 
FSTOMAX     = 9;    %% maximun storage outflow diference     
FSTOMIN     = 10;	%% minimun storage outflow diference
S_STATUS	= 11;	%% storage status 
COST_STO	= 12;	%% storage cost 
COST_OUT	= 13;   %% storage outflow cost
COST_IN     = 14;   %% storage inflow cost
