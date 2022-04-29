function [STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto
%IDX_STO  Defines constants for named column indices to storage matrix.
%   Example:
%   
%   [STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
%    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
%  
%   Some examples of usage, after defining the constants using the line above,
%   are:
%
%    sto(5, STO_0) = 0;         % set to zero the initial storage of store
%                                 unit 5
%
%   The index, name and meaning of each column of the storage matrix 
%   are given below:
%
%   columns 1-14 must be included in input matrix (in case file)
%    1  STO_NODE        node number        
%    2  STO             storage final level (MMSCF)
%    3  STO_0           storage initial level (MMSCF)
%    4  STOMAX          maximum storage (MMSCF)
%    5  STOMIN          minimum storage (MMSCF)
%    6  FSTO            storage outflow diference (MMSCFD)         
%    7  FSTO_OUT        storage outflow (MMSCFD)         
%    8  FSTO_IN         storage inflow (MMSCFD)      
%    9  FSTOMAX         maximum storage outflow diference (MMSCFD)
%    10 FSTOMIN         minimun storage outflow diference (MMSCFD)
%    11 S_STATUS        storage unit status  
%    12 COST_STO        storage cost ($/MMSCF)
%    13 COST_OUT        storage outflow cost ($/MMSCFD)
%    14 COST_IN         storage inflow cost ($/MMSCFD)
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
