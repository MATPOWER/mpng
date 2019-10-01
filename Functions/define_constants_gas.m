%DEFINE_CONSTANTS_GAS  Defines useful constants for indexing data in the 
%   additional cases for the matpower-gas formulation.
%
%   node:
%      DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP,
%      COST_UNP, GD, NGD
%
%   pipeline:
%      F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O
% 
%   compressor:
%      COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, 
%      BETA, GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C
%
%   well: 
%      WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G
%
%   storage:
%       STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,
%       FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN
% 
%   connection struct
%       GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EF
%
%   DEFINE_CONSTANTS calls IDX_NODE, IDX_PLINE, IDX_WELL, IDX_STO and IDX_CONNECT.
% 
%   See also IDX_NODE, IDX_PLINE, IDX_WELL, IDX_STO and IDX_CONNECT. 
%
%   This script is included for convenience for interactive use or
%   for high-level code where maximum performance is not a concern.

%% define named indices into data matrices
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
[GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect;
