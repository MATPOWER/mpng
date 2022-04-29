function mgc_pu = mgc_PU(mgc_real)
% MGC_PU converts the matrices inside the MGC to per unit (p.u.) data. 
%   MGC_PU = MGC_PU(MGC_REAL)
% 
%   This function gets the whole information of the natural case MGC_REAL
%   and converts it into normalized values in MGC_PU, it includes the 
%   constants and the initial values of the variables.
%
%   The bases must be specified inside the case, however they have the 
%   following default values:
%       pressure base =     500 (psi)       
%       gas flow      =     50  (MMSCF)
%       wbase         =     100 (MVA)
%
%   After several tries we found that an accurate selection of the bases
%   will improve the perfomance of the optimal power and natural gas flow.
%   As the problems might differ according to the topology of the gas
%   network we recommend that the user find its own bases. It is also 
%   recomended that the power base in the gas case agrees with the power
%   base of the power case.
%
%   See also MGC_REAL

%   MPNG: MATPOWER - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

%% define constants
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
%%
mgc_pu = mgc_real;
%% get the bases
if isempty(mgc_pu.pbase) 
    pbase = 500;
    mgc_pu.pbase = pbase;
else
    pbase = mgc_pu.pbase;
end
if isempty(mgc_pu.fbase) 
    fbase = 50;
    mgc_pu.fbase = fbase;
else
    fbase = mgc_pu.fbase;
end
if isempty(mgc_pu.wbase)
    wbase = 100;
    mgc_pu.wbase = wbase;
else 
    wbase = mgc_pu.wbase;
end 
%% ------------------------------ Node data ------------------------------
mgc_pu.node.info(:,PR) = mgc_pu.node.info(:,PR).^2/pbase^2;
mgc_pu.node.info(:,PRMAX) = mgc_pu.node.info(:,PRMAX).^2/pbase^2;
mgc_pu.node.info(:,PRMIN) = mgc_pu.node.info(:,PRMIN).^2/pbase^2;
mgc_pu.node.info(:,OVP) = mgc_pu.node.info(:,OVP).^2/pbase^2;
mgc_pu.node.info(:,UNP) = mgc_pu.node.info(:,UNP).^2/pbase^2;
mgc_pu.node.info(:,GD) = mgc_pu.node.info(:,GD)/fbase;
mgc_pu.node.dem = mgc_pu.node.dem/fbase;
% cost
mgc_pu.node.info(:,COST_OVP) = mgc_pu.node.info(:,COST_OVP)*pbase^2;
mgc_pu.node.info(:,COST_UNP) = mgc_pu.node.info(:,COST_UNP)*pbase^2;
mgc_pu.node.demcost = mgc_pu.node.demcost*fbase;

%% ------------------------------ Well data ------------------------------
mgc_pu.well(:,G) = mgc_pu.well(:,G)/fbase;
mgc_pu.well(:,PW) = mgc_pu.well(:,PW).^2/pbase^2;
mgc_pu.well(:,GMAX) = mgc_pu.well(:,GMAX)/fbase;
mgc_pu.well(:,GMIN) = mgc_pu.well(:,GMIN)/fbase;
% cost
mgc_pu.well(:,COST_G) = mgc_pu.well(:,COST_G)*fbase;

%% ---------------------------- Pipeline data ----------------------------
mgc_pu.pipe(:,FG_O) = mgc_pu.pipe(:,FG_O)/fbase;
mgc_pu.pipe(:,K_O) = mgc_pu.pipe(:,K_O)/(fbase/pbase);
mgc_pu.pipe(:,FMAX_O) = mgc_pu.pipe(:,FMAX_O)/fbase;
mgc_pu.pipe(:,FMIN_O) = mgc_pu.pipe(:,FMIN_O)/fbase;
% cost
mgc_pu.pipe(:,COST_O) = mgc_pu.pipe(:,COST_O)*fbase;


%% ---------------------------- Compressor data ----------------------------
mgc_pu.comp(:,FG_C) = mgc_pu.comp(:,FG_C)/fbase;
mgc_pu.comp(:,PC_C) = mgc_pu.comp(:,PC_C)/wbase;
mgc_pu.comp(:,GC_C) = mgc_pu.comp(:,GC_C)/fbase;
mgc_pu.comp(:,RATIO_C) = mgc_pu.comp(:,RATIO_C).^2;
mgc_pu.comp(:,B_C) = mgc_pu.comp(:,B_C)/(wbase/fbase);
mgc_pu.comp(:,ALPHA) = mgc_pu.comp(:,ALPHA)/fbase;
mgc_pu.comp(:,BETA) = mgc_pu.comp(:,BETA)/(fbase/wbase);
mgc_pu.comp(:,GAMMA) = mgc_pu.comp(:,GAMMA)/(fbase/wbase^2);
mgc_pu.comp(:,FMAX_C) = mgc_pu.comp(:,FMAX_C)/fbase;
% cost
mgc_pu.comp(:,COST_C) = mgc_pu.comp(:,COST_C)*fbase;
%% -------------------------------- Storage --------------------------------
mgc_pu.sto(:,STO_0)      = mgc_pu.sto(:,STO_0)/fbase;
mgc_pu.sto(:,STOMAX)     = mgc_pu.sto(:,STOMAX)/fbase;
mgc_pu.sto(:,STOMIN)     = mgc_pu.sto(:,STOMIN)/fbase;
mgc_pu.sto(:,FSTO)       = mgc_pu.sto(:,FSTO)/fbase;
mgc_pu.sto(:,FSTO_OUT)   = mgc_pu.sto(:,FSTO_OUT)/fbase;
mgc_pu.sto(:,FSTO_IN)    = mgc_pu.sto(:,FSTO_IN)/fbase;
mgc_pu.sto(:,FSTOMAX)    = mgc_pu.sto(:,FSTOMAX)/fbase;
mgc_pu.sto(:,FSTOMIN)    = mgc_pu.sto(:,FSTOMIN)/fbase;
%cost
mgc_pu.sto(:,COST_STO)   = mgc_pu.sto(:,COST_STO)*fbase;
mgc_pu.sto(:,COST_OUT)   = mgc_pu.sto(:,COST_OUT)*fbase;
mgc_pu.sto(:,COST_IN)    = mgc_pu.sto(:,COST_IN)*fbase;
