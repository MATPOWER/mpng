function mgc_real = mgc_REAL(mgc_pu)
%% mgc_REAL converts the natural gas case into conventional units.
%   MGC_REAL = MGC_REAL(MGC_PU)
% 
%   This function gets the whole information of the natural case MGC_REAL
%   and converts it into normalized values in MGC_PU, it includes the 
%   constants and the initial values of the variables.
%
%   This function is used in the final stage of MPNG, where the results are
%   reorganized and printed. The whole information of the natural case 
%   (MPC_PU) is converted into conventional values (MGC_REAL), this 
%   includes the constants and the final values of the variables.
%
%   See also MGC_PU

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
mgc_real = mgc_pu;
pbase = mgc_pu.pbase;
fbase = mgc_pu.fbase;
wbase = mgc_pu.wbase;

%% ------------------------------ Node data ------------------------------
mgc_real.node.info(:,PR) = sqrt(mgc_real.node.info(:,PR).*pbase^2);
mgc_real.node.info(:,PRMAX) = sqrt(mgc_real.node.info(:,PRMAX).*pbase^2);
mgc_real.node.info(:,PRMIN) = sqrt(mgc_real.node.info(:,PRMIN).*pbase^2);
mgc_real.node.info(:,OVP) = sqrt(mgc_real.node.info(:,OVP).*pbase^2);
mgc_real.node.info(:,UNP) = sqrt(mgc_real.node.info(:,UNP).*pbase^2);
mgc_real.node.info(:,GD) = mgc_real.node.info(:,GD)*fbase;
mgc_real.node.dem = mgc_real.node.dem*fbase;
% cost
mgc_real.node.info(:,COST_OVP) = mgc_real.node.info(:,COST_OVP)/pbase;
mgc_real.node.info(:,COST_UNP) = mgc_real.node.info(:,COST_UNP)/pbase;
mgc_real.node.demcost = mgc_real.node.demcost/fbase;

%% ------------------------------ Well data ------------------------------
mgc_real.well(:,G) = mgc_real.well(:,G)*fbase;
mgc_real.well(:,PW) = sqrt(mgc_real.well(:,PW).*pbase^2);
mgc_real.well(:,GMAX) = mgc_real.well(:,GMAX)*fbase;
mgc_real.well(:,GMIN) = mgc_real.well(:,GMIN)*fbase;
% cost
mgc_real.well(:,COST_G) = mgc_real.well(:,COST_G)/fbase;

%% ---------------------------- Pipeline data ----------------------------
mgc_real.pipe(:,FG_O) = mgc_real.pipe(:,FG_O)*fbase;
mgc_real.pipe(:,K_O) = mgc_real.pipe(:,K_O)*(fbase/pbase);
mgc_real.pipe(:,FMAX_O) = mgc_real.pipe(:,FMAX_O)*fbase;
mgc_real.pipe(:,FMIN_O) = mgc_real.pipe(:,FMIN_O)*fbase;
% cost
mgc_real.pipe(:,COST_O) = mgc_real.pipe(:,COST_O)/fbase;


%% ---------------------------- Compressor data ----------------------------
mgc_real.comp(:,FG_C) = mgc_real.comp(:,FG_C)*fbase;
mgc_real.comp(:,PC_C) = mgc_real.comp(:,PC_C)*wbase;
mgc_real.comp(:,GC_C) = mgc_real.comp(:,GC_C)*fbase;
mgc_real.comp(:,B_C) = mgc_real.comp(:,B_C)*(wbase/fbase);
mgc_real.comp(:,ALPHA) = mgc_real.comp(:,ALPHA)*fbase;
mgc_real.comp(:,BETA) = mgc_real.comp(:,BETA)*(fbase/wbase);
mgc_real.comp(:,GAMMA) = mgc_real.comp(:,GAMMA)*(fbase/wbase^2);
mgc_real.comp(:,FMAX_C) = mgc_real.comp(:,FMAX_C)*fbase;
% cost
mgc_real.comp(:,COST_C) = mgc_real.comp(:,COST_C)/fbase;
%% -------------------------------- Storage --------------------------------
mgc_real.sto(:,STO_0)    = mgc_real.sto(:,STO_0)*fbase;
mgc_real.sto(:,STOMAX)   = mgc_real.sto(:,STOMAX)*fbase;
mgc_real.sto(:,STOMIN)   = mgc_real.sto(:,STOMIN)*fbase;
mgc_real.sto(:,FSTO)     = mgc_real.sto(:,FSTO)*fbase;
mgc_real.sto(:,FSTO_OUT) = mgc_real.sto(:,FSTO_OUT)*fbase;
mgc_real.sto(:,FSTO_IN)  = mgc_real.sto(:,FSTO_IN)*fbase;
mgc_real.sto(:,FSTOMAX)  = mgc_real.sto(:,FSTOMAX)*fbase;
mgc_real.sto(:,FSTOMIN)  = mgc_real.sto(:,FSTOMIN)*fbase;
%cost
mgc_real.sto(:,COST_STO) = mgc_real.sto(:,COST_STO)/fbase;
mgc_real.sto(:,COST_OUT) = mgc_real.sto(:,COST_OUT)/fbase;
mgc_real.sto(:,COST_IN)  = mgc_real.sto(:,COST_IN)/fbase;