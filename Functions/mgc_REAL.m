function mgcREAL = mgc_REAL(mgc)
%% mgc_REAL converts the natural gas case into conventional units. 
%
%   This function is used in the final stage of MPNG, where the results are
%   reorganized and printed. The whole information of the natural case is
%   converted into conventional values, this includes the constants and the
%   final values of the variables.
%
%   See also MGC_PU

%   MPNG: Matpower - Natural Gas
%   Copyright (c) 2019 - v0.99alpha
%   Sergio García Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo Sánchez - Universidad Nacional de Colombia - Sede Manizales
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
mgcREAL = mgc;
pbase = mgc.pbase;
fbase = mgc.fbase;
wbase = mgc.wbase;

%% ------------------------------ Node data ------------------------------
mgcREAL.node.info(:,PR) = sqrt(mgcREAL.node.info(:,PR).*pbase^2);
mgcREAL.node.info(:,PRMAX) = sqrt(mgcREAL.node.info(:,PRMAX).*pbase^2);
mgcREAL.node.info(:,PRMIN) = sqrt(mgcREAL.node.info(:,PRMIN).*pbase^2);
mgcREAL.node.info(:,OVP) = sqrt(mgcREAL.node.info(:,OVP).*pbase^2);
mgcREAL.node.info(:,UNP) = sqrt(mgcREAL.node.info(:,UNP).*pbase^2);
mgcREAL.node.info(:,GD) = mgcREAL.node.info(:,GD)*fbase;
mgcREAL.node.dem = mgcREAL.node.dem*fbase;
% cost
mgcREAL.node.info(:,COST_OVP) = mgcREAL.node.info(:,COST_OVP)/pbase;
mgcREAL.node.info(:,COST_UNP) = mgcREAL.node.info(:,COST_UNP)/pbase;
mgcREAL.node.demcost = mgcREAL.node.demcost/fbase;

%% ------------------------------ Well data ------------------------------
mgcREAL.well(:,G) = mgcREAL.well(:,G)*fbase;
mgcREAL.well(:,PW) = sqrt(mgcREAL.well(:,PW).*pbase^2);
mgcREAL.well(:,GMAX) = mgcREAL.well(:,GMAX)*fbase;
mgcREAL.well(:,GMIN) = mgcREAL.well(:,GMIN)*fbase;
% cost
mgcREAL.well(:,COST_G) = mgcREAL.well(:,COST_G)/fbase;

%% ---------------------------- Pipeline data ----------------------------
mgcREAL.pipe(:,FG_O) = mgcREAL.pipe(:,FG_O)*fbase;
mgcREAL.pipe(:,K_O) = mgcREAL.pipe(:,K_O)*(fbase/pbase);
mgcREAL.pipe(:,FMAX_O) = mgcREAL.pipe(:,FMAX_O)*fbase;
mgcREAL.pipe(:,FMIN_O) = mgcREAL.pipe(:,FMIN_O)*fbase;
% cost
mgcREAL.pipe(:,COST_O) = mgcREAL.pipe(:,COST_O)/fbase;


%% ---------------------------- Compressor data ----------------------------
mgcREAL.comp(:,FG_C) = mgcREAL.comp(:,FG_C)*fbase;
mgcREAL.comp(:,PC_C) = mgcREAL.comp(:,PC_C)*wbase;
mgcREAL.comp(:,GC_C) = mgcREAL.comp(:,GC_C)*fbase;
mgcREAL.comp(:,B_C) = mgcREAL.comp(:,B_C)*(wbase/fbase);
mgcREAL.comp(:,ALPHA) = mgcREAL.comp(:,ALPHA)*fbase;
mgcREAL.comp(:,BETA) = mgcREAL.comp(:,BETA)*(fbase/wbase);
mgcREAL.comp(:,GAMMA) = mgcREAL.comp(:,GAMMA)*(fbase/wbase^2);
mgcREAL.comp(:,FMAX_C) = mgcREAL.comp(:,FMAX_C)*fbase;
% cost
mgcREAL.comp(:,COST_C) = mgcREAL.comp(:,COST_C)/fbase;
%% -------------------------------- Storage --------------------------------
mgcREAL.sto(:,STO_0)    = mgcREAL.sto(:,STO_0)*fbase;
mgcREAL.sto(:,STOMAX)   = mgcREAL.sto(:,STOMAX)*fbase;
mgcREAL.sto(:,STOMIN)   = mgcREAL.sto(:,STOMIN)*fbase;
mgcREAL.sto(:,FSTO)     = mgcREAL.sto(:,FSTO)*fbase;
mgcREAL.sto(:,FSTO_OUT) = mgcREAL.sto(:,FSTO_OUT)*fbase;
mgcREAL.sto(:,FSTO_IN)  = mgcREAL.sto(:,FSTO_IN)*fbase;
mgcREAL.sto(:,FSTOMAX)  = mgcREAL.sto(:,FSTOMAX)*fbase;
mgcREAL.sto(:,FSTOMIN)  = mgcREAL.sto(:,FSTOMIN)*fbase;
%cost
mgcREAL.sto(:,COST_STO) = mgcREAL.sto(:,COST_STO)/fbase;
mgcREAL.sto(:,COST_OUT) = mgcREAL.sto(:,COST_OUT)/fbase;
mgcREAL.sto(:,COST_IN)  = mgcREAL.sto(:,COST_IN)/fbase;