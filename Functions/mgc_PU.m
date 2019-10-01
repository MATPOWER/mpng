function mgcPU = mgc_PU(mgc)
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
mgcPU = mgc;
pbase = mgc.pbase;
fbase = mgc.fbase;
wbase = mgc.wbase;

%% ------------------------------ Node data ------------------------------
mgcPU.node.info(:,PR) = mgcPU.node.info(:,PR).^2/pbase^2;
mgcPU.node.info(:,PRMAX) = mgcPU.node.info(:,PRMAX).^2/pbase^2;
mgcPU.node.info(:,PRMIN) = mgcPU.node.info(:,PRMIN).^2/pbase^2;
mgcPU.node.info(:,OVP) = mgcPU.node.info(:,OVP).^2/pbase^2;
mgcPU.node.info(:,UNP) = mgcPU.node.info(:,UNP).^2/pbase^2;
mgcPU.node.info(:,GD) = mgcPU.node.info(:,GD)/fbase;
mgcPU.node.dem = mgcPU.node.dem/fbase;
% cost
mgcPU.node.info(:,COST_OVP) = mgcPU.node.info(:,COST_OVP)*pbase;
mgcPU.node.info(:,COST_UNP) = mgcPU.node.info(:,COST_UNP)*pbase;
mgcPU.node.demcost = mgcPU.node.demcost*fbase;

%% ------------------------------ Well data ------------------------------
mgcPU.well(:,G) = mgcPU.well(:,G)/fbase;
mgcPU.well(:,PW) = mgcPU.well(:,PW).^2/pbase^2;
mgcPU.well(:,GMAX) = mgcPU.well(:,GMAX)/fbase;
mgcPU.well(:,GMIN) = mgcPU.well(:,GMIN)/fbase;
% cost
mgcPU.well(:,COST_G) = mgcPU.well(:,COST_G)*fbase;

%% ---------------------------- Pipeline data ----------------------------
mgcPU.pipe(:,FG_O) = mgcPU.pipe(:,FG_O)/fbase;
mgcPU.pipe(:,K_O) = mgcPU.pipe(:,K_O)/(fbase/pbase);
mgcPU.pipe(:,FMAX_O) = mgcPU.pipe(:,FMAX_O)/fbase;
mgcPU.pipe(:,FMIN_O) = mgcPU.pipe(:,FMIN_O)/fbase;
% cost
mgcPU.pipe(:,COST_O) = mgcPU.pipe(:,COST_O)*fbase;


%% ---------------------------- Compressor data ----------------------------
mgcPU.comp(:,FG_C) = mgcPU.comp(:,FG_C)/fbase;
mgcPU.comp(:,PC_C) = mgcPU.comp(:,PC_C)/wbase;
mgcPU.comp(:,GC_C) = mgcPU.comp(:,GC_C)/fbase;
mgcPU.comp(:,B_C) = mgcPU.comp(:,B_C)/(wbase/fbase);
mgcPU.comp(:,ALPHA) = mgcPU.comp(:,ALPHA)/fbase;
mgcPU.comp(:,BETA) = mgcPU.comp(:,BETA)/(fbase/wbase);
mgcPU.comp(:,GAMMA) = mgcPU.comp(:,GAMMA)/(fbase/wbase^2);
mgcPU.comp(:,FMAX_C) = mgcPU.comp(:,FMAX_C)/fbase;
% cost
mgcPU.comp(:,COST_C) = mgcPU.comp(:,COST_C)*fbase;
%% -------------------------------- Storage --------------------------------
mgcPU.sto(:,STO_0) = mgcPU.sto(:,STO_0)/fbase;
mgcPU.sto(:,STOMAX) = mgcPU.sto(:,STOMAX)/fbase;
mgcPU.sto(:,STOMIN) = mgcPU.sto(:,STOMIN)/fbase;
mgcPU.sto(:,FSTO) = mgcPU.sto(:,FSTO)/fbase;
mgcPU.sto(:,FSTO_OUT) = mgcPU.sto(:,FSTO_OUT)/fbase;
mgcPU.sto(:,FSTO_IN) = mgcPU.sto(:,FSTO_IN)/fbase;
mgcPU.sto(:,FSTOMAX) = mgcPU.sto(:,FSTOMAX)/fbase;
mgcPU.sto(:,FSTOMIN) = mgcPU.sto(:,FSTOMIN)/fbase;
%cost
mgcPU.sto(:,COST_STO) = mgcPU.sto(:,COST_STO)*fbase;
mgcPU.sto(:,COST_OUT) = mgcPU.sto(:,COST_OUT)*fbase;
mgcPU.sto(:,COST_IN) = mgcPU.sto(:,COST_IN)*fbase;
%% This will be analyzed later