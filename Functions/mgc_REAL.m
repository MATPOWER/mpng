function mgcREAL = mgc_REAL(mgc)
%% mgc_REAL converts the natural gas case (which had been previously changed 
%   to PU) into conventional units. This function is used in the final
%   stage of MPEG, where the data is reorganize and printed.

%% define constants
[DEM, WELL, NODE_I, NODE_TYPE, PR, PRMAX, PRMIN, OVP, UNP, COST_OVP, ...
    COST_UNP, GD, NGD] = idx_node;
[WELL_NODE, G, PW, GMAX, GMIN, WELL_STATUS, COST_G] = idx_well;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;

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
%% This will be analyzed later