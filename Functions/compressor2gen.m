function mpgc = compressor2gen(mpgc)
% compressor2gen Creates extra negative generators to model the active power
%   consumed by the compressor working with power.
%       
%   MPGC = COMPRESSOR2GEN(MPGC0);
%   
%   The input case (MPGC0) must contain two additional fields, the gas 
%   case (MPGC0.mgc) and an interconnection case (MPGC0.connect).
% 
%   This function creates negative generators to model the demand generated
%   by such compressors. If the MPGC0 contains the fields 'MPGC0.gennames'
%   and 'MPGC0.genfuel' the function also creates a name and a tag for the
%   new generators.
% 
%   See also MPC2GAS_PREP

%   MPNG: MATPOWER - Natural Gas
%   Copyright (c) 2019-2022 - v0.99beta
%   Sergio García-Marín - Universidad Nacional de Colombia - Sede Manizales
%   Wilson González-Vanegas - Universidad Tecnológica de Pereira
%   Carlos E. Murillo-Sánchez - Universidad Nacional de Colombia - Sede Manizales
% 
%   This file is part of MPNG.
%   Covered by the 3-clause BSD License (see LICENSE file for details).
%%
[PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, ...
    VA, BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN] = idx_bus;
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
[PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST, COST] = idx_cost;
[STO_NODE, STO, STO_0, STOMAX, STOMIN, FSTO, FSTO_OUT, FSTO_IN,...
    FSTOMAX, FSTOMIN, S_STATUS, COST_STO, COST_OUT, COST_IN] = idx_sto;
[F_NODE, T_NODE, FG_O, K_O, DIAM, LNG, FMAX_O, FMIN_O, COST_O] = idx_pipe;
[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect;
%%
mpc = mpgc;
mgc = mpgc.mgc;
connect = mpgc.connect;

buses = mpgc.bus(:,BUS_I);
[ng, cg] = size(mpc.gen);
[~, cgcost] = size(mpc.gencost);
iscomp_p = (mgc.comp(:,TYPE_C) == COMP_P);      % compressors working with power  
idcomp_p = find(iscomp_p);                      % id for power-driven compressors
nc_p = length(idcomp_p);                        % number of compressors-power
buscomp = connect.interc.comp(:,BUS_ID);

if size(connect.interc.comp,1) ~= nc_p ...
        || size(connect.interc.comp,2) ~= 2
   error('compressor2gen: matrix for compressors connected to the power network has wrong dimensions.');  
end
if any(~ismember(buscomp,buses))
    error('compressor2gen: trying to connect a compressor to a non-existing bus.');    
end

gen_comp = zeros(nc_p,cg);
gen_comp(:,GEN_BUS) = buscomp;
gen_comp(:,VG) = 1;
gen_comp(:,MBASE) = mpgc.baseMVA;
gen_comp(:,GEN_STATUS) = 1;
%%
gen_comp(:,PMIN) = -inf;
if any(mgc.comp(idcomp_p,RATIO_C) < 1)                     
    gen_comp(mgc.comp(idcomp_p,RATIO_C) < 1,PMAX) = inf;
end   
gencost_comp = zeros(nc_p,cgcost);
gencost_comp(:,MODEL) = POLYNOMIAL;
gencost_comp(:,NCOST) = 1;
gencost_comp(:,COST:end) = 0;
mpgc.gen = [mpgc.gen; gen_comp];
mpgc.gencost = [mpgc.gencost; gencost_comp];
if isfield(mpgc, 'gennames') && iscell(mpgc.gennames)
    j = size(mpgc.gennames,1);
    for k = 1:nc_p
        bus_comp_name = mpgc.busnames{buscomp(k)};
        mpgc.gennames{j+k} = strcat('compressor_',bus_comp_name);
    end
end
if isfield(mpgc, 'genfuel') && iscell(mpgc.genfuel)
    j = size(mpgc.genfuel,1);
    for k = 1:nc_p
        mpgc.genfuel{j+k} = 'compressor';
    end
end
mpgc.gencomp.N = nc_p;
mpgc.gencomp.id = (ng+1:ng+nc_p)';
mpgc.genid.comp = (ng+1:ng+nc_p)';

