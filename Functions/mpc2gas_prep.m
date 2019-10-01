function mpgc = mpc2gas_prep(mpgc,mp)
% mpc2gas_prep Changes the original matpower case in order to
%   add multiple time scenarios and extra generators to model the
%   compressors and the dispatchable load.

%% Initializing things

[COMP_G, COMP_P, F_NODE, T_NODE, TYPE_C, RATIO_C, B_C, Z_C, ALPHA, BETA, ...
    GAMMA, FMAX_C, COST_C, FG_C, GC_C, PC_C] = idx_comp;
[GEN_ID, MAX_ENER, COMP_ID, BUS_ID, NODE_ID, EFF] = idx_connect;
if nargin == 2
    verbose = mp.verbose;
end
%% create extra generators for power power consuption compressors
iscomp_p = (mpgc.mgc.comp(:,TYPE_C) == COMP_P);         % compressors working with power  
iscomp_p_matrix = ~isempty(mpgc.connect.interc.comp);
idcomp_p = find(iscomp_p);                              % id for power-driven compressors
if verbose
    fprintf(1,'     Checking compressors working with power... ');
end
if any(iscomp_p) && iscomp_p_matrix
    if_comp = (idcomp_p(:) ~= mpgc.connect.interc.comp(:,COMP_ID));
    comp_p = mpgc.connect.interc.comp(:,COMP_ID);
    if any(if_comp) || ~isequal(idcomp_p,comp_p)
        error('mpc2gas_prep: id´s for compressors connected to the power network mismatch.');
    end
    mpgc = compressor2gen(mpgc);
    if verbose
        if length(comp_p) == 1
            fprintf(1, 'There is 1.\n');
        else
            fprintf(1, 'There are %u.\n',length(comp_p));
        end
    end
elseif ~isequal(any(iscomp_p),iscomp_p_matrix)
    error('mpc2gas_prep: power-driven compressors info is not consistent.');
elseif  ~any(iscomp_p) && ~iscomp_p_matrix && verbose
    fprintf(1, 'There are none.\n');
end

%% create multiple power time scenarios
if verbose
    fprintf(1,'     Creating multiple time scenarios...\n');
end
time = mpgc.connect.power.time;
pd = mpgc.connect.power.demands.pd;
qd = mpgc.connect.power.demands.qd;
mpgc = multi_period(mpgc,time,pd,qd);

%% create extra generators for non-supplied demands
if verbose
    fprintf(1,'     Creating dispatchable loads...\n');
end
cost_nsd = mpgc.connect.power.cost;
mpgc = nsd2gen(mpgc,cost_nsd);