%% MPNG - Example 2. Optimal power and natural gas flow. 48 node meshed 
%   natural gas system without power case

%% Define useful constants
define_constants                % power system cosntants
define_constants_gas            % natural gas system and connect struct constants
%% Load cases
mpc = loadcase('case2');        % 2 bus power system 
mgc = ng_case48;                % 48 node looped natural gas system
connect = connect_pg_case2;     % interconnection case, it mainly has to fit with the power case
%% putting all together
mpgc = mpc;                         % initialize MPNG case
mpgc.mgc = mgc;                     % adding gas case
mpgc.connect = connect;             % adding connection struct
%% defining option vector
mpopt = mpoption;                   % initialize option struct
mpopt.opf.ac.solver = 'IPOPT';      % current stable sover
mpopt.ipopt.opts.max_iter = 1e8;    % max iterations
%% running program 
res_base = mpng(mpgc,mpopt);             % running MPNG
%% changes in the contingency
mgc_cont = mgc;
mgc_cont.well(9,GMAX) = 200;        % set max injection of well 9
cont_pipe = 38;                     % id for out of service pipeline    
mgc_cont.pipe(cont_pipe,:) = [];         % take pipeline out of service
%% putting all together
mpgc_cont = mpc;                % initialize MPNG case
mpgc_cont.mgc = mgc_cont;           % adding gas case
mpgc_cont.connect = connect;	% adding connection struct
%% running program 
res_cont = mpng(mpgc_cont,mpopt);	% running MPNG
%% Print results 
%% pressures
p1 = res_base.mgc.node.info(:,PR);
p2 = res_cont.mgc.node.info(:,PR);
x = 1:48;
figure
hold on
grid on
box on
stem(x,[p1 p2],'filled')
axis([0 49 1100 1500])
legend('Base case','Contingency case')
xlabel('Node')
ylabel('Pressure [psia]')
title('Nodal Pressure')
% set(gca, 'XTick',x)  
% set(gca, 'XTick',x, 'FontSize',8)
% xtickangle(90)
%% pipelines
fpipe = res_base.mgc.pipe(:,F_NODE);
tpipe = res_base.mgc.pipe(:,T_NODE);
tag_pipe = cell(length(fpipe),1);
for k = 1:length(fpipe) 
    tag_pipe{k} = strcat(num2str(fpipe(k)),'-',num2str(tpipe(k)));    
end
fgo_base = res_base.mgc.pipe(:,FG_O);
fgo_cont = res_cont.mgc.pipe(:,FG_O);
fgo_cont_a = zeros(43,1);
fgo_cont_a(1:cont_pipe-1) = fgo_cont(1:cont_pipe-1);
fgo_cont_a(cont_pipe) = 0;
fgo_cont_a(cont_pipe+1:43) = fgo_cont(cont_pipe:end);
x = 1:length(fgo_base);
figure
hold on
grid on
box on
bar(x,[fgo_base fgo_cont_a])
axis([0 44 -120 1200])
legend('Base case','Contingency case')
xlabel('Pipeline')
ylabel('Gas Flow [MMSCFD]')
title('Gas Flow in All Pipelines')
set(gca,'xtick',1:length(fpipe))
set(gca,'xticklabel',cellstr(tag_pipe));

xtickangle(90)
%% compressors
fcomp = res_base.mgc.comp(:,F_NODE);
tcomp = res_base.mgc.comp(:,T_NODE);
tag_comp = cell(length(fcomp),1);
for k = 1:length(fcomp) 
    tag_comp{k} = strcat(num2str(fcomp(k)),'-',num2str(tcomp(k)));    
end
fgc_base = res_base.mgc.comp(:,FG_C);
fgc_cont = res_cont.mgc.comp(:,FG_C);
x = 1:length(fcomp);
figure
hold on
grid on
box on
bar(x,[fgc_base fgc_cont])
axis([0 length(fcomp)+1 0 1200])
legend('Base case','Contingency case')
xlabel('Compressor')
ylabel('Gas flow [MMSCFD]')
title('Gas Flow in All Compressors')
set(gca,'xtick',1:length(fcomp))
set(gca,'xticklabel',cellstr(tag_comp));
%% injections
g_base = res_base.mgc.well(:,G);
g_cont = res_cont.mgc.well(:,G);
gmax1 = res_base.mgc.well(:,GMAX);
gmax2 = res_cont.mgc.well(:,GMAX);
x = 1:length(g_base);
figure
hold on
grid on
box on
bar(x,[g_base g_cont])
plot(x-0.15,gmax1,'*b','markersize',10)
plot(x+0.15,gmax2,'*r','markersize',10)
axis([0.5 length(x)+.5 0 1000])
legend('Base case (B.C.)','Contingency case (C.C.)','Max. Injection B.C.','Max. Injection C.C.')
xlabel('Gas Wells')
ylabel('Gas Injection [MMSCFD]')
title('Gas Injections at All Wells')