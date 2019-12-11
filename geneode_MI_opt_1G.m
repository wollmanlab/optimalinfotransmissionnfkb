%% Optimization of gene parameters using Simulated Anneling
clc;
clear all; close all;
rng('shuffle');
%% load biologically relevant parameter sets
num = csvread('/home/alok/Desktop/NF_kB/Plos_Comp_Rev2/Code_For_Github/2017-06-27-shared-v1-genes.csv',1,1);
num0 = [num(:,2), num(:,3), num(:,5), num(:,4), num(:,1), num(:,6)]; 
params0i = num0(16,:); % opt one gene parameter set 
params0i(2) = 0.3;
%
pert = 0.25; % apply 25% perturbation for gene parameters
%%
ts = 5; % five different signals (TNF, LPS, polyIC, CpG, Pam3CSK)
load('/home/alok/Desktop/NF_kB/Plos_Comp_Rev2/Code_For_Github/NFkB_normATS_pert25_sig_med.mat');
%% ode simulation of target gene code
TFt = (0:480);
%
sim_time = 451;
% Default tolarances for ODE solver
RELTOL = 1e-6;
ABSTOL = 1e-8;
ode_opt = odeset('RelTol',RELTOL,'AbsTol',ABSTOL);
%
tsp = [50, 100, 150, 300, 450]; % 5 time points
%
c1 = cell(ts,1);
% p = parpool('local',ts);
parfor signal = 1: ts 
     G1 = zeros(1*5,length(NFkB_time_series_differ_signals{signal}));
     %
     for in = 1:length(NFkB_time_series_differ_signals{signal})
         TF = NFkB_time_series_differ_signals{signal}(:,in)';
         params0 = params_pert(pert, params0i); % log-normal perturbation for gene parameters
         delg0 = (params0(1)/params0(3));
         [t,Gt]= ode15s(@geneInitialize_test, [0 sim_time],delg0,ode_opt,TFt,TF,params0); % use ode15s to solve ode
         % geneInitialize_test defines the gene model
         H1 = interp1(t',Gt(:,1)',tsp); % interpolate the G1 data set at time t (0:1:480)
         G1(:,in) = H1';
     end
        c1{signal} = G1;
end
% delete(p)
%% Mutual information of target gene
MI = Mutual_Information_NFkB_dim_exp(ts,c1); % Calculate mutual information
MI0 = MI;
params0=params0i;
%% set starting value for Simulated Anneling optimization
Tat0 = 0.001;
Tsim = 0.001;
Tcut = 1e-30;
Tred = 0.02;
nsa = 1000;
%% Simulated Anneling programme
T_array = zeros(nsa,1);
MI_array = zeros(nsa,1);
params_array = zeros(nsa,6);
for nmove = 1: nsa
    Tat = Tat0;
    T_array(nmove,:) = Tat; 
    MI_array(nmove,:) = MI0;
    params = params0;
    %
    m = randi([1 6]); % randomly opt one parameter among six gene parameters
    n = randi([1, 2]);
    o = rand;
    delta = 0.05;
    fac = (-1)^n*delta*o;
    params0 = log10(params0);
    params0(m) = params0(m) + fac*params0(m); % perturbed the parameter
    params0 = 10.^params0; 
    params0i = params0;
    %%
    c1 = cell(ts,1);
%     p = parpool('local',ts);
    parfor signal = 1: ts
        G1 = zeros(1*5,length(NFkB_time_series_differ_signals{signal}));
        %
        for in = 1:length(NFkB_time_series_differ_signals{signal})
            TF = NFkB_time_series_differ_signals{signal}(:,in)';
            params0 = params_pert(pert, params0i); % log-normal perturbation for gene parameters
            delg0 = (params0(1)/params0(3));
            [t,Gt]= ode15s(@geneInitialize_test, [0 sim_time],delg0,ode_opt,TFt,TF,params0); % use ode15s to solve ode
            % geneInitialize_test defines the gene model
            H1 = interp1(t',Gt(:,1)',tsp); % interpolate the G1 data set at time t (0:1:480)
            G1(:,in) = H1';
        end
        c1{signal} = G1;
    end
%     delete(p)
%% Mutual information of target genes
    MI = Mutual_Information_NFkB_dim_exp(ts,c1);
%% Metropolish test
    dMI = ((-MI) - (-MI0));
    metrop = exp(-dMI/Tat);
    if dMI <0
       MI0 = MI;
       disp([num2str(nmove),'*1*',num2str(m),'***',num2str(dMI),'***',num2str(MI0)]);
       params_array(nmove,:) = params0;
    elseif metrop > rand
        MI0 = MI;
        disp([num2str(nmove),'*2*',num2str(m),'***',num2str(dMI),'***',num2str(MI0)]);
        params_array(nmove,:) = params0;
    else
        params0 = params;
        params_array(nmove,:) = params0;
        disp('error');
    end
    Tat = Tat - (Tat*Tred);
    if Tat < Tcut 
        Tat = Tsim;
    end
    Tat0 = Tat;
 save('/home/alok/Desktop/NF_kB/Plos_Comp_Rev2/Code_For_Github/MI_opt_exp1G16.mat')
end    
