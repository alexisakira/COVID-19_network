clear
close all
clc;

%% figure formatting

set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultAxesTickLabelInterpreter','latex');
set(0,'DefaultLegendInterpreter', 'latex')
   
set(0,'DefaultTextFontSize', 14)
set(0,'DefaultAxesFontSize', 14)
set(0,'DefaultLineLineWidth',2)

temp = get(gca,'ColorOrder');
c1 = temp(1,:);
c2 = temp(2,:);

close all

%% set parameters
R0 = 3; % reproduction number
gamma = 0.1; % recovery probability
D = 10; % average number of connections (network degree)
beta = gamma*R0/D; % transmission probability
%I = 1000; % number of agents
I = 1000; % use a small number for debugging
p_rewire = 0.1; % rewiring probability for WS network
m0 = 1; % number of initial seeds in BA network
y0 = 0.01; % fraction of initially infected

N1vec = [2 5 10]; % upper bounds on number of meetings
Nspec = length(N1vec)+1; % number of specifications (including benchmark)

N0 = 100; % maximum size of meetings after social distancing
Capacity0 = N0*ones(I,1); % capacity vector during social distancing

T = 200; % number of periods
T1 = 50; % number of days to social distance
time = [1:T];

type = {'ERG','WS','BA'}; % type of networks
typeFull = {'Erd\"os-Renyi-Gilbert','Watts-Strogatz','Barab\''asi-Albert'};
Ntype = length(type);