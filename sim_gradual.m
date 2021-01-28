setParam; % set parameters

N1vec = [2:10];
T1 = 10;
T = 200;
time = 1:T;

rng(1);
I0 = binornd(1,y0,[I 1]);
S0 = 1 - I0;
Infection = binornd(1,beta,[I^2 T]);
Recovery = binornd(1,gamma,[I T]);

xMat = zeros(Ntype,T);
yMat = zeros(Ntype,T);
zMat = zeros(Ntype,T);

%% simulation

for j = 1:Ntype % iterate over network type
    % set parameter of network
    if strcmp(type{j},'ERG')
        param = D/(I-1);
    elseif strcmp(type{j},'WS')
        K = round(D/2);
        param = [K p_rewire];
    elseif strcmp(type{j},'BA')
        m = round((I*D-m0*(m0-1))/(2*(I-m0)));
        param = [m0 m];
    end
    
    tic
    Network = getNetwork(I,param,type{j}); % network matrix
    S_init = S0; % initial value
    I_init = I0;
    for k = 1:length(N1vec)
        Capacity1 = N1vec(k)*ones(I,1); % capacity vector during social distancing
        [St,It,Rt] = SIR_network_sim3(Infection(:,T1*(k-1)+1:T1*k),Recovery(:,T1*(k-1)+1:T1*k),Network,Capacity1,S_init,I_init,T1);
        xMat(j,T1*(k-1)+1:T1*k) = sum(St,1)/I; % fraction susceptible
        yMat(j,T1*(k-1)+1:T1*k) = sum(It,1)/I; % fraction infected
        zMat(j,T1*(k-1)+1:T1*k) = sum(Rt,1)/I; % fraction recovered
        Capacity0 = N0*ones(I,1); % capacity vector after social distancing
        [St,It,Rt] = SIR_network_sim3(Infection(:,T1*(k-1)+1:T1*k),Recovery(:,T1*(k-1)+1:T1*k),Network,Capacity1,S_init,I_init,T1);
        S_init = St(:,end); % update initial value
        I_init = It(:,end);
    end
    [St,It,Rt] = SIR_network_sim3(Infection(:,T1*length(N1vec)+1:end),Recovery(:,T1*length(N1vec)+1:end),Network,Capacity0,S_init,I_init,T-T1*length(N1vec));
    xMat(j,T1*length(N1vec)+1:end) = sum(St,1)/I; % fraction susceptible
    yMat(j,T1*length(N1vec)+1:end) = sum(It,1)/I; % fraction infected
    zMat(j,T1*length(N1vec)+1:end) = sum(Rt,1)/I; % fraction recovered
    toc
end

%% plot results

y = yMat;
ymax = ceil(100*max(max(yMat))) + 1;
    
figure
plot(time,100*y);
xline(T1*length(N1vec));
text(T1*length(N1vec)+1,ymax,'Social distancing ends','VerticalAlignment','top');
ylim([0 ymax])
xlabel('Time')
ylabel('Infection rate (\%)')
title('Gradual relaxation of social distancing')
legend(typeFull{1},typeFull{2},typeFull{3},'Location','best')
    
%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_Gradual','-dpdf')
