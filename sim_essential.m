%% simulate model with essential workers

setParam; % set parameters
p_essential = 0.1; % fraction of essential workers
N_essential = round(p_essential*I); % number of essential workers

rng(1);
ind_essential = datasample([1:I],N_essential,'Replace',false); % index of essential agents

rng(1);
I0 = binornd(1,y0,[I 1]);
S0 = 1 - I0;
Infection = binornd(1,beta,[I^2 T]);
Recovery = binornd(1,gamma,[I T]);
xMat = zeros(Nspec,Ntype,T);
yMat = zeros(Nspec,Ntype,T);
zMat = zeros(Nspec,Ntype,T);

for j = 1:length(type) % iterate over network type
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
    
    Network = getNetwork(I,param,type{j}); % network matri for country 1
    
    tic
    for s = 1:Nspec % iterate over specifications
        if s == 1 % benchmark model
            % solve benchmark model
            Capacity0 = N0*ones(I,1); % capacity vector without social distancing
            [St,It,Rt] = SIR_network_sim3(Infection,Recovery,Network,Capacity0,S0,I0,T);        
        else
            Capacity1 = N1vec(s-1)*ones(I,1); % capacity vector during social distancing
            Capacity1(ind_essential) = N1vec(end); % essential workers get mild restrictions
            [St1,It1,Rt1] = SIR_network_sim3(Infection(:,1:T1),Recovery(:,1:T1),Network,Capacity1,S0,I0,T1);
            [St0,It0,Rt0] = SIR_network_sim3(Infection(:,T1+1:end),Recovery(:,T1+1:end),Network,Capacity0,St1(:,end),It1(:,end),T-T1);
            St = [St1 St0];
            It = [It1 It0];
            Rt = [Rt1 Rt0];
        end
        xMat(s,j,:) = sum(St,1)/I; % fraction susceptible
        yMat(s,j,:) = sum(It,1)/I; % fraction infected
        zMat(s,j,:) = sum(Rt,1)/I; % fraction recovered
    end
    toc
end

clear Infection
clear Recovery
save sim_essential

%% plot results

load sim_essential

for s = 1:Nspec % iterate over Nmax
    
    y = squeeze(yMat(s,:,:));
    ymax = ceil(100*max(max(y))) + 1;
    
    figure
    plot(time,100*y);
    xline(T1);
    text(T1+1,ymax,'Social distancing ends','VerticalAlignment','top');
    ylim([0 ymax])
    xlabel('Time')
    ylabel('Infection rate (\%)')
    if s == 1
        title('Essential workers, $N_{\mathrm{max}}=\infty$')
    else
        title(['Essential workers, $N_{\mathrm{max}}=$' num2str(N1vec(s-1))])
    end
    legend(typeFull{1},typeFull{2},typeFull{3},'Location','East')
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    if s==1
        print(fig,'fig_Essential_Ninf','-dpdf')
    else
        print(fig,['fig_Essential_N' num2str(N1vec(s-1))],'-dpdf')
    end
end


for j = 1:Ntype % iterate over network type
    
    y = squeeze(yMat(:,j,:));
    ymax = ceil(100*max(max(y))) + 1;
    
    figure
    plot(time,100*y);
    xline(T1);
    text(T1+1,ymax,'Social distancing ends','VerticalAlignment','top');
    ylim([0 ymax])
    xlabel('Time')
    ylabel('Infection rate (\%)')
    title(['Essential workers, ' typeFull{j}])
    legend('$N_{\mathrm{max}}=\infty$',['$N_{\mathrm{max}}=$' num2str(N1vec(1))],...
        ['$N_{\mathrm{max}}=$' num2str(N1vec(2))],['$N_{\mathrm{max}}=$' num2str(N1vec(3))],'Location','NE')
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_Essential_' type{j}],'-dpdf')
end
