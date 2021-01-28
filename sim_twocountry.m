%% simulate two country model

setParam; % set parameters
p_connect = 1/(10*I); % probability of connecting to a foreigner
rng(0);
A_foreign = binornd(1,p_connect,[I I]); % adjacency matrix for foreigners

I0 = binornd(1,y0,[2*I 1]);
S0 = 1 - I0;
Infection = binornd(1,beta,[4*I^2 T]);
Recovery = binornd(1,gamma,[2*I T]);
x1Mat = zeros(Nspec,Ntype,T);
y1Mat = zeros(Nspec,Ntype,T);
z1Mat = zeros(Nspec,Ntype,T);
x2Mat = zeros(Nspec,Ntype,T);
y2Mat = zeros(Nspec,Ntype,T);
z2Mat = zeros(Nspec,Ntype,T);

for j = 1:length(type) % iterate over network type
    % set parameter of network
    if strcmp(type{j},'ER')
        param = D/(I-1);
    elseif strcmp(type{j},'WS')
        K = round(D/2);
        param = [K p_rewire];
    elseif strcmp(type{j},'BA')
        m = round((I*D-m0*(m0-1))/(2*(I-m0)));
        param = [m0 m];
    end
    
    Network1 = getNetwork(I,param,type{j}); % network matrix for country 1
    Network2 = getNetwork(I,param,type{j}); % network matrix for country 2
    Network = [Network1 A_foreign;A_foreign' Network2]; % international network matrix
    
    tic
    for s = 1:Nspec % iterate over specifications
        if s == 1 % benchmark model
            % solve benchmark model
            Capacity0 = N0*ones(2*I,1); % capacity vector without social distancing
            [St,It,Rt] = SIR_network_sim3(Infection,Recovery,Network,Capacity0,S0,I0,T);        
        else
            % capacity vector during social distancing
            Capacity1 = [N1vec(s-1)*ones(I,1); N0*ones(I,1)]; % country 2 doesn't do social distancing
            [St1,It1,Rt1] = SIR_network_sim3(Infection(:,1:T1),Recovery(:,1:T1),Network,Capacity1,S0,I0,T1);
            [St0,It0,Rt0] = SIR_network_sim3(Infection(:,T1+1:end),Recovery(:,T1+1:end),Network,Capacity0,St1(:,end),It1(:,end),T-T1);
            St = [St1 St0];
            It = [It1 It0];
            Rt = [Rt1 Rt0];
        end
        x1Mat(s,j,:) = sum(St(1:I,:),1)/I; % fraction susceptible
        y1Mat(s,j,:) = sum(It(1:I,:),1)/I; % fraction infected
        z1Mat(s,j,:) = sum(Rt(1:I,:),1)/I; % fraction recovered
        x2Mat(s,j,:) = sum(St(I+1:end,:),1)/I; % fraction susceptible
        y2Mat(s,j,:) = sum(It(I+1:end,:),1)/I; % fraction infected
        z2Mat(s,j,:) = sum(Rt(I+1:end,:),1)/I; % fraction recovered
    end
    toc
end

clear Infection
clear Recovery
save sim_twocountry

%% plot results

load sim_twocountry
% use below to plot results by Nmax
%{
for s = 1:Nspec % iterate over Nmax
    
    y1 = squeeze(y1Mat(s,:,:)); % infection rate
    y2 = squeeze(y2Mat(s,:,:)); % infection rate
    ymax = ceil(100*max(max([y1;y2]))) + 1;
    
    % Country 1
    figure
    plot(time,100*y1);
    if s > 1
        xline(T1);
        text(T1+1,ymax,'Social distancing ends','VerticalAlignment','top');
    end
    ylim([0 ymax])
    xlabel('Time')
    ylabel('Infection rate (\%)')
    if s == 1
        title('Country 1, $N_{\mathrm{max}}=\infty$')
    else
        title(['Country 1, $N_{\mathrm{max}}=$' num2str(N1vec(s-1))])
    end
    legend(typeFull{1},typeFull{2},typeFull{3},'Location','East')
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    if s==1
        print(fig,'fig_Country1_Ninf','-dpdf')
    else
        print(fig,['fig_Country1_N' num2str(N1vec(s-1))],'-dpdf')
    end
    
    % Country 2
    figure
    plot(time,100*y2);
    if s > 1
        xline(T1);
        text(T1+1,ymax,'Social distancing ends','VerticalAlignment','top');
    end
    ylim([0 ymax])
    xlabel('Time')
    ylabel('Infection rate (\%)')
    if s == 1
        title('Country 2, $N_{\mathrm{max}}=\infty$')
    else
        title(['Country 2, $N_{\mathrm{max}}=$' num2str(N1vec(s-1))])
    end
    legend(typeFull{1},typeFull{2},typeFull{3},'Location','East')
    
     %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    if s==1
        print(fig,'fig_Country2_Ninf','-dpdf')
    else
        print(fig,['fig_Country2_N' num2str(N1vec(s-1))],'-dpdf')
    end
end
%}

% use below to plot results by network type
for j = 1:Ntype % iterate over Nmax
    
    y1 = squeeze(y1Mat(:,j,:)); % infection rate
    y2 = squeeze(y2Mat(:,j,:)); % infection rate
    ymax = ceil(100*max(max([y1;y2]))) + 1;
    
    % Country 1
    figure
    plot(time,100*y1);
    xline(T1);
    text(T1+1,ymax,'Social distancing in country A ends','VerticalAlignment','top');
    ylim([0 ymax])
    xlabel('Time')
    ylabel('Infection rate (\%)')
    title(['Country A, ' typeFull{j}])
    legend('$N_{\mathrm{max}}=\infty$',['$N_{\mathrm{max}}=$' num2str(N1vec(1))],...
        ['$N_{\mathrm{max}}=$' num2str(N1vec(2))],['$N_{\mathrm{max}}=$' num2str(N1vec(3))],...
        'Location','East')
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_CountryA_' type{j}],'-dpdf')
    
    % Country 2
    figure
    plot(time,100*y2);
    xline(T1);
    text(T1+1,ymax,'Social distancing in country A ends','VerticalAlignment','top');
    ylim([0 ymax])
    xlabel('Time')
    ylabel('Infection rate (\%)')
    title(['Country B, ' typeFull{j}])
    legend('$N_{\mathrm{max}}=\infty$',['$N_{\mathrm{max}}=$' num2str(N1vec(1))],...
        ['$N_{\mathrm{max}}=$' num2str(N1vec(2))],['$N_{\mathrm{max}}=$' num2str(N1vec(3))],...
        'Location','East')
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_CountryB_' type{j}],'-dpdf')
end
