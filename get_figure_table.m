%% plot figures for large simulations

setParam % set parameters
N1vec = [1 2 3 4 5 10];
Nspec = length(N1vec)+1;

yMaxMean = zeros(Ntype,Nspec); % mean of peak infection rate
yMaxStd = zeros(Ntype,Nspec); % standard deviation
yMaxMed = zeros(Ntype,Nspec); % median
yMaxMeanImmune = zeros(Ntype,Nspec); % mean of peak infection conditional on herd immunity
yMaxMeanNoimmune = zeros(Ntype,Nspec); % mean conditional on no herd immunity
herdImmunity = zeros(Ntype,Nspec); % fraction of herd immunity
zEnd = zeros(Ntype,Nspec); % fraction of recovered

ind_plot1_4 = [2 3 4 5]; % index used for plotting Nmax = 1-4
ind_plot2_10 = [1 3 6 7]; % index used for plotting Nmax = 2-10

for iter = 1:Ntype % change to Ntype later
    load(['sim_' type{iter} '1000.mat']) % load data
    type = {'ERG','WS','BA'}; % type of networks
    
    N1vec = [1 2 3 4 5 10];
    Nspec = length(N1vec)+1;
    
    % compute summary statistics
    yMaxMat = 100*squeeze(max(yMat,[],3)); % peak
    yMaxMean(iter,:) = round(mean(yMaxMat,2),1);
    yMaxStd(iter,:) = round(std(yMaxMat,0,2),1);
    yMaxMed(iter,:) = round(median(yMaxMat,2),1);
    xMatEnd = squeeze(xMat(:,:,end));
    zMatEnd = squeeze(zMat(:,:,end));
    for s = 1:Nspec
        ind_immune = find(beta*D*xMatEnd(s,:) <= gamma);
        ind_noimmune = find(beta*D*xMatEnd(s,:) > gamma);
        yMaxMeanImmune(iter,s) = round(mean(yMaxMat(s,ind_immune)),1);
        yMaxMeanNoimmune(iter,s) = round(mean(yMaxMat(s,ind_noimmune)),1);
    end
    herdImmunity(iter,:) = round(100*sum(beta*D*xMatEnd <= gamma,2)/Nsim,1);
    zEnd(iter,:) = round(100*mean(zMatEnd,2),1);
    % plot histograms for Nmax = 2-10
    figure; hold on
    for s = 1:length(ind_plot2_10)
        histogram(yMaxMat(ind_plot2_10(s),:),[0:max(max(yMaxMat))],'Normalization','pdf');
    end
    hold off
    title(['1000 simulations with ' typeFull{iter}])
    xlabel('Peak infection rate (\%)')
    ylabel('Probability density')
    legend('$N_{\mathrm{max}}=\infty$',['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot2_10(2)-1))],...
        ['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot2_10(3)-1))],['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot2_10(4)-1))])
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_Hist1000_' type{iter} '2-10'],'-dpdf')
    
    % plot histograms for Nmax = 1-4
    figure; hold on
    for s = 1:length(ind_plot1_4)
        histogram(yMaxMat(ind_plot1_4(s),:),[0:max(max(yMaxMat))],'Normalization','pdf');
    end
    hold off
    title(['1000 simulations with ' typeFull{iter}])
    xlabel('Peak infection rate (\%)')
    ylabel('Probability density')
    legend(['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot1_4(1)-1))],['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot1_4(2)-1))],...
        ['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot1_4(3)-1))],['$N_{\mathrm{max}}=$' num2str(N1vec(ind_plot1_4(4)-1))])
    
    %save figure in pdf format
    fig = gcf;
    fig.PaperPositionMode = 'auto';
    fig_pos = fig.PaperPosition;
    fig.PaperSize = [fig_pos(3) fig_pos(4)];
    print(fig,['fig_Hist1000_' type{iter} '1-4'],'-dpdf')
end

setParam % set parameters
N1vec = [1 2 3 4 5 10];
Nspec = length(N1vec)+1;
Reff_mean = zeros(Ntype,T);
Reff_med = zeros(Ntype,T);

for iter = 1:Ntype % change to Ntype later
    load(['sim_' type{iter} '1000.mat']) % load data
    type = {'ERG','WS','BA'}; % type of networks
    
    N1vec = [1 2 3 4 5 10];
    Nspec = length(N1vec)+1;
    
    IMat = [y0*ones(I,1) max(squeeze(yMat(1,:,:)),1e-9)];
    % fraction infected in benchmark model (make positive by taking max with small number)
    ReffMat = 1 + diff(log(IMat),1,2)/gamma; % effective reproduction number
    Reff_mean(iter,:) = mean(ReffMat,1);
    Reff_med(iter,:) = median(ReffMat,1);
end

T1 = 60;

figure
plot(time,Reff_med)
yline(R0);
text(T1-1,R0,'$\mathcal{R}_0$ (basic reproduction number)','VerticalAlignment','top','HorizontalAlignment','right');
xlim([0 T1])
ylim([0 R0+1])
%title('Median effective reproduction number')
xlabel('Time')
ylabel('Reproduction number')
legend(typeFull{1},typeFull{2},typeFull{3},'Location','NE')

%save figure in pdf format
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,'fig_R0','-dpdf')