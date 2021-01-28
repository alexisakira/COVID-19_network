setParam; % set parameters

N1vec = [1 2 3 4 5 10];
Nspec = length(N1vec)+1;

Nsim = 1000; % number of simulations

type = {'WS'}; % use either 'ER', 'WS', or 'BA'
Ntype = length(type);

%% large simulation

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
    
    xMat = zeros(Nspec,Nsim,T);
    yMat = zeros(Nspec,Nsim,T);
    zMat = zeros(Nspec,Nsim,T);
    ppm = ParforProgressbar(Nsim,'title',['Simulating ' type{j}]);
    parfor i = 1:Nsim % iterate over simulations; parallel computing
        
        rng(i); % set seed for replicability
        % set initial condition
        I0 = binornd(1,y0,[I 1]);
        S0 = 1 - I0;
        % draw infections and recovery; inefficient but necessary for
        % apple-to-apple comparison
        Infection = binornd(1,beta,[I^2 T]);
        Recovery = binornd(1,gamma,[I T]);
        
        Network = getNetwork(I,param,type{j}); % network matrix
        for s = 1:Nspec % iterate over specifications
            %h = waitbar((Nspec*(i-1)+s)/(Nsim*Nspec),['Simulating ' typeFull{j}]);
            %tic
            if s == 1 % benchmark model
                rng(i);
                % solve benchmark model
                [St,It,Rt] = SIR_network_sim3(Infection,Recovery,Network,Capacity0,S0,I0,T);        
            else
                rng(i);
                Capacity1 = N1vec(s-1)*ones(I,1);
                [St1,It1,Rt1] = SIR_network_sim3(Infection(:,1:T1),Recovery(:,1:T1),Network,Capacity1,S0,I0,T1);
                [St0,It0,Rt0] = SIR_network_sim3(Infection(:,T1+1:end),Recovery(:,T1+1:end),Network,Capacity0,St1(:,end),It1(:,end),T-T1);
                St = [St1 St0];
                It = [It1 It0];
                Rt = [Rt1 Rt0];
            end
            xMat(s,i,:) = sum(St,1)/I; % fraction susceptible
            yMat(s,i,:) = sum(It,1)/I; % fraction infected
            zMat(s,i,:) = sum(Rt,1)/I; % fraction recovered
            %toc
            %close(h); % close waitbar
        end
        %clear Infection % clear to save memory
        %clear Recovery
        ppm.increment(); 
    end
    delete(ppm);
    clear Infection
    clear Recovery
    save(['sim_large_' type{j} '.mat']); % save results in .mat file
end
