function [St,It,Rt] = SIR_network_sim3(Infection,Recovery,Network,Capacity,S0,I0,T)

% beta:         transmission probability
% gamma:        recovery probability
% Network:      network matrix
% Capacity:     capacity vector
% S0, I0:       initial condition
% T:            number of periods to simulate

%% some error checking

I = size(Network,1); % number of agents
if (size(Infection,1) ~= I^2)||(size(Infection,2) ~= T)
    error('Infection must be I^2 x T')
end
if (size(Recovery,1) ~= I)||(size(Infection,2) ~= T)
    error('Recovery must be I x T')
end

if size(Network,2) ~= I
    error('Network must be square');
end
temp = (Network==1) + (Network==0);
if sum(temp(:)) ~= I^2
    error('Network must contain only 0 ,1')
end
if any(diag(Network) ~= 0)
    error('Diagonal of Network must be 0')
end
% uncomment if network need not be symmetric
%{
if any(Network ~= Network')
    error('Network must be symmetric')
end
%}

if any(size(Capacity) ~= [I 1])
    error('Capacity must be I x 1');
end
if any(size(S0) ~= [I 1])
    error('S0 must be I x 1');
end
if any(size(I0) ~= [I 1])
    error('I0 must be I x 1');
end

if sum(S0==1) + sum(S0==0) ~= I
    error('S0 must contain only 0, 1')
end
if sum(I0==1) + sum(I0==0) ~= I
    error('I0 must contain only 0, 1')
end
if any(S0.*I0>0)
    error('Agents cannot be both susceptible and infected')
end

R0 = 1 - S0 - I0; % initially recovered agents

kappa = min(sum(Network,2),Capacity); % vector of binding capacity

St = zeros(I,T);
It = zeros(I,T);
Rt = zeros(I,T);

for t=1:T % iterate over time
    if t==1
        susceptible = S0; % agents that are susceptible in previous preiod
        infected = I0; % agents that are infected
        recovered = R0; % agents that are recovered
    else
        susceptible = St(:,t-1);
        infected = It(:,t-1);
        recovered = Rt(:,t-1);
    end
    indicator_Infected = zeros(I); % indicator whether i infects j at time t
    for i=1:I % iterate over (potentially) infected agents
        ind_meet = find(Network(i,:)); % index of agents that i can meet in network
        ind_meet = datasample(ind_meet,kappa(i),'Replace',false); % index of agents that i actually meets
        temp = Infection(I*(i-1)+1:I*i,t); % indicator whether i infects others conditional on i infected
        indicator_Infected(i,ind_meet) = infected(i)*(temp(ind_meet).*susceptible(ind_meet))';
        % transmission occurs if i is infected, i meets j, and coin flip is infect
    end
    indicator_Infected = min(sum(indicator_Infected,1),1); % one gets infected if somebody else infects
    St(((indicator_Infected'==0).*susceptible>0),t) = 1; % susceptible that are not infected remain susceptible
    Rt(:,t) = recovered + Recovery(:,t).*infected;
    It(:,t) = 1 - St(:,t) - Rt(:,t); % accounting
end

end

