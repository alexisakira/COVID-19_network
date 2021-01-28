function AdjacencyMatrix = getNetwork(I,param,type)

if nargin < 3
    type = 'ERG'; % default is Erdos-Renyi-Gilbert network
end

if strcmp(type,'ERG') % Erdos-Renyi-Gilbert network
    p = param(1); % connection probability
    % some error checking
    if (p<0)||(p>1)
        error('p must be between 0 and 1')
    end
    temp = binornd(1,p,[I I]); % random matrix
    temp = triu(temp,1); % upper triangular part
    AdjacencyMatrix = temp + temp'; % make connections equivalence relation
    
elseif strcmp(type,'WS') % Watts-Strogatz network
    K = param(1); % distance to adjacent neighors to connect
    beta = param(2); % rewiring probability
    % some error checking
    if (rem(K,1) ~= 0)||(K<1)
        error('K must be a natural number')
    end
    if (beta<0)||(beta>1)
        error('beta must be between 0 and 1')
    end
    % Below, I just copy WattsStrogatz.m function in Matlab.
    
    % Connect each node to its K next and previous neighbors. This constructs
    % indices for a ring lattice.
    s = repelem((1:I)',1,K);
    t = s + repmat(1:K,I,1);
    t = mod(t-1,I)+1;

    % Rewire the target node of each edge with probability beta
    for source=1:I   
        switchEdge = rand(K, 1) < beta;
    
        newTargets = rand(I, 1);
        newTargets(source) = 0;
        newTargets(s(t==source)) = 0;
        newTargets(t(source, ~switchEdge)) = 0;
    
        [~, ind] = sort(newTargets, 'descend');
        t(source, switchEdge) = ind(1:nnz(switchEdge));
    end
    G = graph(s,t);
    AdjacencyMatrix = full(adjacency(G));
    
elseif strcmp(type,'BA') % Barabasi-Albert network
    m0 = param(1); % initial number of nodes
    if m0 >= I
        AdjacencyMatrix = ones(I) - eye(I);
        return
    end
    m = param(2); % number of edges a new node have
    % some error checking
    if (rem(m0,1) ~= 0)||(m0 < 1)
        error('m0 must be a natural number')
    end
    if (rem(m,1) ~= 0)||(m < 1)
        error('m must be a natural number')
    end
    A = ones(m0) - eye(m0); % initial adjacency matrix
    for j = m0+1:I % loop over new nodes
        Nnodes = size(A,1);
        degree = sum(A,2); % vector of degrees of existing nodes
        prob = degree'/sum(degree); % vector of attachment probabilities
        cum_prob = cumsum(prob); % cumulative probabilities
        temp = rand(m,1); % uniform random numbers
        target = sum(bsxfun(@minus,temp,cum_prob)>0,2)+1; % target to attach
        target = unique(min(target,Nnodes)); % not necessary but make sure well-defined
        Anew = [A zeros(Nnodes,1); zeros(1,Nnodes+1)];
        Anew(target,end) = 1; % create edge
        temp = triu(Anew,1); % upper triangular part
        A = temp + temp'; % make connections equivalence relation
    end
    AdjacencyMatrix = A;
else
    error('type of network must be either ER, WS, or BA')
end

end

