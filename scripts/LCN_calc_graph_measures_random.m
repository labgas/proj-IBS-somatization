function Mrandom = LCN_calc_graph_measures_random(w,nr_randomizations,graph_type)
% LCN_calc_graph_measures 
%
% this routines calculates local and global graph measures for random weighted or
% binary undirected networks with no self connections and with the same number of 
% nodes and weight distribution as the original network w.
%
% Assuming that the brain connectivity toolbox is in your matlab path
%
% INPUT
%   w    = matrix representing a binary or weighted graph (weights should 
%          be between 0 and 1)
%   nr_randomizations = number of randomizations
%   graph_type = 'w' (weighted network) or 'b' (binary network). In the latter 
%          case, all nonzero entries of the matrix w will be set to 
%
%__________________________________________________________________________
%
% author: 	Patrick Dupont and Yu Wang
% date: 	February, 2016
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc_graph_measures.m     0.1           last modified: 2016/02/24

nr_nodes = size(w,1);
distribution = w(tril(ones(nr_nodes,nr_nodes),-1) ~= 0);
% initialize
n                   = length(distribution);

if strcmp(graph_type,'w') && min(w(:))>=0 && max(w(:)) <= 1
    for i = 1:nr_randomizations
        clear w_random neworder
        
        neworder = randperm(n);
        w_random = zeros(nr_nodes,nr_nodes);
        w_random(tril(ones(nr_nodes,nr_nodes),-1) ~= 0) = distribution(neworder);
        w_random = w_random + w_random';
        
        disp(['        calculating random graph measures ' num2str(i) ' of ' num2str(nr_randomizations)]);
        Mrandom(i) = LCN_calc_graph_measures(w_random,'w');
    end
elseif strcmp(graph_type,'b')
    % set all nonzero entries of w to 1
    w(w~=0) = 1;
    distribution = w(tril(ones(nr_nodes,nr_nodes),-1) ~= 0);
    for i = 1:nr_randomizations
        clear w_random neworder
        
        neworder = randperm(n);
        w_random = zeros(nr_nodes,nr_nodes);
        w_random(tril(ones(nr_nodes,nr_nodes),-1) ~= 0) = distribution(neworder);
        w_random = w_random + w_random';
   
        disp(['        calculating random graph measures ' num2str(i) ' of ' num2str(nr_randomizations)]);
        Mrandom(i) = LCN_calc_graph_measures(w_random,'b');        
    end
elseif strcmp(graph_type,'w') && (min(w(:))<0 || max(w(:)) > 1)
    M = [];
    disp('ERROR: weights should be between 0 and 1')
else
   M = [];
   disp(['ERROR: unknown graph_type ' graph_type ' - graph_type should be either ''w'' or ''b'''])
end
end