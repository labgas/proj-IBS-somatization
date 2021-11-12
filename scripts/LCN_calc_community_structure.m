function  moduleconsistency = LCN_calc_community_structure(w)

% calculates the communitystructure for a (weigthed) network w and returns
% the module consistency, i.e. the probability that two nodes belong to the
% same community
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	February, 2015
% history: 	
%__________________________________________________________________________
% @(#)LCN_calc_community_structure.m	0.1       last modified: 2015/02/08

nr_randomizations_modularity = 1000;
nr_nodes = size(w,1);

tmp = zeros(nr_nodes,nr_randomizations_modularity);
for i = 1:nr_randomizations_modularity         
    [tmp(:,i) Q]=modularity_und(w,1);
end
moduleconsistency = zeros(nr_nodes,nr_nodes);
for i = 1:nr_nodes-1
    for j = i+1:nr_nodes
        for k = 1:nr_randomizations_modularity
            if tmp(i,k) == tmp(j,k)
               moduleconsistency(i,j) = moduleconsistency(i,j) + 1;         
            end
        end
        moduleconsistency(j,i) = moduleconsistency(i,j);
    end
end
moduleconsistency = moduleconsistency./nr_randomizations_modularity;

end