function hubscore = LCN_calc_hubscore(NODE_DEGREE,PATH_LENGTH,LOCAL_CLUSTER_COEF,BETWEENNESS_CENTRALITY)
% calculates the hub score 
% the hub score is the sum of dummy values for four criteria.
% We gave a score of 1 or 0 depending on whether or not the node belongs 
% to the top 20% of nodes with 
% 1) the highest node degree (NODE_DEGREE)
% 2) the highest betweenness centrality (BETWEENNESS_CENTRALITY), 
% 3) the lowest local cluster coefficient (limited to nodes with a degree
%    >2) (LOCAL_CLUSTER_COEF)
% 4) the lowest average path length (PATH_LENGTH). 
% Typically, we consider nodes with a hub score >= 2 as hubs.
%
% Only nodes with a node degree of at least 3 (you can change this by 
% adapting the variable min_node_degree) can get a hubscore. The rest
% will have a hubscore of 0.
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	February, 2015
% history: 	March 2015 introduced the variable min_node_degree instead of
%                      hard coded in the calculations
%           June 2018 the min_node_degree is not used in case of weighted
%                     networks
%__________________________________________________________________________
% @(#)LCN_calc_hubscore.m	0.2                   last modified: 2018/06/06

threshold = 20; % in percent
% check if we have a weighted or a binary network (in the latter case the
% node degree is an integer.
if sum(abs(ceil(NODE_DEGREE)-floor(NODE_DEGREE))) == 0 % in this case it is a binary network
   min_node_degree = 2; % min node degree to be taken into account
else
   min_node_degree = 0; 
end
N = length(NODE_DEGREE);
hubscore = zeros(N,1);
index_nodes2 = find(NODE_DEGREE > min_node_degree);
if ~isempty(index_nodes2)
   Nt = ceil(threshold*length(NODE_DEGREE(index_nodes2))/100);
   % determine critical values
   tmp = sort(NODE_DEGREE(index_nodes2),'descend');
   thresh1 = tmp(Nt);
   tmp = sort(PATH_LENGTH(index_nodes2),'ascend');
   thresh2 = tmp(Nt);
   tmp = sort(LOCAL_CLUSTER_COEF(index_nodes2),'ascend');
   thresh3 = tmp(Nt);
   tmp = sort(BETWEENNESS_CENTRALITY(index_nodes2),'descend');
   thresh4 = tmp(Nt);
   hubscore(index_nodes2(NODE_DEGREE(index_nodes2) >= thresh1)) = hubscore(index_nodes2(NODE_DEGREE(index_nodes2) >= thresh1))+1;
   hubscore(index_nodes2(PATH_LENGTH(index_nodes2) <= thresh2)) = hubscore(index_nodes2(PATH_LENGTH(index_nodes2) <= thresh2))+1;
   hubscore(index_nodes2(LOCAL_CLUSTER_COEF(index_nodes2) <= thresh3)) = hubscore(index_nodes2(LOCAL_CLUSTER_COEF(index_nodes2) <= thresh3))+1;
   hubscore(index_nodes2(BETWEENNESS_CENTRALITY(index_nodes2) >= thresh4)) = hubscore(index_nodes2(BETWEENNESS_CENTRALITY(index_nodes2) >= thresh4))+1;
end

end
