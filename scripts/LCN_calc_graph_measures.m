function M = LCN_calc_graph_measures(w,graph_type)
% LCN_calc_graph_measures 
%
% this routines calculates local and global graph measures for weighted or
% binary networks. 
%
% Assuming that the brain connectivity toolbox is in your matlab path
%
% INPUT
%   w    = matrix representing a binary or weighted graph (weights should 
%          be between 0 and 1)
%   graph_type = 'w' (weighted network) or 'b' (binary network). In the latter 
%          case, all nonzero entries of the matrix w will be set to 
%
%__________________________________________________________________________
%
% author: 	Patrick Dupont and Yu Wang
% date: 	February, 2016
% history: 	May 2016: weighted networks, betweenness centrality: correction 
%                     of the input of betweenness_wei. It requires a 
%                     connection-length matrix defined now as 1./w
%           May 2017: the nodal degree is transposed when writing to the 
%                     GTA structure to make it similar to the other measures
%           Feb 2018: correction of bug in calculation of the nodal
%                     characteristic path length
%__________________________________________________________________________
% @(#)LCN_calc_graph_measures.m     0.12          last modified: 2017/05/01

nr_nodes = size(w,1);
if strcmp(graph_type,'w') && min(w(:))>=0 && max(w(:)) <= 1
    % local graph measures
    M.C_nodal      = LCN_clustering_coef_wu(w);
    M.degree_nodal = strengths_und(w)'; 
    M.Eloc_nodal   = LCN_efficiency_wei(w,1);
    M.BC_nodal     = betweenness_wei(1./w)/((nr_nodes-1)*(nr_nodes-2)); % introduced 1./w on 27/5/2016 
     
    D_temp = distance_wei(1./w);
    d1_mean = zeros(1,size(D_temp,1));
    for ii = 1:size(D_temp,1)
        d1_mean(ii) = mean(D_temp(isfinite(D_temp(:,ii)),ii));
    end
    M.lambda_nodal = d1_mean';
   
    % global graph measures
    M.C      = mean(M.C_nodal(:));
    M.E      = LCN_efficiency_wei(w,0);
    M.lambda = charpath(D_temp);
    M.BC     = mean(M.BC_nodal(:));
    M.average_node_degree = mean(M.degree_nodal(:))';
        
    % determine the hubscore
    M.hubscore = LCN_calc_hubscore(M.degree_nodal,M.lambda_nodal,M.C_nodal,M.BC_nodal);

    % determine the community structure
    M.moduleconsistency = LCN_calc_community_structure(w);
elseif strcmp(graph_type,'b')
    % set all nonzero entries of w to 1
    w(w~=0) = 1;
    % local graph measures
    M.C_nodal      = clustering_coef_bu(w);
    M.degree_nodal = degrees_und(w)'; 
    M.Eloc_nodal   = efficiency_bin(w,1);
    M.BC_nodal     = betweenness_bin(w)/((nr_nodes-1)*(nr_nodes-2));
     
    D_temp = distance_bin(w);
    d1_mean = zeros(1,size(D_temp,1));
    for ii = 1:size(D_temp,1)
        d1_mean(ii) = mean(D_temp(isfinite(D_temp(:,ii)),ii));
    end
    M.lambda_nodal = d1_mean';
   
    % global graph measures
    M.C      = mean(M.C_nodal(:));
    M.E      = efficiency_bin(w,0);
    M.lambda = charpath(D_temp);
    M.BC     = mean(M.BC_nodal(:));
    M.average_node_degree = mean(M.degree_nodal(:));
        
    % determine the hubscore
    M.hubscore = LCN_calc_hubscore(M.degree_nodal,M.lambda_nodal,M.C_nodal,M.BC_nodal);

    % determine the community structure
    M.moduleconsistency = LCN_calc_community_structure(w);
elseif strcmp(graph_type,'w') && (min(w(:))<0 || max(w(:)) > 1)
    M = [];
    disp('ERROR: weights should be between 0 and 1')
else
   M = [];
   disp(['ERROR: unknown graph_type ' graph_type ' - graph_type should be either ''w'' or ''b'''])
end
end