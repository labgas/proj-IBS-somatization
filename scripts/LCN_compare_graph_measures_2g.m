function LCN_compare_graph_measures_2g(GTA1,GTA2,GTA1r,GTA2r,pvalue_threshold)

% this routine calculates the significance of the difference in group based
% graph measures based on random permutation labeling. If random graph
% measures are available this will be done for the normalized graph
% measures as well since the value of graph measures for networks of the
% same size (= same number of nodes) but different density can differ as a
% result of density as such. 
%
% INPUT:
%   GTA1 = a 1 x 1 structure containing the fields graph_measures, 
%          graph_measures_random and nodenames for group 1 or condition 1
%   GTA2 = a 1 x 1 structure containing the fields graph_measures and
%          graph_measures_random for group 2 or condition 2
%   GTA1r = a 1 x n structure containing the fields graph_measures and
%          optional graph_measures_random using random permutation
%          labelling of group membership or condition membership (paired
%          case) for a random group 1 or condition 1
%   GTA2r = a 1 x n structure containing the fields graph_measures and
%          optional graph_measures_random using random permutation
%          labelling of group membership or condition membership (paired
%          case) for a random group 2 or condition 2
%   pvalue_threshold = only value with p < pvalue_threshold will be
%                      displayed
%__________________________________________________________________________
%
% author: 	Patrick Dupont
% date: 	May, 2016
% history: 	January 2018 analysis of local graph measures, hubscores and
%                        community structure
%           May 2018 will also work for weigthed networks
%__________________________________________________________________________
% @(#)LCN_compare_graph_measures_2g.m  v0.2 
% last modified: 2018/05/31

min_hubscore = 2;

nodenames = GTA1(1).nodenames;
nr_nodes = length(nodenames);

if isfield(GTA1,'graph_measures_random')
   nr_randomizations1 = size(GTA1(1).graph_measures_random,2);
end
if isfield(GTA2,'graph_measures_random')
   nr_randomizations2 = size(GTA2(1).graph_measures_random,2);
end
nr_randomizations_pl = size(GTA1r,2); % number of permutation labelings

if isfield(GTA1r,'graph_measures_random')
   nr_randomizations_pl_norm = size(GTA1r(1).graph_measures_random,2);
end

fprintf('Results for the comparison of group/condition 1 minus group/condition 2\n')
fprintf('-----------------------------------------------------------------------\n');
fprintf('Graph measure \t mean random distribution \t std random distribution \t value1 \t value2 \t actual difference \t p-value \n');

% clustering coefficient
%-----------------------
clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
clear M1 M2 M1n M2n M_dif12 Mn_dif12
distribution1 = zeros(1,nr_randomizations1);
distribution2 = zeros(1,nr_randomizations2);
for i = 1:nr_randomizations1
    distribution1(i) = GTA1(1).graph_measures_random(i).C; 
end
for i = 1:nr_randomizations2
    distribution2(i) = GTA2(1).graph_measures_random(i).C; 
end
M_dif12  = GTA1(1).graph_measures.C - GTA2(1).graph_measures.C;
M1n      = GTA1(1).graph_measures.C./mean(distribution1);
M2n      = GTA2(1).graph_measures.C./mean(distribution2);
Mn_dif12 = M1n - M2n;
distribution_dif12 = zeros(1,nr_randomizations_pl);
for i = 1:nr_randomizations_pl
    distribution_dif12(i) = GTA1r(i).graph_measures.C - GTA2r(i).graph_measures.C;
end
clear p_value_left p_value_right
[p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
p = min(p_value_left,p_value_right);
% if p < pvalue_threshold
   fprintf('clustering coefficient \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.C,GTA2(1).graph_measures.C,M_dif12,p);
% end
if isfield(GTA1r,'graph_measures_random')
   distribution_n_dif12 = zeros(1,nr_randomizations_pl);
   for i = 1:nr_randomizations_pl
       tmp1 = zeros(1,nr_randomizations_pl_norm);
       tmp2 = zeros(1,nr_randomizations_pl_norm);
       for j = 1:nr_randomizations_pl_norm
           tmp1(j) = GTA1r(i).graph_measures_random(j).C;
       end 
       for j = 1:nr_randomizations_pl_norm
           tmp2(j) = GTA2r(i).graph_measures_random(j).C;
       end
       distribution_n_dif12(i) = GTA1r(i).graph_measures.C./mean(tmp1) - GTA2r(i).graph_measures.C./mean(tmp2);
   end
   clear p_value_left p_value_right
   [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
   p = min(p_value_left,p_value_right);
%    if p < pvalue_threshold
      fprintf('normalized clustering coefficient \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
%    end
end

% characteristic path length
%-----------------------
clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
clear M1 M2 M1n M2n M_dif12 Mn_dif12
distribution1 = zeros(1,nr_randomizations1);
distribution2 = zeros(1,nr_randomizations2);
for i = 1:nr_randomizations1
    distribution1(i) = GTA1(1).graph_measures_random(i).lambda; 
end
for i = 1:nr_randomizations2
    distribution2(i) = GTA2(1).graph_measures_random(i).lambda; 
end
M_dif12  = GTA1(1).graph_measures.lambda - GTA2(1).graph_measures.lambda;
M1n      = GTA1(1).graph_measures.lambda./mean(distribution1);
M2n      = GTA2(1).graph_measures.lambda./mean(distribution2);
Mn_dif12 = M1n - M2n;
distribution_dif12 = zeros(1,nr_randomizations_pl);
for i = 1:nr_randomizations_pl
    distribution_dif12(i) = GTA1r(i).graph_measures.lambda - GTA2r(i).graph_measures.lambda;
end
clear p_value_left p_value_right
[p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
p = min(p_value_left,p_value_right);
% if p < pvalue_threshold
  fprintf('characteristic path length \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.lambda,GTA2(1).graph_measures.lambda,M_dif12,p);
% end
if isfield(GTA1r,'graph_measures_random')
   distribution_n_dif12 = zeros(1,nr_randomizations_pl);
   for i = 1:nr_randomizations_pl
       tmp1 = zeros(1,nr_randomizations_pl_norm);
       tmp2 = zeros(1,nr_randomizations_pl_norm);
       for j = 1:nr_randomizations_pl_norm
           tmp1(j) = GTA1r(i).graph_measures_random(j).lambda;
       end 
       for j = 1:nr_randomizations_pl_norm
           tmp2(j) = GTA2r(i).graph_measures_random(j).lambda;
       end
       distribution_n_dif12(i) = GTA1r(i).graph_measures.lambda./mean(tmp1) - GTA2r(i).graph_measures.lambda./mean(tmp2);
   end
   clear p_value_left p_value_right
   [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
   p = min(p_value_left,p_value_right);
%    if p < pvalue_threshold
      fprintf('normalized characteristic path length \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
%    end
end

% global efficiency
%-----------------------
clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
clear M1 M2 M1n M2n M_dif12 Mn_dif12
distribution1 = zeros(1,nr_randomizations1);
distribution2 = zeros(1,nr_randomizations2);
for i = 1:nr_randomizations1
    distribution1(i) = GTA1(1).graph_measures_random(i).E; 
end
for i = 1:nr_randomizations2
    distribution2(i) = GTA2(1).graph_measures_random(i).E; 
end
M_dif12  = GTA1(1).graph_measures.E - GTA2(1).graph_measures.E;
M1n      = GTA1(1).graph_measures.E./mean(distribution1);
M2n      = GTA2(1).graph_measures.E./mean(distribution2);
Mn_dif12 = M1n - M2n;
distribution_dif12 = zeros(1,nr_randomizations_pl);
for i = 1:nr_randomizations_pl
    distribution_dif12(i) = GTA1r(i).graph_measures.E - GTA2r(i).graph_measures.E;
end
clear p_value_left p_value_right
[p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
p = min(p_value_left,p_value_right);
% if p < pvalue_threshold
   fprintf('global efficiency \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.E,GTA2(1).graph_measures.E,M_dif12,p);
% end
if isfield(GTA1r,'graph_measures_random')
   distribution_n_dif12 = zeros(1,nr_randomizations_pl);
   for i = 1:nr_randomizations_pl
       tmp1 = zeros(1,nr_randomizations_pl_norm);
       tmp2 = zeros(1,nr_randomizations_pl_norm);
       for j = 1:nr_randomizations_pl_norm
           tmp1(j) = GTA1r(i).graph_measures_random(j).E;
       end 
       for j = 1:nr_randomizations_pl_norm
           tmp2(j) = GTA2r(i).graph_measures_random(j).E;
       end
       distribution_n_dif12(i) = GTA1r(i).graph_measures.E./mean(tmp1) - GTA2r(i).graph_measures.E./mean(tmp2);
   end
   clear p_value_left p_value_right
   [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
   p = min(p_value_left,p_value_right);
%    if p < pvalue_threshold
      fprintf('normalized global efficiency \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
%    end
end

% betweenness centrality
%-----------------------
clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
clear M1 M2 M1n M2n M_dif12 Mn_dif12
distribution1 = zeros(1,nr_randomizations1);
distribution2 = zeros(1,nr_randomizations2);
for i = 1:nr_randomizations1
    distribution1(i) = GTA1(1).graph_measures_random(i).BC; 
end
for i = 1:nr_randomizations2
    distribution2(i) = GTA2(1).graph_measures_random(i).BC; 
end
M_dif12  = GTA1(1).graph_measures.BC - GTA2(1).graph_measures.BC;
M1n      = GTA1(1).graph_measures.BC./mean(distribution1);
M2n      = GTA2(1).graph_measures.BC./mean(distribution2);
Mn_dif12 = M1n - M2n;
distribution_dif12 = zeros(1,nr_randomizations_pl);
for i = 1:nr_randomizations_pl
    distribution_dif12(i) = GTA1r(i).graph_measures.BC - GTA2r(i).graph_measures.BC;
end
clear p_value_left p_value_right
[p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
p = min(p_value_left,p_value_right);
if p < pvalue_threshold
   fprintf('betweenness centrality \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.BC,GTA2(1).graph_measures.BC,M_dif12,p);
end
if isfield(GTA1r,'graph_measures_random')
   distribution_n_dif12 = zeros(1,nr_randomizations_pl);
   for i = 1:nr_randomizations_pl
       tmp1 = zeros(1,nr_randomizations_pl_norm);
       tmp2 = zeros(1,nr_randomizations_pl_norm);
       for j = 1:nr_randomizations_pl_norm
           tmp1(j) = GTA1r(i).graph_measures_random(j).BC;
       end 
       for j = 1:nr_randomizations_pl_norm
           tmp2(j) = GTA2r(i).graph_measures_random(j).BC;
       end
       distribution_n_dif12(i) = GTA1r(i).graph_measures.BC./mean(tmp1) - GTA2r(i).graph_measures.BC./mean(tmp2);
   end
   clear p_value_left p_value_right
   [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
   p = min(p_value_left,p_value_right);
%    if p < pvalue_threshold
      fprintf('normalized betweenness centrality \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
%    end
end

% NODAL GRAPH MEASURES
% by definition the average node degree is the same for equivalent
% random graphs because of the same weight/connection distribution
fprintf('------------------\n');
fprintf('\n\n'); 
fprintf('Local graph measures \n');
fprintf('---------------------\n');
% node degree
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).degree_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).degree_nodal(node); 
    end
    M_dif12  = GTA1(1).graph_measures.degree_nodal(node) - GTA2(1).graph_measures.degree_nodal(node);
    M1n      = GTA1(1).graph_measures.degree_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.degree_nodal(node)./mean(distribution2);
%     Mn_dif12 = M1n - M2n;
    distribution_dif12 = zeros(1,nr_randomizations_pl);
    for i = 1:nr_randomizations_pl
        distribution_dif12(i) = GTA1r(i).graph_measures.degree_nodal(node) - GTA2r(i).graph_measures.degree_nodal(node);
    end
    clear p_value_left p_value_right
    [p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
    p = min(p_value_left,p_value_right);
    if p < pvalue_threshold && M_dif12~=0
       fprintf('node degree \t %s \t %E \t %E \t %4.2f \t %4.2f \t %4.2f \t %5.4f \n',nodename,mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.degree_nodal(node),GTA2(1).graph_measures.degree_nodal(node),M_dif12,p);
    end
end
% for node = 1:nr_nodes
%     if isfield(GTA1r,'graph_measures_random')
%        distribution_n_dif12 = zeros(1,nr_randomizations_pl);
%        for i = 1:nr_randomizations_pl
%            tmp1 = zeros(1,nr_randomizations_pl_norm);
%            tmp2 = zeros(1,nr_randomizations_pl_norm);
%            for j = 1:nr_randomizations_pl_norm
%                tmp1(j) = GTA1r(i).graph_measures_random(j).degree_nodal(node);
%            end 
%            for j = 1:nr_randomizations_pl_norm
%                tmp2(j) = GTA2r(i).graph_measures_random(j).degree_nodal(node);
%            end
%            distribution_n_dif12(i) = GTA1r(i).graph_measures.degree_nodal(node)./mean(tmp1) - GTA2r(i).graph_measures.degree_nodal(node)./mean(tmp2);
%        end
%        clear p_value_left p_value_right
%        [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
%        p = min(p_value_left,p_value_right);
% %        if p < pvalue_threshold
%           fprintf('normalized node degree \t %s \t %E \t %E \t %E \t %4.2f \t %4.2f \t %5.4f \n',nodename,mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
% %        end
%     end
% end

% nodal clustering coefficient
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).C_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).C_nodal(node); 
    end
    M_dif12  = GTA1(1).graph_measures.C_nodal(node) - GTA2(1).graph_measures.C_nodal(node);
    M1n      = GTA1(1).graph_measures.C_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.C_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    distribution_dif12 = zeros(1,nr_randomizations_pl);
    for i = 1:nr_randomizations_pl
        distribution_dif12(i) = GTA1r(i).graph_measures.C_nodal(node) - GTA2r(i).graph_measures.C_nodal(node);
    end
    clear p_value_left p_value_right
    [p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
    p = min(p_value_left,p_value_right);
    if p < pvalue_threshold && M_dif12~=0
       fprintf('nodal clustering coefficient \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.C_nodal(node),GTA2(1).graph_measures.C_nodal(node),M_dif12,p);
    end
end
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).C_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).C_nodal(node); 
    end
    M1n      = GTA1(1).graph_measures.C_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.C_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    if isfield(GTA1r,'graph_measures_random')
       distribution_n_dif12 = zeros(1,nr_randomizations_pl);
       for i = 1:nr_randomizations_pl
           tmp1 = zeros(1,nr_randomizations_pl_norm);
           tmp2 = zeros(1,nr_randomizations_pl_norm);
           for j = 1:nr_randomizations_pl_norm
               tmp1(j) = GTA1r(i).graph_measures_random(j).C_nodal(node);
           end 
           for j = 1:nr_randomizations_pl_norm
               tmp2(j) = GTA2r(i).graph_measures_random(j).C_nodal(node);
           end
           distribution_n_dif12(i) = GTA1r(i).graph_measures.C_nodal(node)./mean(tmp1) - GTA2r(i).graph_measures.C_nodal(node)./mean(tmp2);
       end
       clear p_value_left p_value_right
       [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
       p = min(p_value_left,p_value_right);
       if p < pvalue_threshold && Mn_dif12~=0
          fprintf('normalized clustering coefficient \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
       end
    end
end

% nodal path length
%------------------
for node = 1:nr_nodes
    nodename = char(nodenames{node});    
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).lambda_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).lambda_nodal(node); 
    end
    M_dif12  = GTA1(1).graph_measures.lambda_nodal(node) - GTA2(1).graph_measures.lambda_nodal(node);
    M1n      = GTA1(1).graph_measures.lambda_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.lambda_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    distribution_dif12 = zeros(1,nr_randomizations_pl);
    for i = 1:nr_randomizations_pl
        distribution_dif12(i) = GTA1r(i).graph_measures.lambda_nodal(node) - GTA2r(i).graph_measures.lambda_nodal(node);
    end
    clear p_value_left p_value_right
    [p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);  
    p = min(p_value_left,p_value_right);
    if p < pvalue_threshold && M_dif12~=0
       fprintf('nodal path length \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.lambda_nodal(node),GTA2(1).graph_measures.lambda_nodal(node),M_dif12,p);
    end
end
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).lambda_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).lambda_nodal(node); 
    end
    M1n      = GTA1(1).graph_measures.lambda_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.lambda_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    if isfield(GTA1r,'graph_measures_random')
       distribution_n_dif12 = zeros(1,nr_randomizations_pl);
       for i = 1:nr_randomizations_pl
           tmp1 = zeros(1,nr_randomizations_pl_norm);
           tmp2 = zeros(1,nr_randomizations_pl_norm);
           for j = 1:nr_randomizations_pl_norm
               tmp1(j) = GTA1r(i).graph_measures_random(j).lambda_nodal(node);
           end 
           for j = 1:nr_randomizations_pl_norm
               tmp2(j) = GTA2r(i).graph_measures_random(j).lambda_nodal(node);
           end
           distribution_n_dif12(i) = GTA1r(i).graph_measures.lambda_nodal(node)./mean(tmp1) - GTA2r(i).graph_measures.lambda_nodal(node)./mean(tmp2);
       end
       clear p_value_left p_value_right
       [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
       p = min(p_value_left,p_value_right);
       if p < pvalue_threshold && Mn_dif12~=0
          fprintf('normalized nodal path length \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
       end
    end
end

% local efficiency
%-----------------------
for node = 1:nr_nodes
    nodename = char(nodenames{node});   
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).Eloc_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).Eloc_nodal(node); 
    end
    M_dif12  = GTA1(1).graph_measures.Eloc_nodal(node) - GTA2(1).graph_measures.Eloc_nodal(node);
    M1n      = GTA1(1).graph_measures.Eloc_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.Eloc_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    distribution_dif12 = zeros(1,nr_randomizations_pl);
    for i = 1:nr_randomizations_pl
        distribution_dif12(i) = GTA1r(i).graph_measures.Eloc_nodal(node) - GTA2r(i).graph_measures.Eloc_nodal(node);
    end
    clear p_value_left p_value_right
    [p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
    p = min(p_value_left,p_value_right);
    if p < pvalue_threshold && M_dif12~=0
       fprintf('local efficiency \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.Eloc_nodal(node),GTA2(1).graph_measures.Eloc_nodal(node),M_dif12,p);
    end
end
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).Eloc_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).Eloc_nodal(node); 
    end
    M1n      = GTA1(1).graph_measures.Eloc_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.Eloc_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    if isfield(GTA1r,'graph_measures_random')
       distribution_n_dif12 = zeros(1,nr_randomizations_pl);
       for i = 1:nr_randomizations_pl
           tmp1 = zeros(1,nr_randomizations_pl_norm);
           tmp2 = zeros(1,nr_randomizations_pl_norm);
           for j = 1:nr_randomizations_pl_norm
               tmp1(j) = GTA1r(i).graph_measures_random(j).Eloc_nodal(node);
           end 
           for j = 1:nr_randomizations_pl_norm
               tmp2(j) = GTA2r(i).graph_measures_random(j).Eloc_nodal(node);
           end
           distribution_n_dif12(i) = GTA1r(i).graph_measures.Eloc_nodal(node)./mean(tmp1) - GTA2r(i).graph_measures.Eloc_nodal(node)./mean(tmp2);
       end
       clear p_value_left p_value_right
       [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
       p = min(p_value_left,p_value_right);
       if p < pvalue_threshold && Mn_dif12~=0
          fprintf('normalized local efficiency \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
       end
    end
end

% nodal betweenness centrality
%-----------------------
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).BC_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).BC_nodal(node); 
    end 
    M_dif12  = GTA1(1).graph_measures.BC_nodal(node) - GTA2(1).graph_measures.BC_nodal(node);
    M1n      = GTA1(1).graph_measures.BC_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.BC_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    distribution_dif12 = zeros(1,nr_randomizations_pl);
    for i = 1:nr_randomizations_pl
        distribution_dif12(i) = GTA1r(i).graph_measures.BC_nodal(node) - GTA2r(i).graph_measures.BC_nodal(node);
    end
    clear p_value_left p_value_right
    [p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
    p = min(p_value_left,p_value_right);
    if p < pvalue_threshold && M_dif12~=0
       fprintf('nodal betweenness centrality \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.BC_nodal(node),GTA2(1).graph_measures.BC_nodal(node),M_dif12,p);
    end
end
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    distribution1 = zeros(1,nr_randomizations1);
    distribution2 = zeros(1,nr_randomizations2);
    for i = 1:nr_randomizations1
        distribution1(i) = GTA1(1).graph_measures_random(i).BC_nodal(node); 
    end
    for i = 1:nr_randomizations2
        distribution2(i) = GTA2(1).graph_measures_random(i).BC_nodal(node); 
    end 
    M1n      = GTA1(1).graph_measures.BC_nodal(node)./mean(distribution1);
    M2n      = GTA2(1).graph_measures.BC_nodal(node)./mean(distribution2);
    Mn_dif12 = M1n - M2n;
    if isfield(GTA1r,'graph_measures_random')
       distribution_n_dif12 = zeros(1,nr_randomizations_pl);
       for i = 1:nr_randomizations_pl
           tmp1 = zeros(1,nr_randomizations_pl_norm);
           tmp2 = zeros(1,nr_randomizations_pl_norm);
           for j = 1:nr_randomizations_pl_norm
               tmp1(j) = GTA1r(i).graph_measures_random(j).BC_nodal(node);
           end 
           for j = 1:nr_randomizations_pl_norm
               tmp2(j) = GTA2r(i).graph_measures_random(j).BC_nodal(node);
           end
           distribution_n_dif12(i) = GTA1r(i).graph_measures.BC_nodal(node)./mean(tmp1) - GTA2r(i).graph_measures.BC_nodal(node)./mean(tmp2);
       end
       clear p_value_left p_value_right
       [p_value_left,p_value_right] = LCN_p_from_distribution(Mn_dif12,distribution_n_dif12,0);
       p = min(p_value_left,p_value_right);
       if p < pvalue_threshold && Mn_dif12~=0
          fprintf('normalized nodal betweenness centrality \t %s \t %E \t %E \t %E \t %E \t  %E \t %5.4f \n',nodename,mean(distribution_n_dif12),std(distribution_n_dif12),M1n,M2n,Mn_dif12,p);
       end
    end
end

% DIFFERENCES IN HUBSCORE
% show the hubs in both groups
fprintf('hubs of group 1\n');
for node = 1:nr_nodes
    if GTA1(1).graph_measures.hubscore(node) >= min_hubscore
       fprintf('hubscore \t %s \t %i \n',char(nodenames{node}),GTA1(1).graph_measures.hubscore(node));
    end
end
fprintf('\n');
fprintf('hubs of group 2\n');
for node = 1:nr_nodes
    if GTA2(1).graph_measures.hubscore(node) >= min_hubscore
       fprintf('hubscore \t %s \t %i \n',char(nodenames{node}),GTA2(1).graph_measures.hubscore(node));
    end
end
fprintf('\n');
fprintf('comparison of hubs of group 1 versus group 2\n');
for node = 1:nr_nodes
    nodename = char(nodenames{node});
    clear distribution1 distribution2 distribution1n distribution2n distribution_dif12 distribution_n_dif12
    clear M1 M2 M1n M2n M_dif12 Mn_dif12
    M_dif12  = GTA1(1).graph_measures.hubscore(node) - GTA2(1).graph_measures.hubscore(node);
    distribution_dif12 = zeros(1,nr_randomizations_pl);
    for i = 1:nr_randomizations_pl
        distribution_dif12(i) = GTA1r(i).graph_measures.hubscore(node) - GTA2r(i).graph_measures.hubscore(node);
    end
    clear p_value_left p_value_right
    [p_value_left,p_value_right] = LCN_p_from_distribution(M_dif12,distribution_dif12,0);
    p = min(p_value_left,p_value_right);
    if p < pvalue_threshold && M_dif12~=0
       fprintf('hubscore \t %s \t %4.2f \t %4.2f \t %i \t %i \t  %i \t %5.4f \n',nodename,mean(distribution_dif12),std(distribution_dif12),GTA1(1).graph_measures.hubscore(node),GTA2(1).graph_measures.hubscore(node),M_dif12,p);
    end
end

% DIFFERENCES IN COMMUNITIES/MODULES
module1  = GTA1(1).graph_measures.moduleconsistency;
module1r = zeros(nr_nodes,nr_nodes);
for i = 1:nr_randomizations_pl
    module1r = module1r + GTA1r(i).graph_measures.moduleconsistency;
end
module1r = module1r./nr_randomizations_pl;
% in random networks, each pair of nodes has an equal probability of being
% in the same module
values = LCN_get_values_tril(module1r);
threshold = mean(values) + 2.*std(values);

module1bin = module1 > threshold;
module1f = (LCN_calc_community_structure(module1bin)>0.95);

assigned1 = zeros(nr_nodes,1);
module1a = zeros(nr_nodes,1);
nr_modules1 = 0;
while min(assigned1) == 0
    clear tmp tmp2
    nr_modules1 = nr_modules1 + 1;
    tmp = find(assigned1 == 0);
    tmp2 = find(module1f(tmp(1),:) == 1);
    module1a(tmp(1)) = nr_modules1;
    module1a(tmp2) = nr_modules1;
    assigned1(tmp(1)) = 1;
    assigned1(tmp2) = 1;
end
if isfield(GTA1(1),'density')
   fprintf('Community structure group 1, density = %4.2f %%\n',GTA1(1).density)
end
for i = 1:max(nr_modules1)
    clear tmp
    fprintf('module %i\n',i);
    tmp = find(module1a == i);
    for j = 1:length(tmp)
        fprintf('\t %s\n', char(nodenames(tmp(j))));
    end
end

module2  = GTA2(1).graph_measures.moduleconsistency;
module2r = zeros(nr_nodes,nr_nodes);
for i = 1:nr_randomizations_pl
    module2r = module2r + GTA2r(i).graph_measures.moduleconsistency;
end
module2r = module2r./nr_randomizations_pl;
% in random networks, each pair of nodes has an equal probability of being
% in the same module
values = LCN_get_values_tril(module2r);
threshold = mean(values) + 2.*std(values);

module2bin = module2 > threshold;
module2f = (LCN_calc_community_structure(module2bin)>0.95);

assigned2 = zeros(nr_nodes,1);
module2a = zeros(nr_nodes,1);
nr_modules2 = 0;
while min(assigned2) == 0
    clear tmp tmp2
    nr_modules2 = nr_modules2 + 1;
    tmp = find(assigned2 == 0);
    tmp2 = find(module2f(tmp(1),:) == 1);
    module2a(tmp(1)) = nr_modules2;
    module2a(tmp2) = nr_modules2;
    assigned2(tmp(1)) = 1;
    assigned2(tmp2) = 1;
end
if isfield(GTA2(1),'density')
   fprintf('Community structure group 2, density = %4.2f %%\n',GTA2(1).density)
end
for i = 1:max(nr_modules2)
    clear tmp
    fprintf('module %i\n',i);
    tmp = find(module2a == i);
    for j = 1:length(tmp)
        fprintf('\t %s\n', char(nodenames(tmp(j))));
    end
end

% % test if these proportions are statistically significantly different
% z_module1 = (module1 - module1r)./sqrt(module1.*(1-module1) + module1r.*(1-module1r)./nr_randomizations_pl);
% % get p-values
% p_module1 = LCN_Normal_z_to_p(z_module1,1);
% 
% % test if these proportions are statistically significantly different
% z_module2 = (module2 - module2r)./sqrt(module2.*(1-module2) + module2r.*(1-module2r)./nr_randomizations_pl);
% % get p-values
% p_module2 = LCN_Normal_z_to_p(z_module2,1);
% 
% % test difference between both groups
% % test if these proportions are statistically significantly different
% z_module12 = (module1 - module2)./sqrt(module1.*(1-module1) + module2.*(1-module2));
% % get p-values
% p_module12 = LCN_Normal_z_to_p(z_module12,1);
end