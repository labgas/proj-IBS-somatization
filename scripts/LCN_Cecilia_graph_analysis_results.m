clear

% GTAfile = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\GTA_correlations_90_regions_raw_new';
GTAfile = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\GTA_correlations_90_regions_corTGMV';
load(GTAfile)

% ANALYZE RESULTS
%++++++++++++++++
average_node_degree_HC        = GTA(1).RESULTS.HC.graph_measures.average_node_degree;
clustering_coefficient_HC     = GTA(1).RESULTS.HC.graph_measures.C;
characteristic_path_length_HC = GTA(1).RESULTS.HC.graph_measures.lambda;
global_efficiency_HC          = GTA(1).RESULTS.HC.graph_measures.E;
betweenness_centrality_HC     = GTA(1).RESULTS.HC.graph_measures.BC;

average_node_degree_IBS_low        = GTA(1).RESULTS.IBS_low.graph_measures.average_node_degree;
clustering_coefficient_IBS_low     = GTA(1).RESULTS.IBS_low.graph_measures.C;
characteristic_path_length_IBS_low = GTA(1).RESULTS.IBS_low.graph_measures.lambda;
global_efficiency_IBS_low          = GTA(1).RESULTS.IBS_low.graph_measures.E;
betweenness_centrality_IBS_low     = GTA(1).RESULTS.IBS_low.graph_measures.BC;

average_node_degree_IBS_high        = GTA(1).RESULTS.IBS_high.graph_measures.average_node_degree;
clustering_coefficient_IBS_high     = GTA(1).RESULTS.IBS_high.graph_measures.C;
characteristic_path_length_IBS_high = GTA(1).RESULTS.IBS_high.graph_measures.lambda;
global_efficiency_IBS_high          = GTA(1).RESULTS.IBS_high.graph_measures.E;
betweenness_centrality_IBS_high     = GTA(1).RESULTS.IBS_high.graph_measures.BC;

nr_randomizations_HC       = size(GTA(1).RESULTS.HC.graph_measures_random,2);
nr_randomizations_IBS_low  = size(GTA(1).RESULTS.IBS_low.graph_measures_random,2);
nr_randomizations_IBS_high = size(GTA(1).RESULTS.IBS_high.graph_measures_random,2);
         
% clustering coefficient
distribution_HC_C = zeros(1,nr_randomizations_HC);
for i = 1:nr_randomizations_HC
    distribution_HC_C(i) = GTA(1).RESULTS.HC.graph_measures_random(i).C; 
end
% characteristic path length
distribution_HC_lambda = zeros(1,nr_randomizations_HC);
for i = 1:nr_randomizations_HC
    distribution_HC_lambda(i) = GTA(1).RESULTS.HC.graph_measures_random(i).lambda; 
end
% global efficiency
distribution_HC_E = zeros(1,nr_randomizations_HC);
for i = 1:nr_randomizations_HC
    distribution_HC_E(i) = GTA(1).RESULTS.HC.graph_measures_random(i).E; 
end
% betweenness centrality
distribution_HC_BC = zeros(1,nr_randomizations_HC);
for i = 1:nr_randomizations_HC
    distribution_HC_BC(i) = GTA(1).RESULTS.HC.graph_measures_random(i).BC; 
end
norm_clustering_coefficient_HC     = GTA(1).RESULTS.HC.graph_measures.C./mean(distribution_HC_C);
norm_characteristic_path_length_HC = GTA(1).RESULTS.HC.graph_measures.lambda./mean(distribution_HC_lambda);
norm_global_efficiency_HC          = GTA(1).RESULTS.HC.graph_measures.E./mean(distribution_HC_E);
norm_betweenness_centrality_HC     = GTA(1).RESULTS.HC.graph_measures.BC./mean(distribution_HC_BC);
% clustering coefficient
distribution_IBS_low_C = zeros(1,nr_randomizations_IBS_low);
for i = 1:nr_randomizations_IBS_low
    distribution_IBS_low_C(i) = GTA(1).RESULTS.IBS_low.graph_measures_random(i).C; 
end
% characteristic path length
distribution_IBS_low_lambda = zeros(1,nr_randomizations_IBS_low);
for i = 1:nr_randomizations_IBS_low
    distribution_IBS_low_lambda(i) = GTA(1).RESULTS.IBS_low.graph_measures_random(i).lambda; 
end
% global efficiency
distribution_IBS_low_E = zeros(1,nr_randomizations_IBS_low);
for i = 1:nr_randomizations_IBS_low
    distribution_IBS_low_E(i) = GTA(1).RESULTS.IBS_low.graph_measures_random(i).E; 
end
% betweenness centrality
distribution_IBS_low_BC = zeros(1,nr_randomizations_IBS_low);
for i = 1:nr_randomizations_IBS_low
    distribution_IBS_low_BC(i) = GTA(1).RESULTS.IBS_low.graph_measures_random(i).BC; 
end
norm_clustering_coefficient_IBS_low     = GTA(1).RESULTS.IBS_low.graph_measures.C./mean(distribution_IBS_low_C);
norm_characteristic_path_length_IBS_low = GTA(1).RESULTS.IBS_low.graph_measures.lambda./mean(distribution_IBS_low_lambda);
norm_global_efficiency_IBS_low          = GTA(1).RESULTS.IBS_low.graph_measures.E./mean(distribution_IBS_low_E);
norm_betweenness_centrality_IBS_low     = GTA(1).RESULTS.IBS_low.graph_measures.BC./mean(distribution_IBS_low_BC);
% clustering coefficient
distribution_IBS_high_C = zeros(1,nr_randomizations_IBS_high);
for i = 1:nr_randomizations_IBS_high
    distribution_IBS_high_C(i) = GTA(1).RESULTS.IBS_high.graph_measures_random(i).C; 
end
% characteristic path length
distribution_IBS_high_lambda = zeros(1,nr_randomizations_IBS_high);
for i = 1:nr_randomizations_IBS_high
    distribution_IBS_high_lambda(i) = GTA(1).RESULTS.IBS_high.graph_measures_random(i).lambda; 
end
% global efficiency
distribution_IBS_high_E = zeros(1,nr_randomizations_IBS_high);
for i = 1:nr_randomizations_IBS_high
    distribution_IBS_high_E(i) = GTA(1).RESULTS.IBS_high.graph_measures_random(i).E; 
end
% betweenness centrality
distribution_IBS_high_BC = zeros(1,nr_randomizations_IBS_high);
for i = 1:nr_randomizations_IBS_high
    distribution_IBS_high_BC(i) = GTA(1).RESULTS.IBS_high.graph_measures_random(i).BC; 
end
norm_clustering_coefficient_IBS_high     = GTA(1).RESULTS.IBS_high.graph_measures.C./mean(distribution_IBS_high_C);
norm_characteristic_path_length_IBS_high = GTA(1).RESULTS.IBS_high.graph_measures.lambda./mean(distribution_IBS_high_lambda);
norm_global_efficiency_IBS_high          = GTA(1).RESULTS.IBS_high.graph_measures.E./mean(distribution_IBS_high_E);
norm_betweenness_centrality_IBS_high     = GTA(1).RESULTS.IBS_high.graph_measures.BC./mean(distribution_IBS_high_BC);

% compare IBS_low vs IBS_high
% ----------------------------
clear GTA1 GTA2 GTA1r GTA2r
GTA1 = GTA.RESULTS.IBS_low;
GTA2 = GTA.RESULTS.IBS_high;
GTA1r = GTA.RESULTS.IBS_low_IBS_high.random_network1;
GTA2r = GTA.RESULTS.IBS_low_IBS_high.random_network2;
GTA1.nodenames = GTA.nodenames;
LCN_compare_graph_measures_2g(GTA1,GTA2,GTA1r,GTA2r,pvalue_threshold)

% compare HC vs IBS_low
%----------------------------
clear GTA1 GTA2 GTA1r GTA2r
GTA1 = GTA.RESULTS.HC;
GTA2 = GTA.RESULTS.IBS_low;
GTA1r = GTA.RESULTS.HC_IBS_low.random_network1;
GTA2r = GTA.RESULTS.HC_IBS_low.random_network2;
GTA1.nodenames = GTA.nodenames;
LCN_compare_graph_measures_2g(GTA1,GTA2,GTA1r,GTA2r,pvalue_threshold)

% compare HC vs IBS_high
%----------------------------
clear GTA1 GTA2 GTA1r GTA2r
GTA1 = GTA.RESULTS.HC;
GTA2 = GTA.RESULTS.IBS_high;
GTA1r = GTA.RESULTS.HC_IBS_high.random_network1;
GTA2r = GTA.RESULTS.HC_IBS_high.random_network2;
GTA1.nodenames = GTA.nodenames;
LCN_compare_graph_measures_2g(GTA1,GTA2,GTA1r,GTA2r,pvalue_threshold)

