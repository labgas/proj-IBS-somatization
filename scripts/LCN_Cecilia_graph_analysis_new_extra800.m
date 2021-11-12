clear

% load the connectivity matrices and results obtained using the script 
% LCN_Cecilia_structural_connectivity_new.m
% resultsfile = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\results_correlations_90_regions_raw_new';
resultsfile = 'J:\GBW-0264_TARGID-Brain-Gut-Axis\LUKAS\GRAPH_CECILIA\results_correlations_90_regions_raw_new';
% resultsfile2 = 'J:\GBW-0264_TARGID-Brain-Gut-Axis\LUKAS\GRAPH_CECILIA\results_correlations_90_regions_corTGMV';
% 
nr_permutations    = 800;
nr_randomizations1 = 10;  % number of randomizations used to normalize the 
%                           graph measures by dividing by the average value 
%                           obtained in this number of equivalent random 
%                           networks 
% 
% % do not change below this line
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
[datadir,filename,~] = fileparts(resultsfile);
cd(datadir)
load(resultsfile)
% 
% calculate weighted graph measures
%---------------------------------------------------------
% network1: positive correlations
% calculate weights
GTA(1).HC.Zpos = zeros(size(RESULTS.HC.Zcorr));
GTA(1).HC.Zpos(RESULTS.HC.Zcorr>0) = RESULTS.HC.Zcorr(RESULTS.HC.Zcorr>0);
GTA(1).HC.wpos = LCN_calc_weights_network(GTA(1).HC.Zpos,1); 

GTA(1).IBS_low.Zpos = zeros(size(RESULTS.IBS_low.Zcorr));
GTA(1).IBS_low.Zpos(RESULTS.IBS_low.Zcorr>0) = RESULTS.IBS_low.Zcorr(RESULTS.IBS_low.Zcorr>0);
GTA(1).IBS_low.wpos = LCN_calc_weights_network(GTA(1).IBS_low.Zpos,1); 

GTA(1).IBS_high.Zpos = zeros(size(RESULTS.IBS_high.Zcorr));
GTA(1).IBS_high.Zpos(RESULTS.IBS_high.Zcorr>0) = RESULTS.IBS_high.Zcorr(RESULTS.IBS_high.Zcorr>0);
GTA(1).IBS_high.wpos = LCN_calc_weights_network(GTA(1).IBS_high.Zpos,1); 

GTA(1).nodenames = RESULTS.name_nodes;
GTA(1).strategy  = {'g'; 'w'; 'pos'};

GTA(1).RESULTS.HC.network = GTA(1).HC.wpos;
GTA(1).RESULTS.IBS_low.network = GTA(1).IBS_low.wpos;
GTA(1).RESULTS.IBS_high.network = GTA(1).IBS_high.wpos;
                      
% determine the graph measures of group HC
clear M Mrandom
fprintf('\t calculating graph measures HC \n');
tic
M = LCN_calc_graph_measures(GTA(1).HC.wpos,'w');
toc
fprintf('\t calculating random graph measures HC \n');
tic
Mrandom = LCN_calc_graph_measures_random(GTA(1).HC.wpos,nr_randomizations1,'w');
toc
GTA(1).RESULTS.HC.graph_measures = M;
GTA(1).RESULTS.HC.graph_measures_random = Mrandom;

% determine the graph measures of group IBS_low
clear M Mrandom
fprintf('\t calculating graph measures IBS_low \n');    
tic
M = LCN_calc_graph_measures(GTA(1).IBS_low.wpos,'w');
toc
fprintf('\t calculating random graph measures IBS_low \n');
tic
Mrandom = LCN_calc_graph_measures_random(GTA(1).IBS_low.wpos,nr_randomizations1,'w');
toc
GTA(1).RESULTS.IBS_low.graph_measures = M;
GTA(1).RESULTS.IBS_low.graph_measures_random = Mrandom;

% determine the graph measures of group IBS_high
clear M Mrandom
fprintf('\t calculating graph measures IBS_high \n');
tic
M = LCN_calc_graph_measures(GTA(1).IBS_high.wpos,'w');
toc
fprintf('\t calculating random graph measures IBS_high \n');
tic
Mrandom = LCN_calc_graph_measures_random(GTA(1).IBS_high.wpos,nr_randomizations1,'w');
toc
GTA(1).RESULTS.IBS_high.graph_measures = M;
GTA(1).RESULTS.IBS_high.graph_measures_random = Mrandom;
save workspace_correlations_90_regions_raw_new_extra800 -v7.3
% 
for rand_i = 1:nr_permutations
    % compare IBS_low and IBS_high
    %-----------------------------    
    Zcorr1_tmp = RESULTS.IBS_low_IBS_high.Zcorr1_rand(:,:,rand_i);
    Zpos1_rand = zeros(size(Zcorr1_tmp));
    Zpos1_rand(Zcorr1_tmp>0) = Zcorr1_tmp(Zcorr1_tmp>0);
    wpos1_rand = LCN_calc_weights_network(Zpos1_rand,1); 
    Zcorr2_tmp = RESULTS.IBS_low_IBS_high.Zcorr2_rand(:,:,rand_i);
    Zpos2_rand = zeros(size(Zcorr2_tmp));
    Zpos2_rand(Zcorr2_tmp>0) = Zcorr2_tmp(Zcorr2_tmp>0);
    wpos2_rand = LCN_calc_weights_network(Zpos2_rand,1); 
           
    clear M Mrandom
    % clear diagonal
    wpos1_rand(1:size(wpos1_rand,1)+1:end) = 0;  %clear diagonal
    wpos2_rand(1:size(wpos2_rand,1)+1:end) = 0;  %clear diagonal
    GTA(1).RESULTS.IBS_low_IBS_high.random_network1(rand_i).network = wpos1_rand;
    GTA(1).RESULTS.IBS_low_IBS_high.random_network2(rand_i).network = wpos2_rand;
                     
    % determine the graph measures of random group1
    clear M Mrandom
    fprintf('\t calculating graph measures of random group1 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos1_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos1_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.IBS_low_IBS_high.random_network1(rand_i).graph_measures = M;
    GTA(1).RESULTS.IBS_low_IBS_high.random_network1(rand_i).graph_measures_random = Mrandom;

    % determine the graph measures of random group2
    clear M Mrandom
    fprintf('\t calculating graph measures of random group2 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos2_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos2_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.IBS_low_IBS_high.random_network2(rand_i).graph_measures = M;
    GTA(1).RESULTS.IBS_low_IBS_high.random_network2(rand_i).graph_measures_random = Mrandom;
end
save workspace_correlations_90_regions_raw_new_extra800 -v7.3
% 
for rand_i = 1:nr_permutations
% for rand_i = 1:200
    % compare HC and IBS_low
    %-----------------------------    
    Zcorr1_tmp = RESULTS.HC_IBS_low.Zcorr1_rand(:,:,rand_i);
    Zpos1_rand = zeros(size(Zcorr1_tmp));
    Zpos1_rand(Zcorr1_tmp>0) = Zcorr1_tmp(Zcorr1_tmp>0);
    wpos1_rand = LCN_calc_weights_network(Zpos1_rand,1); 
    Zcorr2_tmp = RESULTS.HC_IBS_low.Zcorr2_rand(:,:,rand_i);
    Zpos2_rand = zeros(size(Zcorr2_tmp));
    Zpos2_rand(Zcorr2_tmp>0) = Zcorr2_tmp(Zcorr2_tmp>0);
    wpos2_rand = LCN_calc_weights_network(Zpos2_rand,1); 
           
    clear M Mrandom        
    % clear diagonal
    wpos1_rand(1:size(wpos1_rand,1)+1:end) = 0;  %clear diagonal
    wpos2_rand(1:size(wpos2_rand,1)+1:end) = 0;  %clear diagonal
    GTA(1).RESULTS.HC_IBS_low.random_network1(rand_i).network = wpos1_rand;
    GTA(1).RESULTS.HC_IBS_low.random_network2(rand_i).network = wpos2_rand;
                     
    % determine the graph measures of random group1
    clear M Mrandom
    fprintf('\t calculating graph measures of random group1 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos1_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos1_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_low.random_network1(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_low.random_network1(rand_i).graph_measures_random = Mrandom;

    % determine the graph measures of random group2
    clear M Mrandom
    fprintf('\t calculating graph measures of random group2 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos2_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos2_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_low.random_network2(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_low.random_network2(rand_i).graph_measures_random = Mrandom;
end
save workspace_correlations_90_regions_raw_new_extra800  -v7.3
%
for rand_i = 1:nr_permutations
% for rand_i = 1:200
    % compare HC and IBS_high
    %-----------------------------    
    Zcorr1_tmp = RESULTS.HC_IBS_high.Zcorr1_rand(:,:,rand_i);
    Zpos1_rand = zeros(size(Zcorr1_tmp));
    Zpos1_rand(Zcorr1_tmp>0) = Zcorr1_tmp(Zcorr1_tmp>0);
    wpos1_rand = LCN_calc_weights_network(Zpos1_rand,1); 
    Zcorr2_tmp = RESULTS.HC_IBS_high.Zcorr2_rand(:,:,rand_i);
    Zpos2_rand = zeros(size(Zcorr2_tmp));
    Zpos2_rand(Zcorr2_tmp>0) = Zcorr2_tmp(Zcorr2_tmp>0);
    wpos2_rand = LCN_calc_weights_network(Zpos2_rand,1); 
           
    clear M Mrandom
            
    % clear diagonal
    wpos1_rand(1:size(wpos1_rand,1)+1:end) = 0;  %clear diagonal
    wpos2_rand(1:size(wpos2_rand,1)+1:end) = 0;  %clear diagonal
    GTA(1).RESULTS.HC_IBS_high.random_network1(rand_i).network = wpos1_rand;
    GTA(1).RESULTS.HC_IBS_high.random_network2(rand_i).network = wpos2_rand;
                     
    % determine the graph measures of random group1
    clear M Mrandom
    fprintf('\t calculating graph measures of random group1 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos1_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos1_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_high.random_network1(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_high.random_network1(rand_i).graph_measures_random = Mrandom;

    % determine the graph measures of random group2
    clear M Mrandom
    fprintf('\t calculating graph measures of random group2 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos2_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos2_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_high.random_network2(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_high.random_network2(rand_i).graph_measures_random = Mrandom;
end
save workspace_correlations_90_regions_raw_new_extra800  -v7.3

% PART 2 corrected values

clear

resultsfile = 'J:\GBW-0264_TARGID-Brain-Gut-Axis\LUKAS\GRAPH_CECILIA\results_correlations_90_regions_corTGMV';
% 
nr_permutations    = 800;
nr_randomizations1 = 10;  % number of randomizations used to normalize the 
%                           graph measures by dividing by the average value 
%                           obtained in this number of equivalent random 
%                           networks 
% 
% % do not change below this line
% %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 
[datadir,filename,~] = fileparts(resultsfile);
cd(datadir)
load(resultsfile)
% 
% calculate weighted graph measures
%---------------------------------------------------------
% network1: positive correlations
% calculate weights
GTA(1).HC.Zpos = zeros(size(RESULTS.HC.Zcorr));
GTA(1).HC.Zpos(RESULTS.HC.Zcorr>0) = RESULTS.HC.Zcorr(RESULTS.HC.Zcorr>0);
GTA(1).HC.wpos = LCN_calc_weights_network(GTA(1).HC.Zpos,1); 

GTA(1).IBS_low.Zpos = zeros(size(RESULTS.IBS_low.Zcorr));
GTA(1).IBS_low.Zpos(RESULTS.IBS_low.Zcorr>0) = RESULTS.IBS_low.Zcorr(RESULTS.IBS_low.Zcorr>0);
GTA(1).IBS_low.wpos = LCN_calc_weights_network(GTA(1).IBS_low.Zpos,1); 

GTA(1).IBS_high.Zpos = zeros(size(RESULTS.IBS_high.Zcorr));
GTA(1).IBS_high.Zpos(RESULTS.IBS_high.Zcorr>0) = RESULTS.IBS_high.Zcorr(RESULTS.IBS_high.Zcorr>0);
GTA(1).IBS_high.wpos = LCN_calc_weights_network(GTA(1).IBS_high.Zpos,1); 

GTA(1).nodenames = RESULTS.name_nodes;
GTA(1).strategy  = {'g'; 'w'; 'pos'};

GTA(1).RESULTS.HC.network = GTA(1).HC.wpos;
GTA(1).RESULTS.IBS_low.network = GTA(1).IBS_low.wpos;
GTA(1).RESULTS.IBS_high.network = GTA(1).IBS_high.wpos;
                      
% determine the graph measures of group HC
clear M Mrandom
fprintf('\t calculating graph measures HC \n');
tic
M = LCN_calc_graph_measures(GTA(1).HC.wpos,'w');
toc
fprintf('\t calculating random graph measures HC \n');
tic
Mrandom = LCN_calc_graph_measures_random(GTA(1).HC.wpos,nr_randomizations1,'w');
toc
GTA(1).RESULTS.HC.graph_measures = M;
GTA(1).RESULTS.HC.graph_measures_random = Mrandom;

% determine the graph measures of group IBS_low
clear M Mrandom
fprintf('\t calculating graph measures IBS_low \n');    
tic
M = LCN_calc_graph_measures(GTA(1).IBS_low.wpos,'w');
toc
fprintf('\t calculating random graph measures IBS_low \n');
tic
Mrandom = LCN_calc_graph_measures_random(GTA(1).IBS_low.wpos,nr_randomizations1,'w');
toc
GTA(1).RESULTS.IBS_low.graph_measures = M;
GTA(1).RESULTS.IBS_low.graph_measures_random = Mrandom;

% determine the graph measures of group IBS_high
clear M Mrandom
fprintf('\t calculating graph measures IBS_high \n');
tic
M = LCN_calc_graph_measures(GTA(1).IBS_high.wpos,'w');
toc
fprintf('\t calculating random graph measures IBS_high \n');
tic
Mrandom = LCN_calc_graph_measures_random(GTA(1).IBS_high.wpos,nr_randomizations1,'w');
toc
GTA(1).RESULTS.IBS_high.graph_measures = M;
GTA(1).RESULTS.IBS_high.graph_measures_random = Mrandom;
save workspace_correlations_90_regions_corTGMV_new_extra800 -v7.3

for rand_i = 1:nr_permutations
    % compare IBS_low and IBS_high
    %-----------------------------    
    Zcorr1_tmp = RESULTS.IBS_low_IBS_high.Zcorr1_rand(:,:,rand_i);
    Zpos1_rand = zeros(size(Zcorr1_tmp));
    Zpos1_rand(Zcorr1_tmp>0) = Zcorr1_tmp(Zcorr1_tmp>0);
    wpos1_rand = LCN_calc_weights_network(Zpos1_rand,1); 
    Zcorr2_tmp = RESULTS.IBS_low_IBS_high.Zcorr2_rand(:,:,rand_i);
    Zpos2_rand = zeros(size(Zcorr2_tmp));
    Zpos2_rand(Zcorr2_tmp>0) = Zcorr2_tmp(Zcorr2_tmp>0);
    wpos2_rand = LCN_calc_weights_network(Zpos2_rand,1); 
           
    clear M Mrandom
    % clear diagonal
    wpos1_rand(1:size(wpos1_rand,1)+1:end) = 0;  %clear diagonal
    wpos2_rand(1:size(wpos2_rand,1)+1:end) = 0;  %clear diagonal
    GTA(1).RESULTS.IBS_low_IBS_high.random_network1(rand_i).network = wpos1_rand;
    GTA(1).RESULTS.IBS_low_IBS_high.random_network2(rand_i).network = wpos2_rand;
                     
    % determine the graph measures of random group1
    clear M Mrandom
    fprintf('\t calculating graph measures of random group1 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos1_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos1_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.IBS_low_IBS_high.random_network1(rand_i).graph_measures = M;
    GTA(1).RESULTS.IBS_low_IBS_high.random_network1(rand_i).graph_measures_random = Mrandom;

    % determine the graph measures of random group2
    clear M Mrandom
    fprintf('\t calculating graph measures of random group2 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos2_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos2_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.IBS_low_IBS_high.random_network2(rand_i).graph_measures = M;
    GTA(1).RESULTS.IBS_low_IBS_high.random_network2(rand_i).graph_measures_random = Mrandom;
end
save workspace_correlations_90_regions_corTGMV_new_extra800 -v7.3
% 
for rand_i = 1:nr_permutations
% for rand_i = 1:200
    % compare HC and IBS_low
    %-----------------------------    
    Zcorr1_tmp = RESULTS.HC_IBS_low.Zcorr1_rand(:,:,rand_i);
    Zpos1_rand = zeros(size(Zcorr1_tmp));
    Zpos1_rand(Zcorr1_tmp>0) = Zcorr1_tmp(Zcorr1_tmp>0);
    wpos1_rand = LCN_calc_weights_network(Zpos1_rand,1); 
    Zcorr2_tmp = RESULTS.HC_IBS_low.Zcorr2_rand(:,:,rand_i);
    Zpos2_rand = zeros(size(Zcorr2_tmp));
    Zpos2_rand(Zcorr2_tmp>0) = Zcorr2_tmp(Zcorr2_tmp>0);
    wpos2_rand = LCN_calc_weights_network(Zpos2_rand,1); 
           
    clear M Mrandom        
    % clear diagonal
    wpos1_rand(1:size(wpos1_rand,1)+1:end) = 0;  %clear diagonal
    wpos2_rand(1:size(wpos2_rand,1)+1:end) = 0;  %clear diagonal
    GTA(1).RESULTS.HC_IBS_low.random_network1(rand_i).network = wpos1_rand;
    GTA(1).RESULTS.HC_IBS_low.random_network2(rand_i).network = wpos2_rand;
                     
    % determine the graph measures of random group1
    clear M Mrandom
    fprintf('\t calculating graph measures of random group1 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos1_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos1_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_low.random_network1(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_low.random_network1(rand_i).graph_measures_random = Mrandom;

    % determine the graph measures of random group2
    clear M Mrandom
    fprintf('\t calculating graph measures of random group2 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos2_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos2_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_low.random_network2(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_low.random_network2(rand_i).graph_measures_random = Mrandom;
end
save workspace_correlations_90_regions_corTGMV_new_extra800  -v7.3
%  
for rand_i = 1:nr_permutations
% for rand_i = 1:200
    % compare HC and IBS_high
    %-----------------------------    
    Zcorr1_tmp = RESULTS.HC_IBS_high.Zcorr1_rand(:,:,rand_i);
    Zpos1_rand = zeros(size(Zcorr1_tmp));
    Zpos1_rand(Zcorr1_tmp>0) = Zcorr1_tmp(Zcorr1_tmp>0);
    wpos1_rand = LCN_calc_weights_network(Zpos1_rand,1); 
    Zcorr2_tmp = RESULTS.HC_IBS_high.Zcorr2_rand(:,:,rand_i);
    Zpos2_rand = zeros(size(Zcorr2_tmp));
    Zpos2_rand(Zcorr2_tmp>0) = Zcorr2_tmp(Zcorr2_tmp>0);
    wpos2_rand = LCN_calc_weights_network(Zpos2_rand,1); 
           
    clear M Mrandom
            
    % clear diagonal
    wpos1_rand(1:size(wpos1_rand,1)+1:end) = 0;  %clear diagonal
    wpos2_rand(1:size(wpos2_rand,1)+1:end) = 0;  %clear diagonal
    GTA(1).RESULTS.HC_IBS_high.random_network1(rand_i).network = wpos1_rand;
    GTA(1).RESULTS.HC_IBS_high.random_network2(rand_i).network = wpos2_rand;
                     
    % determine the graph measures of random group1
    clear M Mrandom
    fprintf('\t calculating graph measures of random group1 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos1_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos1_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_high.random_network1(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_high.random_network1(rand_i).graph_measures_random = Mrandom;

    % determine the graph measures of random group2
    clear M Mrandom
    fprintf('\t calculating graph measures of random group2 - randomization %i of %i\n',rand_i,nr_permutations)
    M = LCN_calc_graph_measures(wpos2_rand,'w');
    Mrandom = LCN_calc_graph_measures_random(wpos2_rand,nr_randomizations1,'w');
    GTA(1).RESULTS.HC_IBS_high.random_network2(rand_i).graph_measures = M;
    GTA(1).RESULTS.HC_IBS_high.random_network2(rand_i).graph_measures_random = Mrandom;
end
save workspace_correlations_90_regions_corTGMV_new_extra800  -v7.3

