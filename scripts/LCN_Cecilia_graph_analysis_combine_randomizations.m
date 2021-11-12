clear

resultsfile1 = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\workspace5_correlations_90_regions_raw_new.mat';
resultsfile2 = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\workspace_correlations_90_regions_raw_new_extra800.mat';

load(resultsfile1)
GTA1 = GTA;
load(resultsfile2)
GTA2 = GTA;

GTA = GTA1;
for i = 1:800
    GTA.RESULTS.IBS_low_IBS_high.random_network1(i+200) = GTA2.RESULTS.IBS_low_IBS_high.random_network1(i);
    GTA.RESULTS.IBS_low_IBS_high.random_network2(i+200) = GTA2.RESULTS.IBS_low_IBS_high.random_network2(i);
    GTA.RESULTS.HC_IBS_low.random_network1(i+200) = GTA2.RESULTS.HC_IBS_low.random_network1(i);
    GTA.RESULTS.HC_IBS_low.random_network2(i+200) = GTA2.RESULTS.HC_IBS_low.random_network2(i);
    GTA.RESULTS.HC_IBS_high.random_network1(i+200) = GTA2.RESULTS.HC_IBS_high.random_network1(i);
    GTA.RESULTS.HC_IBS_high.random_network2(i+200) = GTA2.RESULTS.HC_IBS_high.random_network2(i);
end
save -v7.3 GTA_correlations_90_regions_raw_new GTA

resultsfile1 = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\workspace5_correlations_90_regions_corTGMV.mat';
resultsfile2 = 'C:\DATA\P9_NEUROGASTRO\Cecilia\PROJECT2_IBS\workspace_correlations_90_regions_corTGMV_new_extra800.mat';

load(resultsfile1)
GTA1 = GTA;
load(resultsfile2)
GTA2 = GTA;

GTA = GTA1;
for i = 1:800
    GTA.RESULTS.IBS_low_IBS_high.random_network1(i+200) = GTA2.RESULTS.IBS_low_IBS_high.random_network1(i);
    GTA.RESULTS.IBS_low_IBS_high.random_network2(i+200) = GTA2.RESULTS.IBS_low_IBS_high.random_network2(i);
    GTA.RESULTS.HC_IBS_low.random_network1(i+200) = GTA2.RESULTS.HC_IBS_low.random_network1(i);
    GTA.RESULTS.HC_IBS_low.random_network2(i+200) = GTA2.RESULTS.HC_IBS_low.random_network2(i);
    GTA.RESULTS.HC_IBS_high.random_network1(i+200) = GTA2.RESULTS.HC_IBS_high.random_network1(i);
    GTA.RESULTS.HC_IBS_high.random_network2(i+200) = GTA2.RESULTS.HC_IBS_high.random_network2(i);
end
save -v7.3 GTA_correlations_90_regions_corTGMV GTA