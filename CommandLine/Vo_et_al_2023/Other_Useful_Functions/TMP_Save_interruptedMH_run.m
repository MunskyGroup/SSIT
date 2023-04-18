clear all
load FixedOnRate_PDOon0_300_Prior1x_0_18.mat

load TMPmh_1_FixedOnRate_PDOon0_300_Prior1x_0_18.mat
J = value~=0;
chainResults{1}.mhSamples = squeeze(smpl(J,1,:));
chainResults{1}.mhValue = value(J);

load TMPmh_2_FixedOnRate_PDOon0_300_Prior1x_0_18.mat
J = value~=0;
chainResults{2}.mhSamples = squeeze(smpl(J,1,:));
chainResults{2}.mhValue = value(J);

load TMPmh_3_FixedOnRate_PDOon0_300_Prior1x_0_18.mat
J = value~=0;
chainResults{3}.mhSamples = squeeze(smpl(J,1,:));
chainResults{3}.mhValue = value(J);

load TMPmh_4_FixedOnRate_PDOon0_300_Prior1x_0_18.mat
J = value~=0;
chainResults{4}.mhSamples = squeeze(smpl(J,1,:));
chainResults{4}.mhValue = value(J);

load TMPmh_5_FixedOnRate_PDOon0_300_Prior1x_0_18.mat
J = value~=0;
chainResults{5}.mhSamples = squeeze(smpl(J,1,:));
chainResults{5}.mhValue = value(J);

save FixedOnRate_PDOon0_300_Prior1x_0_18 FIMZero FIMZeroLog ModelZero chainResults covLog fimResults lTrue2FISH lTrue2MCP fimResultsOne
