%% Plot Comparison
load('IterativeExperimentResults_DUSP1_real_FIMopt.mat');
FIMopt = FIMOptNextExptSaved;
covLogMHFIMopt = covLogMH;
covMHFIMopt = covMH;
exptDesignsFIMopt = exptDesigns;
parmFIMopt = parametersFound;


load('IterativeExperimentResults_DUSP1_real_Intuition.mat');
FIMIntuition = FIMOptNextExptSaved;
covLogMHIntuition = covLogMH;
covMHIntuition = covMH;
exptDesignsIntuition = exptDesigns;
parmIntuition = parametersFound;

load('IterativeExperimentResults_DUSP1_real_random.mat');
FIMRandom = FIMOptNextExptSaved;
covLogMHRandom = covLogMH;
covMHRandom = covMH;
exptDesignsRandom = exptDesigns;
parmRandom = parametersFound;

load('IterativeExperimentResults_DUSP1_real_uniform.mat');
FIMuniform = FIMOptNextExptSaved;
covLogMHuniform = covLogMH;
covMHuniform = covMH;
exptDesignsuniform = exptDesigns;
parmuniform = parametersFound;

%% FIM
fimAvgOpt = zeros(8,8,5);
for i = 1:5
    FIM1 = FIMopt{i};
    A = zeros(size(FIM1{1}));
    for j = 1:10
        A = A +FIM1{j};
    end
    fimAvgOpt(:,:,i) = A/10;
    detFIMopt(i) = det(inv(fimAvgOpt(1:4,1:4,i)));

end

fimAvgIntuition = zeros(8,8,5);
for i = 1:5
    FIM1 = FIMIntuition{i};
    A = zeros(size(FIM1{1}));
    for j = 1:10
        A = A +FIM1{j};
    end
    fimAvgIntuition(:,:,i) = A/10;
    detFIMIntuition(i) = det(inv(fimAvgIntuition(1:4,1:4,i)));
end

fimAvgRandom = zeros(8,8,5);
for i = 1:5
    FIM1 = FIMRandom{i};
    A = zeros(size(FIM1{1}));
    for j = 1:10
        A = A +FIM1{j};
    end
    fimAvgRandom(:,:,i) = A/10;
    detFIMRandom(i) = det(inv(fimAvgRandom(1:4,1:4,i)));
end

fimAvgUniform = zeros(8,8,5);
for i = 1:5
    FIM1 = FIMuniform{i};
    A = zeros(size(FIM1{1}));
    for j = 1:10
        A = A +FIM1{j};
    end
    fimAvgUniform(:,:,i) = A/10;
    detFIMUniform(i) = det(inv(fimAvgUniform(1:4,1:4,i)));

end
%% MH COV
for i = 1:5
    detMHopt(i) =  det(covLogMHFIMopt{i});
end

for i = 1:5
    detMHIntuition(i) =  det(covLogMHIntuition{i});
end

for i = 1:5
    detMHRandom(i) =  det(covLogMHRandom{i});
end

for i = 1:5
    detMHuniform(i) =  det(covLogMHuniform{i});
end
%% Comparing FIM vs Experiment Design
experiment = 1:5;
figure
plot(experiment,detFIMopt,'b',...
    experiment,detFIMIntuition,'r:',...
    experiment,detFIMRandom,'k--', ...
    experiment,detFIMUniform,'g.-',LineWidth=2);
set(gca, 'YScale', 'log')
xlabel('Experiment Round');
ylabel('|FIM^-^1|');
title('|FIM^-^1| vs Experiment Design')
legend('Optimized FIM','Intuition','Random','Uniform')


%Comparing MH cov vs Experiment Design
figure
plot(experiment,detMHopt,'b',experiment,detMHIntuition,'r:',...
    experiment,detMHRandom,'k--',experiment,detMHuniform,'g.-',LineWidth=2);
xlabel('Experiment Round');
ylabel('Log_1_0(MH Cov)');
title('|COV_M_H| vs Experiment Design')
legend('OPtimized FIM','Intuition','Random','Uniform')
