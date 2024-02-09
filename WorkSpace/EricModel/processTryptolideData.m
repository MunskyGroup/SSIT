
X = readtable('DUSP1_GR_dataframes/DUSP1_100nM_Dex_5uM_TPL_R1.csv');

J = (X.conc==100|X.time_index==0);

trptTimes = unique(X.t_TPL(~isnan(X.t_TPL)))

for k = 1:length(trptTimes)
    JJ = J&(X.time_index==0|...
        (X.time_index<=trptTimes(k))&isnan(X.t_TPL))|...
        ((X.time_index>=trptTimes(k))&(X.t_TPL==trptTimes(k)));
    X.(['tryptCond',num2str(k)]) = zeros(size(X,1),1);
    X.(['tryptCond',num2str(k)])(JJ) = k;
end

writetable(X,'DUSP1_GR_dataframes/DUSP1_100nM_Dex_5uM_TPL_R1_Brian.csv');



