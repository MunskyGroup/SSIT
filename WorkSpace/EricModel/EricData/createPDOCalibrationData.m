TMP = readtable('Data010224_Gated_On_Nuc.csv');
% TMP = readtable('EricDataJan23_2024/DUSP1_predict_flitered_data_75min_ConcSweep_030624.csv');
TMP = TMP(~isnan(TMP.RNA_DUSP1_nuc),:);
TMP = TMP(strcmp(TMP.Condition,'DUSP1_timesweep'),:);
p = 0;
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int.*(TMP.Nuc_Area.^p);
maxVal = max(TMP.Nuc_DUSP1_avg_int_tot);
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int_tot/maxVal;
TMP.Nuc_DUSP1_avg_int_tot = round(400*TMP.Nuc_DUSP1_avg_int_tot);
writetable(TMP,'pdoCalibrationData_EricIntensity.csv')

%%
TMP = readtable('Data010224_Gated_On_Nuc.csv');
TMP = TMP(~isnan(TMP.RNA_DUSP1_nuc),:);
TMP = TMP(strcmp(TMP.Condition,'DUSP1_timesweep'),:);

lastval=NaN;
for i = 1:size(TMP,1)
    if TMP.Dex_Conc(i)==0
        if isnan(lastval)
            lastval = i;
        end
    elseif ~isnan(lastval)
        TMP.Dex_Conc(lastval:i) = TMP.Dex_Conc(i);
        lastval = NaN;
    end
end

TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int.*(TMP.Nuc_Area.^p);
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int_tot/maxVal;
TMP.Nuc_DUSP1_avg_int_tot = round(400*TMP.Nuc_DUSP1_avg_int_tot);

totNucRNA = TMP.RNA_DUSP1_nuc;
totNucRNA(~isnan(TMP.DUSP1_ts_size_0)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_0)) + TMP.DUSP1_ts_size_0(~isnan(TMP.DUSP1_ts_size_0));
totNucRNA(~isnan(TMP.DUSP1_ts_size_1)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_1)) + TMP.DUSP1_ts_size_1(~isnan(TMP.DUSP1_ts_size_1));
totNucRNA(~isnan(TMP.DUSP1_ts_size_2)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_2)) + TMP.DUSP1_ts_size_2(~isnan(TMP.DUSP1_ts_size_2));
totNucRNA(~isnan(TMP.DUSP1_ts_size_3)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_3)) + TMP.DUSP1_ts_size_3(~isnan(TMP.DUSP1_ts_size_3));
TMP.totalNucRNA = totNucRNA;
writetable(TMP,'pdoCalibrationData_EricIntensity_DexSweeps.csv')

%%
TMP = readtable('Data010224_Gated_On_Nuc.csv');
TMP = TMP(~isnan(TMP.RNA_DUSP1_nuc),:);
TMP = TMP((TMP.Time_index == 0|TMP.Time_index == 75),:)
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int.*(TMP.Nuc_Area.^p);
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int_tot/maxVal;
TMP.Nuc_DUSP1_avg_int_tot = round(400*TMP.Nuc_DUSP1_avg_int_tot);

totNucRNA = TMP.RNA_DUSP1_nuc;
totNucRNA(~isnan(TMP.DUSP1_ts_size_0)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_0)) + TMP.DUSP1_ts_size_0(~isnan(TMP.DUSP1_ts_size_0));
totNucRNA(~isnan(TMP.DUSP1_ts_size_1)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_1)) + TMP.DUSP1_ts_size_1(~isnan(TMP.DUSP1_ts_size_1));
totNucRNA(~isnan(TMP.DUSP1_ts_size_2)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_2)) + TMP.DUSP1_ts_size_2(~isnan(TMP.DUSP1_ts_size_2));
totNucRNA(~isnan(TMP.DUSP1_ts_size_3)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_3)) + TMP.DUSP1_ts_size_3(~isnan(TMP.DUSP1_ts_size_3));
TMP.totalNucRNA = totNucRNA;

writetable(TMP,'pdoCalibrationData_EricIntensity_ConcSweeps.csv')

%%
TMP = readtable('Data010224_Gated_On_Nuc.csv');
TMP = TMP(TMP.Dex_Conc==10000,:);
TMP = TMP((TMP.Time_index == 0|TMP.Time_index == 75),:);
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int.*(TMP.Nuc_Area.^p);
TMP.Nuc_DUSP1_avg_int_tot = TMP.Nuc_DUSP1_avg_int_tot/maxVal;
TMP.Nuc_DUSP1_avg_int_tot = round(400*TMP.Nuc_DUSP1_avg_int_tot);

totNucRNA = TMP.RNA_DUSP1_nuc;
totNucRNA(~isnan(TMP.DUSP1_ts_size_0)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_0)) + TMP.DUSP1_ts_size_0(~isnan(TMP.DUSP1_ts_size_0));
totNucRNA(~isnan(TMP.DUSP1_ts_size_1)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_1)) + TMP.DUSP1_ts_size_1(~isnan(TMP.DUSP1_ts_size_1));
totNucRNA(~isnan(TMP.DUSP1_ts_size_2)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_2)) + TMP.DUSP1_ts_size_2(~isnan(TMP.DUSP1_ts_size_2));
totNucRNA(~isnan(TMP.DUSP1_ts_size_3)) = totNucRNA(~isnan(TMP.DUSP1_ts_size_3)) + TMP.DUSP1_ts_size_3(~isnan(TMP.DUSP1_ts_size_3));
TMP.totalNucRNA = totNucRNA;

writetable(TMP,'pdoCalibrationData_EricIntensity_ConcHigh.csv')
