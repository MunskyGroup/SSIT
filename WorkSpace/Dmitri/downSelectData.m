CC = readtable('all_data_2017.csv');
TMP = table;
TMP.sca = floor(CC.sca/20);
TMP.time = CC.time;

cols = {'plant_date','seed_trtmt','spray_trtmt'};
for i = 1:length(cols)
    TMP.(cols{i}) = CC.(cols{i});
end

writetable(TMP,'AphidsNewFormatB.csv')