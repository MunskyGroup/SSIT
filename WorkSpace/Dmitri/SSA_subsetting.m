T = readtable('ssa_out.csv');
Tsize = size(T);
rows = 1:Tsize(1, 1);

subset_rows = (mod(rows, 30) == 0);
sum(subset_rows)
subset_T = T(subset_rows, :);
writetable(subset_T, 'ssa_1000.csv');

subset_rows = (mod(rows, 15) == 0);
sum(subset_rows)
subset_T = T(subset_rows, :);
writetable(subset_T, 'ssa_2000.csv');

subset_rows = (mod(rows, 6) == 0);
sum(subset_rows)
subset_T = T(subset_rows, :);
writetable(subset_T, 'ssa_5000.csv');

subset_rows = (mod(rows, 3) == 0);
sum(subset_rows)
subset_T = T(subset_rows, :);
writetable(subset_T, 'ssa_10000.csv');

subset_rows = (mod(rows, 2) == 0);
sum(subset_rows)
subset_T = T(subset_rows, :);
writetable(subset_T, 'ssa_15000.csv');

subset_rows = (mod(rows, 3) > 0); % This will collect 2/3 values
sum(subset_rows)
subset_T = T(subset_rows, :);
writetable(subset_T, 'ssa_20000.csv');

writetable(T, 'ssa_30000.csv'); % This is the entire table