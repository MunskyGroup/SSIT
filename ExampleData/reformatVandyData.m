%% reformatNuertData
% In this script, we reformat the data from Gregor Neuert's lab so that it
% can be loaded into the SSIT.

FileList = {'Result_Exp1_rep1_RNA_CY5_total'
'Result_Exp1_rep2_RNA_CY5_total'
'Result_Exp1_rep2_RNA_TMR_total'
'Result_Exp2_rep1_RNA_CY5_total'
'Result_Exp2_rep1_RNA_TMR_total'
'Result_Exp2_rep2_RNA_CY5_total'
'Result_Exp2_rep2_RNA_TMR_total'
'Result_Exp2_rep3_RNA_CY5_total'};

for i=1:length(FileList)
    FN = ['../ExampleData/NeuertData/',FileList{i}];
    X = importdata([FN,'.csv']);
    timearray = [0,1,2,4,6,8,10,15,20,25,30,35,40,45,50,55];
    D = [];
    for i = 1:16
        time = timearray(i);
        d = X.data(:,i);
        d = d(~isnan(d));
        D = [D;time*ones(size(d)),d];
    end
    T = array2table(D,"VariableNames",{'time','mRNA'});
    writetable(T,[FN,'_FORMATTED.csv'])
end