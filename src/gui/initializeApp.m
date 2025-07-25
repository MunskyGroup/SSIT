function initializeApp(app) 
%% Load an initial model
AA = dir('Models/');
v2 ={AA.name};
I = zeros(1,length(v2),'logical');
for i=1:length(v2)
    if AA(i).name(1) ~= '.'
        I(i)=1;
    end
end
v2 = v2(I);
app.ModelUsePresetExampleTypeDropDown.Items = v2;
app.ModelUsePresetExampleTypeDropDown.Value = 'gene_expression';
