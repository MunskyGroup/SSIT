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
if ~isempty(v2)
    if max(strcmp(v2,'gene_expression'))
        app.ModelUsePresetExampleTypeDropDown.Value = 'gene_expression';
    else
        app.ModelUsePresetExampleTypeDropDown.Value = v2{1};
    end
end
app.SSITModel = SSIT('Empty');
if ~isempty(app.ModelDropDown.Items)
    if max(strcmp(app.ModelDropDown.Items,'M00_Poisson_Process.mat'))
        app.ModelDropDown.Value = 'M00_Poisson_Process.mat';
    else
        app.ModelUsePresetExampleTypeDropDown.Value = app.ModelDropDown.Items{1};
    end
end
