function makeGRPlots(combinedModel,GRpars,splitReps)
arguments
    combinedModel
    GRpars
    splitReps = false
end

GR_Data = 'RonData062025/GR_ALL_gated_with_CytoArea_and_normGR_Feb2825_03.csv';

GRfitCases = {'1','1',101,'GR Fit (1nM Dex)';...
    '10','10',102,'GR Fit (10nM Dex)';...
    '100','100',103,'GR Fit (100nM Dex)'};

combinedGRModel = combinedModel.updateModels(GRpars,false);
nMods = length(combinedGRModel.SSITModels);
ModelGroup = cell(nMods,1);
for i=1:nMods
    if splitReps
        for rep = {'A','B','C'}
            %  Update parameters in original models.
            ModelGroup{i} = combinedGRModel.SSITModels{i};
            ModelGroup{i} = ModelGroup{i}.loadData(GR_Data,...
                {'nucGR','normGRnuc';'cytGR','normGRcyt'},...
                {[],[], ...
                ['(TAB.dex_conc==',GRfitCases{i,1},...
                '|TAB.dex_conc==0)&TAB.time~=20' ...
                '&TAB.time~=40&TAB.time~=60' ...
                '&TAB.time~=90&TAB.time~=150' ...
                '&strcmp(TAB.replica,''',rep{1},''')']});
            ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
            ModelGroup{i}.makeFitPlot([],1,[(i-1)*4+1:i*4],true,'STD');
        end
    else
        %  Update parameters in original models.
        ModelGroup{i} = combinedGRModel.SSITModels{i};
        ModelGroup{i} = ModelGroup{i}.loadData(GR_Data,...
            {'nucGR','normGRnuc';'cytGR','normGRcyt'},...
            {[],[], ...
            ['(TAB.dex_conc==',GRfitCases{i,1},...
            '|TAB.dex_conc==0)&TAB.time~=20' ...
            '&TAB.time~=40&TAB.time~=60' ...
            '&TAB.time~=90&TAB.time~=150']});
        ModelGroup{i}.tSpan = sort(unique([ModelGroup{i}.tSpan,linspace(0,180,30)]));
        ModelGroup{i}.makeFitPlot([],1,[(i-1)*4+1:i*4],true,'STD');
    end
end
end