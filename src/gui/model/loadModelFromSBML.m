function loadModelFromSBML(app)

app.ModelInputTable.Data = {};
app.ModelReactionTable.Data = {};

% Let user select SBML file
[FILENAME,PATHNAME] = uigetfile('*.xml','Select your SBML Model','MultiSelect','on');
sbmlobj = sbmlimport([PATHNAME,FILENAME]);
nR = length(sbmlobj.Reactions);
nS = length(sbmlobj.Species);

% check if 3 or fewer species.
if nS>3
    error('SSIT not yet configured for more than 3 species.')
elseif nS==0
    error('SBML File has zero species!')
elseif nR==0
    error('SBML File has zero reactions!')
end


[S,objSpecies,objReactions]= getstoichmatrix(sbmlobj);

app.ModelReactionTable.Data = {'1','','','-','n'}; 
for i = 1:nR
    app.ModelReactionTable.Data(i,1:5) = {num2str(i),'','','-','n'}; 
    for j = 1:nS
        if S(j,i)<0
            if isempty(app.ModelReactionTable.Data{i,2})
                app.ModelReactionTable.Data{i,2} = ['x',num2str(j),'(',num2str(-S(j,i)),')'];
            else
                app.ModelReactionTable.Data{i,2} = [app.ModelReactionTable.Data{i,2},',x',num2str(j),'(',num2str(-S(j,i)),')'];
            end
        elseif S(j,i)>0
            if isempty(app.ModelReactionTable.Data{i,3})
                app.ModelReactionTable.Data{i,3} = ['x',num2str(j),'(',num2str(S(j,i)),')'];
            else
                app.ModelReactionTable.Data{i,3} = [app.ModelReactionTable.Data{i,3},',x',num2str(j),'(',num2str(S(j,i)),')'];
            end
        end
    end
    
    txt = sbmlobj.Parameter(i).Name;
    for k=1:size(sbmlobj.Reactions(i).Reactants,1)
        txt = [txt,'*x',sbmlobj.Reactions(i).Reactants(k).Name(2:end)];
    end
    app.ModelReactionTable.Data{i,4} = txt;
    
    app.ModelParameterTable.Data{i,1} = sbmlobj.Parameter(i).Name;
    app.ModelParameterTable.Data{i,2} = sbmlobj.Parameter(i).Value;
        
end
app.ReactionsTabOutputs.parameters = app.ModelParameterTable.Data;
app.ReactionsTabOutputs.presetParameters = app.ModelParameterTable.Data(:,2);
app.ReactionsTabOutputs.inputs ={};

IC = zeros(1,3); IC(1:nS) = round(1000*[sbmlobj.Species.Value]);
app.ReactionsTabOutputs.initialCondition = ['[',num2str(IC),']'];

app.NameTable.Data = {};
app.NameTable.Data(1:3,1) = {'x1';'x2';'x3'};
app.NameTable.Data(1:3,2) = {'x1';'x2';'x3'};
app.NameTable.Data(1:nS,2) = {sbmlobj.Species.Name};

end
