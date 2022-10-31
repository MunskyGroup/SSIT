function loadModelFromSBML(app)

% Let user select SBML file
[FILENAME,PATHNAME] = uigetfile('*.xml','Select your SBML Model','MultiSelect','on');
sbmlobj = sbmlimport([PATHNAME,FILENAME]);
nR = length(sbmlobj.Reactions);
nS = length(sbmlobj.Species);

% check if 3 or fewer species.
if nS>3
    error('SSIT not yet configured for more than 3 species.')
end

[S,objSpecies,objReactions]= getstoichmatrix(sbmlobj);

for i = 1:nR
    app.ModelReactionTable.Data{i,:} = {num2str(i),'','','-','n'}; 
    for j = 1:nS
        if S(j,i)<0
            if isempty(app.ModelReactionTable.Data{i,2})
                app.ModelReactionTable.Data{i,2} = ['x',num2str(j),'(',num2str(-S(j,i)),')'];
            else
                app.ModelReactionTable.Data{i,2} = [app.ModelReactionTable.Data{i,2},',x',num2str(j),'(',num2str(-S(j,i)),')'];
            end
        elseif S(j,i)>0
            if isempty(app.ModelReactionTable.Data{i,3})
                app.ModelReactionTable.Data{i,3} = ['x',num2str(j),'(',num2str(-S(j,i)),')'];
            else
                app.ModelReactionTable.Data{i,3} = [app.ModelReactionTable.Data{i,3},',x',num2str(j),'(',num2str(-S(j,i)),')'];
            end
        end
    end
end
