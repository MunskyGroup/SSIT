function [] = addReactionToReactionsTab(app)
% This script adds a new reaction in a new row to the reactions table in
% the reactions tab of the GUI. Through this, the user can then edit the
% defaults to correspond to the reactions in their model.

%%
species = app.SSITModel.species(:,1);
nSp = size(species,1);        
x1 = species{randi(nSp)};
x2 = species{randi(nSp)};
x3 = species{randi(nSp)};

N = size(app.ModelReactionTable.Data,1);          % This evaluates the current number of rows of the Reaction Table
TMP = app.ModelReactionTable.Data;                % Assigns a variable to represent the Reactions Table
TMP{N+1,1} = ['R',num2str(N+1)];              % Adds a new row and places the next number in the sequence as the entry
if strcmp(x1,x2)
    TMP{N+1,2} = [x1,'(2)'];             % Places the value in the second column of the new row
    TMP{N+1,4} = ['k1_*',x2,'*(',x1,'-1)/2'];  % Places the value in the fourth column of the new row
else
    TMP{N+1,2} = [x1,'(2),',x2,'(1)'];             % Places the value in the second column of the new row
    TMP{N+1,4} = ['k1_*',x2,'*(',x1,'-1)*',x2,'/2'];  % Places the value in the fourth column of the new row
end
TMP{N+1,3} = [x3,'(1)'];                   % Places the value in the third column of the new row

prompt = {'Reactants (e.g., "x1(1),x2(2)"):','Products (e.g., "x3(2)"):','Propensity (e.g., k1*x1*x2*(x2-1)/2)'};
dlgtitle = 'New Reaction';
dims = [1 50];  % [height width] of input field
defaultInput = {TMP{N+1,2},TMP{N+1,3},TMP{N+1,4}};
newRxn = inputdlg(prompt, dlgtitle, dims, defaultInput);
% rxn = ssit.parest.propsStorage;
% rxn.props = struct('Reactants',TMP{N+1,2},...
% 'Products',TMP{N+1,3},...
% 'Propensity',TMP{N+1,4});
% propEdtr = ssit.parest.PropEditor(rxn,'props');
% uiwait(propEdtr.UIFigure); % Wait until the app window is closed
% propEdtr.SSIT_GUI
% 
TMP{N+1,2} = newRxn{1};
TMP{N+1,3} = newRxn{2};
TMP{N+1,4} = newRxn{3};

app.ModelReactionTable.Data = TMP;                % Re-assigns the updated values to the Reaction Table
app.ModelhasnotbeenupdatedLabel.Text = 'Model has not been updated.'; % This updates the text by the Update Model Button to notify the user that the model parameters have not been updated

app.DeleteReactionDropDown.Items{end+1} = TMP{N+1,1};
end