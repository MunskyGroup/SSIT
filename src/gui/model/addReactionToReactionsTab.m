function [] = addReactionToReactionsTab(app)
% This script adds a new reaction in a new row to the reactions table in
% the reactions tab of the GUI. Through this, the user can then edit the
% defaults to correspond to the reactions in their model.

            N = size(app.ModelReactionTable.Data,1);          % This evaluates the current number of rows of the Reaction Table
            TMP = app.ModelReactionTable.Data;                % Assigns a variable to represent the Reactions Table
            TMP{N+1,1} = num2str(N+1);              % Adds a new row and places the next number in the sequence as the entry
            TMP{N+1,2} = 'x1(2),x3(1)';             % Places the value in the second column of the new row
            TMP{N+1,3} = 'x2(1)';                   % Places the value in the third column of the new row
            TMP{N+1,4} = 'I1_*k1_*x1*(x1-1)*x3/2';  % Places the value in the fourth column of the new row
            TMP{N+1,5} = 'n';                       % Places a no to keep the new row in the new row
            app.ModelReactionTable.Data = TMP;                % Re-assigns the updated values to the Reaction Table
            app.ModelhasnotbeenupdatedLabel.Text = 'Model has not been updated.'; % This updates the text by the Update Model Button to notify the user that the model parameters have not been updated
end