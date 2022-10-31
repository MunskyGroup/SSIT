function SenResults = exportSensResults(app)
% The results of the Sensitivity for marginal distributions have been
% exported to the workspace as global variables.  To access and plot these
% you can run the codes:
%     global SenResults
%     stairs(SenResults.sensmdist{2,1,5}); 
%       % makes a marginal plot of the sensitivity to parameter 2, for 
%       % species 1 at time 5.
 
if nargout==0
    global SenResults
end
SenResults = [];
SenResults.T_array = unique(eval(app.SensPrintTimesEditField.Value));
Nt = length(SenResults.T_array);
Np = length(app.ReactionsTabOutputs.parameters);
Nd = app.SensFspTabOutputs.solutions.data{1}.p.dim;

for ip=Np:-1:1
    for it=Nt:-1:1
        for id=Nd:-1:1
            INDS = setdiff([1:Nd],id);
            if ~isempty(INDS)
                SenResults.sensmdist{ip,id,it} = double(app.SensFspTabOutputs.solutions.data{it}.S(ip).sumOver(INDS).data);
            else
                SenResults.sensmdist{ip,id,it} = double(app.SensFspTabOutputs.solutions.data{it}.S(ip).data);
            end
        end
    end
end

if nargout==0
    help exportSensResults
end
