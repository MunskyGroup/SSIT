function app=compareFspMomentsAndOde(app)
% function to find the difference between the fsp mean and the ode mean and
% the difference between their variences

T_array = eval(app.FspPrintTimesField.Value);
odeMean=app.FspTabOutputs.odeSolutions;

if isempty(odeMean)
    diff=[]
    odeError=[]
    
else
    for it = 1:length(T_array)
        mdist = ssit.fsp.marginals(app.FspTabOutputs.solutions{it}.states, app.FspTabOutputs.solutions{it}.p);
        for j=1:3
            fspMean(it,j) = [0:length(mdist{j})-1]*mdist{j};
        end
    end
    
    % find the error between the mean fsp solution and the ode solution
    diff=fspMean-odeMean
    
    Error = diff./odeMean;
    J = odeMean==0;
    odeError=sum(abs(Error(~J)))./sum(~J)
end
