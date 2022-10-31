function FSPResultsOut = exportFSPResults(app)
% The results of the FSP for marginal and joint distributions have been
% exported to the workspace as global variables.  To access and plot these
% you can run the codes:
%     global FSPResults
%     plot(FSPResults.T_array,FSPResults.Means(:,1));
%       % plots the mean of species 1 versus time.
%     errorbar(FSPResults.T_array,FSPResults.Means(:,1),sqrt(FSPResults.Var(:,1)));
%       % plots the mean +/- 1 STD of species 1 versus time.
%     plot(FSPResults.Marginals{5}{1});
%       % makes a marginal plot of the pdf for species 1 at time 5.
%     plot(FSPResults.Marginals{5}{2});
%       % makes a marginal plot of the pdf for species 2 at time 5.
%     plot(FSPResults.Marginals{5}{3});
%       % makes a marginal plot of the pdf for species 3 at time 5.
%     contourf(FSPResults.Joints{5}{2,3});
%       % makes a contour plot of the pdf for species 2 and 3 at time 5
%     contourf(FSPResults.Joints{5}{1,3});
%       % makes a contour plot of the pdf for species 1 and 3 at time 5
%     contourf(FSPResults.Joints{5}{1,2});
%       % makes a contour plot of the pdf for species 1 and 2 at time 5

if nargout==0
    global FSPResults
end
FSPResults = [];
Nd = app.FspTabOutputs.solutions{end}.p.dim;
FSPResults.T_array = eval(app.FspPrintTimesField.Value);
FSPResults.Means = zeros(length(FSPResults.T_array),Nd);
FSPResults.Var = zeros(length(FSPResults.T_array),Nd);
for it = length(FSPResults.T_array):-1:1
    for i=1:Nd
        if Nd>1
            for j=i+1:Nd
                k = setdiff([1:Nd],[i,j]);
                if Nd>2
                    joints{i,j} = double(app.FspTabOutputs.solutions{it}.p.sumOver(k).data);
                else
                    joints{i,j} = double(app.FspTabOutputs.solutions{it}.p.data);
                end
            end
            sumInds = setdiff([1:Nd],i);
            marginals{i} = double(app.FspTabOutputs.solutions{it}.p.sumOver(sumInds).data);
        else
            joints = [];
            marginals{i} = double(app.FspTabOutputs.solutions{it}.p.data);
        end
        FSPResults.Means(it,i) = [0:length(marginals{i})-1]*marginals{i};
        SqMeans(it,i) = [0:length(marginals{i})-1].^2*marginals{i};
        FSPResults.Var(it,i) = SqMeans(it,i) - FSPResults.Means(it,i)^2;
    end
    FSPResults.Joints{it} = joints;
    FSPResults.Marginals{it} = marginals;
end
if nargout==0
    help exportFSPResults
else
    FSPResultsOut=FSPResults;
end
