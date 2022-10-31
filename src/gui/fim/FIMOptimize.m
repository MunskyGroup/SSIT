function FIMOptimize(app)
%% Define metric for optimization
switch app.FIMMetricorParameterDropDown.Value
    case 'Determinant'
        met = @(A)-det(A);
    case 'Smallest Eigenvalue'
        met = @(A)-min(eig(A));
    case 'Trace'  
        met = @(A)-trace(A);
    otherwise
        k = find(strcmp(app.SensParDropDown.Items,app.FIMMetricorParameterDropDown.Value));
        ek = zeros(1,length(app.FIMTabOutputs.FIMMatrices{1}));ek(k) = 1;
        met = @(A)(ek*inv(A)*ek');
end
%%
NT = size(app.FIMTabOutputs.FIMMatrices(:,1),1);
if isempty(app.FIMTabOutputs.NcOptimized)
    Nc = zeros(1,NT); 
    Nc(1)=app.CellsEditField.Value;
else
    Nc = app.FIMTabOutputs.NcOptimized;
end
 
Converged = 0; 
while Converged==0
    Converged = 1;
    for i = 1:NT
        while Nc(i)>0
            Ncp = Nc;
            Ncp(i) = Ncp(i)-1;
            k = findBestMove(app.FIMTabOutputs.FIMMatrices(:,1),Ncp,met);
            if k==i
                break
            end
            Nc = Ncp; Nc(k)=Nc(k)+1;
            Converged = 0;
        end
    end
    plot(Nc)
    drawnow
end
app.FIMTabOutputs.NcOptimized = Nc;
for i=1:length(Nc)
    app.FIMTabOutputs.CellsPerTimePoint.props.(['t',num2str(i)])(2) = ...
        [Nc(i)];
end
end

function FIM = totalFim(fims,Nc)
    FIM = 0*fims{1};
    for i = 1:length(fims)
        FIM = FIM+Nc(i)*fims{i};
    end
end

function k = findBestMove(fims,Ncp,met)
obj = zeros(1,length(Ncp));  
FIM0 = totalFim(fims,Ncp);
for i = 1:length(Ncp)
     FIM = FIM0+fims{i};
     obj(i) = met(FIM);
end
[~,k] = min(obj);
end