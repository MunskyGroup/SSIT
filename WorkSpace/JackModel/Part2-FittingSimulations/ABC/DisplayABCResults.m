
load("ABC_Results.mat")

parNames = {'kon', 'koff', 'w', 'kex', 'kr', ...
    'D1', 'D2', 'D3', 'gam1', 'gam2', 'gam3'};

params = vertcat(Parameters(:));
params = cell2mat(params);
dist = reshape(Distances.', 1, []);
accept = reshape(Accepted.', 1, []);


figure;
nParams = size(params,2);
for i = 2:nParams
    for j = 1:i-1
        subplot(nParams-1, nParams-1, (i-2)*(nParams-1)+j);

        scatter(params(:, j), params(:,i), 40, dist, 'filled');

        if j == 1
            ylabel(['Param ', parNames(i)]);
        end

        if i == nParams
            xlabel(['Param ', parNames(j)]);
        end

        if i<nParams
            set(gca, 'XTickLabel', []);
        end

        if j>1
            set(gca, 'XTickLabel', []);
        end
    end
end






















