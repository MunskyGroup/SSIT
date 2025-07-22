function [Distances, Accepted, Parameters] = MHABC(nChains, nSamples, trueResults, Simulator, Proposal, x0, distanceMetric)

    % Use cells for Parameters to avoid 3D slicing issues in parfor
    Parameters = cell(nChains,1);
    Distances = zeros(nChains, nSamples + 1);
    Accepted = zeros(nChains, nSamples + 1);
    
    % Random U values
    U = rand(nChains, nSamples);
    
    parfor c = 1:nChains
        % Local storage for each chain
        localParams = zeros(nSamples + 1, length(x0));
        localAccepted = zeros(1, nSamples + 1);
        localDistances = zeros(1, nSamples + 1);
        
        % Initialize
        x = x0;  % same initial guess for all chains (or make this unique)
        localParams(1, :) = x;
        results_x = Simulator(x);
        d_best = distanceMetric(results_x, trueResults);
        localDistances(1) = d_best;
        
        numAccepted = 0;
        for s = 1:nSamples
            x_prop = Proposal(x);
            
            results_y = Simulator(x_prop);
            d_new = distanceMetric(results_y, trueResults);
            localDistances(s + 1) = d_new;
            rho = d_new-d_best;
            acc = U(c, s) <= exp(rho);
            localAccepted(s + 1) = acc;
            
            if acc
                x = x_prop;
                d_best = d_new;
                numAccepted = numAccepted + 1;
            end
            
            localParams(s + 1, :) = x_prop;
            fprintf('Chain: %d, Sim: %d \n', c, s)
            fprintf('acceptance rate: %1.2f \n', numAccepted/s)
        end
        
        % Store this chain's results
        Parameters{c} = localParams;
        Accepted(c, :) = localAccepted;
        Distances(c, :) = localDistances;
    end
end