function [Distances, Accepted, Parameters] = MHABC(nChains, nSamples, trueResults, Simulator, Priors, x0, distanceMetric)
    Distances = zeros(nChains, nSamples+1);
    Accepted = zeros(nChains, nSamples+1);
    Parameters = zeros(nChains, nSamples+1, length(x0));

    U = rand(nChains, nSamples);

    for c = 1:nChains
        Parameters(c, 1, :) = x0;
        x = Parameters(c, 1, :);
        results_x = Simulator(Parameters(c, 1, :)); 
        d_best = distanceMetric(results_x, trueResults);
        Distances(c, 1) = d_best;
        
        for s = 1:nSamples
            Parameters(c, s+1, :) = Priors(x);
            results_y = Simulator(Parameters(c, s+1, :));
            d_new = distanceMetric(results_y, trueResults);
            Distances(c, s+1) = d_new;
            rho = d_new/d_best;
    
            Ui = U(c,s);
            acc = Ui<=min(rho,1);

            if acc
                d_best = d_new;
                x = Parameters(c, s+1, :);
            end
        end
    end
end
