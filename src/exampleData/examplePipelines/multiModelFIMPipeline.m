function [results, model] = multiModelFIMPipeline(model, Args)
    J = Args.param_of_interest_index;
    model.parameters = Args.pars;

    % (2) Solve FSP for model
    model.solutionScheme = 'FSP';    % Set solution scheme to FSP.
    [FSPsoln,model.fspOptions.bounds] = model.solve;  % Solve the FSP analysis
    
    if Args.makePlots
        % Plot the results from the FSP analysis
        fig21 = figure(21);clf; set(fig21,'Name','Marginal Distributions, gene_on');
        fig22 = figure(22);clf; set(fig22,'Name','Marginal Distributions, gene_off');
        fig23 = figure(23);clf; set(fig23,'Name','Marginal Distributions, RNA');
        fig24 = figure(24);clf; set(fig24,'Name','Marginal Distributions, Distant RNA');
        model.makePlot(FSPsoln,'marginals',[1:4:21],false,[fig21,fig22,fig23, fig24],{'r-','linewidth',2})
    end

    % (3) Solve FSP Sensitivity
    model.solutionScheme = 'fspSens'; % Set solutions scheme to FSP Sensitivity
    [sensSoln,bounds] = model.solve(FSPsoln.stateSpace);  % Solve the sensitivity problem
    
    % (4) Compute FIM using FSP Sensitivity Results
    fimResults = model.computeFIM(sensSoln.sens); % Compute the FIM for full observations and no distortion.
    cellCounts = 100*ones(size(model.tSpan));  % Number of cells in each experiment.
    [FIMTotal,mleCovEstimate,fimMetrics] = model.evaluateExperiment(fimResults,cellCounts);

    A = cell2mat(FIMTotal);
    C = inv(A);
     if isnan(J)
        % If J is NaN, use the entire matrix
        B = A;
        cov_unknown = C;
    else
        % Otherwise, use the specified submatrix
        B = A(J, J);
        cov_unknown = C(J, J);
     end
     cov_known = inv(B);

    known_determinate = det(cov_known);
    unknown_determinate = det(cov_unknown);
    determinates = table(known_determinate, unknown_determinate);
    
    results.determinates = determinates;


end


