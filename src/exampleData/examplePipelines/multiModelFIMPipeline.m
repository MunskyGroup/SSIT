function [results, model] = multiModelFIMPipeline(model, Args)
    J = Args.param_of_interest_index;
    model.parameters = Args.pars;

    % (2) Solve FSP for model
    model.solutionScheme = 'FSP';    % Set solution scheme to FSP.
    [FSPsoln,model.fspOptions.bounds] = model.solve;  % Solve the FSP analysis
    
    if Args.makePlots
        % Plot the results from the FSP analysis
        % for n=model.species
        % fig1 = figure(21);clf; set(fig1,'Name','Marginal Distributions, gene_on');
        % fig2 = figure(22);clf; set(fig2,'Name','Marginal Distributions, gene_off');
        % fig3 = figure(23);clf; set(fig3,'Name','Marginal Distributions, RNA');
        % fig4 = figure(24);clf; set(fig4,'Name','Marginal Distributions, Distant RNA');
        % end
        model.makePlot(FSPsoln,'marginals',[1:4:21],false,[],{'r-','linewidth',2})
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
        B = A([1:J-1, J+1:end], [1:J-1, J+1:end]);
        cov_unknown = C([1:J-1, J+1:end], [1:J-1, J+1:end]);
     end
     cov_known = inv(B);

    known_determinate = det(cov_known);
    unknown_determinate = det(cov_unknown);
    determinates = table(known_determinate, unknown_determinate);
    
    results.determinates = determinates;
    results.params = model.parameters;


end


