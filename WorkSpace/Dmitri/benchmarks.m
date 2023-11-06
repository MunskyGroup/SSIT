clear all
%% Random Matrix
M = 5;
Ns = 6;
Nsim = 5;
Nar = floor(logspace(2,3,Ns));
timeExpm =nan*ones(Nsim,Ns);
timeExpokit =nan*ones(Nsim,Ns);
diff_expm_expokit =nan*ones(Nsim,Ns);
diff_expokit_ode23s=nan*ones(Nsim,Ns);
diff_expm_ode23s=nan*ones(Nsim,Ns);
timeODE23s=nan*ones(Nsim,Ns);
t = 1;
for n = length(Nar):-1:1
    N = Nar(n);

    for j = 1:Nsim
        A = rand(N,N);
        thr = min(maxk(A,5));
        A(A<thr)=0;
        A = sparse(A);
        A = A - diag(sum(A));
        P0 = rand(N,1);
        P0 = P0/sum(P0);

        %% Expm calculation
        tic
        expAt_P = expm(A*t)*P0;
        timeExpm(j,n) = toc;

        %% Expokit Calculation
        m = 15;
        tryAgain=1;
        %     fspTol = fspErrorCondition.fspTol;
        %     nSinks = fspErrorCondition.nSinks;
        fspErrorCondition.tInit = 0;
        tic
        while tryAgain==1
            [~, ~, ~, tExport, solutionsNow, ~, tryAgain, te, PfExpokit] = ssit.fsp_ode_solvers.mexpv_modified_2(t, ...
                A, P0, 1e-8, m,...
                [], [0,t], 1e-3,[], 0, fspErrorCondition);
            if tryAgain==0;break;end
            if m>300
                warning('Expokit expansion truncated at 300');
                [~, ~, ~, tExport, solutionsNow, ~, tryAgain, te, PfExpokit] = ssit.fsp_ode_solvers.mexpv_modified_2(tOut(end), jac, initSolution, fspTol/1e5, m,...
                    [], tOut, fspTol, [length(initSolution)-nSinks+1:length(initSolution)], tStart, fspErrorCondition);
            end
            m=m+5;
        end
        timeExpokit(j,n) = toc;

        %% ODE Calculation
        ode_opts = odeset('Jacobian', A, 'Vectorized','on','JPattern',A~=0,'relTol',1e-8, 'absTol', 1e-10);
        rhs = @(t,x)A*x;
        tic
        [tExport, yout] =  ode15s(rhs, [0,t/2,t], P0);
        timeODE23s(j,n) = toc;

        diff_expm_expokit(j,n) = sum(abs(PfExpokit-expAt_P));
        diff_expm_ode23s(j,n) = sum(abs(expAt_P-yout(end,:)'));
        diff_expokit_ode23s(j,n) = sum(abs(PfExpokit-yout(end,:)'));
    end
    %% Comparisons
    % Compare time
    figure(2)
    subplot(2,1,1);
    loglog(Nar,mean(timeExpm),'-s',Nar,mean(timeExpokit),'-o',Nar,mean(timeODE23s),'-^');
    legend('expm','expokit','ode23s')
    xlabel('Number of states')
    ylabel('Computational time')
    set(gca,'fontsize',15)

    % Compare results
    subplot(2,1,2);
    plot(Nar,mean(diff_expm_expokit),'-o',Nar,mean(diff_expm_ode23s),'-s',Nar,mean(diff_expokit_ode23s),'-^')
    legend('diff_expm_expokit','diff_expm_ode23s','diff_expokit_ode23s')
end

