function [X_array] = runSingleSsa(x0,S,W,T_array,isTimeVarying,signalUpdateRate,parameters)
% Start the simulation.
t = min(T_array);   % initial time of simulation.
x = x0;
iprint = 1;  % The next time at which to record results.
Nsp = size(S,1);  %Number of species.
Nt = length(T_array);
X_array = zeros(Nsp,Nt);
props = [];
if isTimeVarying
        S = [zeros(Nsp,1),S];
        props(1) = signalUpdateRate;
        jt = 1;
else
    jt=0;
end

isAFun = isa(W{1}, 'function_handle');
    
while t<max(T_array)

    %% Choose time of reaction
    if ~isAFun
        %% Command Line Code Approach
        WT = W{1}.hybridFactorVector(t,parameters);
        for i=length(W):-1:1
            if ~W{i}.isTimeDependent||W{i}.isFactorizable
                props(i+jt) = WT(i)*W{i}.stateDependentFactor(x,parameters);     % evaluate the propensity functions at the current state.
            else
                props(i+jt) = W{i}.jointDependentFactor(t,x,parameters);     % evaluate the propensity functions at the current state.
            end
        end
    else
        %% GUI approach
        for i=length(W):-1:1
            props(i+jt) = W{i}(x,t);     % evaluate the propensity functions at the current state.
        end
    end
    w0 = sum(props);  % sum the propensity functions (inverse of ave. waiting time).
    tau = -1/w0*log(rand); % The time until the next reaction.

    %% update time
    t = t+tau;  % The time of the next reaction.
    
    while iprint<=Nt&&t>T_array(iprint)
        X_array(:,iprint) = x;
        iprint=iprint+1;
    end  
        
    if t>max(T_array)
        break;
    end
    
    %% Choose which reaction.
    r2 = w0*rand;    % this is a uniform r.v. betweeen zero and w0;
    j = 1;
    while sum(props(1:j))<r2
        j=j+1;
    end
    % At this point j is the chosen reaction.
    
    %% Update state
    x = x + S(:,j);
end

