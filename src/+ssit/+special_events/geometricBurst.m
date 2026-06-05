function A = geometricBurst(states,burstParameters,prop_val,args,numConstraints)
arguments
    states
    burstParameters
    prop_val
    args
    numConstraints
end
% This function generates an infitesimal generator for a geometric burst
% reaction, where the molecules of each species are generated acording to
% the burst size.
%
% Inputs:
%   states -- list of all states I(nSpecies x nStates)
%   burstParameters -- burst size for each species
%   prop_val -- reaction propensity for burst timing
%   args -- unused but kept to maintain format with other special events
[nSpecies,nStates] = size(states);
A = zeros(nStates,nStates);
B = zeros(nSpecies,nStates);
uniqueStates = cell(1,nSpecies);
uniqueStatesInds = cell(1,nSpecies);
for i = 1:nSpecies
    [uniqueStates{i},uniqueStatesInds{i}] = unique(states(i,:));
end
for j = 1:nStates
    vi = geopdf(states-states(:,j), repmat(1./(burstParameters+1),1,nStates));
    for i = 1:nSpecies
        statesConsistent = min(states([1:i-1,i+1:end],:)==states([1:i-1,i+1:end],j),[],1);
        B(i,j) = 1 - sum(vi(i,statesConsistent));
    end
    v = prod(vi,1);
    v(abs(v)<=1e-12) = 0;
    v(j) = v(j) - 1;

    A(:,j) = prop_val(j)*v;
end
A = sparse(A);
A(nStates+nSpecies+1:nStates+nSpecies*2,:) = prop_val(j)*B;
A(nStates+numConstraints,nStates+numConstraints) = 0;


