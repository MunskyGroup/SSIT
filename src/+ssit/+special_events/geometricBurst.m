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
for j = 1:nStates
    vi = geopdf(states-states(:,j), repmat(1./(burstParameters+1),1,nStates));
    B(:,j) = prop_val(j)*(1-sum(vi,2));
    v = prod(vi,1);
    v(abs(v)<=1e-12) = 0;
    v(j) = v(j) - 1;
    A(:,j) = prop_val(j)*v;
end
A = sparse(A);
A(nStates+nSpecies+1:nStates+nSpecies*2,:) = B;
A(nStates+numConstraints,nStates+numConstraints) = 0;


