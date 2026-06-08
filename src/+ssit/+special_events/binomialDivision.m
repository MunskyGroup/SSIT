function A = binomialDivision(states,binomialParameters,prop_val,args,numConstraints)
arguments
    states
    binomialParameters
    prop_val
    args
    numConstraints
end
% This function generates an infitesimal generator for a division reaction,
% where the molecules of each species are distributed to daughter cells
% with a binomial distribution. 
%
% Inputs:
%   states -- list of all states I(nSpecies x nStates)
%   binomialParameters -- proabilities for daughter to retain each molecule
%   prop_val -- reaction propensity for division timing
%   args: Structure containing arguments
%       args.lineage -- type of division: 
%               'primary' -- follows only specified daughter
%               'bothCoupled' -- follows both daughters, all species are
%                   coupled, meaning that the binomial random numbers are
%                   dependent on each other.
%               'bothUncoupled' -- follows both daughters, all species are
%                   uncoupled, meaning that binomial random numbers are
%                   independent.
%   numConstraints -- number of FSP constraints.
[nSpecies,nStates] = size(states);
A = zeros(nStates+numConstraints);
for j = 1:nStates
    switch args.lineage
        case 'primary'
            v = prod(binopdf(states, repmat(states(:,j),1,nStates), repmat(binomialParameters,1,nStates)),1);
        case 'bothCoupled'
            v = 1/2*(prod(binopdf(states, repmat(states(:,j),1,nStates), repmat(binomialParameters,1,nStates)),1)+...
                prod(binopdf(states, repmat(states(:,j),1,nStates), 1-repmat(binomialParameters,1,nStates)),1));
        case 'bothUncoupled'
            v = 0;
            nPossibilities = 2^nSpecies;
            for iP = 1:nPossibilities
                bin = logical(dec2bin(iP-1,nSpecies)');
                q = binomialParameters;
                q(~bin) = 1 - binomialParameters(~bin);
                v = v + 1/nPossibilities*prod(binopdf(states, repmat(states(:,j),1,nStates), repmat(q,1,nStates)),1);
            end
    end
    v(abs(v)<=1e-12) = 0;
    % v(j) = v(j) - 1;
    % A(:,j) = prop_val(j)*v;

    if isfield(args,'FixedTime')&&args.FixedTime
        A(1:nStates,j) = v;
        % A(nStates+nSpecies+1:nStates+nSpecies*2,j) = B;
    else
        v(j) = v(j) - 1;
        A(1:nStates,j) = prop_val(j)*v;
        % A(nStates+nSpecies+1:nStates+nSpecies*2,j) = prop_val(j)*B;
    end



end
A = sparse(A);
% TODO - this approach may not allow for proper downward FSP expansion as
% transitions into sinks are being set to zero.
% A(end+numConstraints,end+numConstraints) = 0;

