function sol = poolSSABatches(solA, npeA, solB, npeB)
% ssit.ssa.poolSSABatches - Combine two SSA solutions (from SSIT.solve with
% ssaOptions.nSimsPerExpt = npeA and npeB, Nexp = 1) into one solution with
% npeA+npeB trajectories per block, keeping the block layout that
% SSIT.computeLossFunctionSSA expects with enforceIndependence = true.
%
%   sol = ssit.ssa.poolSSABatches(solA, npeA, solB, npeB)
%
% SSA trajectories are stored as trajs(species, time, k), with k running
% over nBlocks consecutive blocks of nSimsPerExpt trajectories (one block
% per time point, so different time points can use independent
% trajectories). Block b of the result is block b of solA followed by
% block b of solB. Fields 'trajs' and (if present) 'trajsDistorted' are
% pooled; all other fields are taken from solB. If solA is empty, solB is
% returned unchanged. For solutions without blocks (enforceIndependence =
% false), pass npeA and npeB = the numbers of trajectories; the result is
% then a plain concatenation.
if isempty(solA) || npeA == 0
    sol = solB;
    return
end
sol = solB;
for f = {'trajs','trajsDistorted'}
    name = f{1};
    if ~isfield(solB, name) || isempty(solB.(name))
        continue
    end
    A = solA.(name);
    B = solB.(name);
    [nSp, nT, kA] = size(A);
    kB = size(B, 3);
    nBlocks = kA/npeA;
    if nBlocks ~= round(nBlocks) || kB ~= nBlocks*npeB
        error('SSIT:poolSSABatches:layout', ...
            ['Cannot pool ''%s'': %d and %d trajectories do not split into ', ...
            'the same number of blocks of %d and %d.'], name, kA, kB, npeA, npeB);
    end
    A = reshape(A, nSp, nT, npeA, nBlocks);
    B = reshape(B, nSp, nT, npeB, nBlocks);
    sol.(name) = reshape(cat(3, A, B), nSp, nT, (npeA+npeB)*nBlocks);
end
end
