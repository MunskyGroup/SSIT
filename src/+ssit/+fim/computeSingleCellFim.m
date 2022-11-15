function F = computeSingleCellFim(px, Sx, distortionOperator)
%Compute the Fisher Information Matrix associated with a single-cell
%distorted observation.
%
% Parameters
% ----------
%   px: :mat:class:`~+ssit.@FspVector.FspVector`
%       Distribution of single-cell responses.
%
%   Sx: 1-D array of :mat:class:`~+ssit.@FspVector.FspVector`
%       Partial derivatives of `px` with respect to parameters.
%
%   distortionOperator: :mat:class:`~+ssit.+pdo.@AbstractDistortionOperator.AbstractDistortionOperator`
%       Distortion operator. Could be left empty (`[]`) if no distortion.
arguments
    px ssit.FspVector;
    Sx (:,1) ssit.FspVector;
    distortionOperator;
end

parCountModel = length(Sx);
if ~isempty(distortionOperator)
    parCountPDO = size(distortionOperator.dCdLam,2);
else
    parCountPDO = 0;
end

F = zeros(parCountModel+parCountPDO, parCountModel+parCountPDO);

if ~isempty(distortionOperator)
    py = distortionOperator.computeObservationDist(px);

    Sy = ssit.FspVector.empty();
    for iPar = 1:parCountModel+parCountPDO
        if iPar<=parCountModel
            Sy(iPar) = distortionOperator.computeObservationDistDiff(px, Sx(iPar), iPar);
        else
            Sy(iPar) = distortionOperator.computeDiffPdoPx(px, [], iPar-parCountModel);
        end            
    end
else
    py = px;
    Sy = Sx;
end

for iPar = 1:parCountModel+parCountPDO
    for jPar = iPar:parCountModel+parCountPDO
%         si = Sy(iPar).data.values;
%         sj = Sy(jPar).data.values;
%         p = py.data.values;
        si = double(Sy(iPar).data);
        sj = double(Sy(jPar).data);
        p = double(py.data);
        si=padarray(si,max(0,size(p)-size(si)),'post');
        sj=padarray(sj,max(0,size(p)-size(sj)),'post');
        if isempty(si)||isempty(sj)
            F(iPar, jPar) = 0;
        else
            indsForSum = p>1e-9;
            F(iPar, jPar) = sum(si(indsForSum).*sj(indsForSum)./p(indsForSum));
        end
    end
    for jPar = 1:iPar-1
        F(iPar, jPar) = F(jPar, iPar);
    end
end
end

function isValid = unionType(x, varargs)
isValid = false;
for arg = varargs
    if isa(x, arg)
        isValid = true;
        return;
    end
end
end