function [PHI,PHI_inv] = radiaBasisPhi(states,nRed,oth,SampleD)
%% Compute PHI
% Here I will use log-distributed RBFs, where each RBF is a multivariate
% gaussian with variance proportional to the position of the RBF.
nSpecies = size(states,1);
nStates = size(states,2);

% Define the centers and widths of all radial bases.
% for i=nSpecies:-1:1
%     center{i} = unique(floor(logspace(-1,log10(nmax(i)),Nred(i))'));
%     width{i} = [1;center{i}(2:end)-center{i}(1:end-1)];
%     m(i,1) = length(center{i});
% end
centers = states(:,unique([1,randperm(nStates,nRed)]));
width = 0*centers;
for i=1:nRed+1
    [~,j] = min(sum(centers(:,i)-centers(:,[1:nRed+1]~=i)).^2); 
    if j<i
        width(:,i) = max(1,abs(centers(:,i)-centers(:,j)));
    else
        width(:,i) = max(1,abs(centers(:,i)-centers(:,j+1)));
    end
end
% for i=1:m(2)
%     centers((i-1)*m(1)+1:i*m(1),:) = [centerx,centery(i)*ones(m(1),1)];
%     width((i-1)*m(1)+1:i*m(1),:) =  [widthx,widthy(i)*ones(m(1),1)];
% end

% X = repmat([0:nmax(1)]',nmax(2)+1,1);
% Y = reshape(repmat([0:nmax(2)],nmax(1)+1,1),prod(nmax+1),1);

% Nphi = length(center(:,1));
PHI = ones(nStates,nRed);
for i=1:nRed
    for j = 1:nSpecies
        PHI(:,i) = PHI(:,i).*exp(-(states(j,:)'-centers(j,i)).^2./max(1,width(j,i)^2));
    end
end
% whos
% W = SampleD*PHI;
% PHI = PHI*diag(W+1e-4);
% whos

PHI(PHI<1e-7)=0;
PHI = sparse(PHI);

if oth==1
    PHI = orth(full(PHI));    
end

if nargout==2
    % PHI_inv = (PHI'*PHI)^-1*PHI';
    % PHI_inv(abs(PHI_inv)<1e-6)=0;
    % PHI_inv=sparse(PHI_inv);
    if oth==1
        PHI_inv = PHI';
    else
        PHI_inv = (PHI'*PHI)^-1*PHI';
    end
end