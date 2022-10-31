function [phi,phi_inv,redOutputs] = getTransformMatrices(redType,n,fspSoln)
redOutputs=[];
switch redType
    case 'Eigen Decomposition'
        [phi,~] = eigs(fspSoln.A_total,n,'largestreal');
        phi = orth(phi);
        phi_inv = phi';
    case 'Linear State Lumping'
        nStates = size(fspSoln.states,2);
        spmax=max(fspSoln.states,[],2);
        nSpecies = size(fspSoln.states,1);
        for i = nSpecies:-1:1
            bins{i} = unique(floor(linspace(1,spmax(i)+1,n+1)))-1;
        end

        phi_map = zeros(nStates,nSpecies);
        for j=1:nStates
            for i=1:nSpecies
                phi_map(j,i) = find(fspSoln.states(i,j)<=bins{i},1,"first");
            end
        end
        binns = max(phi_map);
        cprod = [1,cumprod(binns(1:end-1))]';
        phi_inds = (phi_map-1)*cprod+1;

        phi = zeros(nStates,binns(end)*cprod(end));
        phi(sub2ind(size(phi),(1:nStates)',phi_inds))=1;
        phi = orth(phi);
        phi_inv = phi';
    case 'Logarithmic State Lumping'
        nStates = size(fspSoln.states,2);
        spmax=max(fspSoln.states,[],2);
        for i = length(spmax):-1:1
            m=n;
            bins{i} = unique(floor(logspace(0,log10(spmax(i)+1),m)));
            while length(bins{i})<min(n+1,spmax(i)+1)
                m=ceil(m*1.1);
                bins{i} = unique(floor(logspace(0,log10(spmax(i)+1),m)));
            end
        end

        phi_map = zeros(nStates,3);
        for j=1:nStates
            for i=1:3
                phi_map(j,i) = find(fspSoln.states(i,j)<=bins{i},1,"first");
            end
        end
        binns = max(phi_map);
        cprod = [1,cumprod(binns(1:end-1))]';
        phi_inds = (phi_map-1)*cprod+1;

        phi = zeros(nStates,cprod(end));
        phi(sub2ind(size(phi),(1:nStates)',phi_inds))=1;
        phi = orth(phi);
        phi_inv = phi';

    case 'Balanced Model Truncation (HSV)'
        nStates = size(fspSoln.states,2);
        sys = ss(fspSoln.A_total,fspSoln.P0,eye(nStates),[]);
        [~,redOutputs.info] = balred(sys);
        phi = [];phi_inv = [];

    case 'Proper Orthogonal Decomposition'
        [phi,~,~] = svds(fspSoln.fullSolutionsNow',n);
%         phi = orth([fspSoln.P0,phi]);
        phi_inv = phi';

    case 'Dynamic Mode Decomposition'
        V1 = fspSoln.fullSolutionsNow(1:end-1,:)';
        V2 = fspSoln.fullSolutionsNow(2:end,:)';
        [Ured,Sigred,Wred] = svds(V1,n);
        S = Ured'*V2*Wred*(Sigred^-1);
        [y,~] = eig(S);
        phi = Ured*real(y);
        phi_inv = phi';
    case 'Radial Basis Functions'
        [phi,phi_inv] = ssit.fsp_model_reduction.radiaBasisPhi(fspSoln.states,n,0,1);
end