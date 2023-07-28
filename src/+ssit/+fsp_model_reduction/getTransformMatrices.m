function [phi,phi_inv,redOutputs] = getTransformMatrices(redType,n,fspSoln)
redOutputs=[];

if isfield(fspSoln,'fsp')
    nStates = size(fspSoln.A_total,1);
    nTimes = length(fspSoln.fsp);
    Solns = zeros(nStates,nTimes);
    for i=1:nTimes
        Solns(1:length(fspSoln.fsp{i}.p.data.vals),i) = fspSoln.fsp{i}.p.data.vals;
    end
elseif isfield(fspSoln,'fullSolutionsNow')
    Solns = fspSoln.fullSolutionsNow';
    fspSoln.stateSpace.states = fspSoln.states;
end

switch redType
    case 'No Transform'
        phi = eye(size(fspSoln.A_total,1));
        phi_inv = phi';
    case 'Eigen Decomposition'
        [phi,~] = eigs(fspSoln.A_total,n,0);
        phi_inv = phi';
    case 'Eigen Decomposition Initial'
        [phi,D] = eigs(fspSoln.A_total,n,0);
        [~,I] = sort(real(diag(D)),'descend');
        phi = phi(:,I(1:n));
        phi = orth([Solns(:,1),phi]);
        phi_inv = phi';
    case 'Linear State Lumping'
        nStates = size(fspSoln.stateSpace.states,2);
        spmax=max(fspSoln.stateSpace.states,[],2);
        nSpecies = size(fspSoln.stateSpace.states,1);
        for i = nSpecies:-1:1
            bins{i} = unique(floor(linspace(1,spmax(i)+1,n+1)))-1;
        end

        phi_map = zeros(nStates,nSpecies);
        for j=1:nStates
            for i=1:nSpecies
                phi_map(j,i) = find(fspSoln.stateSpace.states(i,j)<=bins{i},1,"first");
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
        nStates = size(fspSoln.stateSpace.states,2);
        spmax=max(fspSoln.stateSpace.states,[],2);
        nSpecies = size(fspSoln.stateSpace.states,1);
        for i = nSpecies:-1:1
            m=n;
            bins{i} = [0,unique(ceil(logspace(0,log10(spmax(i)),m)))];
            while length(bins{i})<min(n+1,spmax(i)+1)
                m=ceil(m*1.1);
                bins{i} = [0,unique(ceil(logspace(0,log10(spmax(i)),m)))];
            end
        end

        phi_map = zeros(nStates,nSpecies);
        for j=1:nStates
            for i=1:nSpecies
                phi_map(j,i) = find(fspSoln.stateSpace.states(i,j)<=bins{i},1,"first");
            end
        end
        binns = max(phi_map);
        cprod = [1,cumprod(binns(1:end-1))]';
        phi_inds = (phi_map-1)*cprod+1;

        phi = sparse(nStates,binns(end)*cprod(end));
        phi(sub2ind(size(phi),(1:nStates)',phi_inds))=1;
%         phi = sparse(orth(full(phi)));
%         phi_inv = phi';
        phi_inv = phi';
        phi = phi./sum(phi,1);

    case 'Balanced Model Truncation (HSV)'
        nStates = size(fspSoln.states,2);
        sys = ss(fspSoln.A_total,fspSoln.P0,eye(nStates),[]);
        [~,redOutputs.info] = balred(sys);
        phi = [];phi_inv = [];

    case 'Proper Orthogonal Decomposition'

         [phi,D,~] = svds(Solns,n,"largest","Tolerance",1e-18);
%         [phi,D,~] = svds(Solns,n,0);
%         [~,I] = sort(real(diag(D)),'descend');
%         phi = phi(:,I(1:n));
        %         [phi,~,~] = svds(fspSoln.fullSolutionsNow',n);
        phi = orth([Solns(:,1),phi]);
        phi_inv = phi';

    case 'POD 2nd'

         [phi,D,~] = svds([Solns,Solns(:,2:end)-Solns(:,1:end-1)],n,"largest","Tolerance",1e-18);
%         [phi,D,~] = svds(Solns,n,0);
%         [~,I] = sort(real(diag(D)),'descend');
%         phi = phi(:,I(1:n));
        %         [phi,~,~] = svds(fspSoln.fullSolutionsNow',n);
        phi = orth([Solns(:,1),phi]);
        phi_inv = phi';

        
    case 'Dynamic Mode Decomposition'

        V1 = Solns(:,1:end-1);
        V2 = Solns(:,2:end);
        [Ured,Sigred,Wred] = svds(V1,n);
        S = Ured'*V2*Wred*(Sigred^-1);
        [y,~] = eig(S);
        phi = Ured*real(y);
        phi_inv = phi';
    case 'Radial Basis Functions'
        [phi,phi_inv] = ssit.fsp_model_reduction.radiaBasisPhi(fspSoln.states,n,0,1);
end