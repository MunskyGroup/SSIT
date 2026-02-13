function [phi,phi_inv,phiScale,phiPlot,redOutputs] = getTransformMatrices(redOptions,fspSoln,phi)
arguments
    redOptions
    fspSoln
    phi =[];
end
redType = redOptions.reductionType;
redOrder = redOptions.reductionOrder;

redOutputs=[];
phiScale = [];
phiPlot = [];

SolnNeeded = {'Proper Orthogonal Decomposition','POD 2nd','POD Update',...
    'Dynamic Mode Decomposition','Eigen Decomposition Initial',...
    'POD','DMD','POD2','ED'};
if isfield(fspSoln,'fsp')
    nStates = size(fspSoln.A_total,1);
    nTimes = length(fspSoln.fsp);
    
    % Sort the FSP solution into the right order for the statespace.
    if max(contains(SolnNeeded,redType))
        Solns = zeros(nStates,nTimes);
        % inds = state2key(fspSoln.fsp{end}.p.data.subs'-1);
        % inds2 = fspSoln.stateSpace.state2indMap(inds)
        % inds2 = zeros(1,length(inds));
        % for iv=1:length(inds)
        %     try
        %         inds2(iv) = fspSoln.stateSpace.state2indMap(inds(iv));
        %     catch
        %         1+1
        %     end
        % end
        for i=1:nTimes
            % try
                inds = state2key(fspSoln.fsp{i}.p.data.subs'-1);
                inds2 = fspSoln.stateSpace.state2indMap(inds);
                Solns(inds2(1:length(fspSoln.fsp{i}.p.data.vals)),i) = fspSoln.fsp{i}.p.data.vals;
            % catch
            %     1+1
            % end
        end  
    end
elseif isfield(fspSoln,'fullSolutionsNow')
    Solns = fspSoln.fullSolutionsNow';
    fspSoln.stateSpace.states = fspSoln.states;
end

switch redType
    case 'No Transform'
        %'No Transform' --
        %   No model reduction is used.  This is meant only for testing.
        phi = eye(size(fspSoln.A_total,1));
        phi_inv = phi';
    case {'Linear State Lumping','LNSL'}
        %'Linear State Lumping' (LNSL) -- 
        %   State space is divided into linearly distributed bins.
        %   Then within each bin, the proability distribution is assumed to
        %   be constant. The number of bins for each species is defined by
        %   redOrder.
        nStates = size(fspSoln.stateSpace.states,2);
        spmax=max(fspSoln.stateSpace.states,[],2);
        nSpecies = size(fspSoln.stateSpace.states,1);
        for i = nSpecies:-1:1
            bins{i} = unique(floor(linspace(1,spmax(i)+1,redOrder+1)))-1;
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
    case {'Logarithmic State Lumping','LGSL'}
        %'Logarithmic State Lumping' -- 
        %   State space is divided into logarithmically distributed bins.
        %   Then within each bin, the proability distribution is assumed to
        %   be constant. The number of bins for each species is defined by
        %   redOrder.
        nStates = size(fspSoln.stateSpace.states,2);
        spmax=max(fspSoln.stateSpace.states,[],2);
        nSpecies = size(fspSoln.stateSpace.states,1);

        % define bins
        for i = nSpecies:-1:1
            m=redOrder;
            bins{i} = [0,unique(ceil(logspace(0,log10(spmax(i)),m)))];
            while length(bins{i})<min(redOrder+1,spmax(i)+1)
                m=ceil(m*1.1);
                bins{i} = [0,unique(ceil(logspace(0,log10(spmax(i)),m)))];
            end
        end

        % map each state to its corresponding bin
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

        % Remove un-needed lumped states that are not represented in the FSP state
        % space.
        nonEmptyColumns = sum(phi,1)~=0;
        phi = phi(:,nonEmptyColumns);

        phi_inv = phi';
        phi = phi./sum(phi,1);

    case {'Proper Orthogonal Decomposition','POD'}
        %'Proper Orthogonal Decomposition' --
        %   The infintesimal generator is projected onto the orthonorml
        %   basis spanned by a previous solution of the full master
        %   equation at discrete time points.  this is done by perfoming
        %   SVD on the solutions and choosing the output space
        %   corresponding to the 'redOrder' largest singular values.
        if size(Solns,2)>40
            J = unique(round(linspace(1,size(Solns,2),40)));
            Solns = Solns(:,J);
        end
        [phi,D,~] = svds(Solns,redOrder,"largest","Tolerance",1e-18);
        J = find(cumsum(diag(D))<=0.9999*sum(diag(D)),1,"last");
        phi = orth([Solns(:,1),phi(:,1:J)]);
        phi_inv = phi';
    case 'POD Update'
        %'POD Update' --
        %   The infintesimal generator is UPDATED onto the orthonorml
        %   basis spanned by (1) current POD projection, and (2) an NEW
        %   solution of the full master equation at discrete time points.
        if size(Solns,2)>40
            J = unique(round(linspace(1,size(Solns,2),40)));
            Solns = Solns(:,J);
        end
        [phi,D,~] = svds([phi,Solns],redOrder,"largest","Tolerance",1e-18);
        J = find(cumsum(diag(D))<=0.9999*sum(diag(D)),1,"last");
        phi = orth([Solns(:,1),phi(:,1:J)]);
        phi_inv = phi';
    case {'Log Lump QSSA','LGQSSA'}
        %'Log Lump QSSA' -- 
        %   State space is divided into logarithmically distributed bins.
        %   Then within each bin, a quasi-steady-state-assumption is made
        %   to estimate the distribution of proability mass within that
        %   bin. The number of bins for each species is defined by
        %   redOrder.
        nStates = size(fspSoln.stateSpace.states,2);
        spmax=max(fspSoln.stateSpace.states,[],2);
        nSpecies = size(fspSoln.stateSpace.states,1);

        % define bins
        for i = nSpecies:-1:1
            m=redOrder;
            bins{i} = [0,unique(ceil(logspace(0,log10(spmax(i)+1),m)))];
            while length(bins{i})<min(redOrder+1,spmax(i)+1)
                m=ceil(m*1.1);
                bins{i} = [0,unique(ceil(logspace(0,log10(spmax(i)+1),m)))];
            end
        end

        % create map from state to corresponding bins
        phi_map = zeros(nStates,nSpecies);
        for j=1:nStates
            for i=1:nSpecies
                phi_map(j,i) = find(fspSoln.stateSpace.states(i,j)<=bins{i},1,"first");
            end
        end
        binns = max(phi_map);
        cprod = [1,cumprod(binns(1:end-1))]';
        phi_inds = (phi_map-1)*cprod+1;

        % find shape function assuming QSS within each bin
        phiVals = 0*phi_inds;
        nLump = max(phi_inds);
        for iLump = 1:nLump
            J = phi_inds==iLump;
            Alump = fspSoln.A_total(J,J);
            if size(Alump,2)>=2
                Alump = Alump - diag(sum(Alump));
                [qssa,~] = eigs(Alump,1,'largestreal');
                qssa = max(1e-6,qssa/sum(qssa));
            else
                qssa = 1;
            end
            phiVals(phi_inds==iLump) = qssa;
        end
        phi = sparse((1:length(phi_inds)),phi_inds,phiVals,length(phi_inds),max(phi_inds));
        
        % Remove un-needed lumped states that are not represented in the FSP state
        % space.
        nonEmptyColumns = sum(phi,1)~=0;
        phi = phi(:,nonEmptyColumns);
        phi_inv = 1.0*(phi'>=1e-12);

        % plotting shape (constant within bin)
        phiPlot = phi>=1e-12;

        % find bin centers
        centers = cell(nSpecies,1);
        for iSpec = 1:nSpecies
            centers{iSpec} = bins{iSpec}(1:end) + ([bins{iSpec}(2:end),bins{iSpec}(end)+1]-[bins{iSpec}(1:end)])/2 - 0.5;
        end

        % find distance between mesh centers
        uniqueBins = unique(phi_map,'rows');
        [~,J] = sort((uniqueBins-1)*cprod+1,'ascend');
        uniqueBins = uniqueBins(J,:);
        % centersGrid = zeros(nLump,nSpecies);
        % TODO - Check this change.
        centersGrid = zeros(size(uniqueBins,1),nSpecies);

        for iSpec = 1:nSpecies
            centersGrid(:,iSpec) = centers{iSpec}(uniqueBins(:,iSpec));
        end
        phiScale = (1./squareform(pdist(centersGrid)));
        phiScale(isinf(phiScale))=0;

    case {'Eigen Decomposition','ED'}
        %'Eigen Decomposition' --
        %   The infintesimal generator is projected onto the eigenvectors
        %   corresponding to its least negative eigenvalues. The number of
        %   eigenvectors for this projection is defined by 'redOrder'.  
        [phi,~] = eigs(fspSoln.A_total,redOrder,0);
        phi_inv = phi';
    case 'Balanced Model Truncation (HSV)'
        %'Balanced Model Truncation (HSV)' -- 
        %   The projection is chosen according to the Henkel Singular
        %   Values of the linear systems (A, B, C, []), where A is the
        %   infinitesimal generator, B=P0 is the initial probability
        %   distribution, and C=eye is the identity matrix (corresponding
        %   to perfect observations).
        %
        error('TODO - Balanced Truncation Projection has been deactivated.')
        % TODO - this might work much better when C is chosen as PDO that
        % acounts for the loss of observable species.
        % P0 = zeros(size(fspSoln.A_total,1),1);
        % P0(:) = fspSoln.fsp{1}.p.data;
        % nStates = size(fspSoln.stateSpace.states,2);
        % sys = ss(fspSoln.A_total,P0,eye(nStates),[]);
        % [~,redOutputs.info] = balred(sys);
        % phi = [];phi_inv = [];

    case 'QSSA'
        %'QSSA' -- 
        %   Species specified in 'redOptions.qssaSpecies' are assumed to be
        %   in quasi-steady-state assuming that all other species are held
        %   constant. Reaction rates for remaining species are projected
        %   onto a lower dimensional state space. 
        nStates = size(fspSoln.stateSpace.states,2);
        spmax=max(fspSoln.stateSpace.states,[],2);
        nSpecies = size(fspSoln.stateSpace.states,1);

        qssaSpecies = redOptions.qssaSpecies;

        % define bins
        for i = nSpecies:-1:1
            if ismember(qssaSpecies,i)
                bins{i} = [0,spmax(i)+1];
            else
                bins{i} = [0:spmax(i)];
            end
        end

        % create map from state to corresponding bins
        phi_map = zeros(nStates,nSpecies);
        for j=1:nStates
            for i=1:nSpecies
                phi_map(j,i) = find(fspSoln.stateSpace.states(i,j)<=bins{i},1,"first");
            end
        end
        binns = max(phi_map);
        cprod = [1,cumprod(binns(1:end-1))]';
        phi_inds = (phi_map-1)*cprod+1;

        % find shape function assuming QSS within each bin
        phiVals = 0*phi_inds;
        nLump = max(phi_inds);
        for iLump = 1:nLump
            J = phi_inds==iLump;
            Alump = fspSoln.A_total(J,J);
            if size(Alump,2)>=2
                Alump = Alump - diag(sum(Alump));
                [qssa,~] = eigs(Alump,1,'largestreal');
                qssa = max(1e-6,qssa/sum(qssa));
            else
                qssa = 1;
            end
            phiVals(phi_inds==iLump) = qssa;
        end
        phi = sparse((1:length(phi_inds)),phi_inds,phiVals,length(phi_inds),max(phi_inds));
        phi_inv = 1.0*(phi'>=1e-12);

        % plotting shape 
        phiPlot = phi;

    case {'POD 2nd','POD2'}
        %'POD 2nd' --
        %   The infintesimal generator is projected onto the orthonorml
        %   basis spanned by a previous solution of the full master
        %   equation at discrete time points AND the approximate time
        %   derivative of these solutions. This is done by perfoming
        %   SVD on the solutions and choosing the output space
        %   corresponding to the 'redOrder' largest singular values.
        %
        % Ds = Solns(:,2:end)-Solns(:,1:end-1);
        % Ds = Ds./sum(abs(Ds));
        % [phi,~,~] = svds([Solns,fspSoln.A_total*Solns],redOrder,"largest","Tolerance",1e-18);
        [phi,~,~] = svds([Solns,fspSoln.A_total*Solns],redOrder,"largest","Tolerance",1e-18);
        phi = orth([Solns(:,1),phi]);
        phi_inv = phi';

    case {'Dynamic Mode Decomposition','DMD'}
        %'Dynamic Mode Decomposition' --
        V1 = Solns(:,1:end-1);
        V2 = Solns(:,2:end);
        [Ured,Sigred,Wred] = svds(V1,redOrder);
        S = Ured'*V2*Wred*(Sigred^-1);
        [y,~] = eig(S);
        phi = Ured*real(y);
        phi_inv = phi';
    case {'Radial Basis Functions','RBF'}
        error('TODO - Radial Basis Function projection not yet functional')
        % %'Radial Basis Functions' --
        % %  MAster equation is projected onto a set of radial basis
        % %  functions that are randomly placed throughout the statespace.
        % [phi,phi_inv] = ssit.fsp_model_reduction.radiaBasisPhi(fspSoln.stateSpace.states,redOrder,0,1);
end
end
function keys =  state2key( states )
% Hash function to convert a N-dimensional integer vector into a unique
% string "i1  i2  i3 ..."

% N = size(states, 2);
% keys = cell(1,N);
% 
% for n = 1:N
%     str = num2str(states(:,n)');    
%     keys{n} = str;
% end
N = size(states, 2);
keys = mat2cell(uint64(states)',ones(1,N))';

end