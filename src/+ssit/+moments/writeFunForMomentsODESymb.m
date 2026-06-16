function [jacCreated] = writeFunForMomentsODESymb(S,wString,xString,parString, ...
    momentOdeFileName,includeSecondMom,inputExpressions,jacobianFileName)
%% function writeFunForMomentsODE(S,wString,x,momentOdeFileName,includeSecondMom)
% Creates ODE for solving 1st and 2nd moments using an Assumption of
% Gaussian Moments.
% Uses the stoichiometry, propensity functions and species to write out an
% ODE in matrix form that groups together constant, linear and higher
% order terms.
% Inputs:
%    S- stochiometry of reactions. This is a (n,m) vector, where 'n' is the
%       number of species and 'm' is the number of reactions. s(i,j) is the
%       stoichiometry for the i^th speacies in the j^th reaction.
%    w- propensity functions for the reactions. This is a (m,1) vector,
%       where 'm' is the number of reactions. Each element is a string
%       defining the propensity in terms of parameters and species
%       populationds (e.g. 'k1*mRNA', where 'k1' is a parameter and 'mRNA'
%       is a species.
%    x- Cell array containing species involved. This is a (n,1) cell array,
%    where 'n' is the number of species. Each element is a string for the
%    given species (e.g., 'mRNA').
%    parameterNames-- Cell array containing parameters involved. This is a
%    (m,1) cell array, where 'm' is the number of species. Each element is
%    a string for the given parameter (e.g., 'k1')

%    momentOdeFileName- desired name of ode function created for the
%                       moments
% Outputs: 
%    'momentOdeFileName'.m file with moment ode
%
% Example 1 (3 species, 6 rxn system)
% S = [1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]'; 
% syms x1 x2 x3 real
% w = [5; 2*x1; 4; 1*x2; 8; 0.5*x3];
% x = [x1;x2;x3];
% momentOdeFileName = 'tmpMomentOdeRHS';
%
% Example 2 (2 species,4 rxn system)
% S = [1 0; -1 0;0 1;0 -1]'; 
% syms x1 x2 real
% w = [50; x1; 12; 0.5*x2];
% x = [x1;x2];
% momentOdeFileName = 'tmpMomentOdeRHS';
%
% Example 3 (3 species,6 rxn system)
% S = [1 0 0; -1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]'; 
% syms x1 x2 x3 real
% w = [50; x1; 12; 0.5*x2; 3; 0.2*x3];
% x = [x1;x2;x3];
% momentOdeFileName = 'tmpMomentOdeRHS';
arguments
    S
    wString
    xString
    parString
    momentOdeFileName
    includeSecondMom = true
    inputExpressions =[];
    jacobianFileName = [momentOdeFileName,'_jac']
end

nS = size(S,1);   % number of species
nR = size(S,2);   % number of reactions
nP = length(parString); % number of parameters
nI = size(inputExpressions,1);

t = sym('t','real');
%% Create symbolic variables for species
x = sym('x',[nS,1],'real');

%% Create symbolic variables for parameters
ParameterX = sym('ParameterX',[1,nP],'real');

%% Change species names to x1, x2, ... in propensity functions.
specialFuns = {'max','min'}; nSF = 2;
containsSpecialFuns = false; 
for im = 1:nR
    % Change inputs to complete expressions in propensities.
    for iI = 1:nI
        wString{im} = regexprep(wString{im},['\<',inputExpressions{iI,1},'\>'],['(',inputExpressions{iI,2},')']);
    end
    % Change species to symbolic variables in propensities.
    for in = 1:nS
        wString{im} = regexprep(wString{im},['\<',xString{in},'\>'],['x',num2str(in)]);
    end
    % Change parameters to symbolic variables in propensities.
    for ip = 1:nP
        wString{im} = regexprep(wString{im},['\<',parString{ip},'\>'],['ParameterX',num2str(ip)]);
    end

    logs = {'<','>','='};
    k = 0;
    for il = 1:3
        j = strfind(wString{im},logs{il});
        while ~isempty(j)
            j1 = strfind(wString{im}(1:j(1)-1),'(');
            j1 = j1(end);
            nV = -1;j2=j1+1;
            while nV<0
                j2=j2+1;
                nV = -1 - sum(wString{im}(j1+1:j2)=='(') + sum(wString{im}(j1+1:j2)==')');
            end

            % j2 = strfind(wString{im}(j(1)+1:end),')');
            k=k+1;
            subStrings(k,1:2) = {['$',num2str(k)],['piecewise(',wString{im}(j1:j2),',1,0)']};
            wString{im} = strrep(wString{im},wString{im}(j1:j2),['$',num2str(k)]);                
            
            j = strfind(wString{im},logs{il});
        end

    end
    for ik = 1:k
        wString{im} = strrep(wString{im},['$',num2str(ik)],subStrings{ik,2});
    end

    % Change if function contains special functions that would make
    % differentiation incorrect.
    for isp = 1:nSF
        if ~isempty(regexp(wString{im},['\<',specialFuns{isp},'\>'],'once'))
            containsSpecialFuns = true;
        end
    end

    w(im,1) = str2sym(wString{im});
end

%%
nS = size(S,1);   % number of species
np = 3*ones(1,nS); % max order of polynomial for each species. The default is 3,
%                   but this will be expanded if needed.

if ~includeSecondMom
    %%  Create symbolic variables for all 1st and 2nd order moments.
    for i=1:nS % loops through each species
        % eval(['syms v_',num2str(i),' real']); % creates string of symbolic variables v1,v2...vn
        % eval(['Mu(',num2str(i),',1)=v_',num2str(i),';']); % creates string of means
        % eval(['v(',num2str(i),',1)=v_',num2str(i),';']); % creates string of all first moments
    end
    %% Define RHS of mean and 2nd moment ODEs in terms of monomials.
    RHS = S*w; % rhs of ode that describes the means
    %% Finalize RHS and saves as a matlab function.
    try
        matlabFunction(RHS,'Vars',{t,x,ParameterX},'File',momentOdeFileName,'Sparse',true); % save moment equation as a matlab function
    catch
        matlabFunction(RHS,'Vars',{t,x,ParameterX},'File',momentOdeFileName,'Sparse',false); % save moment equation as a matlab function
    end

    % matlabFunctionSSIT(RHS,{t,x,ParameterX},momentOdeFileName); % save moment equation as a matlab function

    try 
        if ~containsSpecialFuns
            jac = jacobian(RHS,x);
            try
                matlabFunction(jac,'Vars',{t,x,ParameterX},'File',jacobianFileName,'Sparse',true); % save moment equation as a matlab function
            catch
                matlabFunction(jac,'Vars',{t,x,ParameterX},'File',jacobianFileName,'Sparse',false); % save moment equation as a matlab function
            end
            % matlabFunctionSSIT(jac,{t,x,ParameterX},jacobianFileName); % save moment equation as a matlab function
            jacCreated = true;
        else
            if exist([jacobianFileName,'.m'],'file')
                delete([jacobianFileName,'.m'])
            end
            jacCreated = false;
        end
    catch
        delete(jacobianFileName)
        jacCreated = false;
    end
else
    %%  Create symbolic variables for all 1st and 2nd order moments.
    for i=1:nS % loops through each species
        eval(['syms v_',num2str(i),' real']); % creates string of symbolic variables v1,v2...vn
        eval(['Mu(',num2str(i),',1)=v_',num2str(i),';']); % creates string of means
        eval(['v(',num2str(i),',1)=v_',num2str(i),';']); % creates string of all first moments
    end
    k=nS;
    for i=1:nS
        for j=i:nS
            k = k+1;
            eval(['syms v_',num2str(k),' real']); % creates string of symoblic variables v1,v2...vn corresponding to second moments
            eval(['Sig(',num2str(i),',',num2str(j),',1)=v_',num2str(k),'-',... % creates string for diag and upper right of covariance matrix
                'v_',num2str(i),'*','v_',num2str(j),';']);
            eval(['Sig(',num2str(j),',',num2str(i),',1)=v_',num2str(k),'-',... % creates string for lower left of covariance matrix
                'v_',num2str(i),'*','v_',num2str(j),';']);
            eval(['v(',num2str(k),',1)=v_',num2str(k),';']); % creates string of means, variances and covariances
        end
    end

    %% Define RHS of mean and 2nd moment ODEs in terms of monomials.
    dxdt = S*w; % rhs of ode that describes the means
    dSigtdt = (S*w)*x' + x*(w'*S') + S*diag(w)*S'; % rhs of ode that describes the covariance matrix
    for i = nS:-1:1 % loop through last species to first species
        [C,T] = coeffs(dxdt(i),x,'all'); % gives the coefficients of the mean ode polynomial with respect to x where t are the terms the coefficient is multiplied by
        if nS==1; C=C'; end % if there is one species transpose C
        C(:) = C(end:-1:1); % creates a vector from c reversing the numbers
        T(:) = T(end:-1:1); % creates a vector from t reversing the numbers
        meanCoeffs{i} = C; % gives coefficients for the mean ode
        for j = nS:-1:i % loop through last species to first species
            [C,T] = coeffs(dSigtdt(i,j),x,'all'); % gives the coefficients of the covariance ode polynomial with respect to x where t are the terms the coefficient is multiplied by
            if nS==1; C=C'; end % if there is one species transpose C
            C(:) = C(end:-1:1); % creates a vector from c reversing the numbers
            T(:) = T(end:-1:1); % creates a vector from t reversing the numbers
            secMomCoeffs{i,j} = C; % gives coefficients for the covariance ode
            for kc=1:nS; sC(kc) = size(C,kc); end % gives size of C for each species
            np = max(np,sC); % finds max order of polynomial (+1)
        end
    end

    %% Form RHS of moment equation as linear equation A * [moments]
    npcell = num2cell(np); % make max order of polynomial a cell
    k=0;
    A = zeros(nS+nS*(nS+1)/2,prod(np),'sym'); % create base matrix to fill in
    for i = 1:nS % loop through each species
        C = meanCoeffs{i}; % call on coefficients from mean ode
        for kc=1:nS; sC(kc) = size(C,kc); end % gives size of C for each species
        if min(sC-np)<0 % if there are less coeff. than max order of polynomial
            C(npcell{:})=0; % then coeff. of max order terms are zero
        end
        J = find(C~=0); % gives vector of coefficients
        A(i,J) = C(J'); % place coefficients in A matrix
        for j = i:nS % loop through each species
            k = k+1;
            C = secMomCoeffs{i,j}; % call on coeff. from covariance ode
            for kc=1:nS; sC(kc) = size(C,kc); end % gives size of C for each species
            if nS>1&&min(sC-np)<0 % if its not species one and there are less coeff. than max order of polynomial
                C(npcell{:})=0; % then coeff. of max order terms are zero
            end
            J = find(C~=0); % gives vector of coefficients
            A(k+nS,J) = C(J'); % add in coefficients in A matrix
        end
    end

    %% Extract columns of A corresponding to 0th, 1st and 2nd, and higher order moments: A0 + A1*[v12] + A2*[vh]
    txt = '[';
    for i=1:nS
        txt = [txt,'ind(:,',num2str(i),'),']; % create three groups for A0,A1 and A2
    end
    if length(np)>1
        txt = [txt(1:end-1),']=ind2sub(np,[1:prod(np)]'');'];
    else
        txt = [txt(1:end-1),']=ind2sub([1,np],[1:prod(np)]'');'];
    end
    eval(txt) % create variable ind as matrix with 3 columns and the max order of the polynomials
    J0 = sum(ind-1,2)==0; % creates vector of ones and zeros with ones when sum of row is 0
    J1 = sum(ind-1,2)==1; % creates vector of ones and zeros with ones when sum of row is 1
    J2 = find(sum(ind-1,2)==2); % find row numbers of when sum is equal to two
    [~,JJ] = sortrows(ind(J2,:)); % pulls out rows given by J2 and reorders them
    J2 = J2(JJ(end:-1:1)); % gives row numbers of ind corresponding to new order from JJ
    Jhigh = sum(ind-1,2)>2; % gives vector of ones and zeros with ones when the sum of the row is greater than 2
    A0 = A(:,J0); % pulls out columns in A corresponding to when J0=1
    A1 = [A(:,J1),A(:,J2)]; % pulls out columns of A corresponding to when J1=1 and column numbers in J2
    A2 = A(:,Jhigh); % pulls out columns of A corresponding to when Jhigh=1

    %% Exclude zero columns from A2;
    inds_high = ind(Jhigh,:); % pulls out rows of ind present in Jhigh
    K = sum(abs(A2))~=0; % find where column sums are equal to zero in A2
    A2 = A2(:,K); % remove columns with zeros
    inds_high = inds_high(K,:); % removes corresponding rows in inds_high from ones removed in A2

    %% Perform moment closure to write higher moments in terms of 1st and 2nd
    try
        [nu] = getUncenteredMoments(inds_high-1,Mu,Sig); % Approximated 3rd and higher moments in terms of the 1st and 2nd moments.
        RHS = A0 + A1*v + A2*nu; % write ode with higher order approximation
    catch
        nu=[]; % if there are no higher order terms
        RHS = A0 + A1*v; % write ode with lower order terms only
    end
    %% Finalize RHS and saves as a matlab function.
    try
        matlabFunction(RHS,'Vars',{t,v,ParameterX},'File',momentOdeFileName,'Sparse',true); % save moment equation as a matlab function
    catch
        matlabFunction(RHS,'Vars',{t,v,ParameterX},'File',momentOdeFileName,'Sparse',true); % save moment equation as a matlab function
    end
    % matlabFunctionSSIT(RHS,{t,v,ParameterX},momentOdeFileName); % save moment equation as a matlab function
    
    try
        if ~containsSpecialFuns
            jac = jacobian(RHS,v);
            try
                matlabFunction(jac,'Vars',{t,v,ParameterX},'File',jacobianFileName,'Sparse',true); % save moment equation as a matlab function
            catch
                matlabFunction(jac,'Vars',{t,v,ParameterX},'File',jacobianFileName,'Sparse',true); % save moment equation as a matlab function
            end
            % matlabFunctionSSIT(jac,{t,v,ParameterX},jacobianFileName); % save moment equation as a matlab function
            jacCreated = true;
        else
            delete([jacobianFileName,'.m'])
            jacCreated = false;
        end
    catch
        delete([jacobianFileName,'.m'])
        jacCreated = false;
    end
end
