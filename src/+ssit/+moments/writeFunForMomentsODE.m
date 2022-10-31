function writeFunForMomentsODE(S,w,x,momentOdeFileName)
%% function writeFunForMomentsODE(S,w,x,momentOdeFileName)
% creates ode for solving moments using the Moment Generating Function for
% Gaussian Distributions.
% Uses the stoichiometry, propensity functions and species to write out an
% ode in matrix form that groups together constant, linear and higher
% order terms.
% Inputs:
%    S- stochiometry of reactions. This is a (n,m) vector, where 'n' is the
%       number of species and 'm' is the number of reactions. s(i,j) is the
%       stoichiometry for the i^th speacies in the j^th reaction.
%    w- propensity functions for the reactions. This is a (m,1) vector,
%       where 'm' is the number of reactions.
%    x- species involved. This is a (n,1) vectorm where 'n' is the number
%       of species
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

%%
n = size(S,1);   % number of species
np = 3*ones(1,n); % max order of polynomial (+1)

%%  Create symbolic variables for all 1st and 2nd order moments.
for i=1:n % loops through each species
    eval(['syms v_',num2str(i),' real']); % creates string of symbolic variables v1,v2...vn
    eval(['Mu(',num2str(i),',1)=v_',num2str(i),';']); % creates string of means
    eval(['v(',num2str(i),',1)=v_',num2str(i),';']); % creates string of all first moments
end
k=n;
for i=1:n
    for j=i:n
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
for i = n:-1:1 % loop through last species to first species
    [C,T] = coeffs(dxdt(i),x,'all'); % gives the coefficients of the mean ode polynomial with respect to x where t are the terms the coefficient is multiplied by
    if n==1; C=C'; end % if there is one species transpose C
    C(:) = C(end:-1:1); % creates a vector from c reversing the numbers
    T(:) = T(end:-1:1); % creates a vector from t reversing the numbers
    meanCoeffs{i} = C; % gives coefficients for the mean ode
    for j = n:-1:i % loop through last species to first species
        [C,T] = coeffs(dSigtdt(i,j),x,'all'); % gives the coefficients of the covariance ode polynomial with respect to x where t are the terms the coefficient is multiplied by
        if n==1; C=C'; end % if there is one species transpose C
        C(:) = C(end:-1:1); % creates a vector from c reversing the numbers
        T(:) = T(end:-1:1); % creates a vector from t reversing the numbers
        secMomCoeffs{i,j} = C; % gives coefficients for the covariance ode
        for kc=1:n; sC(kc) = size(C,kc); end % gives size of C for each species
        np = max(np,sC); % finds max order of polynomial (+1)
    end
end

%% Form RHS of moment equation as linear equation A * [moments]
npcell = num2cell(np); % make max order of polynomial a cell
k=0;
A = zeros(n+n*(n+1)/2,prod(np)); % create base matrix to fill in
for i = 1:n % loop through each species
    C = meanCoeffs{i}; % call on coefficients from mean ode
    for kc=1:n; sC(kc) = size(C,kc); end % gives size of C for each species
    if min(sC-np)<0 % if there are less coeff. than max order of polynomial
        C(npcell{:})=0; % then coeff. of max order terms are zero
    end
    J = find(C); % gives vector of coefficients
    A(i,J) = C(J'); % place coefficients in A matrix
    for j = i:n % loop through each species
        k = k+1;
        C = secMomCoeffs{i,j}; % call on coeff. from covariance ode
        for kc=1:n; sC(kc) = size(C,kc); end % gives size of C for each species
        if n>1&&min(sC-np)<0 % if its not species one and there are less coeff. than max order of polynomial
            C(npcell{:})=0; % then coeff. of max order terms are zero
        end
        J = find(C); % gives vector of coefficients
        A(k+n,J) = C(J'); % add in coefficients in A matrix
    end
end

%% Extract columns of A corresponding to 0th, 1st and 2nd, and higher order moments: A0 + A1*[v12] + A2*[vh]
txt = '[';
for i=1:n
    txt = [txt,'ind(:,',num2str(i),'),']; % create three groups for A0,A1 and A2
end
txt = [txt(1:end-1),']=ind2sub(np,[1:prod(np)]'');'];
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
G = matlabFunction(RHS,'Vars',{[v]},'File',momentOdeFileName); % save moment equation as a matlab function
