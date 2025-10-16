function [momOut] = getUncenteredMoments(s,Mu,Sig)
%% function [Mom] = getUncenteredMoments(s,Mu,Sig)
% Moment Generating Function for Gaussian Distributions
% https://en.wikipedia.org/wiki/Multivariate_normal_distribution
% This function uses gaussian moment closure to write all moments of a
% multivariate gaussian distribution in terms of known means (Mu) and
% covariances (Sig).
% Inputs:
%   s - order of requested moments.  This is an (m,n) vector, where 'n' is
%       the number variables and 'm' is the number of different moments 
%       requested by the user. s(i,j) is the order of the j^th variable in
%       the i^th requested moment.
%   Mu - vector of means (numerical)
%   Sig - covariance matrix (numerical)
% Outputs:
%   Mom(i) = E{x1^s(i,1)*...*xn^s(i,n)}; % Vector of outputs in terms of
%   means and covariances.
% Example 1 (numerial results):
%   s = [1 2];
%   MU = [1;3]; SIGMA = [1 0.1;0.1 2];
%   [mom_mgf] = getUncenteredMoments(s,MU,SIGMA)
%   R = mvnrnd(MU,SIGMA,10000);
%   mom_samples = sum(R(:,1).*R(:,2).^2)/size(R,1)
%
% Example 2 (symbolic results):
%   s = [1 0 0;0,1,0;2,0,1;1,1,1;0,2,0];
%   [mom_mgf] = getUncenteredMoments(s)

%% Create nx1 Mean vector (Mu), nxn Covariance matrix (Sig), and nx1 variable (t)
n = size(s,2); % number of variables;
if nargin==1 % if mu and sig are not inputed
    for i=1:n % loop through each variable
        eval(['syms mu',num2str(i),' real']); % create symbolic variable mu when values not given
        eval(['Mu(',num2str(i),',1)=mu',num2str(i),';']); 
        eval(['syms t',num2str(i)]); % create symbolic variable t when final equation in symbolic
        eval(['t(',num2str(i),',1)=t',num2str(i),';']);
        for j=i:n
            eval(['syms sig',num2str(i),'_',num2str(j),' real']); % create symoblic variable sig when values not given
            eval(['Sig(',num2str(i),',',num2str(j),',1)=sig',num2str(i),'_',num2str(j),';']);
            eval(['Sig(',num2str(j),',',num2str(i),',1)=sig',num2str(i),'_',num2str(j),';']);
        end
    end
else  % Only need symbolic expressions for the variables 't1' to 'tn'.
    for i=1:n % loop through each variable
        eval(['syms t',num2str(i)]); % create symoblic variable t for solving equation
        eval(['t(',num2str(i),',1)=t',num2str(i),';']);
    end
end

%% Create Moment Generating Function
MGF = exp(Mu'*t + 1/2*t'*Sig*t);

m = size(s,1); % number of moments requested by user
for j = 1:m % loop through for each moment
    %% Compute uncentered moment
    mom = MGF; % call on moment generating function
    for i=1:n % loop through each variable
        eval(['mom = diff(mom,t',num2str(i),',',num2str(s(j,i)),');']) % takes derivative of moment generating function with respect to t
        eval(['mom = subs(mom,t',num2str(i),',0);']) % sub values into equation
    end
    
    if nargin>1&&~isa(Mu,'sym') % used if mu and sig values are given
        momOut(j) = eval(mom);  % Returns a numerical value
    else % if numerical values are not given for mu and sig
        eval(['momOut(',num2str(j),',1) = subs(mom,t',num2str(i),',0);']); % gives a symoblic answer
    end
end
